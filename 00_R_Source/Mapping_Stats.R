####-
##-====================================================
##-========================  MAPPING STATS ========================
##-====================================================
####-


#' Calculate Mapping Statistics
#'
#' Calculate genome coverage, exon coverage and more.
#'
#' IMPORTANT: All input vectors must have the same order, i.e. bamFiles[i] must
#' correspond to numReads[i] (or numReads[2*i-1] for paired-end reads).
#' \cr\cr
#' Coverage: The coverage of a reference sequence is defined as the average
#' number of times each base has been mapped.
#' \cr\cr
#' Estimation: Given either (i) the absolute number of mapped nucleotides or
#' (ii) the average length of mapped reads, the estimated coverage is defined
#' as
#' \cr\cr
#' #nt_mapped / #nt_reference
#' \cr\cr
#' or
#' \cr\cr
#' (#nr_reads * avr_len) / #nt_reference
#' \cr\cr
#' The genome coverage is calculated based on (i) and the exon coverage is
#' calculated based on (ii).
#'
#' @export
#' @param bamFiles Character(n). Path to one or more BAM files.
#' @param fqFiles Character(m). Path to raw FASTQ files.
#'   For 'paired=TRUE', this vector must be given with two
#'   files per BAM file: m = 2*n.
#' @param anno Character(1). Path to annotation file. GTF or GFF format.
#'   Exon information is acquired from this file.
#' @param samples Character(n). Sample names used as row names. Defaults to NA.
#' @param numReads Numeric(m). The original number of reads before trimming
#'   and mapping.\cr
#'   If NA, this will be calculated from 'fqFiles', but can
#'   take a considerable amount of time (for large files).
#'   m = 2*n if 'paired=TRUE'.
#' @param numTrimmed Numeric(m). Optional. The number of reads per file, which
#'   remained after trimming. Defaults to NA.
#' @param numNonrRNA Numeric(m). Optional. The number of reads per file, which
#'   remained after rRNA removal. Defaults to NA.
#' @param libSizes Numeric(n). Number of reads per file, which are mapping in
#'   EXONS only. Only required if 'precise = FALSE'. After counting, libsizes is
#'   equal to colSums(counts) with genes in rows and samples in columns.
#'   Defaults to NA.
#' @param featureType Consider only this feature. Usually "exon" or "gene". If
#'   NULL or NA, every feature is considered. Defaults to 'exon'.
#' @param paired Logical(1). If FALSE, assume single-end reads in BAM files. If
#'   TRUE, assume paired-end reads in BAM files. Pay attention to how 'fqFiles'
#'   must be set to correctly work with paired-end BAM files.
#' @param precise Logical(1). If FALSE, exon and genome coverage are estimated.
#'   If TRUE, a precise calculation is used to determine exon coverage. Precise
#'   requires indexed BAM files. Index files must be in the same directory as
#'   the corresponding BAM files. Defaults to FALSE.
#' @param remove.na If TRUE, unused columns are removed. Else, they are kept.
#' @param autoIndex Applicable only if 'precise = TRUE': Logical indicating
#'   whether BAM files should be indexed automatically. If the BAM files are
#'   unsorted, they are sorted before indexing. Defaults to TRUE. If set to
#'   FALSE, execution terminates if index files cannot be found.
#' @param allowMulti Logical. If True, secondary alignments are not considered.
#'   Defaults to FALSE.
#' @param cpus Integer(1) > 0. Maximum number of CPU cores (threads) to use.
#' @return A data.frame.
#' @examples
#' \dontrun{
#' # for single-end FASTQ & BAM files
#' bamFiles <- "f1.bam"
#' fqFiles <- "f1.fastq"
#' anno <- "anno.gtf"
#' libSizes <- c(1000)
#' calc_mapping_stats(bamFiles, fqFiles, anno, libSizes)
#'
#' # for paired-end FASTQ & BAM files
#' bamFiles <- "f1.bam"
#' fqFiles <- c("f1_1.fastq", "f1_2.fastq")
#' anno <- "anno.gtf"
#' libSizes <- c(1000)
#' calc_mapping_stats(bamFiles, fqFiles, anno, libSizes, paired = TRUE)
#' }
#' calc_mapping_stats
calc_mapping_stats <- function
(
    bamFiles,
    fqFiles,
    anno,
    samples = NA,
    numReads = NA,
    numTrimmed = NA,
    numNonrRNA = NA,
    libSizes = NA,
    featureType = "exon",
    paired = FALSE,
    precise = FALSE,
    remove.na = FALSE,
    autoIndex = TRUE,
    allowMulti = FALSE,
    cpus = 1
){
    # check input args
    if (FALSE %in% file.exists(bamFiles) || !file.exists(anno))
        stop("mapping stats: BAM or annotation file(s) do(es) not exist!")
    if ((length(fqFiles) == 1 && is.na(fqFiles)) || length(fqFiles) != length(bamFiles) * (1 + paired))
        stop("mapping stats: fqFiles is NA or has invalid length!")
    if (!is.na(samples) && length(samples) != length(bamFiles) * (1 + paired))
        stop("mapping stats: samples has invalid length!")
    # NOTE: halve the values during scaling!
    if (!is.na(numReads) && (!is.numeric(numReads) || length(numReads) != length(bamFiles) * (1 + paired)))
        stop("mapping stats: numReads is NA, or not numberic, or has invalid length!")
    if (!is.na(numTrimmed) && (!is.numeric(numTrimmed) || length(numTrimmed) != length(bamFiles) * (1 + paired)))
        stop("mapping stats: numTrimmed must be NA, numeric or has invalid length!")
    if (!is.na(numNonrRNA) && (!is.numeric(numNonrRNA) || length(numNonrRNA) != length(bamFiles) * (1 + paired)))
        stop("mapping stats: numNonrRNA must be NA, numeric or has invalid length!")
    if (precise == FALSE && ( (length(libSizes) == 1 && is.na(libSizes))  ||  length(libSizes) != length(bamFiles) ) )
        stop("mapping stats: unprecise mode requires valid libSizes, but set to NA or invalid length!")
    if (length(numReads) == 1 && is.na(numReads))
      numReads <- number_reads_fastq(fqFiles, cpus = cpus)

    if (precise == TRUE && !.is_indexed(bamFiles)) {
        if (autoIndex == TRUE) {
            message("One or more files are not indexed. Auto-Index enabled...")
            bamFiles <- sortBAMs(bamFiles)
        } else {
            stop(
                "Cannot calculate exact exon coverage without indexed BAM files",
                "and autoIndex is FALSE."
            )
        }
    }

    # prepare table
    mapping_stats <- c(
        "reads_raw",
        "reads_trimmed",
        "percentage_reads_trimmed",
        "reads_non_rRNA",
        "percentage_non_rRNA_reads",
        "reads_mapped",
        "reads_unmapped",
        "percentage_reads_mapped",
        "percentage_total_read_loss",
        "genome_coverage",
        "reads_mapping_in_exons",
        "percentage_reads_mapping_in_exons",
        "exon_coverage"
    )

    mapping_stats <- data.frame(matrix(
        nrow = length(fqFiles),
        ncol = length(mapping_stats),
        dimnames = list(basename(fqFiles), mapping_stats)
    ))

    # add trivial data
    prior <- numReads # HACK: temp value of previously, used value for percentage calculations
    mapping_stats$reads_raw <- numReads
    if (length(numTrimmed) > 1 || !is.na(numTrimmed)) {
        # TODO: numTrimmed refers to "remaining after trimming". Isn't that confusing?
        mapping_stats$reads_trimmed <- numTrimmed
        mapping_stats$percentage_reads_trimmed <- (1 - numTrimmed / prior) * 100
        prior <- numTrimmed
    }
    if (length(numNonrRNA) > 1 || !is.na(numNonrRNA)) {
        mapping_stats$reads_non_rRNA <- numNonrRNA
        mapping_stats$percentage_non_rRNA_reads <- (numNonrRNA / prior) * 100
        prior <- numNonrRNA
    }

    # calculate data
    if (precise) {
        map_cov <- mapping_coverage(bamFiles, allowMulti = allowMulti, paired = paired, cpus = cpus)
        exon_cov <- exon_coverage(
            bamFiles,
            anno = anno,
            featureType = featureType,
            paired = paired,
            cpus = cpus
        )

        if (paired) {
            exon_cov$nrReads <- rep(exon_cov$nrReads / 4, each = 2)
            exon_cov$cov     <- rep(exon_cov$cov / 4, each = 2)
            map_cov$records  <- rep(map_cov$records / 2, each = 2)
            map_cov$cov      <- rep(map_cov$cov / 2, each = 2)
        }

        # NOTE: unlike estimation, this automatically excludes multi mapped reads!
        mapping_stats$reads_mapped               <- map_cov$records
        mapping_stats$reads_unmapped             <- pmax(0, prior - map_cov$records)
        mapping_stats$percentage_reads_mapped    <- (map_cov$records / prior) * 100
        mapping_stats$percentage_total_read_loss <-
            pmax(0, (1 - map_cov$records / numReads) * 100)

        mapping_stats$genome_coverage                   <- map_cov$cov
        mapping_stats$reads_mapping_in_exons            <- exon_cov$nrReads
        mapping_stats$percentage_reads_mapping_in_exons <-
            (exon_cov$nrReads / prior) * 100
        mapping_stats$exon_coverage                     <- exon_cov$cov
    } else {
        # estimation
        emc <- estimate_mapping_coverage(bamFiles, allowMulti = allowMulti, cpus = cpus, workers = 5)
        # NOTE: featureCounts uses FRAGMENTS for paired mode. Therefore, each pair produces only 1 count!
        if (paired) emc$readLen <- emc$readLen * 2
        est_exon_cov <- estimate_exon_coverage(
            libSizes = libSizes,
            anno = anno,
            featureType = featureType,
            readLength = emc$readLen,
            allowMulti = allowMulti
        )
        if (paired) {
            emc$nrReads      <- rep(emc$nrReads / 2, each = 2)
            emc$cov          <- rep(emc$cov / 2, each = 2)
            est_exon_cov$cov <- rep(est_exon_cov$cov / 2, each = 2)
            libSizes         <- rep(libSizes, each = 2)
        }

        mapping_stats$reads_mapped               <- emc$nrReads
        mapping_stats$reads_unmapped             <- pmax(0, prior - emc$nrReads)
        mapping_stats$percentage_reads_mapped    <- (emc$nrReads / prior) * 100
        mapping_stats$percentage_total_read_loss <-
            pmax(0, (1 - emc$nrReads / numReads) * 100)

        mapping_stats$genome_coverage <- emc$cov
        mapping_stats$reads_mapping_in_exons <- libSizes
        mapping_stats$percentage_reads_mapping_in_exons <- (libSizes / prior) * 100
        mapping_stats$exon_coverage <- est_exon_cov$cov
    }

    if (remove.na) {
        naCols <- is.na(mapping_stats[1,])
        mapping_stats <- mapping_stats[, !naCols]
    }

    # add file column
    mapping_stats_df <- cbind(fastq_files = fqFiles, mapping_stats)

    # use samples as row names, if available
    if (!is.na(samples)) {
        rownames(mapping_stats_df) = samples
    }

    return(mapping_stats_df)
}


# NOTE: handles empty BAM files
# Averge Read Length Of Mapped Reads
#
# Calculate the average read length of mapped reads using SAMtools.
#
# @rdname read_length_mapped
# @param files BAM files to analyze.
# @param print.nr.reads Logical indicating whether to print absolute number
#   of reads found.
# @param allowMulti Logical. If True, secondary alignments are not considered.
#   Defaults to FALSE.
# @inheritParams parallelizationParameters
# @return If print.nr.reads = TRUE:\cr
# Keyword list. Average read length and number of reads per file (print.nr.reads = TRUE).\cr
# Else:\cr
# Vector. Average read length per file (print.nr.reads = FALSE).
.read_length_mapped <- function
(
    files,
    print.nr.reads = TRUE,
    allowMulti = FALSE,
    cpus = 1,
    workers = 5
){
    driver <- function(f) {
        # see function mapping_coverage()?
        exclude_flag <- if (allowMulti) 4 else 260

        call <- paste(
            "samtools view",
            "-F", exclude_flag,
            f,
            "|",
            "awk '{sum+=length($10)} END {print (NR > 0 ? sum/NR : 0); print NR}'"
        )

        tmp <- as.numeric(system(call, intern = TRUE))
        return(tmp)
    }

    allocCPUS <- allocate_cpus(cpus, length(files), workers)
    res <- unlist(
        parallel::mclapply(files, driver, mc.cores = allocCPUS$in_parallel)
    )

    if (print.nr.reads) {
        return(list(
            readLength = res[seq(1, length(res), by=2)],
            nrReads   = res[seq(2, length(res), by=2)]
        ))
    } else {
      return(res[seq(1, length(res), by=2)])
    }
}


#' Get Genome Size From BAM File
#'
#' Read a BAM file and return the length of, e.g., the genome sequence, which
#' was used for mapping and creating this BAM file.
#'
#' The result of mapping reads against a reference sequence, e.g. a genome, is
#' usually stored in BAM files. The header section of such BAM files from common
#' mapping tools contains the length of each genome. These header lines start
#' with 'SQ' and the length is encoded as 'LN:xxxxx', where 'x' is a digit.
#'
#' @export
#' @param bamFile Path to BAM file to get genome size from.
#' @param per_chrom Logical indicating whether to return chromosome lengths
#'   instead of the total genome length. Defaults to FALSE.
#' @return Either Numeric(1) or Numeric(n).
#' \cr\cr
#' Numeric(1) if per_chrom == FALSE
#' \cr
#' Numeric(n) if per_chrom == TRUE
#' @examples
#' bam <- system.file("extdata", "untreated3_chr4.bam", package = "Geo2RNAseq")
#' genome_size_from_bam(bam)
#' genome_size_from_bam(bam, per_chrom = TRUE)
genome_size_from_bam <- function
(
    bamFile,
    per_chrom = FALSE
){
    # Hack: Genome names may also contain ':'. To catch this issue, is made 'cut' always use the LAST field.
    #       This is achieved by 1. reversing input, 2. cutting out the last field based on ':' and 3. reversing the resulting numbers again
    # call <- paste(
    #   "samtools view -H",
    #   bamfile,
    #   "| grep \"^@SQ\" | rev | cut -f 1 -d \":\" | rev"
    # )
    # genomeLengths <- as.numeric(system(call, intern = TRUE))
    # return(sum(genomeLengths))
    # NOTE: Returns total genome length or 0 for empty BAM files and if 'SQ' lines are missing.
    if (per_chrom)
        return(Rsamtools::scanBamHeader(bamFile[[1]])[[1]]$targets)
    else
        return(sum(as.numeric(Rsamtools::scanBamHeader(bamFile[[1]])[[1]]$targets)))
}


#' Get Chromosome Sizes From Genome
#'
#' Calculate the number of nucleotides per chromosome from a multi-FASTA file.
#'
#' @export
#' @param genome Path to FASTA file, containing one or many (chromosome)
#'   sequences.
#' @return Named Numeric(m).
#' @examples chromosome_sizes_from_genome(system.file("extdata", "dm3_chr4.fa", package = "Geo2RNAseq"))
chromosome_sizes_from_genome <- function
(
    genome
){
    if (!file.exists(genome))
        stop(paste("Could not find genome file:", genome))
    else if (length(genome) > 1)
        stop("Too many inputs. Supply only one genome.")
    else {
        # awk script: Look for header sections. At each header, start a counter 'c'. Else, count letters
        # call <- paste(
        #   "cat",
        #   genome,
        #   "|",
        #   "awk '$0 ~ \">\" {print c; c=0;} $0 !~ \">\" {c += length($0)} END {print c}'"
        # )
        # sizes <- as.numeric(system(call, intern = TRUE)[-1]) # the first is always NA
        g <- Biostrings::readDNAStringSet(genome)
        sizes <- g@ranges@width
        chr <- sapply(strsplit(names(g), split = " \\(.*\\)"), function(x) x[[1]])
        names(sizes) <- chr
        return(sizes)
    }
}



#' Calculate Genome Size From Multi-FASTA File
#'
#' Calculate the number of nucleotides for an entire genome from a multi-FASTA
#' file.
#'
#' @export
#' @param genome Path to FASTA file containing a genome (sequence), which is
#'   often composed of multiple chromosome sequences.
#' @return Numeric(m).
#' @examples genome_size_from_fasta(system.file("extdata", "dm3_chr4.fa", package = "Geo2RNAseq"))
genome_size_from_fasta <- function
(
    genome
){
    if (!file.exists(genome))
        stop(paste("Could not find genome file:", genome))
    else if (length(genome) > 1)
        stop("Too many inputs. Supply only one genome.")
    else {
        #as.numeric(system(paste("grep -v '^>'", genome, "| tr -d '\n' | wc -m"), intern = TRUE))
        g <- Biostrings::readDNAStringSet(genome)
        return(sum(g@ranges@width))
    }
}


#' Calculate Mapping Coverage
#'
#' Calculate mapping coverage using Rsamtools and GenomicAlignments. The
#' coverage is calculated for each BAM file separately.
#'
#' @export
#' @param bamFiles Mapping files in BAM format.
#' @param paired Logical indicating whether the mapping process, which result is
#'   stored in the BAM files, was based on paired-end or single-end reads.
#'   Defaults to FALSE, which means based on single-end reads.
#' @param allowMulti Logical. If True, secondary alignments are not considered.
#'   Defaults to FALSE.
#' @inheritParams parallelizationParameters
#' @return A list with keyword arguments. These are:\cr
#' \describe{
#'   \item{\strong{cov}}{Average mapping coverage across the genome, per BAM file.}
#'   \item{\strong{records}}{Number of mapped reads found, per BAM file.}
#' }
#' @examples
#' mapping_coverage(system.file("extdata", "untreated3_chr4.bam", package = "Geo2RNAseq"))
#' mapping_coverage(system.file("extdata", "untreated3_chr4.bam", package = "Geo2RNAseq"), paired = TRUE)
mapping_coverage <- function
(
    bamFiles,
    paired = FALSE,
    allowMulti = FALSE,
    cpus = 1,
    workers = 5
){
    driver <- function(bam) {
        # do not count unmapped reads
        countSeconds <- if (allowMulti) NA else FALSE
        flag <- Rsamtools::scanBamFlag(isPaired = paired,
                                       isSecondaryAlignment = countSeconds,
                                       isUnmappedQuery = FALSE)
        what <- Rsamtools::ScanBamParam(flag = flag)
        countObj <- Rsamtools::countBam(bam, param = what)

        #TODO: get average & peaks using GenomicAlignments, but values are off!
        # get average per chromosome, then average of average
        #covObj <- GenomicAlignments::coverage(bam, param = what)
        #max <- max(unlist(lapply(covObj, function(chr) max(chr))))
        #avg <- unlist(lapply(covObj, function(chr) mean(chr)))
        genome_size <- genome_size_from_bam(bam)
        cov <- countObj$nucleotides / genome_size
        return(list(cov = cov, nr = countObj$records))
    }

    allocCPUS <- allocate_cpus(cpus, length(bamFiles), workers)
    res <- parallel::mclapply(bamFiles, driver, mc.cores = allocCPUS$in_parallel)

    return(list(
        cov     = sapply(res, function(x) x[["cov"]]),
        records = sapply(res, function(x) x[["nr"]])
    ))
}


#' Estimate Mapping Coverage
#'
#' Estimate mapping coverage based on genome length, average read length and
#' number of mapped reads calculated from BAM files.
#'
#' @export
#' @param files Mapping files in BAM format.
#' @param genomeSize Size of genome used. If 0, the genome size is determined
#'   based on SQ lines in the header region of BAM files. Defaults to 0.
#' @param allowMulti Logical. If True, secondary alignments are not considered.
#'   Defaults to FALSE.
#' @inheritParams parallelizationParameters
#' @return A named list with the following keyword arguments:\cr
#' \describe{
#'   \item{cov    }{Estimated genome coverage per sample.}
#'   \item{nrReads}{Number of reads found in BAM.}
#'   \item{readLen}{Average length of mapped reads.}
#' }
#' @examples estimate_mapping_coverage(system.file("extdata", "untreated3_chr4.bam", package = "Geo2RNAseq"))
estimate_mapping_coverage <- function
(
    files,
    genomeSize = 0,
    allowMulti = FALSE,
    cpus = 1,
    workers = 5
){
    if (genomeSize <= 0)  genomesize <- genome_size_from_bam(files[1])
    len_n_num <- .read_length_mapped(
                    files = files,
                    print.nr.reads = TRUE,
                    allowMulti = allowMulti,
                    cpus = cpus,
                    workers = workers
                 )
    coverage <- len_n_num$nrReads * len_n_num$readLength / genomesize

    return(list(
        cov = coverage,
        nrReads = len_n_num$nrReads,
        readLen = len_n_num$readLength
    ))
}


#' Exon Size
#'
#' Calculate exon size from a feature (annotation) file.
#'
#' @export
#' @param anno Feature file in GFF or GTF format.
#' @param featureType Consider only this feature. Usually "exon" or "gene".
#'   If NULL or NA, every feature is considered. Defaults to "exon".
#' @return Integer(1). Total size of exome.
#' @examples
#' \dontrun{
#' library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
#' exbygene <- GenomicFeatures::exonsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene, "gene")
#' # convert GRangesList to GRanges object
#' exbygene <- unlist(exbygene)
#' rtracklayer::export(exbygene, "dm3.chr4.gtf", format = "gtf")
#' exon_size("dm3.chr4.gtf", featureType = "sequence_feature")
#' }
#' Geo2RNAseq::exon_size
exon_size <- function
(
    anno,
    featureType = "exon"
){
    ranges <- rtracklayer::import(anno)
    ranges <- ranges[ranges$type == featureType]

    ## WRONG
    ## counts regions from overlapping exons or transcript variants multiple times
    # exon_size <- sum(ranges@ranges@width)

    ## WRONG
    ## without multiple counting, but ignores sequence names, e.g. scaffolds,
    ## hence reduces ranges over the entire genome, which gives weird results
    # exon_size <- sum(IRanges::reduce(ranges@ranges)@width) # total length without double counting anything

    ## RIGHT
    ## reduce ranges per scaffold, sum up all reduced ranges of a scaffold,
    ## sum up the sums of all scaffolds
    exon_size <- sum(sapply(
        unique(as.character(ranges@seqnames)),
        function(x) {
            sum(IRanges::reduce(
                ranges@ranges[as.character(ranges@seqnames) == x])@width,
                ignore.strand = T
            )
        }
    ))

    if (exon_size <= 0)
        warning(
            "Exon size of zero. You probably used the wrong 'featureType'. ",
            "Check the 3rd column in the supplied annotation!"
        )
    return(exon_size)
}


#' Estimate Exon Coverage
#'
#' Estimate exon coverage based on the average read length, the number of
#' counted reads per sample (libSizes) and an annotation file. If the average
#' read length is given, the coverage is estimated instantly.
#'
#' @export
#' @param libSizes Integer(n). Counted reads per sample with feature type "exon".
#'   Usually the return value of 'colSums(counts)', with genes in rows and
#'   samples in columns.
#' @param anno Feature file in GFF or GTF format.
#' @param bamFiles Character(n). Path to one or more BAM files. Only needed if
#'   readLength is NA. Else, readLength is calculated from the first given BAM
#'   file. Defaults to NA.
#' @param readLength Numeric(n). Average length of mapped reads. Ignores
#'   'bamFiles' if set. Defaults to NA.
#' @param featureType  Consider only this feature. Usually "exon" or "gene".
#'   If NULL or NA, every feature is considered. Defaults to "exon".
#' @param allowMulti Logical. If True, secondary alignments are not considered.
#'   Defaults to FALSE.
#' @return A named list with the following keyword arguments:\cr
#' \describe{
#'   \item{cov     }{Estimated exon coverage per sample.}
#'   \item{readLen }{Average length of mapped reads.}
#'   \item{exonSize}{Sum of lengths of all exons.}
#' }
#' @examples
#' \dontrun{
#' libSizes <- c(1000, 2000)
#' anno <- "anno.gtf"
#' bamFiles <- c("f1.bam", "f2.bam")
#' read_length <- c(50, 50)
#'
#' estimate_exon_coverage(libSizes, anno, bamFiles)
#' estimate_exon_coverage(libSizes, anno, readLength = read_length)
#' }
#' Geo2RNAseq::estimate_exon_coverage
estimate_exon_coverage <- function
(
    libSizes,
    anno,
    bamFiles = NA,
    readLength = NA,
    featureType = "exon",
    allowMulti = FALSE
){
    exon_size <- exon_size(anno, featureType)
    if (!is.list(readLength) && !is.na(readLength) && length(readLength) > 0 && !is.numeric(readLength))
        readLength <- .read_length_mapped(bamFiles, print.nr.reads = TRUE, allowMulti = allowMulti)$readLength
    coverage <- libSizes * readLength / exon_size

    return(list(
        cov = coverage,
        readLen = readLength,
        exonSize = exon_size
    ))
}



#' Calculate Exon Coverage
#'
#' Calculate the exon coverage per BAM file using Rsamtools.
#'
#' @export
#' @param bamFiles BAM file of mapped reads. Must be sorted and indexed.
#' @param anno Feature file in GFF or GTF format.
#' @param featureType Consider only this feature. Usually "exon" or "gene". If
#'   NULL or NA, every feature is considered.
#' @param autoIndex If TRUE and index files are missing, BAM files are sorted
#'   and indexed first.
#' @param paired Logical indicating the use of paired-end mapped reads.
#' @inheritParams parallelizationParameters
#' @return A list with keyword arguments:\cr
#' \describe{
#'   \item{cov     }{Estimated genome coverage per sample.}
#'   \item{nrReads }{Number of reads found in the BAM files.}
#'   \item{exonsize}{Total number of nucleotides of exons.}
#' }
#' @examples
#' \dontrun{
#' anno <- "anno.gtf"
#' bamFiles <- c("f1.bam", "f2.bam")
#' exon_coverage(bamFiles, anno)
#' }
#' Geo2RNAseq::exon_coverage
exon_coverage <- function
(
    bamFiles,
    anno,
    featureType = "exon",
    autoIndex = TRUE,
    paired = FALSE,
    cpus = 1,
    workers = 5
){
    driver <- function(bam) {
        flag <- Rsamtools::scanBamFlag(
            isPaired = paired,
            isSecondaryAlignment = FALSE,
            isUnmappedQuery = FALSE
        )
        what <- Rsamtools::ScanBamParam(flag = flag, which = ranges)

        countObj <- Rsamtools::countBam(bam, param = what)
        cov <- sum(countObj$nucleotides)/exon_size

        # TODO: 'coverage' uses total genome size. We only need 'ranges' ~> how to?
        # get average per chromosome, then average of average
        #covObj <- GenomicAlignments::coverage(bam, param = what)
        #max <- max(unlist(lapply(covObj, function(chr) max(chr))))
        #cov <- mean(unlist(lapply(covObj, function(chr) mean(chr))))

        return(list(cov = cov, nr = sum(countObj$records)))
    }

    if (!.is_indexed(bamFiles)) {
        if (autoIndex == TRUE) {
            message("One or more files are not indexed. Auto-Index enabled...")
            bamFiles <- sortBAMs(bamFiles)
        } else {
            stop(
                "Cannot calculate exact exon coverage without indexed BAM files",
                "and auto_index is FALSE."
            )
        }
    }

    if (!file.exists(anno))
        stop("Annotation file does not exist: ", shQuote(anno))

    ranges <- rtracklayer::import(anno)
    ranges <- ranges[ranges$type == featureType]
    exon_size <- exon_size(anno, featureType)

    allocCPUS <- allocate_cpus(cpus, length(bamFiles), workers)
    res <- parallel::mclapply(bamFiles, driver, mc.cores = allocCPUS$in_parallel)

    return(list(
        cov      = sapply(res, function(x) x[["cov"]]),
        nrReads  = sapply(res, function(x) x[["nr"]]),
        #max     = sapply(res, function(x) x[["max"]]),
        exonsize = exon_size
    ))
}
