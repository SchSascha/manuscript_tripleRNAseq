#####-
##-====================================================
##-==========================  DUAL  ===========================
##-====================================================
#####-



#' Stats For Combined Dual RNA-seq Mapping
#'
#' After mapping against a combined dual RNA-seq genome, calculate the number of
#' reads mapping to the two genomes, e.g. host and pathogen, respectively.
#'
#' @export
#' @param files Character(n). BAM files resulting from mapping against a
#'   combined two-species genome.
#' @param chromosomes Character(n). All chromosome names of one species of the
#'   two, usually from the genome FASTA file headers of this species.
#' @param outDir Output directory for temporary files. Will be created if it
#'   does not exist.
#' @inheritParams parallelizationParameters
#'
#' @return A list with the number of reads from species 1 and species 2, e.g.
#'   host and pathogen, per BAM file (sample).
#' @examples
#' \dontrun{
#'   number_combinedMapping <- stats_combinedMapping(
#'       files       = c("s1.sort.bam", "s2.sort.bam"),
#'       chromosomes = c("HS-1", "HS-2", "HS-X"),
#'       cpus        = 10
#'    )
#'    names(number_combinedMapping) <- basename(bam_files)
#'    number_combinedMapping <- unlist(number_combinedMapping)
#'    head(number_combinedMapping)
#' }
stats_dualMapping <- function
(
    files,
    chrom_sp1,
    outDir = "tmp",
    cpus = 1
){
    driver <- function(i, files) {
        fnames <- basename(files[i])
        inames <- paste0(files[i], ".bai")

        # BAM files needs to be indexed
        if (!file.exists(inames)) {
            warning("BAM files need to be sorted and indexed. Trying to index. If this fail, use 'sortBAMS' on input file list!")
            system(paste("samtools", "index", files[i]))
        }

        # Use sammtools idxstats to get number of mapped reads
        system(paste("samtools",
                    "idxstats",
                    files[i],
                    ">",
                    file.path(outDir, paste0(fnames, ".csv"))))
        tmp_stats <- read.csv(file.path(outDir, paste0(fnames, ".csv")),
                            header = FALSE,
                            sep = "\t")

        # sum up entries for each set of chromosomes
        reads_sp1 <- sum(tmp_stats[tmp_stats[, 1] %in% chrom_sp1, 3])
        reads_sp2 <- sum(tmp_stats[, 3]) - reads_sp1
        return(list(reads_sp1 = reads_sp1, reads_sp2 = reads_sp2))
    }

    res <- parallel::mclapply(1:length(files), driver, files, mc.cores = cpus)
    return(res)
    ### List
}



#' Stats For Combined Dual RNA-seq Mapping
#'
#' After mapping against a combined triple RNA-seq genome, calculate the number of
#' reads mapping to the three genomes, e.g. host, pathogen 1, pathogen 2, respectively.
#'
#' @export
#' @param files Character(n). BAM files resulting from mapping against a
#'   combined three-species genome.
#' @param chrom_sp1 Character(n). All chromosome names of one species of the
#'   three, usually from the genome FASTA file headers of this species.
#' @param chrom_sp2 Character(n). All chromosome names of another species of
#'   the three, usually from the genome FASTA file headers of this species.
#' @inheritParams parallelizationParameters
#' @inheritParams stats_dualMapping
#'
#' @return A list with the number of reads from species 1 and species 2, e.g.
#'   host and pathogen, per BAM file (sample).
#' @examples
#' \dontrun{
#'   number_combinedMapping <- stats_combinedMapping(
#'       files       = c("s1.sort.bam", "s2.sort.bam"),
#'       chrom_sp1   = c("HS-1", "HS-2", "HS-X"),
#'       chrom_sp2   = c("AFU-1", "AFU-2"),
#'       cpus        = 10
#'    )
#'    names(number_combinedMapping) <- basename(bam_files)
#'    number_combinedMapping <- unlist(number_combinedMapping)
#'    head(number_combinedMapping)
#' }
stats_tripleMapping <- function
(
  files,
  chrom_sp1,
  chrom_sp2,
  outDir = "tmp",
  cpus = 1
){
  driver <- function(i, files) {
    fnames <- basename(files[i])
    inames <- paste0(files[i], ".bai")

    # BAM files needs to be indexed
    if (!file.exists(inames)) {
      warning("BAM files need to be sorted and indexed. Trying to index. If this fail, use 'sortBAMS' on input file list!")
      system(paste("samtools", "index", files[i]))
    }

    # Use sammtools idxstats to get number of mapped reads
    system(paste("samtools",
                 "idxstats",
                 files[i],
                 ">",
                 file.path(outDir, paste0(fnames, ".csv"))))
    tmp_stats <- read.csv(file.path(outDir, paste0(fnames, ".csv")),
                          header = FALSE,
                          sep = "\t")

    # sum up entries for each set of chromosomes
    reads_sp1 <- sum(tmp_stats[tmp_stats[, 1] %in% chrom_sp1, 3])
    reads_sp2 <- sum(tmp_stats[tmp_stats[, 1] %in% chrom_sp2, 3])
    reads_sp3 <- sum(tmp_stats[, 3]) - sum(reads_sp1 + reads_sp2)
    return(list(reads_sp1 = reads_sp1, reads_sp2 = reads_sp2, reads_sp3 = reads_sp3))
  }

  dir.create(outDir, recursive = T, showWarnings = F)
  res <- parallel::mclapply(1:length(files), driver, files, mc.cores = cpus)
  return(res)
}





#' Calculate Mapping Stats for Dual RNA-seq Mapping
#'
#' After mapping against a combined dual RNA-seq genome, calculate the number of
#' reads mapping to the two genomes, e.g. host and pathogen, respectively.
#'
#' @export
#' @param number_reads Numeric(n). Number of quality controlled reads per sample.
#' @param raw_fastq_files Character(n). Paths to original FASTQ reads files.
#' @param name_sp1 Character(1). Name of species 1.
#' @param name_sp2 Character(1). Name of species 2.
#' @param genome_sp1 Character(1). Full path to genome FASTA reference for
#'   first species.
#' @param genome_sp2 Same for second species.
#' @param anno_sp1 Character(1). Full path to annotation for first species.
#' @param anno_sp2 Same for second species.
#' @param counts_sp1 Numeric(n x m). Number of reads per sample per gene for
#'   species one. This usually is acquired by running featureCounts first,
#'   and subsetting the genes into each species second.
#' @param counts_sp2 Same for second species.
#' @param featureType Annotation feature type to use for exon entries. This
#'   reference to the 3rd column of a GTF foramtted annotation file. Defaults
#'   to 'exon'.
#' @return A list with the number of reads from species 1 and species 2, e.g.
#'   host and pathogen, per BAM file (sample).
#' @examples
#' \dontrun{
#'   number_combinedMapping <- stats_combinedMapping(
#'       files       = c("s1.sort.bam", "s2.sort.bam"),
#'       chrom_sp1   = c("HS-1", "HS-2", "HS-X"),
#'       cpus        = 10
#'    )
#'    names(number_combinedMapping) <- basename(bam_files)
#'    number_combinedMapping <- unlist(number_combinedMapping)
#'    head(number_combinedMapping)
#' }
calc_dual_mapping_stats <- function(number_reads,
                                    raw_fastq_files,
                                    name_sp1, name_sp2,
                                    genome_sp1, genome_sp2,
                                    anno_sp1,   anno_sp2,
                                    counts_sp1, counts_sp2,
                                    featureType = "exon") {

    if (T %in% !file.exists(c(genome_sp1, genome_sp2)))
        stop("One or more genome files were not found on defined path!")
    if (T %in% !file.exists(c(anno_sp1, anno_sp2)))
        stop("One or more annotation files were not found on defined path!")

    # get chromosome names
    chromosomeNames_sp1 <- system(
        paste("perl", "-ne", "'if (/^>(\\S+)/) { print \"$1\n\" }'", genome_sp1),
        intern = TRUE
    )
    chromosomeNames_sp2 <- system(
        paste("perl", "-ne", "'if (/^>(\\S+)/) { print \"$1\n\" }'", genome_sp2),
        intern = TRUE
    )

    if (length(chromosomeNames_sp1) == 0 || length(chromosomeNames_sp2) == 0)
        stop("One of the supplied FASTA genomes is wrong. No FASTA header were found.")

    # get chromosome and exon sizes
    exonsize_sp1   <- exon_size(anno_sp1, featureType)
    exonsize_sp2   <- exon_size(anno_sp2, featureType)
    if (exonsize_sp1 == 0 || exonsize_sp2 == 0)
        stop("One of the given exon files has zero length! Did you use the correct 'featureType'?")

    genomesize_sp1 <- sum(as.numeric(chromosome_sizes_from_genome(genome_sp1)))
    genomesize_sp2 <- sum(as.numeric(chromosome_sizes_from_genome(genome_sp2)))
    if (genomesize_sp1 == 0 || genomesize_sp2 == 0)
        stop("One of the given genome files has zero length!")

    number_exon_sp1 <- colSums(counts_sp1)
    number_exon_sp2 <- colSums(counts_sp2)


    # remaining chromosomes are sp2
    number_combinedMapping <- stats_dualMapping(
        files       = bam_files,
        chrom_sp1   = chromosomeNames_sp1,
        cpus        = 10
    )
    names(number_combinedMapping) <- basename(bam_files)
    number_combinedMapping <- unlist(number_combinedMapping)

    # atm, it only handles single RNA-seq
    mapping_stats <- c(
        "reads mapped sp1",
        "reads mapped sp2",
        "percentage reads mapped sp1",
        "percentage reads mapped sp2",
        "genome coverage sp1",
        "genome coverage sp2",
        "reads mapping in exons sp1",
        "reads mapping in exons sp2",
        "percentage reads mapping in exons sp1",
        "percentage reads mapping in exons sp2",
        "exon coverage sp1",
        "exon coverage sp2"
    )

    # substitute all whitespaces in colnames with underscores
    mapping_stats <- sapply(mapping_stats, function(x) {
        gsub(" ", "_", x)
    }, USE.NAMES = FALSE)

    mapping_stats <- data.frame(matrix(
        nrow = length(samples),
        ncol = length(mapping_stats),
        dimnames = list(samples, mapping_stats)
    ))

    # mapping to genome in general
    mapping_stats$reads_mapped_sp1 <- number_combinedMapping[grep("reads_sp1$", names(number_combinedMapping))]
    mapping_stats$reads_mapped_sp2 <- number_combinedMapping[grep("reads_sp2$", names(number_combinedMapping))]

    mapping_stats$percentage_reads_mapped_sp1 <- mapping_stats$reads_mapped_sp1 / number_reads * 100
    mapping_stats$percentage_reads_mapped_sp2 <- mapping_stats$reads_mapped_sp2 / number_reads * 100

    readlength <- .read_length_mapped(bam_files, cpus = MAX_CPUS)$readLength

    mapping_stats$genome_coverage_sp1 <- mapping_stats$reads_mapped_sp1 * readlength / genomesize_sp1
    mapping_stats$genome_coverage_sp2 <- mapping_stats$reads_mapped_sp2 * readlength / genomesize_sp2

    # mapping to exons only
    mapping_stats$reads_mapping_in_exons_sp1 <- number_exon_sp1
    mapping_stats$reads_mapping_in_exons_sp2 <- number_exon_sp2
    mapping_stats$exon_coverage_sp1 <- number_exon_sp1 * readlength / exonsize_sp1
    mapping_stats$exon_coverage_sp2 <- number_exon_sp2 * readlength / exonsize_sp2

    mapping_stats$percentage_reads_mapping_in_exons_sp1 <- number_exon_sp1 / number_reads * 100
    mapping_stats$percentage_reads_mapping_in_exons_sp2 <- number_exon_sp2 / number_reads * 100

    # replace sp1/2 in colnames with actual names of the species
    colnames(mapping_stats) <- gsub("sp1", name_sp1, colnames(mapping_stats))
    colnames(mapping_stats) <- gsub("sp2", name_sp2, colnames(mapping_stats))

    # add column with FASTQ file names to mapping stats
    mapping_stats_df <- cbind(
        "fastq_files" = basename(raw_fastq_files),
        mapping_stats
    )

    return(mapping_stats_df)
}




#' Calculate Mapping Stats for Triple RNA-seq Mapping
#'
#' After mapping against a combined triple RNA-seq genome, calculate the number of
#' reads mapping to the three genomes, e.g. host, pathogen 1, pathogen 2, respectively.
#'
#' @export
#' @param name_sp3 Name of species 3.
#' @param genome_sp3 Same for third species.
#' @param anno_sp3 Same for third species.
#' @param counts_sp3 Same for third species.
#' @inheritParams calc_dual_mapping_stats
#' @return A list with the number of reads from species 1, 2 and 3, e.g.
#'   host, pathogen 1 and pathogen 2, per BAM file (sample).
#' @examples
#' \dontrun{
#'   number_combinedMapping <- stats_combinedMapping(
#'       files       = c("s1.sort.bam", "s2.sort.bam"),
#'       chrom_sp1   = c("HS-1", "HS-2", "HS-X"),
#'       chrom_sp2   = c("AFU-1", "AFU-2"),
#'       cpus        = 10
#'    )
#'    names(number_combinedMapping) <- basename(bam_files)
#'    number_combinedMapping <- unlist(number_combinedMapping)
#'    head(number_combinedMapping)
#' }
calc_triple_mapping_stats <- function(number_reads,
                                      raw_fastq_files,
                                      name_sp1, name_sp2, name_sp3,
                                      genome_sp1, genome_sp2, genome_sp3,
                                      anno_sp1,   anno_sp2,   anno_sp3,
                                      counts_sp1, counts_sp2, counts_sp3,
                                      featureType = "exon") {

    if (T %in% !file.exists(c(genome_sp1, genome_sp2, genome_sp3)))
        stop("One or more genome files were not found on defined path!")
    if (T %in% !file.exists(c(anno_sp1, anno_sp2, anno_sp3)))
        stop("One or more annotation files were not found on defined path!")

    # get chromosome names
    chromosomeNames_sp1 <- system(
        paste("perl", "-ne", "'if (/^>(\\S+)/) { print \"$1\n\" }'", genome_sp1),
        intern = TRUE
    )
    chromosomeNames_sp2 <- system(
        paste("perl", "-ne", "'if (/^>(\\S+)/) { print \"$1\n\" }'", genome_sp2),
        intern = TRUE
    )
    chromosomeNames_sp3 <- system(
        paste("perl", "-ne", "'if (/^>(\\S+)/) { print \"$1\n\" }'", genome_sp3),
        intern = TRUE
    )

    if (length(chromosomeNames_sp1) == 0 || length(chromosomeNames_sp2) == 0 || length(chromosomeNames_sp3) == 0)
        stop("One of the supplied FASTA genomes is wrong. No FASTA header were found.")

    # get chromosome and exon sizes
    exonsize_sp1   <- exon_size(anno_sp1, featureType)
    exonsize_sp2   <- exon_size(anno_sp2, featureType)
    exonsize_sp3   <- exon_size(anno_sp3, featureType)
    if (exonsize_sp1 == 0 || exonsize_sp2 == 0 || exonsize_sp3 == 0)
        stop("One of the given exon files has zero length! Did you use the correct 'featureType'?")

    genomesize_sp1 <- sum(as.numeric(chromosome_sizes_from_genome(genome_sp1)))
    genomesize_sp2 <- sum(as.numeric(chromosome_sizes_from_genome(genome_sp2)))
    genomesize_sp3 <- sum(as.numeric(chromosome_sizes_from_genome(genome_sp3)))
    if (genomesize_sp1 == 0 || genomesize_sp2 == 0 || genomesize_sp3 == 0)
        stop("One of the given genome files has zero length!")

    number_exon_sp1 <- colSums(counts_sp1)
    number_exon_sp2 <- colSums(counts_sp2)
    number_exon_sp3 <- colSums(counts_sp3)


    # remaining chromosomes are sp2
    number_combinedMapping <- stats_tripleMapping(
        files       = bam_files,
        chrom_sp1   = chromosomeNames_sp1,
        chrom_sp2   = chromosomeNames_sp2,
        cpus        = 10
    )
    names(number_combinedMapping) <- basename(bam_files)
    number_combinedMapping <- unlist(number_combinedMapping)


    # precise mapping stats requires BAM sorting.
    # If TRUE, Bioconductor functions are used to determine mapping stats.
    precise <- FALSE
    if (!exists("number_trimmed")) {
        warning("Variable 'number_trimmed' undefined. Setting it to NA.")
        number_trimmed <- NA
    }
    if (!exists("number_nonrRNA")) {
        warning("Variable 'number_nonrRNA' undefined. Setting it to NA.")
        number_nonrRNA <- NA
    }


    # atm, it only handles single RNA-seq
    mapping_stats <- c(
        "reads mapped sp1",
        "reads mapped sp2",
        "reads mapped sp3",
        "percentage reads mapped sp1",
        "percentage reads mapped sp2",
        "percentage reads mapped sp3",
        "genome coverage sp1",
        "genome coverage sp2",
        "genome coverage sp3",
        "reads mapping in exons sp1",
        "reads mapping in exons sp2",
        "reads mapping in exons sp3",
        "percentage reads mapping in exons sp1",
        "percentage reads mapping in exons sp2",
        "percentage reads mapping in exons sp3",
        "exon coverage sp1",
        "exon coverage sp2",
        "exon coverage sp3"
    )

    # substitute all whitespaces in colnames with underscores
    mapping_stats <- sapply(mapping_stats, function(x) {
        gsub(" ", "_", x)
    }, USE.NAMES = FALSE)

    mapping_stats <- data.frame(matrix(
        nrow = length(samples),
        ncol = length(mapping_stats),
        dimnames = list(samples, mapping_stats)
    ))

    # mapping to genome in general
    mapping_stats$reads_mapped_sp1 <- number_combinedMapping[grep("reads_sp1$", names(number_combinedMapping))]
    mapping_stats$reads_mapped_sp2 <- number_combinedMapping[grep("reads_sp2$", names(number_combinedMapping))]
    mapping_stats$reads_mapped_sp3 <- number_combinedMapping[grep("reads_sp3$", names(number_combinedMapping))]

    mapping_stats$percentage_reads_mapped_sp1 <- mapping_stats$reads_mapped_sp1 / number_reads * 100
    mapping_stats$percentage_reads_mapped_sp2 <- mapping_stats$reads_mapped_sp2 / number_reads * 100
    mapping_stats$percentage_reads_mapped_sp3 <- mapping_stats$reads_mapped_sp3 / number_reads * 100

    readlength <- .read_length_mapped(bam_files, cpus = MAX_CPUS)$readLength

    mapping_stats$genome_coverage_sp1 <- mapping_stats$reads_mapped_sp1 * readlength / genomesize_sp1
    mapping_stats$genome_coverage_sp2 <- mapping_stats$reads_mapped_sp2 * readlength / genomesize_sp2
    mapping_stats$genome_coverage_sp3 <- mapping_stats$reads_mapped_sp3 * readlength / genomesize_sp3

    # mapping to exons only
    mapping_stats$reads_mapping_in_exons_sp1 <- number_exon_sp1
    mapping_stats$reads_mapping_in_exons_sp2 <- number_exon_sp2
    mapping_stats$reads_mapping_in_exons_sp3 <- number_exon_sp3
    mapping_stats$exon_coverage_sp1 <- number_exon_sp1 * readlength / exonsize_sp1
    mapping_stats$exon_coverage_sp2 <- number_exon_sp2 * readlength / exonsize_sp2
    mapping_stats$exon_coverage_sp3 <- number_exon_sp3 * readlength / exonsize_sp3

    mapping_stats$percentage_reads_mapping_in_exons_sp1 <- number_exon_sp1 / number_reads * 100
    mapping_stats$percentage_reads_mapping_in_exons_sp2 <- number_exon_sp2 / number_reads * 100
    mapping_stats$percentage_reads_mapping_in_exons_sp3 <- number_exon_sp3 / number_reads * 100


    # replace sp1/2 in colnames with actual names of the species
    colnames(mapping_stats) <- gsub("sp1", name_sp1, colnames(mapping_stats))
    colnames(mapping_stats) <- gsub("sp2", name_sp2, colnames(mapping_stats))
    colnames(mapping_stats) <- gsub("sp3", name_sp3, colnames(mapping_stats))

    # add column with FASTQ file names to mapping stats
    mapping_stats_df <- cbind(
        "fastq_files" = basename(raw_fastq_files),
        mapping_stats
    )

    return(mapping_stats_df)
}












