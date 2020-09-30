#####-
##-====================================================
##-==========================  COUNTING ===========================
##-====================================================
#####-

#' Read FeatureCounts Files
#'
#' Read 'featureCounts' output files and create 3 new files: a counting table,
#' an annotation-like table and a counting summary. The annotation-like table
#' contains read lengths as additional field.
#'
#' @export
#' @seealso \code{\link{run_featureCounts}}
#' @param featureCountsOutDir Path to directory directory of featureCounts. Only
#'   works for output from our package function \code{\link{run_featureCounts}}.
#' @param countFile Path to count file. NULL by default. Overwrites
#'   featureCountsOutDir if set.
#' @param annotFile Path to annotation file. NULL by default. Overwrites
#'   featureCountsOutDir if set.
#' @param summaryFile Path to summary file. NULL by default. Overwrites
#'   featureCountsOutDir if set.
#' @return A list with keyword arguments: \cr
#' \describe{
#'   \item{countFile}{Path to count file.}
#'   \item{sumFile  }{Path to summary file.}
#'   \item{counts   }{Matrix with genes in rows and samples in columns.}
#'   \item{summary  }{Matrix summarizing the counting result. Samples
#'                    in rows, various counting values in columns.}
#'   \item{genes    }{Vector of counted genes.}
#'   \item{anno     }{Contains additional information (like 'genes',
#'                             'Length') found in the annotation file supplied
#'                             to featureCounts.}
#'   \item{tool     }{The name of the tool.}
#' }
#'
#' @examples
#' \dontrun{
#' read.featureCounts.files("path/to/counting/output")
#' }
#' read.featureCounts.files
read.featureCounts.files <- function(
    featureCountsOutDir,
    countFile = NULL,
    annotFile = NULL,
    summaryFile = NULL
) {
    # get file paths
    if (!is.null(featureCountsOutDir) && file.exists(featureCountsOutDir)) {
        if (is.null(countFile))
            countFile   <- file.path(featureCountsOutDir, "fc.counts.tsv")
        if (is.null(summaryFile))
            summaryFile <- file.path(featureCountsOutDir, "fc.summary.tsv")
        if (is.null(annotFile))
            annotFile   <- file.path(featureCountsOutDir, "fc.annotation.tsv")
    } # else: use user given files

    # read and process data
    counts <- as.matrix(read.table(
        countFile,
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE,
        sep = "\t"
    ))
    summary <- as.matrix(read.table(
        summaryFile,
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE,
        sep = "\t"
    ))
    anno <- read.table(
        annotFile,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
    )
    genes <- as.character(anno$GeneID) # same for all count entries

    return(
        list(
            countFile = countFile,
            sumFile   = summaryFile,
            counts    = counts,
            summary   = summary,
            genes     = genes,
            anno      = anno,
            tool      = "featureCounts"
        )
    )
}



# TODO OPT: parallel mode for paired. Execution corrupts with current code!
#' Count Reads Per Feature Using FeatureCounts
#'
#' Count the number of reads mapping to each genomic feature, given one or more
#' mapping files and a set of genomic features.
#'
#' @export
#' @references Yang Liao, Gordon K Smyth and Wei Shi. featureCounts: an
#'   efficient general-purpose program for assigning sequence reads to genomic
#'   features. Bioinformatics, 30(7):923-30, 2014.
#' @seealso \code{\link{read.featureCounts.files}}
#' @param files SAM or BAM files.
#' @param annotation GFF or GTF file of genomic features. \strong{GTF format is
#'   highly recommended}.
#' @param outDir Path to output directory.
#' @param isGTF Set to FALSE, if annotation is not a GTF. TRUE by default, which
#'   means assuming GTF format by default.
#' @param IDtype If GTF annotation: Specify the attribute type used to group
#'   features (e.g. exons) into meta-features (e.g. genes). "gene_id" by
#'   default. This attribute type is usually the gene identifier.
#' @param featureType If GTF annotation: Only rows of this feature type will be
#'   included for read counting. "exon" by default.
#' @param strandSpecific Indicate if strand-specific read counting should be
#'   performed. Acceptable values: 0 (unstranded), 1 (stranded) and 2
#'   (reversely stranded). 0 by default. For paired-end reads, the strand of the
#'   first read defines the strand of the whole fragment.
#' @param isPairedEnd Logical(1) indicating whether the reads are paired-end.
#'   FeatureCounts is \strong{significantly slower} when using this option.
#'   FALSE by default.
#' @param useMetaFeatures If TRUE, summarization will be performed on the
#'   meta-feature level (e.g. genes) instead of on the feature level (e.g.
#'   exons). TRUE by default.
#' @param allowMultiOverlap If TRUE, reads are allowed to be assigned to more
#'   than one feature. FALSE by default.
#' @param countMultiMapped If TRUE, multi-mapped reads will be counted. FALSE by
#'   default.
#' @inheritParams parallelizationParameters
#' @return A list with keyword arguments: \cr
#' \describe{
#'   \item{countFile}{Path to count file.}
#'   \item{sumFile  }{Path to summary file.}
#'   \item{counts   }{Matrix with genes in rows and samples in columns.}
#'   \item{summary  }{Matrix summarizing the counting result. Samples
#'                    in rows, various counting values in columns.}
#'   \item{genes    }{Vector of counted genes.}
#'   \item{anno     }{Contains additional information (like 'genes',
#'                    'Length') found in the annotation file supplied
#'                    to featureCounts.}
#'   \item{calls    }{Vector of all calls to featureCounts.}
#'   \item{version  }{Version of Rsubread (featureCounts' R package).}
#'   \item{tool     }{The name of the tool.}
#' }
#'
#' @examples
#' \dontrun{
#'   run_featureCounts(
#'     files = c("file1.fastq", "file2.fastq"),
#'     annotation = "anno.gtf",
#'     IDtype = "gene_id",
#'     featureType = "exon",
#'     cpus = 20,
#'     workers = 10,
#'   )
#' }
#' run_featureCounts
run_featureCounts <- function(
    files,
    annotation,
    outDir = "counting",
    isGTF = TRUE,
    IDtype = "gene_id",
    featureType = "exon",
    strandSpecific = "0",
    isPairedEnd = FALSE,
    useMetaFeatures   = TRUE,
    allowMultiOverlap = FALSE,
    countMultiMapped  = FALSE,
    cpus = 1,
    workers = 10
) {
    driver <- function(i) {
        # HACK: featureCounts prints a lot of text. Remove that.
        invisible(capture.output(
            res <- Rsubread::featureCounts(
                files                  = files[i],
                annot.ext              = annotation,
                isGTFAnnotationFile    = isGTF,
                GTF.attrType           = IDtype,
                GTF.featureType        = featureType,
                strandSpecific         = strandSpecific,
                isPairedEnd            = isPairedEnd,
                requireBothEndsMapped  = TRUE,
                useMetaFeatures        = useMetaFeatures,
                allowMultiOverlap      = allowMultiOverlap,
                countMultiMappingReads = countMultiMapped,
                nthreads               = allocCPUS$for_call[i]
            )
        ))

        inputName <- tools::file_path_sans_ext(basename(files[i]))
        outAnnotFile <- file.path(outDir, paste("fc", "annotation", "tsv", sep = "."))
        if (i == 1)
            write.table(res$annotation,
                        outAnnotFile,
                        sep = "\t",
                        quote = F,
                        row.names = F)
        names(res$stat)[[2]] <- paste0(inputName, ".bam")
        sumFile <- file.path(splitPath, paste0(inputName, ".summary"))
        write.table(
            res$stat,
            sumFile,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE
        )

        colnames(res$counts) <- c("counts")
        # we need the call for statistics
        call <- paste0(
            "featureCounts(files=", shQuote(files[i]),
            ", annot.ext=", shQuote(annotation),
            ", isGTFAnnotationFile=", isGTF,
            ", GTF.attrType=", shQuote(IDtype),
            ", GTF.featureType=", shQuote(featureType),
            ", strandSpecific=", strandSpecific,
            ", isPairedEnd=", isPairedEnd,
            ", requireBothEndsMapped=TRUE",
            ", useMetaFeatures=", useMetaFeatures,
            ", allowMultiOverlap=", allowMultiOverlap,
            ", countMultiMappingReads=", countMultiMapped,
            ", nthreads=", allocCPUS$for_call[i],
            ")"
        )
        return(
            list(
                counts     = res$counts,
                annotation = res$annotation,
                summary    = res$stat,
                call       = call
            )
        )
    }

    # check if input files exist
    if (FALSE %in% file.exists(asPairVector(files))) {
        f <- asPairVector(files)
        stop(paste(
            "Files do not exist:\n",
            paste(shQuote(f[!file.exists(f)]),
            collapse = "\n ")
        ))
    }

    message("featureCounts running... Paired: ", isPairedEnd)
    version <- as.character(packageVersion("Rsubread"))
    splitPath <- file.path(outDir, "split")
    dir.create(splitPath, recursive = TRUE, showWarnings = FALSE)

    # NOTE: In paired mode, only one file at a time! BAM files are sorting in
    #       working directory. Parallel instances will cause errors.
    allocCPUS <- allocate_cpus(cpus, length(files), if (isPairedEnd) 1 else workers)
    bio_par <- BiocParallel::MulticoreParam(
        workers = allocCPUS$in_parallel,
        tasks = length(files),
        stop.on.error = TRUE,
        progressbar = TRUE
    )
    allRes <- tryCatch({
        BiocParallel::bplapply(1:length(files), driver, BPPARAM = bio_par)
    }, error = identity)
    .stop_at_error(allRes, from = match.call()[[1]])

    # combine summaries and counts to a matrix, respectively
    inputNames <- basename(files)
    genes      <- allRes[[1]]$annotation$GeneID

    summary <- sapply(allRes, function(x) {x[["summary"]][[2]]} )
    rownames(summary) <- allRes[[1]]$summary[[1]]
    colnames(summary) <- inputNames
    summary <- t(summary)
    sumfile <- file.path(outDir, "fc.summary.tsv")
    write.table(summary, sumfile, sep = "\t", col.names = NA, quote = F)

    counts   <- sapply(allRes, function(x) x[["counts"]])
    rownames(counts) <- genes
    colnames(counts) <- inputNames
    countfile <- file.path(outDir, "fc.counts.tsv")
    write.table(counts, countfile, sep = "\t", col.names = NA, quote = F)

    return(
        list(
            countFile = countfile,
            sumFile   = sumfile,
            counts    = counts,
            summary   = summary,
            genes     = genes,
            anno      = allRes[[1]]$annotation,
            calls     = sapply(allRes, function(x) {x[["call"]]} ),
            version   = version,
            tool      = "featureCounts"
        )
    )
}
