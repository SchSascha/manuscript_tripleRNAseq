####-
##-====================================================
##-==========================  QC ===========================
##-====================================================
####-


#' Make Read Statistics Using FastQC
#'
#' Wrapper function for calling FastQC, a tool for measuring the quality of
#' RNA-seq reads.
#'
#' FastQC expects FASTA quality files (FASTQ) as input. FastQC's output is an
#' extensive HTML quality report for each input file, i.e. RNA-seq sample.
#'
#' @export
#' @param files List of paths to input files.
#' @param outDir Path to output directory.
#' @param extend If TRUE, output file names will be extended by '.fastq'.
#'   Example: the file 'A.html' will be renamed to 'A.fastq.html'. This is
#'   important to create decently separated MultiQC reports.
#' @param fastqc_exec Optional. Path to FastQC binary.
#' @inheritParams parallelizationParameters
#' @return A list with keyword arguments:\cr
#' \describe{
#'   \item{calls  }{Character(n). System calls to FastQC.}
#'   \item{version}{Character(1). Version string of FastQC.}
#'   \item{tool   }{The name of the tool.}
#' }
#' @references Andrews S. (2017). FastQC: a quality control tool for high
#'   throughput sequence data.
#'   \cr
#'   Available online at:
#'   \url{http://www.bioinformatics.babraham.ac.uk/projects/fastqc}
#' @seealso FASTQ format on Wikipedia
#'   \url{https://en.wikipedia.org/wiki/FASTQ_format} or the original publication:
#'   \cr
#'   Cock, Peter JA, et al. "The Sanger FASTQ file format for sequences with
#'   quality scores, and the Solexa/Illumina FASTQ variants." Nucleic acids
#'   research 38.6 (2009): 1767-1771.
#' @examples
#' \dontrun{
#' files <- c("f1.fastq", "f2.fastq")
#' run_FastQC(files = files, outDir = "fastq", extend = FALSE)
#'
#' # after some FASTQ processing
#' files <- c("f1.trimo.fastq", "f2.trimo.fastq")
#' run_FastQC(files = files, outDir = "fastq", extend = TRUE)
#' }
#' run_FastQC
run_FastQC <- function
(
    files,
    outDir,
    cpus = 1,
    workers = 10,
    extend = FALSE,
    fastqc_exec = ""
){
    driver <- function(i){
        call <- paste(
            fastqc_exec,
            files[i],
            "--outdir" , outDir,
            "--threads", allocCPUS$for_call[i]
        )

        system(call)

        if (extend) {
            out_file <- file.path(outDir, sub("\\.(fastq|fq)(\\.gz)?$", "", basename(files[i])))
            out_file_zip <-  paste0(out_file, "_fastqc.zip")
            out_file_html <- paste0(out_file, "_fastqc.html")
            file.rename(out_file_zip, paste0(out_file, ".fastq_fastqc.zip"))
            file.rename(out_file_html, paste0(out_file, ".fastq_fastqc.html"))
        }

        return(call)
    }

    # check if input files exist
    if (FALSE %in% file.exists(asPairVector(files))) {
        f <- asPairVector(files)
        stop(paste("Files do not exist:\n", paste(shQuote(f[!file.exists(f)]), collapse = "\n ")))
    }

    # look for FastQC executable
    if (fastqc_exec == "" || is.null(fastqc_exec) || is.na(fastqc_exec))
        fastqc_exec <- .get_executable("fastqc", sysVar = "FASTQC_EXEC")
    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

    allocCPUS <- allocate_cpus(cpus, length(files), workers)
    bio_par <- BiocParallel::MulticoreParam(
        workers <- allocCPUS$in_parallel,
        tasks <- length(files)
    )
    allRes <- tryCatch({
        BiocParallel::bplapply(1:length(files), driver, BPPARAM = bio_par)
    }, error = identity)
    .stop_at_error(allRes, from = match.call()[[1]])

    version <- system(paste(fastqc_exec, "--version"), intern=TRUE)
    version <- gsub("\\w+ ", "", version)
    return(list(calls = unlist(allRes),
                version = version,
                tool = "fastqc"))
}


#' Get MultiQC Configuration File
#'
#' Return path to default MultiQC configuration file as used by GEO2RNAseq.
#'
#' @export
#' @return Character(1). Path to configuration file.
#' @examples
#' get_MultiQC_config()
get_MultiQC_config <- function() {
    return(system.file(file.path("extdata", "multiqc_config.yaml"), package = "Geo2RNAseq"))
}



#' Run MultiQC
#'
#' Wrapper function for calling the MultiQC, a tool that aggregates results
#' from many bioinformatic tools and compiles them into a single interactive
#' HTML report.
#'
#' Collecting and showing the results/reports of specific tools can be (de-)
#' activated by (de-) activating their corresponding MultiQC modules.
#'
#' @export
#' @param dir Execution path. By default, the current working directory.
#' @param outDir Output path for MultiQC files. Defaults to 'dir'.
#' @param tools Optional. Character(n). If defined, only these MultiQC modules
#'   will be used.
#' @param config Optional. Path to MultiQC config. You may also call
#'   \code{\link{get_MultiQC_config}} to get the default configuration file.
#' @param exclude Optional. Exclude these modules/tools.
#' @param force Logical indicating whether previous results of MultiQC
#'   will be overwritten. Defaults to TRUE.
#' @references Ewels, Philip, et al. "MultiQC: summarize analysis results for
#'   multiple tools and samples in a single report." Bioinformatics 32.19
#'   (2016): 3047-3048.
#' @seealso MultiQC is available online at: \url{http://multiqc.info/}.
#' @return Character(1). The system call to MultiQC.
#' @examples
#' \dontrun{
#' run_MultiQC()
#' }
#' run_MultiQC
run_MultiQC <- function
(
    dir = getwd(),
    outDir = dir,
    tools = NULL,
    config = NULL,
    exclude = c("bowtie2"),
    force = TRUE
) {
    tools_param <- ""
    if (!is.null(tools))  tools_param <- paste("-m", tools, collapse = " ")

    exclude_param <- ""
    if (!is.null(exclude))  exclude_param <- paste("-e", exclude, collapse = " ")

    config_param <- ""
    if (!is.null(config)) {
        if (file.exists(config))
            config_param <- paste("-c", shQuote(config))
        else
            stop(paste("Config file does not exist:", shQuote(config)))
    }

    # find executable
    multiqc_exec <- .get_executable("multiqc", sysVar = "MULTIQC_EXEC")

    call <- paste(
        multiqc_exec,
        if (force) "--force" else "",
        "-o",
        shQuote(outDir),
        exclude_param,
        tools_param,
        config_param,
        shQuote(dir)
    )

    # MultiQC output is on stderr, but the progress bar is still on stdout
    # --> ignoring stderr gives no output except the progress bar, which is nice
    system(call, ignore.stderr = TRUE)
    writeLines(paste(
        "See MultiQC report at",
        shQuote(file.path(dir, "multiqc_report.html"))
    ))

    return(call)
}



#' Make Flagstats
#'
#' Create BAM stats using the 'flagstat' tool from SAMtools.
#'
#' @export
#' @param files Path to BAM files.
#' @param outDir Output directory for flagstat files.
#' @inheritParams parallelizationParameters
#' @return Character(n). Paths to flagstat files created. NA, if SAMtools binary
#'   cannot be found.
#' @examples
#' \dontrun{
#' files <- c("f1.bam", "f2.bam")
#' make_flagstats(files = files, outDir = "flagStats")
#' }
make_flagstats <- function
(
    files,
    outDir,
    cpus = 1
){
    driver <- function(file) {
        out_file <- file.path(outDir, paste0(basename(file), ".flagstat"))
        call <- paste(sam_exec,
            "flagstat",
            file,
            ">",
            out_file)
        system(call)
        return(out_file)
    }

    sam_exec <- tryCatch({
        .get_executable("samtools", sysVar = "SAMTOOLS_EXEC")
    }, error = function(e) {
        message("Cannot create flagstats from SAMtools: cannot find executable on your system.")
        return(NA)
    })

    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
    allocCPUS <- allocate_cpus(cpus, length(files))
    allRes <- parallel::mclapply(
        files,
        driver,
        mc.cores = allocCPUS$in_parallel
    )
    return(unlist(allRes))
}
