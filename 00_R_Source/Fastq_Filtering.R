####-
##-====================================================
##-==========================  TRIMMING ===========================
##-====================================================
####-

#' Locked Arguments For Trimmomatic
#'
#' Show all locked arguments for Trimmomatic.
#'
#' These arguments are used by GEO2RNAseq's \code{\link{run_Trimmomatic}}
#' wrapper function for the tool Trimmomatic. They cannot be accessed directly
#' as arguments to the wrapper function. Instead, they are set automatically.
#'
#' @export
#' @param paired Logical determining whether paired-end reads are used.
#'   Defaults to FALSE.
#' @return Character(n).
#' @examples args_Trimmomatic_locked(FALSE)
args_Trimmomatic_locked <- function
(
    paired=FALSE
){
    if (paired) {
        return(c("trimmomatic", "PE", "-threads",
            "ILLUMINACLIP", "LEADING", "TRAILING", "SLIDINGWINDOW", "MINLEN"
        ))
    } else {
        return(c("trimmomatic", "SE", "-threads",
            "ILLUMINACLIP", "LEADING", "TRAILING", "SLIDINGWINDOW", "MINLEN"
        ))
    }
}


#' Arguments For Trimmomatic
#'
#' Show arguments used for Trimmomatic by GEO2RNAseq. They are used for the
#' wrapper function \code{\link{run_Trimmomatic}}. Locked arguments must not be
#' set as additional arguments. Instead, they are set automatically.
#'
#' @export
#' @inheritParams args_Trimmomatic_locked
#' @return Nothing.
#' @examples args_Trimmomatic(FALSE)
args_Trimmomatic <- function
(
    paired=FALSE
){
    writeLines(c(
        "We use only these (locked) arguments (in the same order as shown):",
        paste(args_Trimmomatic_locked(paired), collapse = " ")
    ))
}



#' Run Trimmomatic
#'
#' Wrapper function for Trimmomatic, a tool for quality and adapter trimming of
#' RNA-seq reads.
#'
#' @export
#' @param files Matrix or list of files in FASTQ format. Matrix must have pairs
#'   in rows: first column for forward reads, second column for backward reads.
#' @param outDir Defaults to "./fastq". Output directory for trimmed FASTQ files.
#' @param is.paired Defaults to FALSE. Logical indicating whether reads are paired-end.
#' @param windowsize Defaults to 15. Windows size for sliding window trimming.
#' @param qualcut Defaults to 25. Cut-off quality for sliding window trimming.
#' @param phred Defaults to "-phred33". Encoding of the FASTQ file. Usually
#'   "-phred33" or "-phred64".
#' @param leading Defaults to 3. Remove low quality bases from the beginning.
#'   As long as a base has a value below this threshold the base is removed and
#'   next base will be investigated.
#' @param trailing Defaults to 3. Remove low quality bases from the beginning.
#'   As long as a base has a value below this threshold the base is removed and
#'   the next the base will be investigated.
#' @param minlen Defaults to 30. Reads shorter than this size will be removed
#'   completely.
#' @param adapters Defaults to NA. Optional. Path to adapter FASTA file.
#' @inheritParams parallelizationParameters
#' @param overwrite Defaults to FALSE. Logical indicating if Trimmomatic should
#'   overwrite files which already exist. If FALSE, the path to the existing
#'   file be will returned along with a message.
#' @param use.existing Defaults to TRUE. If set to FALSE, execution will stop
#'   if one or more output files already exist. No trimming is performed
#'   in that case.
#' @param compress If TRUE, output files will be gz compressed. Defaults to
#'   FALSE.
#' @param trimoJar Defaults to the empty string "". Path to the .jar file of
#'   Trimmomatic.
#'
#' @references Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A
#'   flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.
#' @seealso Available online at: \url{http://www.usadellab.org/cms/?page=trimmomatic}\cr
#'   FASTQ format on Wikipedia \url{https://en.wikipedia.org/wiki/FASTQ_format}
#' @return A list with keyword arguments:\cr
#' \describe{
#'   \item{files    }{Paths to output files.\cr
#'                    For is.paired=FALSE, this is a list.\cr
#'                    For is.paired=TRUE, this is a 2 column matrix. Row i is read-pair i.}
#'   \item{calls    }{System calls as strings.}
#'   \item{input    }{Number of reads read (per file).}
#'   \item{surviving}{Number of reads remaining after trimming (per file).}
#'   \item{log      }{Path to log files.}
#'   \item{tool     }{The name of the tool.}
#' }
#' If 'outDir' contains trimmed files, overwrite is FALSE and use.existing
#' is TRUE, only files with suffix 'trimo.fastq' or 'trimo.pe.fastq' will be
#' returned.
#' @examples
#' \dontrun{
#'   files <- c("f1.fastq", "f2.fastq")
#'   run_Trimmomatic(files = files)
#'
#'   files <- c("f1_1.fastq", "f1_2.fastq", "f2_1.fastq", "f2_2.fastq")
#'   run_Trimmomatic(files, paired = TRUE)
#' }
#' run_Trimmomatic
run_Trimmomatic <- function
(
    files,
    outDir = "./fastq",
    is.paired = FALSE,
    windowsize = 15,
    qualcut = 25,
    phred = "-phred33",
    leading = 3,
    trailing = 3,
    minlen = 30,
    adapters = NA,
    cpus = 1,
    workers = 4,
    overwrite=FALSE,
    use.existing=TRUE,
    compress = FALSE,
    trimoJar = ""
){
    if (length(list.files(outDir, pattern = "\\.trimo(\\.pe)*\\.fastq$")) > 0)
        message(paste0("Directory \"", outDir, "\" not empty! overwrite? ", overwrite))
    # check if input files exist
    if (FALSE %in% file.exists(asPairVector(files))) {
        f <- asPairVector(files)
        stop(paste(
            "Files do not exist:\n",
            paste(shQuote(f[!file.exists(f)]), collapse = "\n ")
        ))
    }

    if (!use.existing) {
      outNames <- paste0(sub("(fastq|fq)(\\.gz)?$", "", basename(files)), "trimo.fastq")
      outFiles <- file.path(outDir, outNames)
      if (TRUE %in% file.exists(outFiles))
            stop("One or more files already exist and use.existing is FALSE!")
    }

    writeLines("Trimming ...")
    if (is.paired) {
        trimmingres <- .run_Trimmomatic_paired(
            trimoJar = trimoJar,
            pairedFiles = files,
            outDir = outDir,
            cpus = cpus,
            workers = workers,
            overwrite = overwrite,
            windowsize = windowsize,
            qualcut = qualcut,
            phred = phred,
            leading = leading,
            trailing = trailing,
            minlen = minlen,
            adapters = adapters,
            compress = compress
        )
    } else {
        trimmingres <- .run_Trimmomatic_unpaired(
            trimoJar = trimoJar,
            fastqFiles = files,
            outDir = outDir,
            cpus = cpus,
            workers = workers,
            overwrite = overwrite,
            windowsize = windowsize,
            qualcut = qualcut,
            phred = phred,
            leading = leading,
            trailing = trailing,
            minlen = minlen,
            adapters = adapters,
            compress = compress
        )
    }

    return(trimmingres)
}


#' Get Trimmomatic Adapters
#'
#' Return FASTA file containing the adapter sequences used by Trimmomatic. The
#' location of the adapter file depends on the location of the Trimmomatic jar
#' file.
#'
#' If you changed the default path of the adapters, supply the new path using
#' the 'adapters' argument of the \code{\link{run_Trimmomatic}} wrapper function.
#'
#' @export
#' @inheritParams args_Trimmomatic_locked
#' @param trimoJar Path to trimmomatic jar, if it cannot be found in your
#'   PATH environment.
#' @seealso \code{\link{list_executables}}
#' @return Character(1). Path to adapter file.
get_Trimmomatic_adapters <- function
(
    trimoJar = NA,
    paired = FALSE
) {
    if (is.na(trimoJar) || is.null(trimoJar))
        trimoJar <- ""

    adapters <- if (paired)
        file.path("adapters", "TruSeq3-PE-2.fa")
    else
        file.path("adapters", "TruSeq3-SE.fa")

    if (file.exists(file.path(trimoJar, adapters)))
        return(file.path(trimoJar, adapters))
    else if (file.exists(file.path(dirname(trimoJar), adapters)))
        return(file.path(dirname(trimoJar), adapters))
    else {
        message(
            paste0(
                "Trimmomatic installation does not contain adapter files. Using ",
                "adapters from Trimmomatic version 0.36."
            )
        )
        return(file.path(
            system.file("extdata", package = "Geo2RNAseq"),
            adapters
        ))
    }
}


# Run Trimmomatic in single-end mode
#
# Single-End mode only.
#
# @rdname run_Trimmomatic_unpaired
# @inherit run_Trimmomatic
# @param fastqFiles List of paths to FASTQ files.
.run_Trimmomatic_unpaired <- function
(
    trimoJar,
    fastqFiles,
    outDir = "",
    windowsize = 15,
    qualcut = 25,
    phred = "-phred33",
    leading = 3,
    trailing = 3,
    minlen = 30,
    adapters = NA,
    cpus = 1,
    workers = 4,
    compress = FALSE,
    overwrite = FALSE
){
    driver <- function(i) {
        resDir <- if (outDir == "" || is.na(outDir) || is.null(outDir))
            dirname(fastqFiles[i])
        else
            outDir

        dir.create(resDir, recursive = TRUE, showWarnings = FALSE)

        # NOTE: Trimmomatic can read 'fastq' and also 'fastq.gz' files
        file_name_core <- sub("\\.(fastq|fq)(\\.gz)?'?$", "", basename(fastqFiles[i]))
        end <- if (compress) ".trimo.fq.gz" else ".trimo.fastq"
        trimmed_file <- file.path(resDir, paste0(file_name_core, end))
        log_file     <- file.path(resDir, paste0(file_name_core, ".trimlog"))

        # the file might already exist
        if (file.exists(trimmed_file) && !overwrite) {
            message(
                "File ", shQuote(trimmed_file),
                " already exists. Skipped for trimming."
            )
            # if the log file exists, parse the information from there
            if (file.exists(log_file)) {
                out <- readLines(log_file)
                reads <- grep("Input Reads: ", out, ignore.case = TRUE, value = TRUE)
                input <- as.numeric(sub("Input Reads: (\\d+).*", "\\1", reads))
                surviving <- as.numeric(sub(".*Surviving: (\\d+).*", "\\1", reads))
                return(list(
                    files = trimmed_file,
                    call = paste("SE", out[[2]]),
                    input = input,
                    surviving = surviving,
                    log = log_file
                ))
            } else {
                return(list(
                    files = trimmed_file,
                    call = "NOT USED",
                    input = NA,
                    surviving = NA,
                    log = log_file
                ))
            }
        }

        if (binType == "jar")
            call <- paste("java -jar", trimoJar, "SE")
        else if (binType == "bashSE")
            call <- trimoJar
        else
            call <- paste(trimoJar, "SE")

        # NOTE: Trimmomatic executes operations in order
        call <- paste(
            call,
            "-threads",
            allocCPUS$for_call[i],
            phred,
            fastqFiles[i],
            trimmed_file,
            paste0("ILLUMINACLIP:", adapters, ":2:30:10"),
            paste0("LEADING:", leading),
            paste0("TRAILING:", trailing),
            paste0("SLIDINGWINDOW:", windowsize, ":", qualcut),
            paste0("MINLEN:", minlen),
            "2>&1"
        )

        message(call)

        # capture number of reads after trimming from Trimmomatic output
        out <- system(call, intern = TRUE)
        write(out, file = log_file)
        reads <- grep("Input Reads: ", out, ignore.case = TRUE, value = TRUE)
        input <- as.numeric(sub("Input Reads: (\\d+).*", "\\1", reads))
        surviving <- as.numeric(sub(".*Surviving: (\\d+).*", "\\1", reads))

        return(list(
            file = trimmed_file,
            call = call,
            input = input,
            surviving = surviving,
            log = log_file
        ))
    }

    # Trimmomatic can be found on the system in two ways:
    # (1) a .jar to be executed
    # (2) as BASH script (usually when installed via apt-get install)
    binType <- if (!is.null(trimoJar) && grepl("\\.jar$", trimoJar)) "jar" else "none"
    if (trimoJar == "" || is.null(trimoJar) || is.na(trimoJar)) {
        trimoJar <- tryCatch({
            ex <- .get_executable("trimmomatic", sysVar = "TRIMMOMATIC_EXEC")
            binType <- if (grepl("\\.jar$", ex)) "jar" else "bash"
            ex
        }, error = function(e) {
            NA
        })

        if (is.na(trimoJar)) {
            trimoJar <- tryCatch({
                ex <- .get_executable("TrimmomaticSE", sysVar = "TRIMMOMATIC_EXEC")
                # workaround: TrimmomaticSE is a link to PE,
                # but the PE script checks for the corresponding call
                substr(ex, nchar(ex)-1, nchar(ex)-1) <- "S"
                binType <- "bashSE"
                ex
            }, error = function(e) {
                stop("Cannot find trimmomatic binary! Consider using 'trimoJar' ",
                     "argument or setting 'TRIMMOMATIC_EXEC' as system variable.")
            })
        }
    }

    if (is.na(adapters))  adapters <- get_Trimmomatic_adapters(trimoJar = trimoJar, paired = FALSE)

    allocCPUS <- allocate_cpus(cpus, length(fastqFiles), workers)
    allRes <- parallel::mclapply(1:length(fastqFiles), driver, mc.cores = allocCPUS$in_parallel)

    return(list(
        files     = sapply(allRes, function(x) x[["file"]]),
        calls     = sapply(allRes, function(x) x[["call"]]),
        input     = sapply(allRes, function(x) x[["input"]]),
        surviving = sapply(allRes, function(x) x[["surviving"]]),
        log       = sapply(allRes, function(x) x[["log"]]),
        tool      = "trimmomatic"
    ))
}


# Run Trimmomatic in paired-end mode
#
# Paired-end mode only.
#
# @rdname run_Trimmomatic_paired
# @inherit run_Trimmomatic
# @param pairedFiles Two column matrix or list of paths to FASTQ files.
.run_Trimmomatic_paired <- function
(
    trimoJar,
    pairedFiles,
    outDir = "",
    windowsize = 15,
    qualcut = 25,
    phred = "-phred33",
    leading = 3,
    trailing = 3,
    minlen = 30,
    adapters = NA,
    cpus = 1,
    workers = 2,
    compress = FALSE,
    overwrite=FALSE
) {
    driver <- function(i) {
        resDir = if (outDir == "" || is.na(outDir) || is.null(outDir))
            dirname(pairedFiles[i, 1])
        else
          outDir

        dir.create(resDir, recursive = TRUE, showWarnings = FALSE)

        file_name_core_1 <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(pairedFiles[i, 1]))
        end <- if (compress) ".fq.gz" else ".fastq"
        trimmed_file_PE_1 <- file.path(resDir, paste0(file_name_core_1, ".trimo.pe", end))
        trimmed_file_SE_1 <- file.path(resDir, paste0(file_name_core_1, ".trimo.se", end))
        file_name_core_2  <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(pairedFiles[i, 2]))
        trimmed_file_PE_2 <- file.path(resDir, paste0(file_name_core_2, ".trimo.pe", end))
        trimmed_file_SE_2 <- file.path(resDir, paste0(file_name_core_2, ".trimo.se", end))
        log_file          <- file.path(resDir, paste0(file_name_core_1, ".trimlog"))

        # the files might already exist
        if (file.exists(trimmed_file_PE_1) && file.exists(trimmed_file_PE_2) && !overwrite) {
            message(
                "Files ",
                shQuote(paste0(trimmed_file_PE_1, "; ", trimmed_file_PE_2)),
                " already exist. Skipped for trimming."
            )
            # if the log file exists, parse the information from there
            if (file.exists(log_file)) {
                out <- readLines(log_file)
                reads <- grep("Input Read", out, ignore.case = TRUE, value = TRUE)
                input <- as.numeric(sub("Input Read Pairs: (\\d+).*", "\\1", reads))
                input <- rep(input, each = 2)
                names(input) <- pairedFiles[i,]
                surviving <- as.numeric(sub(".*Both Surviving: (\\d+).*", "\\1", reads))
                surviving <- rep(surviving, each = 2)
                names(surviving) <- pairedFiles[i,]
                return(list(
                    files = c(trimmed_file_PE_1, trimmed_file_PE_2),
                    call = paste("PE", out[[2]]),
                    input = input,
                    surviving = surviving,
                    log = log_file
                ))
            } else {
                return(list(
                    files = c(trimmed_file_PE_1, trimmed_file_PE_2),
                    call = "NOT USED",
                    input = c(NA,NA),
                    surviving = c(NA,NA),
                    log = log_file
                ))
            }
        }

        if (binType == "jar")
            call <- paste("java -jar", trimoJar, "PE")
        else if (binType == "bashPE")
            call <- trimoJar
        else
            call <- paste(trimoJar, "PE")

        # NOTE: Trimmomatic executes operations in order
        call <- paste(
            call,
            "-threads",
            allocCPUS$for_call[i],
            phred,
            pairedFiles[i, 1],
            pairedFiles[i, 2],
            trimmed_file_PE_1,
            trimmed_file_SE_1,
            trimmed_file_PE_2,
            trimmed_file_SE_2,
            paste0("ILLUMINACLIP:", adapters, ":2:30:10"),
            paste0("LEADING:", leading),
            paste0("TRAILING:", trailing),
            paste0("SLIDINGWINDOW:", windowsize, ":", qualcut),
            paste0("MINLEN:", minlen),
            "2>&1"
        )

        # capture number of reads after trimming from Trimmomatic output.
        # NOTE: different to SE version ...
        out <- system(call, intern = TRUE)
        write(out, file <- log_file)
        reads <- grep("Input Read", out, ignore.case = TRUE, value = TRUE)
        input <- as.numeric(sub("Input Read Pairs: (\\d+).*", "\\1", reads))
        input <- rep(input, each = 2)
        names(input) <- pairedFiles[i,]
        surviving <- as.numeric(sub(".*Both Surviving: (\\d+).*", "\\1", reads))
        surviving <- rep(surviving, each = 2)
        names(surviving) <- pairedFiles[i,]

        resfiles <- c(trimmed_file_PE_1, trimmed_file_PE_2)
        return(list(
            files = resfiles,
            call = call,
            input = input,
            surviving = surviving,
            log = log_file
        ))
    }

    # Trimmomatic can be found on the system in two ways:
    # (1) a .jar to be executed
    # (2) as Bash script (usually when installed via apt-get install)
    binType <- if (!is.null(trimoJar) && grepl("\\.jar$", trimoJar)) "jar" else "none"
    if (trimoJar == "" || is.null(trimoJar) || is.na(trimoJar)) {
        trimoJar <- tryCatch({
            ex <- .get_executable("trimmomatic", sysVar = "TRIMMOMATIC_EXEC")
            binType <- if (grepl("\\.jar$", ex)) "jar" else "bash"
            ex
        }, error = function(e) {
            NA
        })

        if (is.na(trimoJar)) {
            trimoJar <- tryCatch({
                ex <- .get_executable("TrimmomaticPE", sysVar = "TRIMMOMATIC_EXEC")
                # workaround: TrimmomaticSE is a link to PE,
                # but the PE script checks for the corresponding call
                substr(ex, nchar(ex)-1, nchar(ex)-1) <- "S"
                binType <- "bashPE"
                ex
            }, error = function(e) {
                stop("Cannot find trimmomatic binary! Consider using 'trimoJar' ",
                     "argument or setting 'TRIMMOMATIC_EXEC' as system variable.")
            })
        }
    }

    if (is.vector(pairedFiles) || is.list(pairedFiles))
        pairedFiles <- asPaired(pairedFiles)
    if (is.na(adapters))
        adapters <- get_Trimmomatic_adapters(trimoJar = trimoJar, paired = TRUE)

    allocCPUS <- allocate_cpus(cpus, nrow(pairedFiles), workers)
    allRes <- parallel::mclapply(1:nrow(pairedFiles), driver, mc.cores = allocCPUS$in_parallel)

    input <- unlist(sapply(allRes, function(x) x[["input"]], simplify=FALSE))
    names(input) <- asPairVector(pairedFiles)
    surviving <- unlist(sapply(allRes, function(x) x[["surviving"]], simplify=FALSE))
    names(surviving) <- asPairVector(pairedFiles)
    return(list(
        files     = t(sapply(allRes, function(x) x[["files"]])),
        calls     = sapply(allRes, function(x) x[["call"]]),
        input     = input,
        surviving = surviving,
        log       = sapply(allRes, function(x) x[["log"]]),
        tool      = "trimmomatic"
    ))
}



####-
##-====================================================
##-========================  rRNA REMOVAL =========================
##-====================================================
####-


#' Run SortMeRNA Indexer
#'
#' Index FASTA sequences, which are used by SortMeRNA as rRNA database. This
#' needs to be done only once.
#'
#' If using SortMeRNA the first time, the reference FASTA sequences must be
#' indexed first, before they can be used. Usually, they are shipped together
#' with SortMeRNA and can be found in the 'rRNA_databases' subdirectory. The
#' index files will be stored in the 'index' subdirectory.
#' \cr\cr
#' Unless defined otherwise, this function will try to find the indexing
#' executable and reference sequence files automatically.
#' \cr\cr
#' See the documentation of \code{list_executable} for further advice of how to
#' set up external tools.
#'
#' @export
#' @param refDir Path to directory containing rRNA reference FASTA files. If NA,
#'   the function searches for this directory based on the detected SortMeRNA
#'   binary.
#' @return Character(n). System calls used for indexing.
#' @examples
#' \dontrun{
#'   run_sortmerna_indexer()
#'   run_sortmerna_indexer("path/to/sortmernaDirectory")
#' }
#' run_sortmerna_indexer
run_sortmerna_indexer <- function
(
    refDir = NA
){
  # find and create the database string
    db <- function(idx) {
        message("Indexing: ", idx)
        flush.console()
        sortmeDir <- dirname(.get_executable("sortmerna", sysVar = "SORTMERNA_EXEC"))
        if (is.na(refDir) || is.null(refDir) || refDir == "")
            refDir <- .find_upward(sortmeDir, "rRNA_databases")
        idxDir <- normalizePath(paste0(refDir, "/../index"))
        dir.create(idxDir, recursive = TRUE, showWarnings = FALSE)
        if (!file.exists(refDir) || !file.exists(idxDir))
            stop(paste0(
                "Cannot find SortMeRNA databases directory. Make sure to",
                "put database in:\n ", refDir
            ))
        ref <- grep(paste0(idx, ".*\\.fasta"), dir(refDir), value = TRUE)[1]
        if (!file.exists(file.path(refDir, ref)))
            stop(paste0(
                "Could not find database FASTA file: ",
                file.path(refDir, ref)
            ))

        return(paste(
            file.path(refDir, ref),
            file.path(idxDir, idx),
            sep = ","
        ))
    }

    merge <- function(db_list) {
        message("Merging SortMeRNA databases for 'fast' mode. See 'run_SortMeRNA' for more information.")
        sortmeDir <- dirname(.get_executable("sortmerna", sysVar = "SORTMERNA_EXEC"))
        if (is.na(refDir) || is.null(refDir) || refDir == "")
            refDir <- .find_upward(sortmeDir, "rRNA_databases")
        target_file <- file.path(refDir, "all_dbs_combined.fasta")
        if (file.exists(target_file))
            file.remove(target_file)
        file.create(target_file)

        for (i in 1:length(db_list)) {
            ref <- grep(paste0(db_list[i], ".*\\.fasta"), dir(refDir), value = TRUE)[1]
            ref <- file.path(refDir, ref)
            if (!file.exists(ref))
                stop(paste0("Could not find database FASTA file: ", ref))
            else
                system(paste("cat", ref, ">>", target_file))
        }
    }

    driver <- function(idx) {
        exec <- .get_executable(
            "indexdb_rna",
            derive = "sortmerna",
            sysVar = "SORTMERNA_EXEC",
            only.dir = TRUE
        )
        call <- paste(exec, "--ref", db(idx))
        system(call)
        return(call)
    }

    dbs <- c(
        "silva-bac-23s", "silva-bac-16s",
        "silva-arc-23s", "silva-arc-16s",
        "silva-euk-28s", "silva-euk-18s",
        "rfam-5.8s"    , "rfam-5s"
    )

    # merge first
    merge(dbs)
    dbs <- c(dbs, "all_dbs_combined")
    # start indexing
    res <- lapply(dbs, driver)
    return(res)
}


# SortMeRNA --> alternative to riboPicker
# TODO PACKAGE: skip merging of paired-end if split files already exist?
# ATTENTION: first, build index with indexdb_rna
# has to be done once for a fresh SortMeRNA install
# NOTE: --num-alignments = 1 leads to same behavior as --feeling-lucky
# NOTE: fastest setup for arguments
# NOTE: concatenating the databases is significantly faster at the cost of sensitivity
# NOTE: paired-end reads must be merged into a single file (to maintain order)
#       or synchronized afterwards

#' Run SortMeRNA
#'
#' Wrapper function for SortMeRNA, a tool for the removal of rRNA sequences from
#' RNA-seq reads (and more).
#'
#' SortMeRNA uses several, different database files. It is significantly faster
#' to merge all databases of interest and use this single 'merged' database,
#' instead of using them one after another. However, this also leads to a small
#' loss of sensitivity. Only the 'fast' mode merges databases. More precisely,
#' it will merge ALL databases shipped with SortMeRNA.
#'
#' @export
#' @references Kopylova E., Noe L. and Touzet H., "SortMeRNA: Fast and accurate
#'   filtering of ribosomal RNAs in metatranscriptomic data", Bioinformatics
#'   (2012), doi: 10.1093/bioinformatics/bts611. Available online at:
#'   \url{http://bioinfo.lifl.fr/RNA/sortmerna/}
#' @param files Input files in FASTQ format.
#' @param outDir Defaults to "fastq". Main directory for output files.
#' @param mode String. Defaults to "fast". There are several options to
#'   influence the speed and sensitivity of SortMeRNA:
#'   \cr
#'   \emph{fast}: check all rRNAs with less sensitivity
#'   \cr
#'   \emph{all}: check all rRNAs with default sensitivity
#'   \cr
#'   \emph{bac}: only bacteria
#'   \cr
#'   \emph{euk}: only eukaryotes
#'   \cr
#'   \emph{arc}: only archea
#' @param paired Logical(1) indicating whether paired-end reads are used.
#'   Defaults to FALSE.
#' @param removeMerged For paired=TRUE only. Logical(1) indicating whether
#'   intermediate merge files should be removed after completion. Defaults to
#'   FALSE. Only for paired data. If TRUE, merged files will be removed after
#'   ALL of them were split into paired-end files again.
#' @param split Logical(1) Defaults to TRUE.
#' @param overwrite Defaults to FALSE. Logical. If TRUE, all existing files
#'   will be overwritten. If FALSE, existing files will be returned.
#' @inheritParams parallelizationParameters
#' @return A data.frame with the following keyword arguments:\cr
#' \describe{
#'   \item{rrna   }{Number of rRNA reads.}
#'   \item{nonrrna}{Number of non-rRNA reads.}
#'   \item{files  }{Single-end: Filtered output files (without rRNA), vector.\cr
#'                  Paired-end: Filtered output files (without rRNA), 2 column matrix.}
#'   \item{files1 }{Paired-end. Filtered output files (without rRNA), forward-reads.}
#'   \item{files2 }{Paired-end. Filtered output files (without rRNA), reverse-reads.}
#'   \item{calls  }{Strings of system calls for SortMeRNA.}
#'   \item{tool   }{The name of the tool.}
#' }
#' @examples
#' \dontrun{
#'   files <- c("f1.fastq", "f2.fastq")
#'   run_SortMeRNA(files = files)
#'
#'   files <- c("f1_1.fastq", "f1_2.fastq", "f2_1.fastq", "f2_2.fastq")
#'   run_SortMeRNA(files = files, paired = TRUE)
#' }
#' run_SortMeRNA
run_SortMeRNA <- function
(
    files,
    outDir = "fastq",
    mode = "fast",
    paired = FALSE,
    removeMerged = FALSE,
    split = TRUE,
    overwrite = FALSE,
    cpus = 1,
    workers = 3
){
    merge_paired <- function(i, pairedFiles) {
        base <- tools::file_path_sans_ext(basename(pairedFiles[i, 1]))
        mergedFile <- file.path(outDir, paste0(base, ".merged.fastq"))
        if (!overwrite && file.exists(mergedFile)) {
            writeLines(paste("Use existing:", mergedFile))
            return(mergedFile)
        } else {
            writeLines(paste("Worker: ", i, "with files:", pairedFiles[i,]))
            sortmeDir <- dirname(.get_executable("sortmerna", sysVar = "SORTMERNA_EXEC"))
            merger <- .find_upward(sortmeDir, file.path("scripts", "merge-paired-reads.sh"))

            # unzip compressed files before merging
            if (grepl("\\.(gz|bz2)$", pairedFiles[i,1])) {
                pairedFiles[i,1] <- gunzip(pairedFiles[i,1])[1]
            } else if (grepl("\\.zip$", pairedFiles[i,1])) {
                pairedFiles[i,1] <- unzip(pairedFiles[i,1])
            }
            if (grepl("\\.(gz|bz2)$", pairedFiles[i,2])) {
                pairedFiles[i,2] <- gunzip(pairedFiles[i,2])[1]
            } else if (grepl("\\.zip$", pairedFiles[i,2])) {
                pairedFiles[i,2] <- unzip(pairedFiles[i,2])
            }

            call <- paste("bash", merger, pairedFiles[i,1], pairedFiles[i,2], mergedFile)
            system(call)
            return(mergedFile)
        }
    }

    split_paired <- function(i, mergedFiles, inFileNames) {
        base1 <- tools::file_path_sans_ext(inFileNames[2*i - 1])
        base2 <- tools::file_path_sans_ext(inFileNames[2*i    ])
        out1 <- file.path(outDir, paste0(base1, ".non_rrna.fastq"))
        out2 <- file.path(outDir, paste0(base2, ".non_rrna.fastq"))
        if (!overwrite && (file.exists(out1) && file.exists(out2))) {
            writeLines(paste("Use existing:", c(out1, out2)))
            return(c(out1, out2))
        } else {
            writeLines(paste("Worker: ", i, "with file:", mergedFiles[i]))
            sortmeDir <- dirname(.get_executable("sortmerna", sysVar = "SORTMERNA_EXEC"))
            unmerger <- .find_upward(sortmeDir, file.path("scripts", "unmerge-paired-reads.sh"))

            call <- paste("bash", unmerger, mergedFiles[i], out1, out2)
            system(call)
            return(c(out1, out2))
        }
    }

    # helper function to compile string for reference and index file
    sortmerna_db <- function(idx) {
        sortmeDir <- dirname(.get_executable("sortmerna", sysVar = "SORTMERNA_EXEC"))
        refDir <- .find_upward(sortmeDir, "rRNA_databases")
        idxDir <- .find_upward(sortmeDir, "index")
        if (!file.exists(refDir) || !file.exists(idxDir))
          # TODO PACKAGE: if database indexer works, update this error message!
          stop(paste0("Cannot find SortMeRNA databases or indices. Make sure to",
                      "put database in:\n ", refDir, "\nAnd index files in:\n ",
                      idxDir, "\n"))
        ref <- grep(paste0(idx, ".*\\.fasta"), dir(refDir), value = TRUE)[1]
        # ref has to be exact file name with extension
        # idx not, must only share basename with ref

        return(paste(
          file.path(refDir, ref),
          file.path(idxDir, idx),
          sep = ","
        ))
    }

    # helper function to get databses based on the chosen 'mode'
    get_db <- function(mode) {
        if (mode == "fast")
            # concatenated DB --> much faster but a little bit less sensitive
            return(sortmerna_db("all_dbs_combined"))
        else {
            dbs  <- c(
                "silva-bac-23s", "silva-bac-16s",
                "silva-arc-23s", "silva-arc-16s",
                "silva-euk-28s", "silva-euk-18s",
                "rfam-5.8s", "rfam-5s"
            )
            if (mode == "all")
                return( paste(lapply(dbs,           sortmerna_db), collapse=":") )
            else if (mode == "bac")
                return( paste(lapply(dbs[c(1:2,8)], sortmerna_db), collapse=":") )
            else if (mode == "arc")
                return( paste(lapply(dbs[c(3:4,8)], sortmerna_db), collapse=":") )
            else if (mode == "euk")
                return( paste(lapply(dbs[c(5:7,8)], sortmerna_db), collapse=":") )
            else
                stop(paste("Unknown mode for SortMeRNA:", mode))
        }
    }


    driver <- function(i) {
        writeLines(paste("\nRunner", i, "..."))
        prefix     <- tools::file_path_sans_ext(basename(files[[i]]))
        norRNAfile <- file.path(outDir, paste0(prefix, ".non_rrna"))
        rRNAfile   <- file.path(subDir, paste0(prefix, ".rrna"))
        log <- paste0(rRNAfile, ".log")

        if (!overwrite && !(FALSE %in% file.exists(c(paste0(norRNAfile, ".fastq"), log)))) {
            writeLines(paste("Use existing:", norRNAfile))
            call <- "NOT USED"
        } else {
            args <- paste(
                "--ref"           , get_db(mode),
                "--reads"         , files[i],
                "--num_alignments", 1,
                "--aligned"       , rRNAfile,
                "--other"         , norRNAfile,
                "--fastx",
                "--log",
                "-m"              , mem,
                "-a"              , allocCPUS$for_call[i]
            )
            # input is a *merged* FASTQ file. If one read is aligned, reject both.
            if (paired) args <- paste(args, "--paired_in")

            # actual SortMeRNA call
            # ignore stdout, but save stderr
            call <- paste(sortmerna_exec, args)
            stderr <- system(call, intern = TRUE, ignore.stdout = TRUE)

            # stop and print error if error occurred
            if (length(stderr) != 0) stop(stderr)
        }

        # get number of rRNA and non-rRNA reads from logfile. Log is always in the same directory as the rRNA file.
        rrna <- as.integer(system(
            paste(
                "perl",
                "-ne",
                "'if (/Total reads passing E-value threshold = (\\d+)/i) { print $1 }'",
                log
            ),
            intern = TRUE
        ))
        nonrrna <- as.integer(system(
            paste(
                "perl",
                "-ne",
                "'if (/Total reads failing E-value threshold = (\\d+)/i) { print $1 }'",
                log
            ),
            intern = TRUE
        ))

        writeLines(paste("Runner", i, "finished."))
        norRNAfile <- paste0(norRNAfile, ".fastq") # SortMeRNA adds that ending itself - but we (may) need the full file path later
        return(list(
            rrna = rrna,
            nonrrna = nonrrna,
            file = norRNAfile,
            call = call
        ))
    }

    # check if input files exist
    if (FALSE %in% file.exists(asPairVector(files))) {
        f <- asPairVector(files)
        stop(paste(
            "Files do not exist:\n",
            paste(shQuote(f[!file.exists(f)]), collapse = "\n ")
        ))
    }

    # find SortMeRNA executable
    sortmerna_exec <- .get_executable("sortmerna", sysVar = "SORTMERNA_EXEC")

    # create subDir for log and rRNA files, if it does not exist
    subDir <- file.path(outDir, "sortmerna")
    dir.create(subDir, recursive = TRUE, showWarnings = FALSE)

    # paired-end files MUST be merged into a single file beforehand.
    # Merger script is comes with SortMeRNA
    if (paired) {
        # remember these to name split files correctly
        inFileNames <- basename(asPairVector(files))

        if (is.vector(files)) files <- asPaired(files)
        # NOTE: merging will only happen if the resulting file does not exist already
        writeLines("Merging paired files...")
        merged <- sapply(1:nrow(files), merge_paired, files)
        writeLines("Paired mode. Using merged files:")
        print(merged)
        files <- merged
    }

    # spawn multiple instances of run_sortmerna() in parallel
    allocCPUS <- allocate_cpus(cpus, length(files), workers)

    mem <- as.integer(.free_mem2() / allocCPUS$in_parallel)
    max_mem <- as.integer(system(
      paste(
        sortmerna_exec,
        "-h",
        "|",
        "perl",
        "-ne",
        "'if (/maximum -m INT is (\\d+)/) { print $1 }'"
      ),
      intern = TRUE
    ))
    if (mem > max_mem) mem <- max_mem

    bio_par <- BiocParallel::MulticoreParam(
        workers = allocCPUS$in_parallel,
        tasks = length(files),
        progressbar = TRUE
    )
    allRes <- tryCatch({
        BiocParallel::bplapply(1:length(files), driver, BPPARAM = bio_par)
    }, error = identity)
    .stop_at_error(allRes, from = match.call()[[1]])

    res <- list(
        rrna    = sapply(allRes, function(x) x[["rrna"]]),
        nonrrna = sapply(allRes, function(x) x[["nonrrna"]]),
        files   = sapply(allRes, function(x) x[["file"]]),
        calls   = sapply(allRes, function(x) x[["call"]])
    )
    # return number of rRNA reads, non-rRNA reads, and file name as data.frame
    # res <- as.data.frame(t(matrix(unlist(allRes), nrow = 4)), stringsAsFactors = FALSE)
    # colnames(res) <- c("rrna", "nonrrna", "files", "calls")

    if (paired && split) {
        writeLines("Splitting ...")
        mergedNonrRNAFiles <- res$files
        unmergedFiles <- t(sapply(
            1:length(mergedNonrRNAFiles),
            split_paired,
            mergedNonrRNAFiles,
            inFileNames
        ))
        res$files <- asPairVector(unmergedFiles)
        # TODO: after splitting, more reads are removed. the log files are "out of date"
        #  ~> for precise counts, count the number of lines in FASTQ files....
        res$rrna    <- rep(as.integer(res$rrna) %/% 2, each = 2)
        res$nonrrna <- rep(as.integer(res$nonrrna) %/% 2, each = 2)
        res$calls   <- rep(res$calls, each = 2)

        if (removeMerged) {
            file.remove(merged[file.exists(merged)])
            file.remove(mergedNonrRNAFiles[file.exists(mergedNonrRNAFiles)])
        }
    }

    res         <- as.data.frame(res)
    res$files   <- as.character(res$files)
    res$rrna    <- as.integer(res$rrna)
    res$nonrrna <- as.integer(res$nonrrna)
    res$calls   <- as.character(res$calls)
    res$version <- system(paste(.get_executable("sortmerna", sysVar = "SORTMERNA_EXEC"), "--version 2>&1"), intern=TRUE)[[2]]
    res$tool    <- "sortmerna"
    return(res)
}


# Free Memory
#
# Return the remaining amount of RAM in mb available on the machine.
#
# @return Integer(1). Remaining RAM in mega bytes.
.free_mem <- function() {
    mem_info <- system("free", intern = TRUE)
    value <- sub(
      ".*Cache:\\s+(\\d+)\\s+(\\d+).*",
      "\\2",
      mem_info[3],
      ignore.case = TRUE
    )
    return(as.integer(as.numeric(value) / 1024)) # kB --> divide by 1024 --> MB
}


# Free Memory 2
#
# Return the remaining amount of RAM in mb available on the machine.
#
# @return Integer(1). Remaining RAM in mega bytes.
.free_mem2 <- function() {
    mem_info <- system("free", intern = TRUE)
    arr <- strsplit(mem_info, "\\s+")
    col <- grep("available", arr[[1]], ignore.case = TRUE)
    if (length(col) == 0) {
        col <- grep("free", arr[[1]], ignore.case = TRUE)
        if (length(col) == 0)
            stop("free_mem2: could not find memory column")
    }
    if (arr[[2]][1] != "Mem:")
        stop("free mem 2: invalid row for memory search?")
    value <- arr[[2]][col]
    # NOTE: do not force usage of swap... ~ 500 MB buffer
    # NOTE: at least 1 MB
    # NOTE: kB --> divide by 1024 --> MB
    return(max(as.integer(as.numeric(value) / 1024)-500, 1))
}
