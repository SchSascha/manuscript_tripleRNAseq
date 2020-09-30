####-
##-====================================================
##-=========================  HELPER ===========================
##-====================================================
####-

#' Parallelization Arguments
#'
#' Dummy function to inherit argument descriptions from.
#'
#' @keywords internal
#' @name parallelizationParameters
#' @param cpus Integer(1) > 0. Maximum number of CPU cores (threads) to use.
#' @param workers Integer(1) > 0. Maximum number of parallel I/O processes. A
#'   fraction of 'cpus' is assigned to each worker. If set to NA, workers is set
#'   to min( length(files), cpus ).
NULL


# Note: The usage could be improved with mutex locks or thread-safe counters.
#       However, BiocParallel has this feature only for R 3.4.1 or later.
# NOTE: synchronicity mutexes don't work with bplapply.

#' Allocate CPUs
#'
#' Determine good parallelization parameters for other functions, which can work
#' on multiple input files in parallel. Calculates the number of parallel calls
#' to those functions and the number of CPUs/threads per call.
#'
#' If running on a system with multiple CPUs or CPU cores available, parallel
#' execution of commands usually speeds up calculations. In the context of
#' RNA-seq data pre-processing, this is especially true if multiple input files,
#' e.g. FASTQ files or BAM files, can be processed independently in a parallel
#' manner.
#' \cr\cr
#' In general, it is faster to process more files (with less threads per file)
#' in parallel, than processing less files with more threads in parallel.
#' \cr\cr
#' The function suggests 2 or more threads per file only if the number of
#' files to process is smaller than the number of available CPUs/cores.
#'
#' @export
#' @param numThreads Integer(1). Maximum number of CPUs/threads available.
#' @param numCalls Integer(1).   Number of independent calls to execute.
#' @param numWorkers Integer(1). Maximum number of running processes. If NA,
#'   numWorkers is set to numCalls. If there are less workers than threads,
#'   more than 1 thread is assigned to some/all workers.
#' @return A list with keyword arguments: \cr
#' \describe{
#'   \item{in_parallel}{Suggested number of parallel calls.}
#'   \item{per_call   }{Suggested number of CPUs/threads per call.}
#'   \item{for_call   }{Specific number of CPUs/threads for each call.
#'     Useful when numThreads is not a multiple of numCalls.}
#' }
#' @examples allocate_cpus(60, 10, 5)
allocate_cpus <- function
(
    numThreads,
    numCalls,
    numWorkers = NA
){
    if (is.list(numCalls))
        numCalls <- length(numCalls)
    if (numThreads <= 0 || numCalls <= 0)
        stop("Number of CPUs and files must be greater zero!")
    if (is.na(numWorkers) || is.null(numWorkers) || !is.numeric(numWorkers))
        # bypass negative numbers
        numWorkers <- max(1, min(numCalls, numThreads))
    else
        # bypass negative numbers
        numWorkers <- max(1, min(numWorkers, numCalls, numThreads))

    in_parallel <- min(numWorkers, numThreads)
    per_call    <- numThreads %/% in_parallel # %/% is integer division

    for_call <- rep(per_call, numCalls)
    if (numThreads > numWorkers && numThreads %% in_parallel > 0)
        for (i in 1:(numThreads %% in_parallel))
            for_call[i] <- for_call[i]+1

    return(list(
        in_parallel = in_parallel,
        per_call = per_call,
        for_call = for_call
    ))
}



#' Perform sanity checks on rownames, colnames and content of design matrix
#'
#' @param designMatrix DEG design matrix
#' @param counts Counts matrix as returned by featureCounts. If NA, column
#'   name check is ignored.
.check_design_matrix <- function(designMatrix, counts = NA) {
    if (length(rownames(designMatrix)) == 0)
        stop("'designMatrix' does not have row names!")

    if (is.matrix(counts) || is.data.frame(counts)) {
        if (length(colnames(counts)) == 0)
            stop("'counts' object does not have column names!")

        if (FALSE %in% (rownames(designMatrix) %in% colnames(counts)))
            stop(paste(
                "One or more samples defined in design matrix, but not in counts!",
                paste(rownames(designMatrix)[!(rownames(designMatrix) %in% colnames(counts))], collapse = "\n"),
                collapse = "\n"
            ))
    } else if (!(is.na(counts) || is.null(counts)))
        stop("Given 'counts' object is neither matrix nor data.frame!")


    if (FALSE %in% (designMatrix %in% c("control", "treatment", "none")))
        stop("Invalid content strings in design matrix! Only 'control', 'treatment' or 'none' are allowed as values.
             Check the vignette for more information using ?Geo2RNAseq.")

    if (FALSE %in% grepl("_[vV][sS]_", colnames(designMatrix)))
        warning(paste("One or more comparisons seems to miss the '_vs_' substring. These are:",
                      paste(grep("_[vV][sS]_", colnames(designMatrix), value = T, invert = T), collapse = "\n"),
                      sep = "\n"))
    return(TRUE)
}



# TODO PACKAGE: make R implementation without system call?
# Number Of Reads In FASTQ Files
#
# Calculate the number of reads per FASTQ file.
#
# @param files Character(n). Paths to FASTQ files.
# @inheritParams parallelizationParameters
# @return Integer(n). Number of reads per file.
number_reads_fastq <- function
(
    files,
    cpus = 1,
    workers = 15
){
    driver <- function(f) {
        call <- paste("wc -l", f)
        r <- system(call, intern = TRUE)
        r <- unlist(strsplit(r, " "))[1]
        num <- as.integer(r) %/% 4
        return(num)
    }

    if (is.matrix(files)) files <- asPairVector(files)
    if (FALSE %in% file.exists(files))
        stop("number_reads: one or more input files do not exist!")

    allocCPUS <- allocate_cpus(cpus, length(files), workers)
    res <- parallel::mclapply(files, driver, mc.cores = allocCPUS$in_parallel)
    res <- as.integer(unlist(res))
    names(res) <- files
    return(res)
}



# Count Lines
#
# Count the number of lines in all files.
# @param files List of files to count in.
# @inheritParams parallelizationParameters
# @return List of number of lines per file.
countLines <- function
(
    files,
    cpus=1
){

    driver <- function(f){
        tmp <- system(paste("wc -l", f), intern=TRUE)
        return(as.integer(strsplit(tmp," ")[[1]][[1]]))
    }

    res <- parallel::mclapply(files, driver, mc.cores = min(cpus, length(files)))
    return(res)
}



#' Synchronize Names Based On Sorted Order
#'
#' Given a list of characters 'names', put them in the same order as the
#' characters in 'order'.
#'
#' Assumption: if 'order' and 'names' are sorted, then they are synchronized.
#' This is used to synchronize the order of file names returned by
#' \code{\link{dir}} or \code{\link{list.files}} with the order of file names
#' given by, e.g., a metadata table (in which the file order might be different).
#' File extensions don't need to match.
#' \cr\cr
#' Example:
#' \cr
#' 'order' and 'names' were sorted. 'order' has been changed
#' afterwards. This change is now applied to 'names'. If 'order =
#' c("b", "a", "c")', 'names' are sorted first and then re-ordered by '2,1,3'.
#'
#' @export
#' @param order Atomic vector(n). Usually character.
#' @param names Atomic vector(n). Usually character.
#' @return Character vector. The re-ordered version of 'names'.
#'
#' @examples
#' names <- as.character(1:5)
#' order <- c("b", "a", "c", "e", "d")
#' print(sort(names))
#' print(sort(order))
#' synchronize_sorted(order, names)
synchronize_sorted <- function
(
    order,
    names
){
    if (!is.atomic(order) || !is.atomic(names) || !is.vector(order) || !is.vector(names))
        stop("Synchronize: input variables must be atomic vectors! Try using 'unlist'.")
    if (length(order) != length(names))
        stop("synchronize: 'order' does not have same length as 'names'!")

    order <- reorder_idx(sort(order), order)
    return(sort(names)[order])
}



#' Archive Pipeline Data
#'
#' Zip important result files.
#'
#' @export
#' @param zipFile Character(1). Zip file to be created. Defaults to "results.zip".
#' @param IDF Character(1). IDF table file to be added to zip file. Defaults to
#'   NA.
#' @param SDRF Character(1). SDRF table file to be added to zip file.  Defaults
#'   to NA.
#' @param driver Character(1). "Driver" (pipeline) R-script file to be added to
#'   zip file. Defaults to NA.
#' @param workspace Character(1). RData workspace file to be added to zip file.
#'   Defaults to NA.
#' @param statsDir Character(1). Directory with stat files, e.g. counts, DEGs,
#'   ... . Defaults to "stats".
#' @param includeXLS Logical(1). Include XLS files? Defaults to TRUE.
#' @param includePDF Logical(1). Include PDF files? Defaults to FALSE.
#' @param includeMappingStats Logical(1). Include mapping stats file? Defaults
#'   to FALSE.
#' @param copyToDir Character(1). If not NA, copy resulting zip file to this
#'   directory. Defaults to NA.
#'
#' @return Nothing.
#' @examples
#' \dontrun{
#'    archive(zipFile = "results.zip")
#' }
archive <- function
(
    zipFile = "results.zip",
    IDF = NA,
    SDRF = NA,
    driver = NA,    # "rnaseq_preprocessing_pipeline.R"
    workspace = NA, # "R_workspace.RData"
    statsDir = "stats",
    includeXLS = TRUE,
    includePDF = FALSE,
    includeMappingStats = FALSE,
    copyToDir = NA
) {
    zipFile <- zipFile
    filesToZip <- c()

    if (!is.na(IDF))       filesToZip <- c(filesToZip, IDF)
    if (!is.na(SDRF))      filesToZip <- c(filesToZip, SDRF)
    if (!is.na(driver))    filesToZip <- c(filesToZip, driver)
    if (!is.na(workspace)) filesToZip <- c(filesToZip, workspace)

    statsPattern <- "(\\.csv$)"
    if (includeXLS) statsPattern <- paste(statsPattern, "(\\.xls$)", sep = "|")
    if (includePDF) statsPattern <- paste(statsPattern, "(\\.pdf$)", sep = "|")

    filesToZip <- c(
        filesToZip,
        file.path(
            statsDir,
            list.files(
                statsDir,
                pattern = statsPattern,
                ignore.case = TRUE,
                recursive = TRUE
            )
        )
    )

    # do not zip mapping_stats files
    if (!includeMappingStats) filesToZip <- grep("mapping_stats", filesToZip, invert = TRUE, value = TRUE)

    # TODO test with not existing file name(s)
    sapply(filesToZip, function(file){
        if (!file.exists(file)) stop(paste(file, "does not exist"))
    })

    zip(zipFile, filesToZip)

    if (!is.na(copyToDir)) {
        if (dir.exists(copyToDir)) {
            success <- file.copy(zipFile, copyToDir)
            # file.copy() returns a logical vector
            # in this case with length == 1, because length(zipFile) == 1
            # --> testing for TRUE/FALSE only works as expected for
            # length(zipFile) == length(success) == 1
            if (!success) {
                warning(paste(
                    "Could not copy", shQuote(zipFile), "to", shQuote(copyToDir)
                ))
            }
        } else {
            warning(paste("Directory", shQuote(copyToDir), "does not exist."))
        }
    }
}



#' Conditions From DesignMatrix
#'
#' Given a design matrix with samples in rows and comparisons of two conditions
#' in columns, assign a condition to each sample.
#'
#' Each comparison must be defined as '<A>_VS_<B>', where <A> and <B> are
#' variable conditions. <A> is considered treatment, <B> is considered control
#' for the corresponding comparison.
#' \cr\cr
#' If multiple, different labels (conditions) are assigned to the same sample,
#' all assignments are reported.
#'
#' @export
#' @param designMatrix Character matrix(n x m). Samples in rows, comparisons in
#' columns.
#' @return Character(n). Condition per sample.
#' @examples
#' dm <- matrix("none", nrow = 4, ncol = 2)
#' colnames(dm) <- c("A_VS_B", "A_VS_C")
#' rownames(dm) <- paste("Sample", 1:nrow(dm))
#' dm[1:2,] <- "treatment"
#' dm[3,"A_VS_B"] <- "control"
#' dm[4,"A_VS_C"] <- "control"
#' conditions_from_design(dm)
conditions_from_design <- function
(
    designMatrix
){
    .check_design_matrix(designMatrix)

    conds <- rep("none", nrow(designMatrix))

    for (i in 1:ncol(designMatrix)) {
        s <- colnames(designMatrix)[i]
        # first group is 'treatment', second group is 'control'
        groups <- unlist(strsplit(s, "_[(VS)|(vs)|(Vs)]{2}_"))
        tre <- designMatrix[, i] == "treatment"
        con <- designMatrix[, i] == "control"
        if (TRUE %in% !(conds[tre] %in% c("none", groups[1])))
            warning(
                "A sample is assigned multiple, different conditions. Overwrite.\n",
                " Keeping: ", paste(conds[tre], collapse = ", "), "\n",
                " for: ", paste(which(tre), collapse = ","), "\n",
                " new: ", paste(groups[1])
            )
        if (TRUE %in% !(conds[con] %in% c("none", groups[2])))
            warning(
                "A sample is assigned multiple, different conditions. Overwrite.\n",
                " Keeping: ", paste(conds[con], collapse = ", "), "\n",
                " for: ", paste(which(con), collapse = ","), "\n",
                " new: ", paste(groups[2])
            )
        conds[tre] <- groups[1]
        conds[con] <- groups[2]
    }

  return(conds)
}



#' Detect High Coverage In Counts
#'
#' Report any gene, which holds more than 10\% of all reads of a sample. Such
#' genes might be coding for possibly unwanted rRNA molecules.
#'
#' rRNA is the most abundant type of RNA in cells. It can be depleted using
#' several wet lab techniques. Sometimes, the depletion is not as successful as
#' desired, resulting in samples "contaminated" with high levels of rRNA reads.
#' \cr\cr
#' Many RNA-seq experiments are designed to measure transcript, i.e. mRNA, levels.
#' Therefore, reads arising from rRNA genes should be removed as much as possible
#' for this kind of experiments and before applying statistical measurements.
#'
#' @export
#' @param counts Matrix or data.frame. Genes in rows, samples in columns.
#' @param libSizes Integer(n). Must be the same size as the
#'   number of columns in 'counts'. Defaults to colSums(counts).
#' @return Nothing.
#' @examples
#' mat <- matrix(10, ncol = 2, nrow = 10)
#' colnames(mat) <- c("Cond_A", "Cond_B")
#' rownames(mat) <- paste("Sample", 1:10)
#' mat[8,1] <- 1000
#' detect_high_coverage(mat)
detect_high_coverage <- function
(
    counts,
    libSizes = colSums(counts)
){
    if (length(libSizes) != ncol(counts))
        stop("Length of 'libSizes' must match ncol(counts)!")

    for (i in 1:ncol(counts)) {
        perc = counts[, i] / libSizes[i] * 100
        names(perc) = rownames(counts)
        if (max(perc) > 10) {
            writeLines(paste(
                "Attention! At least one gene in condition",
                shQuote(colnames(counts)[i]),
                "holds more than 10% off all reads.",
                "This might be rRNA contamination."
            ))
#             writeLines(c(
#                 head(sort(perc, decreasing = TRUE)),
#                 head(names(sort(perc, decreasing = TRUE)))
#             ))
            print(head(sort(perc, decreasing = TRUE)))
        }
    }
}


#' Write Count Table
#'
#' Write the content of a matrix or data.frame object to disk in tabulated CSV
#' or XLS format.
#'
#' @export
#' @param file Prefix path to output file. Appropriate file extension are added
#'   automatically.
#' @param counts Matrix or data.frame.
#' @param as.xls If TRUE, the output format is XLS.
#'   Else, tabulated CSV is used.
#' @param sheetNames In the case of XLS output, sheets are named according to
#'   sheetNames. Names are clipped at 31 characters.
#' @param rnames Character. Row names. See \code{\link{write.table}}. Defaults
#'   to rownames(counts).
#' @param cnames Character. Column names. See \code{\link{write.table}}.
#'   Defaults to colnames(counts).
#' @return Nothing.
#' @examples
#' mat <- matrix(10, ncol = 2, nrow = 10)
#' colnames(mat) <- c("Cond_A", "Cond_B")
#' rownames(mat) <- paste0("gene_", 1:10)
#' write_count_table("count", mat)
write_count_table <- function
(
      file,
      counts,
      as.xls = FALSE,
      sheetNames = "",
      rnames = rownames(counts),
      cnames = colnames(counts)
){
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

    if (as.xls) {
        if (!is.logical(rnames)) {
            if (!is.na(rnames) && !is.null(rnames) && length(rnames) > 0) {
                rownames(counts) <- rnames
                rnames = TRUE
            } else
                rnames = FALSE
        }

        if (!is.logical(cnames)) {
            if (!is.na(cnames) && !is.null(cnames) && length(cnames) > 0) {
                colnames(counts) <- cnames
                cnames = TRUE
            } else
              cnames = FALSE
        }

        # slow for huge tables, e.g. human genome and many conditions
        WriteXLS::WriteXLS(
            as.data.frame(counts),
            ExcelFileName = paste0(file, ".xls"),
            SheetNames = substr(sheetNames, 1, 31), # truncate to at most 31 chars
            AdjWidth = TRUE,
            BoldHeaderRow = TRUE,
            FreezeRow = 1,
            FreezeCol = 1,
            row.names = rnames,
            col.names = cnames
        )
    } else
        write.table(
          counts,
          file = paste0(file, ".csv"),
          sep = "\t",
          row.names = rnames,
          col.names = cnames
        )
}



# Stop At Error
#
# Only for bplapply calls. Terminate if a 'remote-error' occurred and show the
# corresponding error message.
#
# @rdname stop_at_error
# @param res BiocParallel::bplapply return value.
# @param from Name of the function or any identifier that indicates the origin
#'   of res (or the error source).
# @param warn If TRUE, a warning is returned instead of an error.
.stop_at_error = function
(
    res,
    from = NA,
    warn = FALSE
){
    if ("bplist_error" %in% class(res))
        res <- attr(res, "result")

    error_inds = sapply(
        res,
        function(x)
            TRUE %in% (c("remote-error", "remote_error") %in% class(x))
    )

    if (sum(error_inds) > 0) {
        if (!warn) {
            stop(
                paste(res[[which(error_inds == TRUE)]], collapse = "\n"),
                "\n",
                paste("From: ", from)
            )
        } else
            warning(
                paste(res[[which(error_inds == TRUE)]], collapse = "\n"),
                "\n",
                paste("From: ", from)
            )
    }
}


# Resolve Link
#
# Normalize a relative path from an absolute path.
# e.g. "/home/user/test1" & "../test2" -> "/home/user/test2"
#
# @rdname resolve_link
# @param abs_path Absolute path to start from.
# @param rel_path Relative path to resolve.
# @return Character(1).
.resolve_link <- function
(
    abs_path,
    rel_path
){
    relink_path <- normalizePath(file.path(abs_path, rel_path), mustWork = FALSE)
    if (file.exists(rel_path)) {
        return(rel_path)
    } else if (file.exists(relink_path)) {
        return(relink_path)
    } else
        stop(paste("Found a link, but cannot resolve it! ", rel_path))
}


#' Get Executable
#'
#' Find the path to a given executable. If the executable is not in the PATH
#' system variable, try with workarounds based on other information.
#'
#' @param exec Character(1). The name of the executable to find.
#' @param derive Character(1). Name of other executable to use to get exec.
#' @param sysVar System variable name.
#' @param only.dir If TRUE and sysVar is NOT NA, only the directory name of
#'   sysVar is used and appended by 'exec'.
#' @param alt Character(1). Alternative binary name.
#' @return Character(1). Path to the executable.
.get_executable <- function
(
    exec,
    derive = "",
    sysVar = NA,
    only.dir = FALSE,
    alt = ""
){
    # Check if an OS variable has been set. If it did, it might point to a wrong
    # item. Make sure to notify the user!
    if (!is.na(sysVar) && !is.null(sysVar) && Sys.getenv(sysVar) != "") {
        if (file.exists(Sys.getenv(sysVar))) {
            if (only.dir)
                sys_exec <- file.path(dirname(Sys.getenv(sysVar)), exec)
            else
                sys_exec <- Sys.getenv(sysVar)
        } else {
            warning(
                "System variable ",
                sysVar,
                " defined, but file doesn't exist: ",
                shQuote(Sys.getenv(sysVar))
            )
            sys_exec <- Sys.which(exec)[1]
        }
    } else
        # set from PATH binaries
        sys_exec <- Sys.which(exec)[1]

    if (sys_exec == "") {
        if (alt != "")
            return(.get_executable(alt))

        if (derive == "")
            if (alt == "")
                stop(paste("Cannot find:", exec))


        # try getting exec from 'derive' path
        dir_exec <- Sys.which(derive)
        if (dir_exec == "")
            stop(paste("Cannot find ", exec, "or", derive, "!"))

        # may be a link
        link <- Sys.readlink(dir_exec)
        if (link != "")
            dir_exec <- .resolve_link(dirname(dir_exec), link)

        sys_exec <- file.path(dirname(dir_exec), exec)
        if (!file.exists(sys_exec))
            stop(paste("Found", derive, "but cannot find", exec, "!"))
    } else {
        link <- Sys.readlink(sys_exec)
        if (link != "")
            sys_exec <- .resolve_link(dirname(sys_exec), link)
    }

    return(sys_exec)
}



#' Get Program Version
#'
#' Return a string describing the version of a given (external) program/tool.
#'
#' @export
#' @param tool Character(1). Can be either "fastq-dump", "multiqc", "trimmomatic"
#'   , "samtools", "sortmerna", "bowtie2", "tophat2", "hisat2" or "fastqc".
#' @return Character(1).
#' @examples get_program_version("fastqc")
get_program_version <- function
(
    tool
){
    tools <- c("fastq-dump", "multiqc", "trimmomatic", "samtools", "sortmerna",
        "bowtie2", "tophat2", "hisat2", "fastqc")
    sVars <- c("FASTQ_DUMP_EXEC", "MULTIQC_EXEC", "TRIMMOMATIC_EXEC",
        "SAMTOOLS_EXEC", "SORTMERNA_EXEC", "BOWTIE2_EXEC", "TOPHAT2_PATH",
        "HISAT2_EXEC", "FASTQC_EXEC")
    aVars <- c("", "", "TrimmomaticSE", "", "", "", "", "", "")

    idx <- which(tool == tools)
    if (length(idx) <= 0)
        stop(paste(
            "Unknown tool:",
            tool,
            ". Choose one of these: ",
            paste(tools, collapse = " "))
        )

    # Get executable first. If it doesn't exist on the system, the functin will
    # throw an error.
    exec <- tryCatch({
        .get_executable(tools[idx], sysVar = sVars[idx], alt = aVars[idx])
    }, error = function(e) {
        return("Binary not found.")
    })

    if (!file.exists(exec))
        # NOTE: if this return message is changed, also change it in "list_executables"!
        return("Binary not found.")

    # get version string. Acquisition is different for every tool...
    ver <- ""
    if (tool == "fastq-dump") {
        ver <- grep("\\d+$", system(paste(exec, "--version"), intern = TRUE), value = TRUE)
        ver <- paste("fastq-dump :", gdata::last(unlist(strsplit(ver, " "))))
    } else if (tool ==  "multiqc") {
        ver <- system(paste(exec, "--version"), intern=TRUE)
        ver <- gdata::last(unlist(strsplit(ver, split = " ")))
    } else if (tool ==  "trimmomatic") {
        if (grepl("\\.jar$", exec))
            ver <- system(paste("java -jar", exec, "-version"), intern = TRUE)
        else if (grepl("Trimmomatic[PS]E$", exec))
            ver <- "Indeterminable"
        else
            ver <- system(paste(exec, "-version"), intern = TRUE)
    } else if (tool ==  "samtools") {
        suppressWarnings(ver <- grep("Version:", system2(exec, stdout = TRUE, stderr = TRUE), value = TRUE))
        ver <- sub("Version: ", "", ver)
    } else if (tool ==  "sortmerna")
        ver <- system(paste(exec, "--version 2>&1"), intern=TRUE)[[2]]
    else if (tool ==  "bowtie2") {
        ver <- grep("bowtie2-align", system(paste(exec, "--version"), intern = TRUE), value = TRUE)
      ver <- gdata::last(unlist(strsplit(ver, split = " ")))
    } else if (tool ==  "tophat2")
        ver <- system(paste(exec, "--version"), intern = TRUE)
    else if (tool ==  "hisat2") {
        ver <- grep("hisat.*version", system(paste(exec, "--version"), intern = TRUE), value = TRUE)
        ver <- gdata::last(unlist(strsplit(ver, split = " ")))
    } else if (tool ==  "fastqc") {
        ver <- system(paste(exec, "--version"), intern=TRUE)
        ver <- gdata::last(unlist(strsplit(ver, split = " ")))
        ver <- gsub("[a-zA-Z]", "", ver)
    }

    return(ver)
}



#' List Executables
#'
#' List all tools found from within the R session. Most binaries may be found
#' within the directories defined by the PATH system variable. Also show other
#' system variables that can be defined for each tool. Ill defined system
#' variables, e.g. those which point to a non-existing file, are reported with a
#' warning.
#'
#' @export
#' @return Nothing.
#' @examples list_executables()
list_executables <- function
(){
    tools <- c("fastq-dump", "multiqc", "trimmomatic", "samtools", "sortmerna",
        "bowtie2", "tophat2", "hisat2", "fastqc")
    sVars <- c("FASTQ_DUMP_EXEC", "MULTIQC_EXEC", "TRIMMOMATIC_EXEC",
        "SAMTOOLS_EXEC", "SORTMERNA_EXEC", "BOWTIE2_EXEC", "TOPHAT2_PATH",
        "HISAT2_EXEC", "FASTQC_EXEC")
    aVars <- c("", "", "TrimmomaticSE", "", "", "", "", "", "")

    writeLines("\n=== List of installed tools ===")
    for (i in 1:length(tools)) {
        writeLines(paste("-->", tools[i]))
        writeLines(paste("    sys var:", sVars[i], "=", Sys.getenv(sVars[i])))
        # HACK: the function already catches any obvious errors
        ver <- get_program_version(tools[i])
        if (ver == "Binary not found.") {
            writeLines("    *Binary not found on your system.*")
        } else {
            writeLines(c(
                paste("    exec:", shQuote(.get_executable(
                    tools[i],
                    sysVar = sVars[i],
                    alt = aVars[i]
                ))),
                paste("    version:", ver)
            ))
        }
    }
}



#' As Pair Vector
#'
#' Return an ordered vector of paired files given by a 2 column matrix of paired
#' files. Ordered means that the two files of a pair, as given by the matrix,
#' are always "neighbors" in the vector.
#'
#' Consider mat[i,1] and mat[i,2] to be a pair of, e.g. paired-end FASTQ files.
#' 'asPairVector' returns a vector of length 2i, where
#' \cr\cr
#' vector[j] = mat[i,1]
#' \cr
#' and
#' \cr
#' vector[j+1] = mat[i,2]
#' \cr
#' for all i with j = 2i
#' \cr\cr
#' 'asPairVector' is the inverse function of 'asPaired':
#' v == asPairVector(asPaired(v))
#'
#' @export
#' @param  mat Matrix or vector. A vector is returned unchanged.
#' @return Vector. Type is the same as the given matrix.
#' @seealso Inverse function \code{\link{asPaired}}.
#' @examples
#' mat <- matrix("0", ncol = 2, nrow = 5)
#' mat[,1] <- paste0("a", 1:5)
#' mat[,2] <- paste0("b", 1:5)
#' asPairVector(mat)
#'
#' # inverse function
#' asPaired(asPairVector(mat))
asPairVector <- function
(
    mat
){
    if (is.vector(mat) || is.null(mat) || is.na(mat) || nrow(mat) < 1)
        return(mat)
    if (!is.matrix(mat))
        stop(paste(
            "Invalid argument for mat: expected a matrix, but got:",
            typeof(mat)
        ))
    if (ncol(mat) != 2)
        stop(paste(
            "Invalid number of rows for 'asPairVector'. Expected 2, but got: ",
            ncol(mat)
        ))

    # note: optimized for 2 column
    even     <- 1:length(mat) %% 2 == 1
    v        <- vector(typeof(mat), length(mat))
    v[even]  <- mat[,1]
    v[!even] <- mat[,2]
    return(v)
}


#' As Paired Matrix
#'
#' Given an even numbered list of files, return a matrix m, so that m[,1] are
#' the odd and m[,2] are the even files.
#' \cr
#' For paired-end FASTQ files, the first pair would then be m[1,1] and m[1,2].
#'
#' @export
#' @param v A list or vector. For FASTQ files, it should be ordered, so
#'   that:
#'   \enumerate{
#'     \item files with uneven index are forward strand
#'     \item files with even index are reverse strand
#'     \item i and i+1 are a pair, if i is even
#'   }
#' @return A n row, 2 column matrix. Each row index is a pair.
#' @seealso Inverse function \code{\link{asPairVector}}.
#' @examples
#' v <- vector("character", 10)
#' v[1:10 %% 2 == 0] <- paste0("f", 1:5, "_1.fastq")
#' v[1:10 %% 2 == 1] <- paste0("f", 1:5, "_2.fastq")
#' print(v)
#' asPaired(v)
#'
#' # inversed
#' asPairVector(asPaired(v))
asPaired <- function
(
    v
){
    pairedFiles = matrix(NA, nrow=length(v)/2, ncol=2)
    for ( i in seq(2, length(v), 2) )
        pairedFiles[i/2,] = c(v[[i-1]], v[[i]])

    return(pairedFiles)
}


# Find In Upward Directory
#
# Given a directory, try to find the file given by a relative path 'what'.
#' Unless 'what' is found, follows each higher directory until root. If 'what'
#' wasn't found, an error is thrown.
#
# @rdname find_upward
# @param dir Relative or absolute directory to start searching from. If relative,
#   the current working directory is used to normalize the path.
# @param what Relative path to file of interest.
# @return Character(1).
.find_upward <- function
(
    dir,
    what
){
    dir <- normalizePath(dir)
    if (file.exists(file.path(dir, what)))
        return(normalizePath(file.path(dir, what)))
    else {
        p <- tryCatch({
            parts <- unlist(strsplit(dir, split = .Platform$file.sep))
            for (i in (length(parts)-1):1) {
              path <- paste(parts[1:i], collapse = .Platform$file.sep)
              if (file.exists(file.path(path, what)))
                return(normalizePath(file.path(path, what)))
            }
        }, error = function(e) {
            stop(paste("Could not find", what, "in directory tree of", dir, "."))
        })

        return(p)
    }
}



# Is Synchronized Pair
#
# Given the paths of two paired FASTQ files, check if every line i of file1
# corresponds to line i in file2. If this fails for any index, the whole tests
# fails, i.e. both files are desynchronized. Mostly used for testing the
# synchronized files.
#
# @keywords internal
# @rdname is_synchronized_pair
# @param file1 Path to forward read file.
# @param file2 Path to reverse read file.
# @return  A boolean(1). TRUE, if all reads are synchronized, FALSE otherwise.
# @seealso \code{\link{synchronize_fastq}}
.is_synchronized_pair <- function
(
    file1,
    file2
){
    if (FALSE %in% file.exists(c(file1, file2)))
        stop(
            "One or more input files do not exist:\n ",
            paste(c(file1, file2), collapse = "\n ")
       )

    f1_str <- ShortRead::FastqStreamer(file1, n = 1e5)
    f2_str <- ShortRead::FastqStreamer(file2, n = 1e5)
    on.exit(close(f1_str), add = TRUE)
    on.exit(close(f2_str), add = TRUE)

    # filter a BStringSet 'id' and convert it to character vector
    id_to_str <- function(id) {
        splits <- strsplit(as.character(id), split = " ", fixed = TRUE)
        return(vapply(splits, function(s) s[[1]], "character"))
    }
    f1_fq <- ShortRead::yield(f1_str)
    f1_ids <- id_to_str(f1_fq@id)
    f2_fq <- ShortRead::yield(f2_str)
    f2_ids <- id_to_str(f2_fq@id)

    i <- 1
    repeat{
        # terminal condition - no new data fetched in one of the three input files
        if (0 %in% c(length(f1_fq), length(f2_fq)))
            break

        message("Processing: ", i * 1e5)
        i <- i + 1

        if (length(f1_ids) != length(f2_ids) || FALSE %in% (f1_ids == f2_ids))
            return(FALSE)
        f1_fq <- ShortRead::yield(f1_str)
        f1_ids <- id_to_str(f1_fq@id)
        f2_fq <- ShortRead::yield(f2_str)
        f2_ids <- id_to_str(f2_fq@id)
    }

    return(TRUE)
}



# NOTE: robustness check
# 1. empty FASTQ files ~> TEST
# NOTE: measurements
# - bulk processing 1.000 reads   ~ 5.10 seconds for uncompressed out
# - bulk processing 10.000 reads  ~ 1.65 seconds for uncompressed out
# - bulk processing 100.000 reads ~ 1.30 seconds for uncompressed out
# - bulk processing 100.000 reads ~ 4.90 seconds for gz compressed out
# - writing single read: ~ 0.004 seconds

#' Synchronize Paired-End FASTQ Files
#'
#' Remove all unpaired, single-end reads from two (possibly) desynchronized
#' FASTQ files. One of the original, synchronized files and both
#' (de-) synchronized files are required.
#'
#' After processing paired-end read data, the reads in the resulting paired
#' FASTQ files may not be synchronized anymore, i.e. single-end reads were
#' introduced.
#' \cr\cr
#' Reading and writing of FASTQ files is done by using ShortRead.
#' \cr\cr
#' This is a very memory efficient paired-end FASTQ synchronizer with a memory
#' footprint of 100 MB to 300 MB and O(3*|reads|) runtime.
#'
#' @export
#' @param orig Path to original first or second mate FASTQ file.
#' @param file1 Path to first mate FASTQ file.
#' @param file2 Path to second mate FASTQ file.
#' @param outDir Path to output directory. Defaults to the working directory.
#' @param chunkSize The number of reads loaded in memory, per file, at the same
#'    time. Defaults to 1e+05.
#' @param full If TRUE, the identifier line is repeated on the third
#'   line of the FASTQ record. Defaults to FALSE.
#' @param compress If TRUE, output files are gz compressed. Defaults to FALSE.
#' @param checkEqual For testing purposes. If TRUE, output IDs are
#'   checked for equality twice. Defaults to FALSE.
#' @param overwrite Logical indicating whether existing output files
#'   should be overwritten. Defaults to FALSE.
#' @return A list keywords:\cr
#' \describe{
#'   \item{files  }{Vector(2). Paths to output files.}
#'   \item{numSync}{Numeric(1). Number of reads kept.}
#'   \item{loss1  }{Numeric(1). Number of reads removed from file 1.}
#'   \item{loss2  }{Numeric(1). Number of reads removed from file 2.}
#' }
#' @references Morgan M, Anders S, Lawrence M, Aboyoun P, Pages H and Gentleman
#'   R (2009). "ShortRead: a Bioconductor package for input, quality assessment
#'   and exploration of high-throughput sequence data." Bioinformatics, 25, pp.
#'   2607-2608. doi: 10.1093/bioinformatics/btp450, http://dx.doi.org10.1093/bioinformatics/btp450.
#' @examples
#' \dontrun{
#' orig <- "sample1_1.fastq"
#' # after some processing
#' file1 <- "sample1_1.trimo.fastq"
#' file2 <- "sample1_2.trimo.fastq"
#' synchronize_fastq(orig, file1, file2)
#' }
#' synchronize_fastq
synchronize_fastq <- function
(
      orig,
      file1,
      file2,
      outDir = getwd(),
      chunkSize = 1e5,
      full = FALSE,
      compress = FALSE,
      checkEqual = FALSE,
      overwrite = FALSE
){
    if (FALSE %in% file.exists(c(orig, file1, file2)))
        stop(
            "One or more input files do not exist:\n ",
            paste(c(orig, file1, file2), collapse = "\n ")
        )
    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
    fout1 <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(file1))
    fout1 <- file.path(outDir, paste0(fout1, ".sync.fastq"))
    fout2 <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(file2))
    fout2 <- file.path(outDir, paste0(fout2, ".sync.fastq"))
    if (compress) {
        fout1 <- paste0(fout1, ".gz")
        fout2 <- paste0(fout2, ".gz")
    }
    if (file.exists(fout1) || file.exists(fout2)) {
        if (overwrite)
            file.remove(fout1, fout2)
        else
            stop("
                One or both output files already existing and not overwrite:\n ",
                fout1,
                "\n ",
                fout2
            )
    }
    message("Synchronizing paired-end files:\n ", file1, "\n ", file2, "\nBased on:\n ", orig)

    # TODO: form orig file, only the header is needed ~ find a work around to NOT load the whole file!
    # TODO: test speed, RAM requirements, etc.
    or_str <- ShortRead::FastqStreamer(orig , n = chunkSize, verbose = FALSE)
    f1_str <- ShortRead::FastqStreamer(file1, n = chunkSize)
    f2_str <- ShortRead::FastqStreamer(file2, n = chunkSize)
    on.exit(close(or_str), add = TRUE)
    on.exit(close(f1_str), add = TRUE)
    on.exit(close(f2_str), add = TRUE)

    # filter a BStringSet 'id' and convert it to character vector
    id_to_str <- function(id) {
        splits <- strsplit(as.character(id), split = " ", fixed = TRUE)
        return(vapply(splits, function(s) s[[1]], "character"))
    }
    or_fq <- ShortRead::yield(or_str)
    or_ids <- id_to_str(or_fq@id)
    f1_fq <- ShortRead::yield(f1_str)
    f1_ids <- id_to_str(f1_fq@id)
    f2_fq <- ShortRead::yield(f2_str)
    f2_ids <- id_to_str(f2_fq@id)

    # variables for statistics - number of kept and removed reads
    num_sync <- 0
    num_lost_1 <- 0
    num_lost_2 <- 0

    # Process each chunk of reads until one of the three input files is empty
    i <- 1
    repeat{
        # terminal condition - no new data fetched in one of the three input files
        if (0 %in% c(length(or_fq), length(f1_fq), length(f2_fq)))
            break

        message("Processing: ", i * chunkSize)
        i <- i + 1
        in_f1 <- or_ids %in% f1_ids
        in_f2 <- or_ids %in% f2_ids
        # keep only reads that are present in at least one file besides orig
        in_sync <- (in_f1 & in_f2)
        num_in_sync <- sum(in_sync)
        # adapt 'sync' to fit the vector of reads actually present in F1 and F2
        sync_f1 <- in_sync[in_f1]
        sync_f2 <- in_sync[in_f2]

        # only for debugging and validity checks
        if (checkEqual) {
              # JUST FOR TESTING
              i1 <- id_to_str(f1_fq[sync_f1]@id)
              i2 <- id_to_str(f2_fq[sync_f2]@id)
              stopifnot(!(FALSE %in% (i1 == i2)))
        }

        # NOTE: if num_in_sync is zero, nothing is written out. Safe that time...
        if (length(in_sync) > 0 && num_in_sync > 0) {
            ShortRead::writeFastq(f1_fq[sync_f1], fout1, "a", compress = compress, full = full)
            ShortRead::writeFastq(f2_fq[sync_f2], fout2, "a", compress = compress, full = full)
        }

        # just some statistics
        valid   <- in_f1 | in_f2
        num_sync <- num_sync + num_in_sync
        num_lost_1 <- num_lost_1 + sum(!in_f1[valid])
        num_lost_2 <- num_lost_2 + sum(!in_f2[valid])
        or_fq <- ShortRead::yield(or_str)
        or_ids <- id_to_str(or_fq@id)

        # if not all reads from f1 were used, keep them on top of the next list
        if (length(in_sync) < length(f1_fq)) {
            f1_fq <- f1_fq[length(in_sync)+1 : length(f1_fq)]
            f1_fq <- ShortRead::append(f1_fq, ShortRead::yield(f1_str))
        } else {
            f1_fq <- ShortRead::yield(f1_str)
        }
        f1_ids <- id_to_str(f1_fq@id)

        # if not all reads from f2 were used, keep them on top of the next list
        if (length(in_sync) < length(f2_fq)) {
            f2_fq <- f2_fq[length(in_sync)+1 : length(f2_fq)]
            f2_fq <- ShortRead::append(f2_fq, ShortRead::yield(f2_str))
        } else {
            f2_fq <- ShortRead::yield(f2_str)
        }
        f2_ids <- id_to_str(f2_fq@id)
    }
    num_lost_1 <- num_lost_1 + length(f1_fq)
    num_lost_2 <- num_lost_2 + length(f2_fq)

    # iterate remaining reads... can only be in f1 or f2
    while (length(f1_fq <- ShortRead::yield(f1_str))) {
        num_lost_1 <- num_lost_1 + length(f1_fq)
    }
    while (length(f2_fq <- ShortRead::yield(f2_str))) {
        num_lost_2 <- num_lost_2 + length(f2_fq)
    }
    message(
        "Synchronized: ", num_sync,
        "\nRemoved F1: ", num_lost_1,
        "\nRemoved F2: ", num_lost_2
    )

    return(list(
        files = c(fout1, fout2),
        numSync = num_sync,
        loss1 = num_lost_1,
        loss2 = num_lost_2)
    )
}



# alternative unzip function which automatically detects the compression format
# and WITHOUT parallelization (might give no speed improvement anyway)
# --> sometimes (e.g. data from GATC) it is a mixture of gz and bz2
# Unzip Files
#
# Wrapper function. Tries to detect compression format and decompresses
# appropriately.
#
# @rdname unzip
# @param files Character vector or list. Paths to compressed files.
# @param dir Character. Path to output directory. By default, same directory as
#   input files.
# @return Character(n). Relative paths of output files.
.unzip = function
(
    files,
    dir = ""
) {
    out_files = list()
    file_names = basename(files)
    if (!file.exists(dir)) dir.create(dir, recursive = TRUE)

    for (i in 1:length(files)) {
        out = if (dir == "") basename(files[i]) else dir
        if (grepl(".tar.gz$", files[i])) {
            system(paste("tar -zxvf", files[i], "-C", dir))
            out_files[i] = file.path(out, substr(file_names[i], 1, nchar(file_names[i])-7))
        }else if (grepl(".gz$", files[i])) {
            system(paste("pigz", "-d", files[i]))
            out_files[i] = file.path(out, substr(file_names[i], 1, nchar(file_names[i]) - 3))
            file.rename(from = files[i], to = out_files[i]) # move to output directory
        } else if (grepl(".bz2$", files[i])) {
            system(paste("bunzip2",  files[i]))
            out_files[i] = substr(file_names[i], 1, nchar(file_names[i]) - 4)
            file.rename(from = files[i], to = out_files[i])
        } else {
            stop(paste(
                "Unsupported file format for decompression detected:", files[i]
            ))
        }
    }

    return(unlist(out_files))
}



#' Calculate Reorder Column Indices
#'
#' Calculate an index vector 'inds' of indices for 'in_names' so that:\cr
#' in_names[inds] == order_names
#' \cr\cr
#' 'in_names' must be a subset of 'order_names'.
#'
#' @export
#' @param in_names Character vector(n).
#' @param order_names Character vector(n).
#' @return Integer(n). Indices to achieve reordering.
#' @examples
#' unOrder <- c("a", "b", "c", "d")
#' newOrder <- c("b", "c", "a", "d")
#' inds <- reorder_idx(unOrder, newOrder)
#' print(inds)
#' print(unOrder[inds])
reorder_idx <- function
(
    in_names,
    order_names
){
    res = rep(NA,length(order_names))
    for(i in 1:length(order_names))
        res[i] = which(in_names==order_names[i])

    return(res)
}


# NOTE: in TPM paper, average read length is used. However, it can be shown that it is unnecessary
#' Calculate TPM
#'
#' Transcripts per million (TPM) is a read count normalization. TPM fixes
#' certain issues related to the reads per kilobase (RPK) normalization.
#'
#' @export
#' @param countings Matrix of counted reads. Genes in rows, samples in columns.
#' @param genelengths List of gene lengths.
#'   length(genelength) == nrow(countings) must be True.
#' @param libsizes List of library sizes.
#'   length(libsizes) == ncol(countings) must be True.
#' @return Matrix with TPM normalized values.
#' @references Wagner, Guenter P., Koryu Kin, and Vincent J. Lynch. "Measurement
#'   of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among
#'   samples." Theory in biosciences 131.4 (2012): 281-285.
#' @seealso \code{\link{get_rpkm}}
#' @examples
#' genLen <- c(1:10)
#' input <- cbind(c(1:10), 3*c(1:10), 10*c(1:10))
#' libSizes <- colSums(input)
#' get_tpm(input, genLen, libSizes)
get_tpm <- function
(
    countings,
    genelengths,
    libsizes
){
    RPK <- (countings / genelengths)
    # NOTE: min: if colSums is zero, we get NA. This way, we force 0
    return(t(t(RPK) / pmax(colSums(RPK), 1)) * 10^6)
}


#' Calculate RPKM
#'
#' Reads per kilobase per million mapped reads (RPKM) is a read count
#' normalization.
#'
#' @export
#' @inheritParams get_tpm
#' @return Matrix with RPKM normalized values.
#' @references Wagner, Guenter P., Koryu Kin, and Vincent J. Lynch. "Measurement
#'   of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among
#'   samples." Theory in biosciences 131.4 (2012): 281-285.
#' @seealso \code{\link{get_tpm}}
#' @examples
#' genLen <- c(1:10)
#' input <- cbind(c(1:10), 3*c(1:10), 10*c(1:10))
#' libSizes <- colSums(input)
#' get_rpkm(input, genLen, libSizes)
get_rpkm <- function
(
  countings,
  genelengths,
  libsizes
){
  return((t(t(countings)/libsizes)*10^9)/(genelengths))
}


#' Helper function for get_mrn
.mrn_size_factors <- function
(
    counts
){
    if (ncol(counts) == 1)
        stop("Cannot comute MRN size factors: only 1 sample!")

    geomeans <- exp(rowMeans(log(counts)))
    sizefactors <- apply(counts, 2, function(c) median((c / geomeans)[geomeans > 0]))
    return(sizefactors)
}


#' Calculate MRN
#'
#' Median by-Ratio Normalization (MRN) is a read count normalization as used
#' in EBSeq. DESeq2 and edgeR use slightly modified versions of MRN (e.g. TMM).
#'
#' @export
#' @inheritParams get_mrn
#' @return Matrix with MRN normalized values.
#' @references Simon Anders and Wolfgang Huber. Differential expression
#'   analysis for sequence count data. Genome Biology (2010) 11:R106
#' @seealso \code{\link{get_tpm}}
#' @details Counts are normalized based on the geometric mean, so that at most
#'   on halve of genes is up or down-regulated. It has been shown to be more
#'   robust against library size effects. Unlike TMM, all data is considered
#'   (and not just a quantile). MRN shows decent performance when used for log
#'   fold change filtering after DEG analysis when using different tools (see
#'   volcano plots. The data is centered very close to 0).
#' @examples
#' input <- cbind(c(1:10), 3*c(1:10), 10*c(1:10))
#' get_mrn(input, genLen, libSizes)
get_mrn <- function
(
    countings
){
    sizefactors <- .mrn_size_factors(countings)
    stopifnot(length(sizefactors) == ncol(countings))
    mrn = t(t(countings) / sizefactors)
    return(mrn)
}



# Detect Tools
#
# Helper function. Tries to read the tools used for DEG analysis, i.e. the tools
# present in the resulting data.frame.
#
# @rdname detect_tools
# @param degTable Single list element (data.frame) as returned by
#   \code{\link{calculate_DEGs}}.
# @return Character(n).
.detect_tools <- function
(
    degTable
){
    cols <- colnames(degTable)
    candidates <- sapply(cols, function(c) typeof(degTable[[c]]) == "logical")
    tools <- cols[candidates]
    message("Found the following tools:\n", paste(tools, collapse = "\n"))
    hasEntries <- sapply(tools, function(x) TRUE %in% degTable[,x])
    if (length(tools) > sum(hasEntries))
        message("Skip tools with 0 hits:\n", paste(tools[!hasEntries], collapse = "\n"))
    tools <- tools[hasEntries]
    if (length(tools) == 0)
        warning("Could not detect any tools!")
    return(tools)
}


#' Make DEG-List
#'
#' Take a data.frame (any list element returned by \code{\link{calculate_DEGs}})
#' and return a named list of differentially expressed genes (DEGs) reported per
#' tool.
#'
#' This is a helper function for the visualization of overlapping (gene) sets.
#'
#' @export
#' @include DEG_Analysis.R
#' @param degTable Single list element (data.frame) as returned by
#'   \code{\link{calculate_DEGs}}.
#' @param tools Character(n). Names of statistical tools for the calculation of
#'   DEGs. If NA, the method tries to guess the tools based on the given
#'   data.frame. In general, only boolean columns are searched. If not NA, only
#'   the given tools are used to create the list object.
#' @return A list of DEGs per tool.
#' @examples
#' \dontrun{
#' deg_res <- calculate_DEGs("with/valid/input")
#' make_DEG_list(deg_res$DEGs[[1]])
#' }
#' make_DEG_list
make_DEG_list <- function
(
    degTable,
    tools = NA
){
   if (!is.vector(tools) || is.na(tools) || is.null(tools)) {
        tools <- .detect_tools(degTable)
   } else {
        # fix any issue regarding capital letters
        inds <- tolower(colnames(degTable)) %in% tolower(tools)
        tools <- colnames(degTable)[inds]
   }

   # For each tool, look at their top hits (boolean vector).
   # And then grab the corresponding gene names.
   DEGlist <- lapply(tools, function(i) as.character(degTable$id[degTable[, i]]))
   names(DEGlist) <- tools

   return(DEGlist)
}


# Get DEG-List
#
# Wrapper function to safe code. DEGlist (or any named list) is returned
# unchanged. Otherwise, the DEGlist is created using degTable and tools.
#
# @rdname get_DEGlist
# @inheritParams make_DEG_list
# @param DEGlist Named Character list. If supplied, it is returned unchanged.
# @return A list of DEG genes per tool.
.get_DEGlist <- function
(
  degTable = NA,
  tools = NA,
  DEGlist = NA
){
  if (!is.na(DEGlist) && !is.null(DEGlist) && !is.null(names(DEGlist)))
    return(DEGlist)
  else
    return(make_DEG_list(degTable, tools))
}



####-
##-====================================================
##-=========================  CONVERTERS ==========================
##-====================================================
####-


# merge counts for A/B genes of diploid genomes
# only C. albicans assembly A22

# Merge Count Results For Diploid Genomes
#
# Diploid genomes have 2 slightly different versions of the same chromosome.
# They usually have 2 versions of each gene as well. This function takes a
# counting table having at most 2 alleles of the same gene and sums these
# counts of these alleles to a single gene.\cr
# The mapping file describes which alleles should be mapped to the same gene.\cr
# This function was tested to map the diploid A22 genome of Candida albicans
# to the haploid A21 version.
#
# @keywords internal
# @param counts Matrix(n x m).
# @param mapping asdf
# @param allele_regex asdf
# @param gene_length Integer(n).
# @return A list with keywords:\cr
# \describe{
#   \item{counts    }{New counting table for haploid gene set. Counts for
#                     alleles of the same gene are summed.}
#   \item{genelength}{The lengths of the diploid genes used for merge.}
# }
.merge_diploid <- function
(counts, mapping, allele_regex, gene_length) {
    shortGeneID <- unique(sub(allele_regex, "", rownames(counts)))

    res = data.frame(matrix(
        data = NA,
        nrow = length(shortGeneID),
        ncol = ncol(counts)
    ))
    colnames(res) = colnames(counts)

    IDs = c()
    # TODO super slow for loop
    for (j in 1:length(shortGeneID)) {
        longGeneID_idx = which(grepl(shortGeneID[j], rownames(counts), fixed = TRUE) == TRUE)
        select = counts[longGeneID_idx, ]
        # HACK: for just one number, we do not get a matrix, but a vector
        if (length(longGeneID_idx) == 1) {
            IDs = c(IDs, rownames(counts)[longGeneID_idx])
            res[j,] = select
        } else {
            IDs = c(IDs, rownames(select)[1]) # only keep first gene ID (= "_A")

            if (nrow(select) > 1) {
                res[j,] = colSums(select)
            } else {
                res[j,] = select
            }
        }
    }

    counts = res
    rownames(counts) = IDs

    # build hash table for much faster id lookup
    mapping_hash = new.env()
    for (i in 1:nrow(mapping)) {
        mapping_hash[[mapping[i, 1]]] = mapping[i, 2]
    }

    rownames = c()
    for (i in rownames(counts)) {
        rownames = c(rownames, mapping_hash[[i]])
    }
    rownames(counts) = rownames

    # TODO is it correct to only and always use the length of the first gene ID ("_A")?
    gene_length = gene_length[IDs]
    names = c()
    for (i in names(gene_length)) {
        names = c(names, mapping_hash[[i]])
    }
    names(gene_length) = names

    return(list(counts = counts, gene_length = gene_length))
}



#' Convert BAM Files To SAM Files
#'
#' Convert a list of files in BAM format to SAM format using SAMtools.
#'
#' @export
#' @param files List of files in BAM format.
#' @param samDir Output directory for SAM files. If NA, the input file directory
#'   is used. Defaults to NA.
#' @param overwrite Logical. If TRUE, files with the ".sam" suffix are
#'   overwritten. Defaults to FALSE.
#' @inheritParams parallelizationParameters
#' @return Character(n). The paths to the created SAM files.
#' @references Li, Jun, and Robert Tibshirani. "Finding consistent patterns: a
#'   nonparametric approach for identifying differential expression in RNA-Seq
#'   data." Statistical methods in medical research 22.5 (2013): 519-536.
#' @seealso \code{\link{sam_to_bam}}
#' @examples
#' \dontrun{
#' bam_to_sam(c("f1.bam", "f2.bam"))
#' }
#' bam_to_sam
bam_to_sam <- function
(
    files,
    samDir = NA,
    cpus = 1,
    workers = 2,
    overwrite = FALSE
){
    driver <- function(i) {
        if (is.na(samDir))
            samDir = dirname(files[i])
        samName = paste0(tools::file_path_sans_ext(basename(files[i])), ".sam")
        outputSam = file.path(samDir, samName)

        if (file.exists(outputSam))
            message("\nDestination exists: ", outputSam, "\n")
        else
            # asSam add the '.sam' suffix itself
            Rsamtools::asSam(
                file = files[i],
                destination = sub("\\.sam$", "", outputSam)
            )

        return(outputSam)
    }

    allocCPUS <- allocate_cpus(cpus, length(files), workers)
    bio_par <- BiocParallel::MulticoreParam(
        workers = allocCPUS$in_parallel,
        tasks = length(files),
        progressbar = TRUE
    )
    allRes <- BiocParallel::bptry({
        BiocParallel::bplapply(1:length(files), driver, BPPARAM = bio_par)
    })
    .stop_at_error(allRes, from = match.call()[[1]])
    return(unlist(allRes))
}



#' Convert SAM To BAM
#'
#' Convert a list of files in SAM format to BAM format using SAMtools.
#'
#' @export
#' @param files List of files in SAM format.
#' @param bamDir Output directory for BAM files. If NA, input file directory is
#'   used. Defaults to NA.
#' @param sort Logical indicating whether to sort entries of output BAM files.
#' @param inMemory Only has an effect if 'sort' is TRUE. If SAMtools is
#'   available, the SAM files are converted to BAM files and sorted before
#'   written to disk. This is very effective on most systems, but requires a
#'   UNIX operating system. If SAMtools is not available, Rsamtools is used
#'   instead. In this case, the BAM files are first written to disk and
#'   sorted afterwards. Defaults to TRUE.
#' @inheritParams bam_to_sam
#' @return Character(n). The paths to the created BAM files.
#' @references Li, Jun, and Robert Tibshirani. "Finding consistent patterns: a
#'   nonparametric approach for identifying differential expression in RNA-Seq
#'   data." Statistical methods in medical research 22.5 (2013): 519-536.
#' @seealso \code{\link{sam_to_bam}}
#' @examples
#' \dontrun{
#' sam_to_bam(c("f1.sam", "f2.sam"))
#' }
#' sam_to_bam
sam_to_bam <- function
(
    files,
    bamDir = NA,
    sort = TRUE,
    inMemory = TRUE,
    cpus = 1,
    workers = 2,
    overwrite = FALSE
){
    # This driver is used to optimize sorting. Can only be used if SAMtools is
    # installed on the system.
    pipeDriver <- function(i) {
        if (is.na(bamDir))
            bamDir = dirname(files[i])
        bamName = paste0(tools::file_path_sans_ext(basename(files[i])), ".sorted.bam")
        outSortBam = file.path(bamDir, bamName)

        if (!overwrite && file.exists(outSortBam)) {
            message("Destination file already exists: ", shQuote(outSortBam), ". Skipped!")
        } else {
            writeLines(paste("Runner", i, "started ..."))
            call = paste(
                samtools_exec,
                "view -bS",
                files[i],
                "|",
                samtools_exec,
                "sort -",
                sub("\\.bam$", "", outSortBam)
            )
            system(call)
            writeLines(paste("Runner", i, "finished."))
        }

        return(outSortBam)
    }


    # Rsamtools driver.
    rsamDriver <- function(i) {
        if (is.na(bamDir))
            bamDir = dirname(files[i])
        bamName = paste0(tools::file_path_sans_ext(basename(files[i])), ".bam")
        sortBamName = paste0(tools::file_path_sans_ext(basename(files[i])), ".sorted.bam")
        outBam = file.path(bamDir, bamName)
        outSortBam = file.path(bamDir, sortBamName)

        if (sort && !overwrite && file.exists(outSortBam)) {
            message("\nBAM already exists: ", shQuote(outSortBam), ". Skipped!")
            return(outSortBam)
        }

        # convert to BAM first
        if (!overwrite && file.exists(outBam)) {
            message("BAM already exists: ", shQuote(outBam), ". Skipped!")
        } else {
            writeLines(paste("\nRunner:", i, "started..."))
            Rsamtools::asBam(files[i], sub("\\.bam$", "", outBam), indexDestination = FALSE)
        }

        # Sort, if demanded
        if (sort){
            Rsamtools::sortBam(outBam, sub("\\.bam$", "", outSortBam))
            file.remove(outBam)
            writeLines(paste("\nRunner:", i, "finished."))
            return(outSortBam)
        } else {
            writeLines(paste("\nRunner:", i, "finished."))
            return(outBam)
        }
    }

    # prepare parallelization
    allocCPUS <- allocate_cpus(cpus, length(files), workers)
    bio_par <- BiocParallel::MulticoreParam(
        workers = allocCPUS$in_parallel,
        tasks = length(files),
        progressbar = TRUE
    )

    driver <- rsamDriver # function pointer!
    if (sort && inMemory) {
        # in-memory sorting. This mode is inherently UNIX only.
        samtools_exec <- tryCatch({
            .get_executable("samtools", sysVar = "SAMTOOLS_EXEC")
        }, error = function(e) {
            return(NA)
        })

        if (is.na(samtools_exec)) {
            message("Could not find SAMtools on your system. In-memory sorting disabled.")
        } else {
            BiocParallel::bpprogressbar(bio_par) <- FALSE
            driver <- pipeDriver # function pointer!
        }
    }

    allRes <- BiocParallel::bptry({
        BiocParallel::bplapply(1:length(files), driver, BPPARAM = bio_par)
    })
    .stop_at_error(allRes, from = match.call()[[1]])
    return(unlist(allRes))
}


# Is Indexed
#
# Check a list of file paths for existing BAM files in the same directory.
.is_indexed <- function(files) {
    if (!(is.vector(files) || is.list(files)))
        stop("Wrong input argument type. Vector or list expected!")
    inds <- paste0(files, ".bai")
    return(!(FALSE %in% file.exists(inds)))
}



#' Index BAM Files
#'
#' Index a list of BAM files using Rsamtools. Index files are placed in the same
#' path as their input files.
#'
#' @export
#' @param files Path to input BAM files.
#' @inheritParams bam_to_sam
#' @return Character. Paths to new index files.
#' @examples
#' \dontrun{
#' indexBAMs(c("f1.bam", "f2.bam"))
#' }
#' indexBAMs
indexBAMs <- function
(
    files,
    overwrite = FALSE,
    cpus = 1,
    workers = 10
){
    driver <- function(bam) {
        outBam <- paste0(bam, ".bai")
        Rsamtools::indexBam(bam)
        return(outBam)
    }

    allocCPUS <- allocate_cpus(cpus, length(files), workers)
    bio_par <- BiocParallel::MulticoreParam(
        workers = allocCPUS$in_parallel,
        tasks = length(files),
        progressbar = TRUE
    )
    allRes = tryCatch({
        BiocParallel::bplapply(files, driver, BPPARAM = bio_par)
    }, error = identity)
    .stop_at_error(allRes, from = match.call()[[1]])
    return(unlist(allRes))
}



#' Sort BAM Files
#'
#' Sort and index a list of BAM files using Rsamtools.
#'
#' @export
#' @param files Path to input BAM files.
#' @param overwrite If TRUE, existing output files are overwritten.
#'   If FALSE, existing files are skipped. Defaults to FALSE.
#' @param delete.orig If TRUE, original BAM files are deleted afterwards.
#'   Defaults to FALSE.
#' @return Character. Paths to new and sorted/indexed BAM output files.
#' @examples
#' \dontrun{
#' sortBAMs(c("f1.bam", "f2.bam"))
#' }
#' sortBAMs
sortBAMs <- function
(
    files,
    overwrite = FALSE,
    delete.orig = FALSE
){
    res = list()
    for (i in 1:length(files)) {
        # NOTE: .bam is added by sortBAM
        outBam <- paste0(tools::file_path_sans_ext(files[[i]]), ".sorted")
        if (!overwrite && file.exists(paste0(outBam, ".bam"))) {
            message("Output file already exists: ", shQuote(outBam), ". Skipped!")
        } else {
            writeLines(paste("Sorting: ", shQuote(files[[i]])))
            Rsamtools::sortBam(files[[i]], outBam)
            writeLines("Indexing ...")
            Rsamtools::indexBam(paste0(outBam, ".bam"))
        }
        res[[i]] <- paste0(outBam, ".bam")
    }

    if (delete.orig)
        file.remove(files)

    return(unlist(res))
}



#' Filter Mapping Quality
#'
#' Remove alignments with a quality score lower than a given 'quality' value.
#'
#' @export
#' @param files Path to input BAM files.
#' @param overwrite If TRUE, existing output files are overwritten. If FALSE,
#'   existing files are skipped. Defaults to FALSE.
#' @param paired If TRUE, only pairs are retained. Defaults to FALSE.
#' @param quality Minimum quality of an alignment to be retained.
#'   Defaults to 10.
#' @param delete.orig If TRUE, original BAM files are deleted afterwards.
#'   Defaults to FALSE.
#' @return Character(n). Paths to new and quality filtered BAM output files.
#' @examples
#' \dontrun{
#' sortBAMs(c("f1.bam", "f2.bam"))
#' }
#' filterMappingQuality
filterMappingQuality <- function
(
    files,
    overwrite = FALSE,
    paired = FALSE,
    quality = 10,
    delete.orig = FALSE
){
    res = list()
    for (i in 1:length(files)) {
        # NOTE: .bam is added by sortBAM
        outBam <- paste0(tools::file_path_sans_ext(files[[i]]), ".mq", quality, ".bam")
        if (!overwrite && file.exists(outBam)) {
            message("Output file already exists: ", shQuote(outBam), ". Skipped!")
        } else {
            # requires sorted, indexed BAM files
            flag <- Rsamtools::scanBamFlag(isPaired = paired, isUnmappedQuery = FALSE)
            what <- Rsamtools::ScanBamParam(flag = flag, mapqFilter = 10)
            filterObj <- Rsamtools::filterBam(files[i], destination = outBam, param = what)
        }
        res[[i]] <- outBam
    }

    if (delete.orig)
        file.remove(files)

    return(unlist(res))
}
