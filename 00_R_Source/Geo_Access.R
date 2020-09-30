###-
##-====================================================
##-=========================  GEO ACCESS ==========================
##-====================================================
###-



#' Parse DATA.FILE Column
#'
#' Meta data objects usually have a column for input files (e.g. FASTQ)
#' and one for the library strategy (e.g. SINGLE-END or PAIRED-END).
#' Split the Data.File column into single and paired end files.
#' Paired end information should be in the 'Library.Strategy' column.
#'
#' @rdname parse_DATA_FILE
#' @param mat data.frame from (parsed) SDRF or meta table object.
#' @param lib_strat_col Default to "Library.Strategy". Custom column name for
#'   library strategy. Indicate single-end reads by the keyword 'SINGLE' and
#'   paired-end reads by the keyword 'PAIRED'.
#' @param data_file_col Defaults to "Data.File". Custom column name for file paths.
#' @return A list with keyword arguments:\cr
#' \describe{
#'   \item{seFiles  }{Character. List of singe end FASTQ files.}
#'   \item{peFiles  }{Character. List of paired end FASTQ files.}
#'   \item{seOrder  }{re-order indices. If seFiles were sorted, seFiles[order]
#'                    will list them based on the SDRF. Important, when files
#'                    are read from drive using dir() or list.files()}
#'   \item{peOrder  }{re-order indices. If peFiles were sorted, peFiles[order]
#'                    will list them based on the SDRF. Important, when files
#'                    are read from drive using dir() or list.files()}
#'   \item{seIndices}{Row indices of SE files. Only important for mixed mode.}
#'   \item{peIndices}{Row indices of PE files. Only important for mixed mode.}
#'   \item{read_mode}{Character. Sample library strategy. Can be 'single' or
#'                   'paired' if all samples have same strategy. 'mixed'
#'                    otherwise.}
#' }
.parse_DATA_FILE <- function
(
    mat,
    lib_strat_col = "Library.Strategy",
    data_file_col = "Data.File"
){
    # initiate return values
    se_files  = NULL
    pe_files  = NULL
    se_order  = NULL
    pe_order  = NULL
    se_index  = NULL
    pe_index  = NULL
    read_mode = NULL

    if (is.null(mat[[lib_strat_col]])) {
        # try other two columns
        if (!is.null(mat[["Library.Strategy"]])) {
            warning("Wrong library strategy column [", lib_strat_col, "] supplied. Using 'Library.Strategy' instead.")
            lib_strat_col <- "Library.Strategy"
        } else if (!is.null(mat[["Layout"]])) {
            warning("Outdates metadata format detected. Using 'Layout' for instead of [", lib_strat_col, "].")
            lib_strat_col <- "Layout"
        } else
            stop(paste("Could not find library strategy column:", shQuote(lib_strat_col)))
    }

    # detect single, paired or mixed end
    if (TRUE %in% (grepl("paired", mat[[lib_strat_col]], ignore.case = TRUE))) {
        paired_idx <- grep("paired", mat[[lib_strat_col]], ignore.case = TRUE)
        pe_files <- unlist(strsplit(as.character(mat[paired_idx,][[data_file_col]]), ";\\s?"))
        # files must be unique. Also, missing entries are detected this way
        if (length(unique(pe_files)) != length(pe_files))
            stop("Data.File column must contain unique file paths! Did you forget to fill this in?")
        pe_order <- reorder_idx(
            sort(mat[paired_idx,][[data_file_col]]),
            mat[paired_idx,][[data_file_col]]
        )
        pe_index <- paired_idx
        read_mode <- "paired"
    }

    if (TRUE %in% (grepl("single", mat[[lib_strat_col]], ignore.case = TRUE))) {
        single_idx <- grep("single", mat[[lib_strat_col]], ignore.case = TRUE)
        se_files <- as.character(mat[single_idx,][[data_file_col]])
        if (length(unique(se_files)) != length(se_files))
            stop("Data.File column must contain unique file paths! Did you forget to fill this in?")
        se_order <- reorder_idx(
            sort(mat[single_idx,][[data_file_col]]),
            mat[single_idx,][[data_file_col]]
        )
        se_index <- single_idx

        if (is.null(read_mode)) {
            read_mode <- "single"
        } else {
            read_mode <- "mixed"
            message("Single- and Paired-end FASTQ files present. ",
                    "Use mixed pipeline strategy!")
        }
    }

    if (is.null(se_files) && is.null(pe_files))
        warning("Could not find any single-end or paired-end files. Did you forget to fill this in column 'Data.File'?")
    return(list(
        seFiles = se_files,
        peFiles = pe_files,
        seOrder = se_order,
        peOrder = pe_order,
        seIndices = se_index,
        peIndices = pe_index,
        readMode = read_mode
    ))
}



#' Convert SRA To FASTQ
#'
#' Convert SRA to FASTQ files using the NCBI SRA toolkit.
#'
#' This is basically a wrapper function for the tool 'fastq-dump' from the NCBI
#' SRA toolkit.
#'
#' @export
#' @param sraFiles Character(n). Path(s) to SRA file(s).
#' @param delete.sra Logical(1) indicating whether to remove SRA files after
#'   conversion. Defaults to FALSE.
#' @return A list with keyword arguments:\cr
#' \describe{
#'   \item{fastq  }{Paths to output FASTQ files.}
#'   \item{version}{Version string of fastq-dump.}
#' }
#' @examples
#' \dontrun{
#'   files <- c("f1.sra", "f2.sra")
#'   sra_to_fastq(files)
#' }
#' sra_to_fastq
sra_to_fastq <- function
(
    sraFiles,
    delete.sra  = FALSE
){
    # filter invalid files
    inds <- grep(".sra$", sraFiles)
    if (length(inds) != length(sraFiles)) {
        invalid <- grep(".sra$", sraFiles, invert = TRUE, value = TRUE)
        warning(paste("Found invalid files. Ignoring:\n", paste(invalid, collapse = "\n ")), immediate. = TRUE)
        sraFiles <- sraFiles[inds]
    }

    if (length(sraFiles) == 0) {
        writeLines("Nothing to convert.")
        return("")
    }

    # find sra-toolkit "fastq-dump" executable
    sra_exec <- .get_executable("fastq-dump", "FASTQ_DUMP_EXEC")
    sra_version <- grep(
        "\\d+$",
        system(paste(sra_exec, "--version"), intern = TRUE),
        value = TRUE
    )
    sra_version <- gdata::last(unlist(strsplit(sra_version, " ")))

    fastq_files <- paste0(tools::file_path_sans_ext(sraFiles), ".fastq")

    # NOTE: fastq-dump can be called with multiple files, but doing so leads to
    #       system errors.
    logs <- vector(mode = "character", length(sraFiles))
    for (i in 1:length(sraFiles)) {
        if (!file.exists(fastq_files[[i]])) {
            writeLines(paste("Converting: ", sraFiles[[i]]))
            call <- paste(
                sra_exec,
                "--split-3",
                sraFiles[[i]]
            )
            log <- system(call, intern = TRUE)
            logs[i] <- log[[1]]
        } else {
            cat(paste(
                "Skipping:",
                fastq_files[[i]],
                "Already exists.\n",
                sep = "\n"
            ))
        }
    }

    write(logs, "sra_log.txt")
    return(list(fastq = fastq_files, version = sra_version))
}



#' Get GEO Data From GSE
#'
#' Access the Gene Expression Omnibus (GEO) database and acquire raw sequencing
#' data as well as metadata information for a given GSE accession number.
#'
#' Will download SRA files automatically. For large data sets, it is advised to
#' have both, a strong internet connection and a stable R session (e.g. R
#' started within a screen session to resume it if the connection was lost).
#' \cr\cr
#' The SRA toolkit (available at the NCBI website) is required to convert files
#' in SRA format to FASTQ format. If installed, the conversion is handled by
#' this function.
#' \cr\cr
#' If gzip is available, output FASTQ files are compressed to .fastq.gz files.
#'
#' @export
#' @param accession Character. GSE accession number from any RNA-seq experiment
#'   available at GEO.
#' @param outDir Path to directory for non-temporary files (SRA, FASTQ, logs).
#'   Defaults to current working directory.
#' @param compress Logical(1). If TRUE, FASTQ files will be gz compressed. Only
#'   if 'sra_only' is FALSE.
#' @param sra_only Logical(1). If TRUE, downloaded SRA files will not be converted.
#'   Otherwise, 'fastq-dump' (if available) from the SRA toolkit is used to
#'   convert SRA files to FASTQ files.
#' @param rm_temp Logical(1). If TRUE, intermediate files are removed. This
#'   includes SRA files if conversion to FASTQ is enabled.
#' @return A list with the following keyword arguments:
#' \describe{
#'   \item{SDRFfile  }{Path to output metadata file in SDRF format.}
#'   \item{IDFfile   }{Path to output metadata file in IDF format.}
#'   \item{METAfile  }{Path to output metadata file in tabulated CSV format.}
#'   \item{SDRF      }{Metadata data.frame from SDRF file.}
#'   \item{META      }{Metadata data.frame from tabulated CSV file.}
#'   \item{sra_files }{Paths to SRA files. NA if 'rm_temp = TRUE'.}
#'   \item{se_files  }{Paths to converted single-end FASTQ files.}
#'   \item{pe_files  }{Paths to converted paired-end FATQ files.}
#'   \item{se_index  }{Row-indices of single-end FASTQ files. This should be
#'                     used to update metadata data.frames if the user works
#'                     with both, single-and paired end files at the same time.}
#'   \item{pe_index  }{Row-indices of paired-end FASTQ files. This should be
#'                     used to update metadata data.frames if the user works
#'                     with both, single-and paired end files at the same time.}
#' }
#' @examples
#' \dontrun{
#'   res <- getGEOdata("GSE41749")
#' }
#' Geo2RNAseq::getGEOdata
getGEOdata <- function
(
    accession,
    outDir = getwd(),
    compress = TRUE,
    sra_only = FALSE,
    rm_temp = FALSE
){
    # Test if external tools are installed or give advice how to fix/bypass them
    if (sra_only && get_program_version("fastq-dump") == "Binary not found.")
        stop("In order to use 'getGEOdata', please install 'fastq-dump' from the SRA toolkit or set 'sra_only=TRUE'!")

    if (compress) {
        tryCatch({
            .get_executable("gzip")
        }, error = function(e){
            stop("'gzip' was not found on your system. This is required to compress FASTQ files. Install it or set 'compress = FALSE'.")
        })
    }

    # prepare output directories
    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
    outDir <- normalizePath(outDir)
    if (!grepl("/$", outDir))
        outDir <- paste0(outDir, "/")
    gseDir <- file.path(outDir, accession, "")
    dir.create(gseDir, recursive = TRUE, showWarnings = FALSE)

    if (length(list.files(outDir, pattern = "_SDRF.tsv$")) > 0) {
        message("Metadata SDRF file already present. Using existing table. Delete it and call this function again if this was not intended!")
    } else {
        print("Downloading metadata information")
        scripts <- system.file("extdata/BASH", package = "Geo2RNAseq")
        scr1 <- file.path(scripts, "getMetaGSE.sh")
        call <- paste(scr1, accession, outDir)
        # get metadata
        ret  <- tryCatch({
            system(call, intern = TRUE)
        }, error = function(e) {
            stop("Failed to download metadata information. Make sure that the given ",
                 "GSE accession is correct!")
        })

        # convert to SDRF and IDF tables
        scr2 <- file.path(scripts, "IDFSDRF.sh")
        call <- paste(scr2, gseDir, outDir)
        system(call)
    }

    # grab files
    print(paste0("Output tables were downloaded to: ", shQuote(outDir)))
    SDRFfile <- list.files(outDir, pattern = "_SDRF.tsv$", full.names = TRUE)
    IDFfile <- list.files(outDir, pattern = "_IDF.tsv$", full.names = TRUE)
    METAfile <- list.files(outDir, pattern = paste0("_", accession, ".tsv$"), full.names = TRUE)

    # check if it is a re-written meta file
    skip <- if (grepl("^Source", readLines(SDRFfile, n=1))) 1 else 0
    META <- read.csv(METAfile, sep = "\t", header = TRUE, as.is = TRUE, skip = skip)
    SDRF <- read.csv(SDRFfile, sep = "\t", header = TRUE, as.is = TRUE, skip = skip)

    # here comes the hard part: download FASTQ but pay attention to mixed-data sets
    # paired-end files will be downloaded and put into >1< row
    # single-end files will be downloaded and put into >1< row
    # if both types are present, the table must be combined from both

    layout <- rep("none", nrow(META))
    layout[grep("PAIRED", META$Library.Layout, ignore.case = TRUE)] <- "PAIRED"
    layout[grep("SINGLE", META$Library.Layout, ignore.case = TRUE)] <- "SINGLE"
    s_index <- which(layout == "SINGLE")
    p_index <- which(layout == "PAIRED")
    SRR <- META$SRR.Acc
    SRR_s <- META$SRR.Acc[s_index]
    SRR_p <- META$SRR.Acc[p_index]
    dest_files <- file.path(outDir, paste(
        META$SRR.Acc,
        META$GSM,
        META$GSE,
        sep = "_"
    ))
    sra_files <- paste0(dest_files, ".sra")

    if (length(SRR) > 0) {
        # process single-end FASTQ files
        for (i in 1:length(SRR)) {
            srr <- SRR[[i]]
            sra_file <- sra_files[[i]]
            # the first six characters
            f6 <- substr(srr, 1, 6)
            sra_url <- paste0(
                "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/",
                f6,
                "/",
                srr,
                "/",
                srr,
                ".sra"
            )
            # check, if file already exists OR if download was aborted
            header <- RCurl::getURL(sra_url, nobody=1L, header=1L)
            header <- strsplit(header, "\r\n")[[1]][1]
            ftp_size <- as.numeric(strsplit(header, " ")[[1]][2])

            if (!file.exists(sra_file) || file.info(sra_file)$size != ftp_size) {
                print(paste0("Downloading SRA to: ", shQuote(sra_file)))
                download.file(url=sra_url, destfile = sra_file)
            } else {
                print(paste("SRA file", shQuote(sra_file), "already present."))
            }
        }
    }

    # convert to FASTQ format
    fq_se_files <- NA
    fq_pe_files <- NA
    if (!sra_only) {
        print("Converting to FASTQ format...")
        bin <- .get_executable("fastq-dump")

        # single-end SRA files produce just one FASTQ file
        if (length(s_index) > 0) {
            fq_se_files <- paste0(dest_files[s_index], ".fastq")
            if (compress) fq_se_files <- paste0(fq_se_files, ".gz")
            sra_se_files <- paste0(dest_files[s_index], ".sra")
            for (i in 1:length(sra_se_files)) {
                if (file.exists(fq_se_files[i])) {
                    print(paste(
                        "FASTQ file",
                        shQuote(fq_se_files[i]),
                        "already present."
                    ))
                } else {
                    print(paste0("Converting: ", shQuote(sra_se_files[i])))
                    if (compress)
                        call <- paste(
                            bin,
                             "--split-3",
                             "--gzip",
                             "-O",
                             outDir, sra_se_files[i]
                         )
                    else
                        call <- paste(
                            bin,
                            "--split-3",
                            "-O",
                            outDir, sra_se_files[i]
                        )
                    system(call)
                }
            }

            # update both tables
            META$Data.File[s_index] <- fq_se_files
            SDRF$Data.File[s_index] <- fq_se_files
        }

        # paired-end SRA files produce two FASTQ files
        if (length(p_index) > 0) {
            fq_pe_files <- rep(dest_files[p_index], each = 2)
            forw <- seq(1, length(fq_pe_files), by = 2)
            back <- seq(2, length(fq_pe_files), by = 2)
            fq_pe_files[forw] <- paste0(fq_pe_files[forw], "_1.fastq")
            fq_pe_files[back] <- paste0(fq_pe_files[back], "_2.fastq")
            if (compress) fq_pe_files <- paste0(fq_pe_files, ".gz")

            sra_pe_files <- paste0(dest_files[p_index], ".sra")
            for (i in 1:length(sra_pe_files)) {
                if (file.exists(fq_pe_files[forw][i]) && file.exists(fq_pe_files[back][i])) {
                    print(paste(
                        "FASTQ files",
                        shQuote(fq_pe_files[forw][i]),
                        shQuote(fq_pe_files[back][i]),
                        "already present."
                    ))
                } else {
                    print(paste0("Converting: ", shQuote(sra_pe_files[i])))
                    if (compress)
                        call <- paste(
                            bin,
                            "--split-3",
                            "--gzip",
                            "-O",
                            outDir, sra_pe_files[i]
                        )
                    else
                        call <- paste(
                            bin,
                            "--split-3",
                            "-O",
                            outDir,
                            sra_pe_files[i]
                        )
                    system(call)
                }
            }

            # update both tables
            comb_files <- paste(fq_pe_files[forw], fq_pe_files[back], sep = ";")
            META$Data.File[p_index] <- comb_files
            SDRF$Data.File[p_index] <- comb_files
        }
    } else {
        # add SRA to tables
        META$SRA.Files <- sra_files
        SDRF$SRA.Files <- sra_files
    }

    # write out updates meta-tables
    write.table(META, METAfile, quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(SDRF, SDRFfile, quote = FALSE, sep = "\t", row.names = FALSE)

    if (rm_temp && !sra_only) {
        file.remove(sra_files)
        sra_files <- rep(NA, length(sra_files))
    }

    return(list(
        SDRFfile = SDRFfile,
        IDFfile = IDFfile,
        METAfile = METAfile,
        SDRF = SDRF,
        META = META,
        sra_files = sra_files,
        se_files = fq_se_files,
        pe_files = fq_pe_files,
        se_index = s_index,
        pe_index = p_index
    ))
}



#' Write Meta Table
#'
#' Write a metadata data.frame to disk.
#'
#' @export
#' @param table A data.frame. Usually the output of \code{\link{make_meta_table}}.
#' @param as.xls Logical indicating whether to write the meta table as XLS file.
#'   If FALSE, tab separated CSV format is used instead. Defaults to TRUE.
#' @param outDir Path to output directory. Defaults to working directory.
#' @param overwrite Logical indicating whether to overwrite existing files.
#' @return Path to output file.
#' @seealso \code{\link{make_meta_table}}
#' @examples
#' mat <- matrix("fill", nrow = 2, ncol = 3)
#' df <- data.frame(mat)
#' colnames(df) <- c("files", "conditions", "Layout")
#' write_meta_table(df, overwrite = TRUE)
write_meta_table <- function
(
    table,
    as.xls = TRUE,
    outDir = getwd(),
    overwrite = FALSE
){
    out_xls <- file.path(outDir, "meta_data.xls")
    out_csv <- file.path(outDir, "meta_data.csv")
    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

    if (file.exists(out_xls) || file.exists(out_csv)) {
        if (overwrite)
            message("Metadata files already exist. Overwrite enabled...")
        else
            stop("Metadata files already exist. Overwrite is FALSE.")
    }

    if (as.xls) {
        WriteXLS::WriteXLS(
            table,
            ExcelFileName = out_xls,
            SheetNames = "Main",
            AdjWidth = TRUE,
            BoldHeaderRow = TRUE,
            FreezeRow = 1,
            FreezeCol = 1,
            row.names = FALSE
        )
    } else {
        write.table(table, file = out_csv, sep = "\t")
    }

    return(if (as.xls) out_xls else out_csv)
}


#' Make and Write Meta Table
#'
#' Create a simple metadata table based on the output of \code{\link{getGEOdata}}.
#'
#' The returned metadata table is not sufficient to run the entiry GEO2RNAseq
#' pipeline. More information must be provided.
#'
#' @export
#' @param metas NA or 'metas' list object as returned by \code{\link{getGEOdata}}.
#'   If set to NA or NULL, a generic, empty metadata table with 'nrRows' rows
#'   will be created. Defaults to NA.
#' @param files A 'read_files' list object as returned by \code{\link{getGEOdata}}.
#'   Alternativly, a character vector of paths. Conditions and and other sample
#'   information must correspond to the correct files! Defaults to the empty
#'   string "".
#' @param as.xls Logical indicating whether to write the metadata table as XLS
#'   file. If FALSE, tabulated CSV format is used instead. Defaults to TRUE.
#' @param detect.pairs Logical. If TRUE, the method will group paired-end FASTQ
#'   files together, i.e. both files will share the same row. Will only work
#'   with the default naming scheme, i.e. "*_1.fastq" for forward reads and
#'   "*_2.fastq" for reverse reads. Also works for '.gz' and '.bzip' file
#'   extensions.
#' @param nrRows Integer(1) > 0. If 'files = ""' and 'metas = NA', an empty
#'   metadata table with this many rows will be created. Defaults to 1.
#' @param overwrite Logical indicating whether to overwrite existing files.
#' @param outDir Output directory. Defaults to working directory.
#' @return A list with keyword arguments:\cr
#' \describe{
#'   \item{file}{Output file of the metadata.}
#'   \item{meta}{The metadata as data.frame.}
#' }
#' @examples
#' make_meta_table(nrRows = 1)
#'
#' \dontrun{
#'    geo <- getGEOdata("GSE41749")
#'    make_meta_table(metas = geo$metas, files = geo$sra_files, overwrite = TRUE)
#' }
make_meta_table <- function
(
    metas = NA,
    files = "",
    as.xls = TRUE,
    detect.pairs = TRUE,
    nrRows = 1,
    overwrite = FALSE,
    outDir = getwd()
){
    cols <- c("Specimen", "Organism", "Strain", "Date", "Library.Strategy", "Dual", "Data.File")
    row <- c("FILL", "FILL", "NA", "NA", "SINGLE", "FALSE", "FILL")

    if (is.null(metas) || is.na(metas)) {
        # create an empty table
        writeLines("Create empty metadata table")
        if (files != "" && !is.null(files) && !is.na(files)) {
            mat <- matrix(rep(row, each = length(files)), length(files), ncol = length(cols))
            colnames(mat) <- cols
            rownames(mat) <- as.character(1:length(files))
            mat[,"Data.File"] <- files
        } else if (nrRows >= 1) {
            mat <- matrix(rep(row, each = nrRows), nrow = nrRows)
            colnames(mat) <- cols
            rownames(mat) <- as.character(1:nrRows)
            mat[,"Data.File"] <- paste0("FILL", 1:nrRows)
        } else
            stop("Cannot make metadata table with less than 1 rows!")

        mat <- as.data.frame(mat)

    } else if (is.list(metas)) {
        # create table from meta information table
        if (detect.pairs && files != "" && !is.na(files) && !is.null(files)) {
            # detect paired files
            endings <- c(".fastq", ".fastq.gz", ".fastq.bzip")
            pattern <- paste(endings, collapse = "|")
            endless <- sub(pattern, "", files)
            for_idx <- grep("_1$", endless)
            rev_idx <- grep("_2$", endless)
            # get pairs in a robust way
            forward <- sub("_1$", "", endless[for_idx])
            reverse <- sub("_2$", "", endless[rev_idx])

            # helps to find reverse mate. If not found, mark it for removal
            find_mate <- function(x) {
                idx <- grep(x, reverse)
                return(if (length(idx) == 0) 0 else idx)
            }
            pair_map <- sapply(forward, find_mate)

            # remove single mates ~> outlier detection
            if (0 %in% pair_map) {
                single <- pair_map %in% 0
                warning(paste(
                    "One or more mates are single!\n",
                    paste(forward[single], collapse = "\n")
                ))
                # exclude single mates from synchronization
                forward <- forward[!single]
                #pair_map <- pair_map[!single]
            }

            # reorder reverse to fit forward order
            rev_idx <- rev_idx[pair_map]
            files[for_idx] <- paste0(files[for_idx], "; ", files[rev_idx])
            # remove single reverse entries
            files <- files[-rev_idx]
        }

        # the 'metas' object is a list of GSM entries. This results in a matrix.
        mat <- as.data.frame(sapply(metas, function(gsm) unlist(gsm)))
        mat <- base::cbind(mat, Data.File = files)
        # reorder columns to a certain pattern (see 'cols' for order!)
        other <- colnames(mat)[which(!is.element(colnames(mat), cols))]
        reorder <- c(cols, other)
        order_idx <- reorder_idx(colnames(mat), reorder)
        mat <- mat[order_idx]
    } else
        stop("Wrong object for 'metas'. Expected a list of GSM entries.")

    out_file <- write_meta_table(
        mat,
        as.xls = as.xls,
        outDir = outDir,
        overwrite = overwrite
    )

    return(list(file = out_file, table = mat))
}


#' Parse Metadata CSV Files
#'
#' Parse information contained in CSV metadata files.
#'
#' @export
#' @param file Path to meta file (CSV format).
#' @param skip Number of lines to skip in CSV file. Defaults to 1.
#' @param lib_strat_col Column name for library strategy. Defaults to "Library.Strategy".
#' @return
#' List with keyword arguments: \cr
#' \describe{
#'   \item{table       }{Metadata as data.frame.}
#'   \item{nrSamples   }{Integer. Number of samples.}
#'   \item{samples     }{Character. Parsed sample names.}
#'   \item{seFiles     }{Character. List of single-end FASTQ files.}
#'   \item{peFiles     }{Character. List of paired-end FASTQ files.}
#'   \item{seOrder     }{Integer(n). If seFiles were sorted, seFiles[seOrder]
#'                      will return the files in the same order as defined in
#'                      the metadata table. Important if file names were
#'                      acquired using \code{\link{dir}} or
#'                      \code{\link{list.files}}.}
#'   \item{peOrder     }{Integer(n). If peFiles were sorted, peFiles[peOrder]
#'                      will return the files in the same order as defined in
#'                      the metadata table. Important if file names were
#'                      acquired using \code{\link{dir}} or
#'                      \code{\link{list.files}}.}
#'   \item{seIndices   }{Row indices of SE files. Only important for mixed mode.}
#'   \item{peIndices   }{Row indices of PE files. Only important for mixed mode.}
#'   \item{read_mode   }{Character. Sample library strategy. Can be 'single' or
#'                       'paired' if all samples have the same strategy. 'mixed'
#'                       otherwise.}
#'   \item{conditions  }{Character. Parsed conditions.}
#'   \item{designMatrix}{Design Matrix as used for differentially expressed
#'                       gene (DEG) analysis. Each column defines a separate
#'                       test.}
#' }
#' @examples
#' res <- make_meta_table(nrRows = 1, as.xls = FALSE, overwrite = TRUE)
#' parse_meta_csv(res$file, skip = 0)
parse_meta_csv <- function(
    file,
    skip = 1,
    lib_strat_col = "Library.Strategy"
){
    meta <- read.csv(
        file,
        sep = "\t",
        header = TRUE,
        as.is = TRUE,
        row.names = 1,
        skip = skip
    )
    colnames(meta) <- gsub("\\s+", ".", colnames(meta))
    number_samples <- nrow(meta)
    samples <- meta[["Specimen"]]
    conditions <- meta[["Specimen"]]

    # split file information of Data.File into SE and PE file information
    pdf <- .parse_DATA_FILE(as.data.frame(meta), lib_strat_col = lib_strat_col)

    # TODO PACKAGE: generalize that more?
    # design_matrix: which samples to compare and test for differential expression
    design_matrix <- as.matrix(meta[grep("_vs_", colnames(meta), ignore.case = TRUE)])
    rownames(design_matrix) <- meta[, 1]
    colnames(design_matrix) <- gsub("\\W", "_", colnames(design_matrix))

    if (ncol(design_matrix) <= 0)
        warning("Did not find any design matrix information. You have to add this manually to the xls file.")

    return(list(
        table = meta,
        nrSamples = number_samples,
        samples = samples,
        seFiles = pdf$seFiles,
        peFiles = pdf$peFiles,
        seOrder = pdf$seOrder,
        peOrder = pdf$peOrder,
        seIndices = pdf$seIndices,
        peIndices = pdf$peIndices,
        readMode = pdf$readMode,
        conditions = conditions,
        designMatrix = design_matrix
    ))
}


#' Parse Metadata XLS Files
#'
#' Parse information contained in XLS metadata files.
#'
#' @export
#' @inheritParams parse_meta_csv
#' @inherit parse_meta_csv return return
#' @param sheet Integer(1). XLS files can have more than one sheet, but only
#'   one sheet can be loaded with this function. The page with this 'sheet'
#'   index will be used.
#' @examples
#' res <- make_meta_table(nrRows = 1, as.xls = TRUE, overwrite = TRUE)
#' parse_meta_xls(res$file)
parse_meta_xls <- function
(
    file,
    sheet = 1,
    lib_strat_col = "Library.Strategy"
){
    meta <- gdata::read.xls(file, sheet = sheet)
    # HACK: remove all factors.
    meta <- sapply(meta, function(x) as.character(x), simplify = FALSE)
    meta <- as.data.frame(meta, stringsAsFactors = FALSE)
    rownames(meta) <- gsub("\\s+", ".", rownames(meta))
    number_samples <- nrow(meta)
    samples <- meta[["Specimen"]]
    conditions <- meta[["Specimen"]]

    # split file information of Data.File into SE and PE file information
    pdf <- .parse_DATA_FILE(meta, lib_strat_col = lib_strat_col)

    # TODO: generalize that more?
    # design_matrix: which samples to compare and test for differential expression
    design_matrix <- as.matrix(meta[grep("_vs_", colnames(meta), ignore.case = TRUE)])
    rownames(design_matrix) <- meta[, 1]
    colnames(design_matrix) <- gsub("\\W", "_", colnames(design_matrix))

    if (ncol(design_matrix) <= 0)
        warning("Did not find any design matrix information. You have to add this manually to the XLS file.")

    return(list(
        table = meta,
        nrSamples = number_samples,
        samples = samples,
        seFiles = pdf$seFiles,
        peFiles = pdf$peFiles,
        seOrder = pdf$seOrder,
        peOrder = pdf$peOrder,
        seIndices = pdf$seIndices,
        peIndices = pdf$peIndices,
        readMode = pdf$readMode,
        conditions = conditions,
        designMatrix = design_matrix
    ))
}


#' Parse SDRF file
#'
#' Parse information contained in SDRF files.
#'
#' Usually, the first row of a SDRF file is the header (column names) and the
#' second and following rows are the samples. If this is not the case, use
#' 'skip' to tell \code{parse_SDRF} how many rows to skip before the first
#' sample row and 'col.name.row' to tell the header row.
#'
#' @export
#' @param file Path to SDRF file.
#' @param skip Number of lines to skip in SDRF file. Defaults to 1.
#' @param col.name.row Integer. Row index containing all column names. Defaults
#'    to 'skip'.
#' @param comment.char Character indicating lines to ignore. By default, read all.
#' @param sep Field separator of SDRF. Defaults to tab.
#' @inherit parse_meta_csv return
#' @seealso \url{https://wiki.nci.nih.gov/display/TCGA/Sample+and+Data+Relationship+Format}
parse_SDRF <- function
(
    file,
    skip = 1,
    col.name.row = skip,
    comment.char = "",
    sep = "\t"
) {
    if (!grepl(".[tc]sv$", file))
        warning("File ending should be '.csv' or '.tsv'. If you supplied an Excel file, please use the function 'parse_meta_xls' instead.")
    SDRF <- readLines(file)

    # find "Fold Change Assignment" column --> start of design matrix
    # it is one row ABOVE the anticipated header row
    fca_first_col <- grep(
        "Fold Change Assignment",
        unlist(strsplit(SDRF[col.name.row - 1], sep))
    )
    if (length(fca_first_col) == 0) {
        message("Fold Change Assignment column not found. Ignored.")
        fca_first_col <- 1
    }

    SDRF <- read.table(
        text = SDRF,
        sep = sep,
        as.is = TRUE,
        comment.char = comment.char,
        skip = skip,
        # remove comment char and transform into vector --> header
        col.names = unlist(strsplit(sub("#\\s?", "", SDRF[col.name.row]), "\t"))
    )
    SDRF <- SDRF[nchar(SDRF[, 1]) > 0,] # remove lines with no content

    number_samples <- nrow(SDRF)
    samples <- SDRF[["Specimen"]]
    conditions <- SDRF[["Specimen"]]

    # split file information of Data.File into SE and PE file information
    pdf <- .parse_DATA_FILE(SDRF)

    # design_matrix: which samples to compare and test for differential expression
    fca_cols <- grep("_vs_", colnames(SDRF), ignore.case = T)
    design_matrix <- as.matrix(SDRF[fca_cols[fca_cols >= fca_first_col]])
    rownames(design_matrix) <- SDRF[, 1]
    colnames(design_matrix) <- sub("\\.\\d$", "", colnames(design_matrix))
    colnames(design_matrix) <- gsub("\\W", "_", colnames(design_matrix))

    if (ncol(design_matrix) <= 0)
        warning(
            "Did not find any design matrix information. You have to add this manually to the SDRF file."
        )

    return(list(
        table = SDRF,
        nrSamples = number_samples,
        samples = samples,
        seFiles = pdf$seFiles,
        peFiles = pdf$peFiles,
        seOrder = pdf$seOrder,
        peOrder = pdf$peOrder,
        seIndices = pdf$seIndices,
        peIndices = pdf$peIndices,
        readMode = pdf$readMode,
        conditions = conditions,
        designMatrix = design_matrix
    ))
}


#' Update Metadata
#'
#' Updates a metadata data.frame based on the return object of a wrapper
#' function. This includes the wrapper functions for FastQC, Trimmomatic,
#' SortMeRNA, TopHat2, HISAT2, featureCounts and all DEGtools. The function may
#' determine 'ind' and 'paired' automatically.
#'
#' The arguments 'paired' and 'ind' can be determined automatically if, and only
#' if, the dataset is cleary single- or paired-end. In case of both (a mixed
#' dataset), the 'paired' variable must be defined. Based on that, either the
#' single-end or paired-end rows will be overwritten.
#'
#' @export
#' @param metaRes The return object of a parser function from this package.
#' @param toolRes The return object of a wrapper function from this package.
#' @param paired Logical(1). Indicates if the tool was used in paired-end or
#'   single-end mode. If NA, the function tries to determine it based on the
#'   given meta-data object.
#' @param ind Integer(n). Row indices to write the entries of 'toolRes' into.
#'   If NULL, the function tries to determine it based on the given meta-data
#'   object.
#' @param rank Integer(1). The rank of a tools in terms of execution, e.g.
#'   FastQC is usually the first tool. This parameter can be used to overwrite
#'   existing entries or create new ones. It is advised to set this value.
#'   If -1, the next non-existing rank in the metadata is used.
#' @return Updated version of the 'metaRes' input argument.
#' @examples
#' \dontrun{
#'   # update SDRFinfo with new metadata
#'   SDRFinfo <- updateMeta(
#'       metaRes = SDRFinfo,
#'       toolRes = fq_raw_res # result from FastQC wrapper
#'   )
#'   SDRFinfo <- updateMeta(
#'       metaRes = SDRFinfo,
#'       toolRes = trim_res, # result from Trimmomatic wrapper
#'       rank = 2
#'   )
#' }
updateMeta <- function(
    metaRes,
    toolRes,
    paired = NA,
    ind = NULL,
    rank = -1
){
    meta <- metaRes$table
    if (!is.data.frame(meta))
        stop("'meta' is not a data.frame!")
    if (!is.list(toolRes))
        stop("'res' is not a list!")

    # parse tool string
    tool <- tolower(toolRes$tool[[1]])
    if (length(tool) == 0)
        stop("'res' argument does not contain a 'tool' entry!")

    if (is.null(ind)) {
        # define ind
        if (metaRes$readMode == "single") {
            if (!is.na(paired) && paired == TRUE) {
                warning("meta file indicates single-end mode, but paired = TRUE was defined. Using paired-end indices. Make sure this wasn't an accident!")
                ind <- metaRes$peIndices
            } else {
                ind <- metaRes$seIndices
                paired = FALSE
            }
        } else if (metaRes$readMode == "paired") {
            if (!is.na(paired) && paired == FALSE) {
                warning("meta file indicates paired-end mode, but paired = FALSE was defined. Using single-end indices. Make sure this wasn't an accident!")
                ind <- metaRes$seIndices
            } else {
                ind <- metaRes$peIndices
                paired = TRUE
            }
        } else if (metaRes$readMode == "mixed") {
            if (is.na(paired))
                stop("Meta file indicates mixed read mode, but 'paired' is NA. Set 'paired' to either 'TRUE' or 'FALSE'!")
            else if (paired == TRUE)
                ind <- metaRes$peIndices
            else if (paired == FALSE)
                ind <- metaRes$seIndices
            else
                stop("Wrong value for 'paired': value must be boolean!")
        }
    } else {
        # ind was defined by user. Check for mode
        if (is.na(paired) || is.null(paired) || !is.logical(paired)) {
            if (tool != "deg" && metaRes$readMode == "mixed") {
                stop("Row indices were supplied, but mode is read mode is 'Mixed' and paired is 'NA'. Set 'paired' to either 'TRUE' or 'FALSE'!")
            } else if (metaRes$readMode == "single" || tool == "deg")
                paired = FALSE
            else if (metaRes$readMode == "paired")
                paired = TRUE
            stopifnot(is.logical(paired))
        }
    }

    # determine rank - if undefined
    if (!is.numeric(rank) || rank < 0) {
        if (!"Data.Processing.Software" %in% colnames(meta))
            rank <- 0
        else {
            rank <- 1
            while(paste0("Data.Processing.Software.", rank) %in% colnames(meta))
                rank <- rank + 1
        }
        message("Selected rank: ", shQuote(as.character(rank)))
    }

    # determine column index
    if (tool %in% c("fastqc", "sortmerna")) {

        # merge two entries into one column in paired-mode
        col <- if (rank == 0) "Data.Processing.Software" else paste0("Data.Processing.Software.", rank)
        meta[ind, col] <- paste(toolRes$tool, toolRes$version)

        col <- if (rank == 0) "Parameters.And.Values" else paste0("Parameters.And.Values.", rank)
        if (!paired)
            meta[ind, col] <- toolRes$calls
        else
            meta[ind, col] <- paste0(
                toolRes$calls[seq(1, length(toolRes$calls), by = 2)],
                ";",
                toolRes$calls[seq(2, length(toolRes$calls), by = 2)]
            )

    } else if (tool %in% c("trimmomatic", "hisat", "tophat", "featurecounts")) {

        col <- if (rank == 0) "Data.Processing.Software" else paste0("Data.Processing.Software.", rank)
        meta[ind, col] <- paste(toolRes$tool, toolRes$version)
        col <- if (rank == 0) "Parameters.And.Values" else paste0("Parameters.And.Values.", rank)
        meta[ind, col] <- toolRes$calls
    } else if (tool == "deg") {

        col <- if (rank == 0) "Data.Processing.Software" else paste0("Data.Processing.Software.", rank)
        meta[ind, col] <- paste(toolRes$tools, collapse = ",")
        col <- if (rank == 0) "Description" else paste0("Description.", rank)
        meta[ind, col] <- toolRes$desc
    } else
        stop("Unknown tool: ", toolRes$tool)

    newMetaRes <- metaRes
    newMetaRes$table <- meta
    return(newMetaRes)
}
