####-
##-====================================================
##-==========================  MAPPING  ===========================
##-====================================================
####-

#' Make TopHat2 (Bowtie2) genome index
#'
#' Wrapper function for creating Bowtie2 index files, which are used by
#' TopHat2.
#'
#' @export
#' @param genomeFile Path to genome FASTA file.
#' @param customOut Optional. Path to directory were the index files will be
#'   saved. By default, the 'genomeFile' directory is used.
#' @return A list with keyword arguments: \cr
#' \describe{
#'   \item{outpref}{Prefix of index files. This is used as input for Bowtie2 or
#'                  TopHat2.}
#'   \item{call   }{System call used to build the index files.}
#' }
#' @examples
#' \dontrun{
#' make_Tophat_index("path/to/genome.fa")
#' }
#' make_Tophat_index
make_Tophat_index <- function
(
    genomeFile,
    customOut = NULL
){
    build_exec <- .get_executable("bowtie2-build", "bowtie2", sysVar = "BOWTIE2_EXEC")

    writeLines(c("Make bowtie2 index using:", build_exec))
    outpref <- if (is.null(customOut)) tools::file_path_sans_ext(genomeFile) else customOut
    call <- paste(build_exec, genomeFile, outpref)
    system(call)

    return(list(outpref = outpref, call = call))
}


#' List Additional TopHat2 Arguments
#'
#' These arguments are not required to run TopHat2, but their usage usually
#' leads to better results. They are active by default and can be omitted by
#' setting 'no.default' to TRUE.
#'
#' @export
#' @param paired Logical indicating whether paired-end files are used.
#'   Defaults to FALSE.
#' @return Character(1).
#' @seealso \code{\link{run_Tophat}}
#' @examples args_Tophat_default()
args_Tophat_default <- function
(
    paired = FALSE
){
    if (paired)
        return("-g 2 --b2-very-sensitive --no-coverage-search --no-mixed --no-discordant")
    else
        return("-g 2 --b2-very-sensitive --no-coverage-search")
}


#' List TopHat2 Arguments
#'
#' List the arguments used for running TopHat2 by GEO2RNA-seq. More precisely,
#' these arguments are used for the function \code{\link{run_Tophat}}.
#'
#' @export
#' @param paired Logical indicating whether paired-end files are used.
#'   Defaults to FALSE.
#' @return A list with keyword arguments. Keywords are 'default' and 'locked'.
#' @seealso \code{\link{run_Tophat}}
#' @examples args_Tophat_default()
args_Tophat <- function
(
    paired = FALSE
){
    writeLines(c(
        "We need these arguments for input, output and multi-threading:",
        paste(args_Tophat_locked(paired), collapse = " "),
        "And those arguments in addition:",
        paste(args_Tophat_default(paired), collapse = " ")
    ))
    return(list(default = args_Tophat_default(paired), locked = args_Tophat_locked(paired)))
}


#' List Locked TopHat2 Arguments
#'
#' TopHat2 requires certain arguments to define input and output files, as well
#' as certain modes. These arguments are set up by the wrapper function and
#' should not be supplied using the 'addArgs' argument of
#' \code{\link{run_Tophat}}.
#'
#' @export
#' @param paired Logical indicating whether paired-end files are used.
#'   Defaults to FALSE.
#' @return Character(n).
#' @seealso \code{\link{run_Tophat}}
#' @examples args_Tophat_locked()
args_Tophat_locked <- function
(
    paired = FALSE
){
    if (paired)
        return(c("tophat2", "-p", "-o", "-G"))
    else
        return(c("tophat2", "-p", "-o", "-G"))
}


#' Run TopHat2
#'
#' Wrapper function for TopHat2, a tool for precise mapping of reads in FASTQ
#' format to a reference genome in FASTA format.
#'
#' The reference genome must be indexed once before it can be used. See
#' \code{\link{make_Tophat_index}} to index genome reference files.
#' \cr\cr
#' See \code{\link{args_Tophat_default}} for default TopHat2 arguments. It is
#' not allowed to set locked arguments defined by \code{\link{args_Tophat_locked}}
#' via 'addArgs' .
#'
#' @export
#' @param files List of files in FASTQ format.
#' @param index Path to index files INCLUDING the genome prefix.
#' @param outDir General directory for output files. Multiple subdirectory will
#'   be created. Output BAM files will put in '<outDir>/bamfiles/'.
#'   Defaults to "mapping".
#' @param is.paired Logical indicating if paired-end reads are used. Defaults to
#'   FALSE.
#' @param transcriptome Path to transcriptome files as prefix (NO file
#'   extensions!). If supplied, 'anno' is ignored and the given transcriptome is
#'   used instead.
#' @param anno Additional annotation. TopHat2 will try to map to these gene
#'   locations first. Any unmapped read will then be searched against the whole
#'   genome. Defaults to NA.
#' @param addArgs Additional arguments. It is not allowed to set locked
#'   arguments. They are set by the corresponding wrapper function. See details.
#'   Defaults to NA.
#' @param no.default If TRUE, do not use default parameters, see details.
#'   Defaults to FALSE.
#' @inheritParams parallelizationParameters
#' @param overwrite Logical indicating if TopHat2 should overwrite files if
#'   'outDir' already contains some files.
#'   If the directory is not empty and overwrite is FALSE, all BAM files found
#'   in '<outDir>/bamfiles' will be returned.
#' @param use.existing Logical indicating if existing output BAM files
#'   should be returned. Defaults to TRUE. If set to FALSE, execution will stop
#'   if existing BAM files would be overwritten.
#' @param tophat_exec Optional. Path to TopHat2 binary.
#' @references Kim D, Pertea G, Trapnell C, Pimentel H, Kelley R, Salzberg SL.
#'   TopHat2: accurate alignment of transcriptomes in the presence of insertions,
#'   deletions and gene fusions. Genome Biology 2013, 14:R36. Available online
#'   at: \url{https://ccb.jhu.edu/software/tophat/index.shtml}
#' @return A list with keyword arguments: \cr
#' \describe{
#'   \item{files  }{Paths to output BAM files.}
#'   \item{calls  }{System calls to TopHat2.}
#'   \item{version}{TopHat2 version.}
#'   \item{tool   }{The name of the tool.}
#' }
#' @examples
#' \dontrun{
#' files <- c("f1.fastq", "f2.fastq")
#' genome <- "genome.fa"
#' index <- "genome"
#' anno <- "annotation.gtf"
#' run_Tophat(files = files, index = index, anno = anno)
#'
#' files <- c("f1_1.fastq", "f1_2.fastq", "f2_1.fastq", "f2_2.fastq")
#' newArgs <- "-g 3"
#' run_Tophat(files = files, index = index, anno = anno, is.paired = TRUE,
#'            addArgs = newArgs, no.default = TRUE)
#' }
#' run_Tophat
run_Tophat <- function
(
    files,
    index,
    outDir = "./mapping",
    is.paired = FALSE,
    transcriptome = NA,
    anno = NA,
    addArgs = NA,
    no.default = FALSE,
    cpus = 1,
    workers = 5,
    overwrite = FALSE,
    use.existing = TRUE,
    tophat_exec = ""
){
    # check if input files exist
    if (FALSE %in% file.exists(asPairVector(files))) {
        f <- asPairVector(files)
        stop(paste(
            "Files do not exist:\n",
            paste(shQuote(f[!file.exists(f)]), collapse = "\n ")
        ))
    }

    # abort execution if one or more outfiles already exist but should not be used
    if (!use.existing) {
        bamDir <- file.path(outDir, "bamfiles")
        if (is.paired)
            f <- if (is.matrix(files)) files[,1] else asPaired(files)[,1]
        else
            f <- files

        bam_files <- paste0(tools::file_path_sans_ext(basename(f)), ".bam")
        bam_files <- file.path(bamDir, bam_files)
        if (TRUE %in% file.exists(bam_files))
            stop("One or more files already exist and use.existing is FALSE!")
    }

    writeLines("TopHat2 - mapping ...")
    if (is.paired) {
        resmapping <- .run_Tophat_paired(
            pairedFiles = files,
            index = index,
            outDir = outDir,
            transcriptome = transcriptome,
            anno = anno,
            addArgs = addArgs,
            no.default = no.default,
            cpus = cpus,
            workers = workers,
            overwrite = overwrite,
            tophat_exec = tophat_exec
        )
    } else {
        resmapping <- .run_Tophat_unpaired(
            fqFiles = files,
            index = index,
            outDir = outDir,
            transcriptome = transcriptome,
            anno = anno,
            addArgs = addArgs,
            no.default = no.default,
            cpus = cpus,
            workers = workers,
            overwrite = overwrite,
            tophat_exec = tophat_exec
        )
    }

    return(resmapping)
}


# Unpaired (Single-end) TopHat2
#
# This version is for single-end reads only.
#
# @rdname run_Tophat_unpaired
# @inherit run_Tophat
# @param fqFiles List of files in FASTQ format.
.run_Tophat_unpaired <- function
(
    fqFiles,
    index,
    outDir,
    transcriptome = NA,
    anno = NA,
    share = TRUE,
    addArgs = NA,
    no.default = FALSE,
    cpus = 1,
    workers = 5,
    overwrite = FALSE,
    tophat_exec = NA
){
    driver <- function(i) {
        # Note: make a directory for each FASTQ file
        fileName <- basename(tools::file_path_sans_ext(fqFiles[i]))
        tophatOutDir <- file.path(outDir, fileName)
        if (!file.exists(tophatOutDir))
            dir.create(tophatOutDir)
        bamFile <- file.path(bamdir, paste0(fileName, ".bam"))
        if (!overwrite && file.exists(bamFile)) {
            message("File already exists: ", shQuote(bamFile), ". Skipped!")
            return(list(bam = bamFile, call = "NOT USED"))
        }

        # set additional arguments. This includes:
        # 1. checking addArgs for locked arguments
        # 2. dealing with annotation files
        args <- if (!no.default) args_Tophat_default(paired=FALSE) else ""
        if (TRUE %in% (addArgs %in% args_Tophat_locked(paired=FALSE)))
            stop("You supplied one or more locked arguments! See 'args_Tophat_locked' for more information.")
        if (!(is.na(addArgs) || is.null(addArgs)))
            args <- paste(args, addArgs)
        if (!is.na(transcriptome))
            args <- paste(args, paste0("--transcriptome-index=", transcriptome))
        else if (!(is.na(anno) || is.null(anno))) {
          if (!file.exists(anno))
                stop("Given annotation file does not exist:", anno)
          else
                args <- paste(args, "-G", anno)
        }

        # make cmd call
        call = paste(
            tophat_exec,
            "-p",     allocCPUS$for_call[i],# number of threads to use
            args,
            "-o",     tophatOutDir,  # bam files will be put here
            index,            # genome index file
            fqFiles[i],    # input FASTQ file
            ">",
            file.path(outDir, paste0(fileName, ".tophatout.txt")), # just some statistics from tophat
            "2>&1"
        )
        system(call)

        # rename accepted_hits.bam to '<FASTQ file>.bam' and move it to './Mapping/Bamfiles'
        file.rename(from=file.path(tophatOutDir, "accepted_hits.bam"), to=bamFile )
        alignStatFile <- file.path(outDir, paste0(fileName, ".align_summary.txt"))
        file.rename(from=file.path(tophatOutDir, "align_summary.txt"), to=alignStatFile)

        return(list(bam = bamFile, call = call))
    }

    if (FALSE %in% file.exists(fqFiles))
        warning(paste0(
            "One or more files do not exist:\n",
            paste(fqFiles[!file.exists(fqFiles)], collapse="\n")
        ))
    if (tophat_exec == "" || is.null(tophat_exec) || is.na(tophat_exec))
        tophat_exec <- .get_executable("tophat2", "tophat", "TOPHAT2_PATH")
    else if (!file.exists(tophat_exec))
        stop(paste(
            "Cannot find TopHat2 binary: user given file does not exist at",
            "path:", tophat_exec
        ))
    version <- system(paste(tophat_exec, "--version"), intern = TRUE)

    bamdir <- file.path(outDir, "bamfiles")
    dir.create(bamdir, recursive = TRUE, showWarnings = FALSE)

    # check if transcriptome is valid and exists
    if (!is.na(transcriptome)) {
        if (!file.exists(paste0(transcriptome, ".fa"))) {
            if (!file.exists(paste0(tools::file_path_sans_ext(transcriptome), ".fa"))) {
                warning("Transcriptome path to long. Remove file endings for next run!")
                transcriptome <- tools::file_path_sans_ext(transcriptome)
            } else
                stop(paste(
                    "Cannot find FASTA of transcriptome. Supply only the prefix of files!",
                    transcriptome
                ))
        }
    } else if (!is.na(anno) && !is.null(anno) && file.exists(anno) && share) {
        # create transcriptome index beforehand.
        # This will speed up mapping on large genomes!
        transcriptome_dir   <- file.path(outDir, "transcriptome")
        transcriptome <- file.path(
            transcriptome_dir,
            tools::file_path_sans_ext(basename(anno))
        )
        call <- paste(
            tophat_exec,
            "-G",
            anno,
            paste0("--transcriptome-index=", transcriptome),
            "-o",
            transcriptome_dir, # also store transcriptome index logfiles here
            index
        )
        system(call)
    }

    allocCPUS <- allocate_cpus(cpus, length(fqFiles), workers)
    bio_par <- BiocParallel::MulticoreParam(
        workers = allocCPUS$in_parallel,
        tasks = length(fqFiles),
        progressbar = TRUE
    )
    allRes <- tryCatch({
        BiocParallel::bplapply(1:length(fqFiles), driver, BPPARAM = bio_par)
    }, error = identity)
    .stop_at_error(allRes, from = match.call()[[1]])

    return(list(
        files = sapply(allRes, function(x) x[["bam"]]),
        calls = sapply(allRes, function(x) x[["call"]]),
        version = version,
        tool = "tophat"
    ))
}


# Paired-end TopHat2
#
# This is for paired-end reads only.
#
# @rdname run_Tophat_paired
# @inherit run_Tophat
# @param pairedFiles Matrix or list of files in FASTQ format. Matrix must have
#   pairs as rows (first column is for forward reads, second
#   column is for backward reads).
.run_Tophat_paired <- function
(
    pairedFiles,
    index,
    outDir = "./mapping",
    transcriptome = NA,
    anno = NA,
    share = TRUE,
    addArgs = NA,
    no.default = FALSE,
    cpus = 1,
    workers = 5,
    overwrite = FALSE,
    tophat_exec = NA
){
    driver <- function(i) {
        # Note: make a directory for each FASTQ file
        fileName <- basename(tools::file_path_sans_ext(pairedFiles[i,1]))
        tophatOutDir <- file.path(outDir, fileName)
        bamFile <- file.path(bamdir, paste0(fileName, ".bam"))
        if (!file.exists(tophatOutDir))
            dir.create(tophatOutDir)
        if (!overwrite && file.exists(bamFile)) {
            message("File already exists: ", shQuote(bamFile), ". Skipped!")
            return(list(bam = bamFile, call = "NOT USED"))
        }

        # get additional arguments
        args <- if (!no.default) args_Tophat_default(paired=TRUE) else ""
        if (TRUE %in% (addArgs %in% args_Tophat_locked(paired=FALSE)))
            stop("You supplied one or more locked arguments! See 'args_Tophat_locked' for more information.")
        if (!(is.na(addArgs) || is.null(addArgs)))
            args <- paste(args, addArgs)
        if (!(is.na(anno) || is.null(anno))) {
            if (!file.exists(anno))
                stop("Given annotation file does not exist:", anno)
            else {
                args <- paste(args, "-G", anno)
                if (share)
                    args <- paste(
                        args,
                        paste0("--transcriptome-index=", transcriptome)
                    )
            }
        }

        call=paste(
            tophat_exec,
            "-p", allocCPUS$for_call[i], # number of threads to use
            args,
            "-o", tophatOutDir, # BAM files will be put here
            index,              # genome index file
            pairedFiles[i,1],   # 1 input FASTQ file
            pairedFiles[i,2],   # 2 input FASTQ file
            ">",
            file.path(outDir, paste0(fileName, ".tophatout.txt")), # just some statistics from TopHat2
            "2>&1"
        )
        system(call)

        # rename accepted_hits.bam to '<FASTQ file>.bam' and move it to './Mapping/Bamfiles'
        file.rename(from=file.path(tophatOutDir, "accepted_hits.bam"), to=bamFile)
        alignStatFile <- file.path(outDir, paste0(fileName, ".align_summary.txt"))
        file.rename(from=file.path(tophatOutDir, "align_summary.txt"), to=alignStatFile)

        return(list(bam = bamFile, call = call))
    }

    if (FALSE %in% file.exists(pairedFiles))
        warning(paste0(
            "One or more files do not exist:\n",
            paste(pairedFiles[!file.exists(pairedFiles)], collapse="\n")
        ))
    if (tophat_exec == "" || is.null(tophat_exec) || is.na(tophat_exec))
        tophat_exec <- .get_executable("tophat2", "tophat", "TOPHAT2_PATH")
    if (!file.exists(tophat_exec))
        stop(paste(
            "Cannot find TopHat2 binary: user given file does not exist at",
            "path:", tophat_exec
        ))
    version <- system(paste(tophat_exec, "--version"), intern = TRUE)

    if (is.vector(pairedFiles) || is.list(pairedFiles)) pairedFiles <- asPaired(pairedFiles)

    bamdir <- file.path(outDir, "bamfiles")
    dir.create(bamdir, recursive = TRUE, showWarnings = FALSE)

    # check if transcriptome is valid and exists
    if (!is.na(transcriptome)) {
        if (!file.exists(paste0(transcriptome, ".fa"))) {
            if (!file.exists(paste0(tools::file_path_sans_ext(transcriptome), ".fa"))) {
                warning("Transcriptome path to long. Remove file endings for next run!")
                transcriptome <- tools::file_path_sans_ext(transcriptome)
            } else
                stop(paste(
                    "Cannot find FASTA of transcriptome. Supply only the prefix of files!",
                    transcriptome
                ))
        }
    } else if (!is.na(anno) && !is.null(anno) && file.exists(anno) && share) {
        # create transcriptome index beforehand.
        # This will speed up mapping on large genomes!
        transcriptome_dir   <- file.path(outDir, "transcriptome")
        transcriptome <- file.path(
            transcriptome_dir,
            tools::file_path_sans_ext(basename(anno))
        )
        call <- paste(
            tophat_exec,
            "-G",
            anno,
            paste0("--transcriptome-index=", transcriptome),
            "-o",
            transcriptome_dir, # also store transcriptome index logfiles here
            index
        )
        system(call)
    }

    allocCPUS <- allocate_cpus(cpus, nrow(pairedFiles), workers)
    bio_par <- BiocParallel::MulticoreParam(
        workers = allocCPUS$in_parallel,
        tasks = nrow(pairedFiles),
        progressbar = TRUE
    )
    allRes <- tryCatch({
        BiocParallel::bplapply(1:nrow(pairedFiles), driver, BPPARAM = bio_par)
    }, error = identity)
    .stop_at_error(allRes, from = match.call()[[1]])

    return(list(
        files = sapply(allRes, function(x) x[["bam"]]),
        calls = sapply(allRes, function(x) x[["call"]]),
        version = version,
        tool = "tophat"
    ))
}


#' Make HISAT2 genome index
#'
#' Wrapper function for creating HISAT2 index files.
#'
#' If you want to use exon and splice site information, you need to use the
#' two python scripts (see arguments) which can be found in your HISAT2
#' installation directory.
#'
#' @export
#' @inheritParams make_Tophat_index
#' @inheritParams parallelizationParameters
#' @inherit make_Tophat_index return return
#' @param exonFile Tabular files as generated by 'hisat2_extract_exons.py'.
#' @param spliceFile Tabular file as generated by 'hisat2_extract_splice_sites.py'.
#' @examples
#' \dontrun{
#' make_HiSat2_index("path/to/genome.fasta")
#' make_HiSat2_index("path/to/genome.fasta",
#'  exonFile = "exons.tab",
#'  spliceFile = "splice.tab")
#' }
#' make_HiSat2_index
make_HiSat2_index <- function
(
    genomeFile,
    customOut = NULL,
    exonFile = "",
    spliceFile = "",
    cpus = 1
){
    build_exec <- .get_executable("hisat2-build", "hisat2", "HISAT2_EXEC")

    writeLines(c("Make HISAT2 index using:", build_exec))
    args <- paste("-p", cpus)
    if (exonFile != "") {
        if (file.exists(exonFile))
            args <- paste(args, "--exon", exonFile)
        else
            stop(paste("Exon file doesn't exist at path: ", exonFile))
    }

    if (spliceFile != "") {
        if (file.exists(spliceFile))
            args <- paste(args, "--ss", spliceFile)
        else
            stop(paste("Splice site files doesn't exist at path:", spliceFile))
    }

    outpref <- if (is.null(customOut)) tools::file_path_sans_ext(genomeFile) else customOut
    call <- paste(build_exec, args, genomeFile, outpref)
    system(call)

    return(list(outpref = outpref, call = call))
}


#' List Additional HISAT2 Arguments
#'
#' Return the default arguments used to run HISAT2.
#'
#' Input and output arguments are not included. They must not be
#' overwritten directly. See \code{\link{args_HiSat2_locked}} for further
#' information.
#'
#' @export
#' @param paired Show arguments for paired-end or single-end reads.
#' @return Character(1).
#' @examples args_HiSat2_default()
args_HiSat2_default <- function
(
    paired = FALSE
){
    if (paired) {
        return("-k 2 --no-unal --no-mixed --no-discordant")
    } else {
        return("-k 2 --no-unal")
    }
}


#' List HISAT2 Arguments
#'
#' List the arguments used for running HISAT2 by GEO2RNA-seq. More precisely,
#' these arguments are used for the function  \code{\link{run_Hisat2}}.
#'
#' @export
#' @param paired Logical indicating whether paired-end files are used.
#'   Defaults to FALSE.
#' @return A list with keyword arguments. Keywords are: 'default' and 'locked'.
#' @seealso \code{\link{run_Hisat2}}
#' @examples args_HiSat2()
args_HiSat2 <- function
(
    paired = FALSE
){
    writeLines(c(
        "We need these arguments for input, output and multi-threading:",
        paste(args_HiSat2_locked(), collapse = " "),
        "And those arguments in addition:",
        paste(args_HiSat2_default(paired), collapse = " ")
    ))
    return(list(
        default = args_HiSat2_default(paired),
        locked = args_HiSat2_locked()
    ))
}

#' List Locked HISAT2 Arguments
#'
#' HISAT2 requires certain arguments to define input and output files, as well
#' as certain modes. These arguments are set up by the wrapper function and
#' should not be supplied using the 'addArgs' argument of
#  \code{\link{run_Hisat2}}.
#'
#' @export
#' @seealso \code{\link{run_Hisat2}}, \code{\link{args_HiSat2}}
#' @return Character(n).
#' @examples args_HiSat2_locked()
args_HiSat2_locked <- function(){
    return(c("hisat2", "-p", "--summary-file", "-x", "-U"))
}


# NOTES:
# only SAM output
# -U for unpaired; -1 <m1> -2 <m2> for paired;
# for unpaired: --un, --al
# for paired: --no-mixed, --no-discordant, --un-conc & --al-conc

#' Run HISAT2
#'
#' Wrapper function for HISAT2, a very fast tool for aligning reads in FASTQ
#' format to a reference genome in FASTA format.
#'
#' The reference genome must be indexed once before it can be used. See
#' \code{\link{make_HiSat2_index}} to index genome reference files. Exon and
#' splice site information must be generated from outside the R environment.
#' See the description of \code{\link{make_HiSat2_index}} for more details.
#' \cr\cr
#' See \code{\link{args_HiSat2_default}} for default HISAT2 arguments. It is
#' not allowed to set locked arguments defined by \code{\link{args_HiSat2_locked}}
#' via 'addArgs' .
#'
#' @export
#' @inheritParams run_Tophat
#' @param splice Additional HISAT2 specific argument.
#' @param as.bam If TRUE, HISAT2 output will be converted to BAM files.
#'   If SAMtools is installed on your system, this will be done in memory (no
#'   SAM files). Defaults to TRUE.
#' @param keepSam Only has an effect if 'as.bam=TRUE'. If FLASE, intermediate
#'   SAM files are removed. If TRUE, they are kept. Defaults to FALSE.
#' @param hisat_exec Path to HISAT2 binary.
#' @return A list with keywords:\cr
#' \describe{
#'   \item{files       }{Character vector. Path to alignment output files.}
#'   \item{summaryFiles}{Character vector. Paths to output summary files.}
#'   \item{calls       }{Character vector. System calls to HISAT2.}
#'   \item{version     }{Version of HISAT2.}
#'   \item{tool        }{The name of the tool.}
#' }
#' @references Kim, Daehwan, Ben Langmead, and Steven L. Salzberg. "HISAT: a
#'   fast spliced aligner with low memory requirements." Nature methods 12.4
#'   (2015): 357-360.
#' Available online at: \url{http://www.ccb.jhu.edu/software/hisat/index.shtml}
#' @examples
#' \dontrun{
#' files <- c("f1.fastq", "f2.fastq")
#' genome <- "genome.fa"
#' index <- "genome"
#' run_Hisat2(files = files, index = index)
#'
#' files <- c("f1_1.fastq", "f1_2.fastq", "f2_1.fastq", "f2_2.fastq")
#' newArgs <- "-k 3"
#' run_Hisat2(files = files, index = index, is.paired = TRUE,
#'            addArgs = newArgs, no.default = TRUE)
#' }
#' run_Hisat2
run_Hisat2 <- function
(
    files,
    index,
    outDir = "./mapping",
    is.paired = FALSE,
    addArgs = "",
    splice = TRUE,
    no.default = FALSE,
    as.bam = TRUE,
    keepSam = FALSE,
    cpus = 1,
    workers = 5,
    overwrite = FALSE,
    use.existing = TRUE,
    hisat_exec = ""
){
    # check if input files exist
    if (FALSE %in% file.exists(asPairVector(files))) {
        f <- asPairVector(files)
        stop(paste(
            "Files do not exist:\n",
            paste(shQuote(f[!file.exists(f)]), collapse = "\n ")
        ))
    }

    if (!use.existing) {
        bamDir <- file.path(outDir, "bamfiles")
        samDir <- file.path(outDir, "samfiles")
        if (is.paired)
            f <- if (is.matrix(files)) files[,1] else asPaired(files)[,1]
        else
            f <- files

        bam_files <- paste0(tools::file_path_sans_ext(basename(f)), ".bam")
        bam_files <- file.path(bamDir, bam_files)
        if (TRUE %in% file.exists(bam_files))
            stop("One or more files already exist in BAM directory and use.existing is FALSE!")
        sam_files <- paste0(tools::file_path_sans_ext(basename(f)), ".sam")
        sam_files <- file.path(samDir, sam_files)
        if (TRUE %in% file.exists(sam_files))
            stop("One or more files already exist in sam directory and use.existing is FALSE!")
    }

    writeLines("HISAT2 - mapping ...")
    if (is.paired) {
        res <- .run_HiSat2_paired(
            files = files,
            outDir = outDir,
            index = index,
            args = addArgs,
            use_default = !no.default,
            splice = splice,
            as.bam = as.bam,
            keepSam = keepSam,
            cpus = cpus,
            workers = workers,
            overwrite = overwrite,
            hisat_exec = hisat_exec
        )
    } else {
        res <- .run_HiSat2_unpaired(
            files = files,
            outDir = outDir,
            index = index,
            args = addArgs,
            use_default = !no.default,
            splice = splice,
            as.bam = as.bam,
            keepSam = keepSam,
            cpus = cpus,
            workers = workers,
            overwrite = overwrite,
            hisat_exec = hisat_exec
        )
    }
    return(res)
}


# Run HISAT2 Single-End
#
# This wrapper function is for single-end files only.
#
# @rdname run_HiSat2_unpaired
# @inherit run_Hisat2
.run_HiSat2_unpaired <- function
(
    files,
    outDir,
    index,
    args = "",
    use_default=TRUE,
    splice = TRUE,
    as.bam = TRUE,
    keepSam = FALSE,
    cpus = 1,
    workers = 10,
    overwrite = FALSE,
    hisat_exec = NA
) {
    driver <- function(i) {
        fileName   <- tools::file_path_sans_ext(basename(files[i]))
        samOut     <- file.path(samDir, paste0(fileName, ".sam"))
        bamOut     <- file.path(bamDir, paste0(fileName, ".bam"))
        summaryOut <- file.path(outDir, paste0(fileName, ".summary.txt"))

        # SAM or BAM files may already exist. SAM files may just be converted to BAM files
        if (!overwrite && as.bam && file.exists(bamOut)) {
            message("\nBAM file already exists:", shQuote(bamOut), " . Skipped!\n")
            return(list(map = bamOut, summary=summaryOut, call="NOT USED"))
        }

        if (!overwrite && file.exists(samOut)) {
            message("\nSAM file already exists:", shQuote(samOut), " . Skipped!\n")
            if (as.bam) {
                writeLines("\nConverting to BAM...")
                Rsamtools::asBam(samOut, sub("\\.bam$", "", bamOut), indexDestination=FALSE)

                if (!keepSam)
                  file.remove(samOut)
                return(list(map = bamOut, summary=summaryOut, call="NOT USED"))
            } else {
                return(list(map = samOut, summary=summaryOut, call="NOT USED"))
            }
        }

        if (use_default)
            args <- paste(args, args_HiSat2_default())
        if (!splice)
            args <- paste(args, "--no-spliced-alignment")
        if (compareVersion(version, "2.1.0") >= 0)
            args <- paste(args, "--new-summary --summary-file", summaryOut)
        call <- paste(
            hisat_exec,
            args,
            "-p", allocCPUS$for_call[i],
            "-x", index,
            "-U", shQuote(files[i])
            #"-S", samOut
        )
        if (as.bam && !is.na(samtools_exec))
            call <- paste(call, "|", samtools_exec, "view -bS - >", bamOut)
        else
            call <- paste(call, "-S", samOut)
        # run HISAT2
        msg <- system(call, intern = TRUE)
        message("\n", msg, "\n")

        # if SAMtools is not available, use Rsamtools for BAM conversion
        if (as.bam) {
            if (!is.na(samtools_exec)) {
                return(list(map = bamOut, summary=summaryOut, call=call))
            } else {
                writeLines(paste("\nConverting:", samOut))
                Rsamtools::asBam(samOut, sub("\\.bam$", "", bamOut), indexDestination=FALSE)
                if (!keepSam)
                    file.remove(samOut)
                return(list(map = bamOut, summary=summaryOut, call=call))
            }
        } else
            return(list(map = samOut, summary=summaryOut, call=call))
    }


    if (hisat_exec == "" || is.null(hisat_exec) || is.na(hisat_exec))
        hisat_exec <- .get_executable("hisat2", sysVar = "HISAT2_EXEC")
    else if (!file.exists(hisat_exec))
        stop(paste(
            "Cannot find HISAT2 binary: user given file does not exist at",
            "path:", hisat_exec
        ))

    # NOTE: do NOT add 'hisat' as prefix. Version is needed as arguments!
    version <- grep(
        "hisat.*version",
        system(paste(hisat_exec, "--version"), intern = TRUE),
        value = TRUE
    )
    version <- gdata::last(unlist(strsplit(version, split = " ")))

    samtools_exec <- tryCatch({
        .get_executable("samtools", sysVar = "SAMTOOLS_EXEC")
    }, error = function(e) {
        return(NA)
    })

    if (is.na(samtools_exec) && as.bam)
        message("Could not find SAMtools on your system. In-memory BAM conversion disabled.")

    samDir <- file.path(outDir, "samfiles")
    if (!file.exists(samDir))
        dir.create(samDir, recursive = TRUE, showWarnings = FALSE)
    bamDir <- file.path(outDir, "bamfiles")
    if (!file.exists(bamDir))
        dir.create(bamDir, recursive = TRUE, showWarnings = FALSE)

    allocCPUS <- allocate_cpus(cpus, length(files), if (!as.bam) 2 else workers)
    bio_par <- BiocParallel::MulticoreParam(
        workers = allocCPUS$in_parallel,
        tasks = length(files),
        progressbar = TRUE
    )
    allRes <- tryCatch({
        BiocParallel::bplapply(1:length(files), driver, BPPARAM = bio_par)
    }, error = identity)
    .stop_at_error(allRes, from = match.call()[[1]])

    return(list(
        files        = sapply(allRes, function(x) x[["map"]]),
        summaryFiles = sapply(allRes, function(x) x[["summary"]]),
        calls        = sapply(allRes, function(x) x[["call"]]),
        version      = version,
        tool         = "hisat"
    ))
}


# Run Paired HISAT2
#
# This wrapper function is for paired-end files only.
#
# @rdname run_HiSat2_paired
# @inherit run_Hisat2
.run_HiSat2_paired <- function
(
    files,
    outDir,
    index,
    args = "",
    use_default=TRUE,
    splice = TRUE,
    as.bam = TRUE,
    keepSam = FALSE,
    cpus = 1,
    workers = 5,
    overwrite = FALSE,
    hisat_exec = NA
) {
    driver <- function(i) {
        fileName   <- tools::file_path_sans_ext(basename(files[i,1]))
        samOut     <- file.path(samDir, paste0(fileName, ".sam"))
        bamOut     <- file.path(bamDir, paste0(fileName, ".bam"))
        summaryOut <- file.path(outDir, paste0(fileName, ".summary.txt"))

        # SAM or BAM files may already exist. SAM files may just be converted to BAM files
        if (!overwrite && as.bam && file.exists(bamOut)) {
            message("\nBAM file already exists:", shQuote(bamOut), " . Skipped!\n")
            return(list(map = bamOut, summary=summaryOut, call="NOT USED"))
        }

        if (!overwrite && file.exists(samOut)) {
            message("\nSAM file already exists:", shQuote(samOut), " . Skipped!\n")
            if (as.bam) {
                writeLines("\nConverting to BAM...\n")
                Rsamtools::asBam(samOut, sub("\\.bam$", "", bamOut), indexDestination=FALSE)

                if (!keepSam)
                  file.remove(samOut)
                return(list(map = bamOut, summary=summaryOut, call="NOT USED"))
            } else {
                return(list(map = samOut, summary=summaryOut, call="NOT USED"))
            }
        }

        if (use_default)
            args <- paste(args, args_HiSat2_default(paired = TRUE))
        if (!splice)
            args <- paste(args, "--no-spliced-alignment")
        if (compareVersion(version, "2.1.0") >= 0)
            args <- paste(args, "--new-summary --summary-file", summaryOut)

        call <- paste(
            hisat_exec,
            args,
            "-p", allocCPUS$for_call[i],
            "-x", index,
            "-1", shQuote(files[i,1]),
            "-2", shQuote(files[i,2])
            #"-S", samOut
        )

        if (as.bam && !is.na(samtools_exec))
            call <- paste(call, "|", samtools_exec, "view -bS - >", bamOut)
        else
            call <- paste(call, "-S", samOut)
        # run HISAT2
        msg <- system(call, intern = TRUE)
        message("\n", msg, "\n")

        # if SAMtools is not available, use Rsamtools for BAM conversion
        if (as.bam) {
            if (!is.na(samtools_exec)) {
                return(list(map = bamOut, summary=summaryOut, call=call))
            } else {
                Rsamtools::asBam(sub("\\.bam$", "", bamOut), indexDestination=FALSE)
                if (!keepSam)
                    file.remove(samOut)
                return(list(map = bamOut, summary=summaryOut, call=call))
            }
        } else
            return(list(map = samOut, summary=summaryOut, call=call))
    }

    # nothing exists or overwrite == TRUE
    if (hisat_exec == "" || is.null(hisat_exec) || is.na(hisat_exec))
        hisat_exec <- .get_executable("hisat2", sysVar = "HISAT2_PATH")
    else if (!file.exists(hisat_exec))
        stop(paste(
            "Cannot find HISAT2 binary: user given file does not exist at",
            "path:",
            hisat_exec
        ))

    version <- grep(
        "hisat.*version",
        system(paste(hisat_exec, "--version"), intern = TRUE),
        value = TRUE
    )
    version <- gdata::last(unlist(strsplit(version, split = " ")))

    # convert to 2 column matrix, if neccessary
    if (is.vector(files) || is.list(files))
        files <- asPaired(files)

    samtools_exec <- tryCatch({
        .get_executable("samtools", sysVar = "SAMTOOLS_EXEC")
    }, error = function(e) {
        return(NA)
    })

    if (is.na(samtools_exec) && as.bam)
        message("Could not find SAMtools on your system. In-memory BAM conversion disabled.")

    samDir <- file.path(outDir, "samfiles")
    dir.create(samDir, recursive = TRUE, showWarnings = FALSE)
    bamDir <- file.path(outDir, "bamfiles")
    dir.create(bamDir, recursive = TRUE, showWarnings = FALSE)

    allocCPUS <- allocate_cpus(cpus, nrow(files), if (!as.bam) 2 else workers)
    bio_par <- BiocParallel::MulticoreParam(
        workers = allocCPUS$in_parallel,
        tasks = nrow(files),
        progressbar = FALSE
    )
    allRes <- tryCatch({
        BiocParallel::bplapply(1:nrow(files), driver, BPPARAM = bio_par)
    }, error = identity)
    .stop_at_error(allRes, from = match.call()[[1]])

    return(list(
        files        = sapply(allRes, function(x) x[["map"]]),
        summaryFiles = sapply(allRes, function(x) x[["summary"]]),
        calls        = sapply(allRes, function(x) x[["call"]]),
        version      = version,
        tool         = "hisat"
    ))
}
