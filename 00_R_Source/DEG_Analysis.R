#####-
##-====================================================
##-=========================  STATISTICS ==========================
##-====================================================
#####-


#' Create a design matrix based on the given metadata
#'
#' A design matrix describes, in each column, which samples from the control
#' group and which others from the treatment group should be compaired to each
#' other. Each column represents a test. Each test represents a comparison
#' between two conditions. The tests are independent from each other. Usually,
#' one is interested in differentially expressed genes between the two
#' groups/conditions.
#'
#' A design matrix will only contain the value 'treatment', 'control' or 'none'.
#' Column names have the form <A>_VS_<B>, where <A> is a variable string describing
#' the treatent, and <B> is a variable string describing the control group.
#'
#' @export
#' @param meta data.frame. Samples in rows. Samples must be given as row.names
#'   or by a column names 'Sample'. Must contain at least one column describing
#'   condition and may contain another column containing time information.
#'   See examples below.
#' @param condCol Integer or Character. Define the column containing condition
#'   annotation. Defaults to 'Condition'.
#' @param condVal Boolean or Character. If TRUE, all pairwise combinations of
#'   conditions are calculated. Otherwise, supply the name of a condition value
#'   (as contained in 'condCol'-column). In this case, all pairwise comparisons
#'   against the supplied value are calculated. Defaults to 'TRUE'.
#' @param timeCol Integer or Character. Optional. If supplied, condition and
#'   time information can be combined to perform more precise comparisons.
#'   See 'modes' argument. Defaults to NULL.
#' @param timeVal Boolean or Character. Optional. If TRUE, all pairwise combinations
#'   of times are calculated. Otherwise, supply the name of a time value
#'   (as contained in 'timeCol'-column). In this case, all pairwise comparisons
#'   against the supplied value are calculated. Defaults to 'TRUE'.
#' @param mode Character. Defines which comparisons to create. Must be one of either:
#' \describe{
#'   \item{ALL}{        Create all tests based on condition and time.}
#'   \item{COND}{       Only consider condition column.}
#'   \item{CombFixCond}{Use condition and time column. Compare different time points at same condition.}
#'   \item{CombFixTime}{Use condition and time column. Compare different condition at same time.}
#'   \item{Combined}{   Return the results for CombFixCond and CombFixTime.}
#' }
#' @param groupConds list. Can be used to group together values in the condition column.
#'   Must have the form: list( "groupName1" = c(...), "groupName2" = c(...))
#' @param groupTimes list. Same as groupConds, but for time vector
#' @param switchTreatContrl Vector of test indices (from the 'tests' return
#'   value), where treatment and control should be switched: "A_VS_B" becomes
#'   "B_VS_A". An empty vector by default, hence no switching. The order of
#'   conditions in each test will be based on the order of appearance in 'meta'.
#' @return
#' \describe{
#'   \item{dm}{Design Matrix as data.frame}
#'   \item{tests}{List of pairwise tests}
#' }
#' @examples
#'
#'  samples <- paste(rep(c("A", "B", "C"), each = 3), rep(1:3, 3), sep = "-")
#'  meta <- data.frame(Sample = samples, row.names = samples)
#'  meta$Condition <- rep(c("A", "B", "C"), each = 3)
#'  meta$Time <- rep(1:3, 3)
#'
#'  # Only compare conditions - All against All
#'  createDesignMatrix(meta)
#'  # Only compare conditions - All against A
#'  createDesignMatrix(meta, condVal = "A")
#'  # Only compare times - All against All
#'  createDesignMatrix(meta, condCol = "Time")
#'
#'  # group time points together
#'  createDesignMatrix(meta, condCol = "Time", groupConds = list('12h' = c(1, 2)))
#'
#'  # Combine condition and time - Compare Time condition-wise
#'  createDesignMatrix(meta, timeCol = "Time", mode = "CombFixCond")
#'  # Combine condition and time - Compare Conditions time-wise
#'  createDesignMatrix(meta, timeCol = "Time", mode = "CombFixTime")
#'
#'  # even more grouping points together
#'  createDesignMatrix(meta, timeCol = "Time", mode = "CombFixCond", groupTimes = list('12h' = c("1", "2")), groupConds = list('AB' = c("A", "B")))
createDesignMatrix <- function(
    meta,
    condCol = "Condition",
    condVal = TRUE,
    timeCol = NULL,
    timeVal = TRUE,
    mode = "COND",
    groupConds = NULL,
    groupTimes = NULL,
    switchTreatContrl = c()
) {
    ### 0. parse data and fix potential issues

    if (is.integer(condCol) && (condCol < 1 || condCol > ncol(meta)))
        stop("condCol value is integer and out of bounds! Check your number of columns!")
    if (is.character(condCol) && !(condCol %in% colnames(meta)))
        stop("condCol value is character. Given name of condition column does not exist in given meta table!")
    if (is.integer(timeCol) && (timeCol < 1 || timeCol > ncol(meta)))
        stop("timeCol value is integer and out of bounds! Check your number of columns!")
    if (is.character(timeCol) && !(timeCol %in% colnames(meta)))
        stop("timeCol value is character. Given name of time column does not exist in given meta table!")

    if (is.null(rownames(meta)) && !("Sample" %in% colnames(meta)))
        stop("Please add a column named 'Sample' or add row names to the object given as 'meta' argument!")

    if (!("Sample" %in% colnames(meta)))
        meta <- cbind(Sample = row.names(meta), meta)

    if (mode %in% c("ALL", "CombFixCond", "CombFixTime", "Combined") && is.null(timeCol))
        stop("For mode ", shQuote(mode), ", the 'timeCol' argument must be defined!")

    # replace "" strings by something meaningful
    inds <- which(meta[,condCol] == "" | is.null(meta[,condCol]))
    if (length(inds) > 0)
        meta[inds, condCol] <- "none"

    # replace '_' sign by '-'
    meta[,condCol] <- gsub("_", "-", meta[,condCol])
    if (condVal != TRUE)
        condVal <- gsub("_", "-", condVal)
    if (timeVal != TRUE)
        timeVal <- gsub("_", "-", timeVal)
    if (!is.null(timeCol))
        meta[,timeCol] <- gsub("_", "-", meta[,timeCol])

    # group data together, if specified
    if (!is.null(groupConds)) {
        # replace "_" with "-"
        groupConds <- lapply(groupConds, function(x) gsub("_", "-", x))
        for (i in 1:length(groupConds)) {
            rind <- which(meta[,condCol] %in% groupConds[[i]])
            meta[rind, condCol] <- names(groupConds)[[i]]
        }
    }

    if (!is.null(groupTimes)) {
        if (is.null(timeCol))
            stop("Cannot group time points: 'timeCol' is NULL!")
        # replace "_" with "-"
        groupTimes <- lapply(groupTimes, function(x) gsub("_", "-", x))

        for (i in 1:length(groupTimes)) {
            rind <- which(meta[,timeCol] %in% groupTimes[[i]])
            meta[rind, timeCol] <- names(groupTimes)[[i]]
        }
    }

    # replace any string resembling '_vs_' and grab unique values
    meta[,condCol] <- gsub("_[vV][sS]_", "_x_", meta[,condCol])
    conds <- unique(meta[,condCol])
    if (!is.null(timeCol)) {
        meta[,timeCol] <- gsub("_[vV][sS]_", "_x_", meta[,timeCol])
        times <- unique(meta[,timeCol])
        if (length(unique(times)) <= 1)
            stop("There must be more than 1 distinct group! If 'groupTimes' was used, make sure not to group all values together!")
    }

    if (length(unique(conds)) <= 1)
        stop("There must be more than 1 distinct group! If 'groupConds' was used, make sure not to group all values together!")


    ### 1. Create list of pairwise tests to perform
    if(mode == "ALL") {

        # (3) Generate all possible tests
        ct.comb <- unique(paste(meta[,condCol], meta[,timeCol], sep = "_"))
        if (length(ct.comb) >= 45)
            stop("Cannot make combinations with more than 44 distinct feature.")
        cond.comb.mat <- gtools::combinations(length(ct.comb), 2, ct.comb, repeats.allowed = FALSE)
        cond.comb.list <- lapply(1:nrow(cond.comb.mat), function(r) cond.comb.mat[r,])

    } else {

        # Make Condition Combinations - MANDETORY
        if (condVal == TRUE) {
            # generate all condition combinations
            if (length(conds) >= 45)
                stop("Cannot make combinations with more than 44 distinct feature. You have to reduce the number of different conditions!")
            cond.comb.mat <- gtools::combinations(length(conds), 2, conds, repeats.allowed = FALSE)
            cond.comb.list <- lapply(1:nrow(cond.comb.mat), function(r) cond.comb.mat[r,])
        } else {
            # generate combinations relative to one condition
            if (!(condVal %in% conds))
                stop("Argument 'condVal' must be either TRUE or an existing column name!")
            cond.comb.list <- lapply(conds[condVal != conds], function(v) c(condVal, v))
        }

        # Make Time Combinations - OPTIONAL
        if (!is.null(timeCol)) {
            if (timeVal == TRUE) {
                if (length(times) >= 45)
                    stop("Cannot make combinations with more than 44 distinct features. You have to reduce the number of different times!")
                time.comb.mat <- gtools::combinations(length(times), 2, times, repeats.allowed = FALSE)
                time.comb.list <- lapply(1:nrow(time.comb.mat), function(r) time.comb.mat[r,])
            } else {
                if (!(timeVal %in% times))
                    stop("Argument 'timeVal' must be either NULL, TRUE or an existing column name!")
                time.comb.list <- lapply(times[timeVal != times], function(v) c(timeVal, v))
            }
        }

        .makeCombined <- function(mode) {
            if (mode == "CombFixCond") {
                l <- lapply(conds,
                            function(cond) lapply(time.comb.list,
                                                  function(tpair) paste(paste(cond, tpair, sep = "_"),
                                                                        collapse = "_VS_")))
                return(unlist(l, recursive = T))
            } else if (mode == "CombFixTime") {
                l <- lapply(times,
                            function(time) lapply(cond.comb.list,
                                                  function(cpair) paste(paste(cpair, time[1], sep = "_"),
                                                                        collapse = "_VS_")))
                return(unlist(l, recursive = T))
            }
        }

        # Combination sector
        if (mode == "COND") {
            pairwise.tests <- sapply(cond.comb.list, function(cpair) paste(cpair, collapse = "_VS_"))
        } else if (mode == "CombFixCond") {
            pairwise.tests <- .makeCombined(mode)
        } else if (mode == "CombFixTime") {
            pairwise.tests <- .makeCombined(mode)
        } else if (mode == "Combined") {
            pairwise.tests <- c(.makeCombined("CombFixCond"), .makeCombined("CombFixTime"))
        } else
            stop("Unknown mode: ", shQuote(mode), ". See documentation for advice.")

        # switch control and treatment - OPTIONAL
        if (length(switchTreatContrl) > 0)
            pairwise.tests[switchTreatContrl] <- lapply(pairwise.tests[switchTreatContrl], function(pair) paste(rev(unlist(strsplit(pair, "_VS_"))), collapse = "_VS_"))
    }


    # 2. Start of table creation
    l <- lapply(pairwise.tests, function(test) {
        ret <- rep("none", length(meta$Sample))
        vs <- strsplit(test, "_VS_", fixed = T)[[1]]

        if (!is.null(timeCol) && mode != "COND") {
            # split up time condition
            tmp <- strsplit(vs, "_", fixed = T)
            treat <- tmp[[1]]
            contr <- tmp[[2]]
        } else {
            treat <- vs[[1]]
            contr <- vs[[2]]
        }

        if (length(treat) == 1) {
            stopifnot(length(contr) == 1)
            treatInds <- which(meta[, condCol] == treat[1])
            ret[treatInds] <- "treatment"
            contrInds <- which(meta[, condCol] == contr[1])
            ret[contrInds] <- "control"
        } else {
            stopifnot(length(contr) == 2 && length(treat) == 2)
            treatInds <- which(meta[, condCol] == treat[1] & meta[, timeCol] == treat[2])
            ret[treatInds] <- "treatment"
            contrInds <- which(meta[, condCol] == contr[1] & meta[, timeCol] == contr[2])
            ret[contrInds] <- "control"
        }

        if (!("control" %in% unique(ret) && "treatment" %in% unique(ret)))
            warning("Dropping comparison: ", shQuote(test), " - No Compatible Samples found.")
        return(ret)
    })
    names(l) <- unlist(pairwise.tests)
    designMatrix <- do.call(cbind, l)
    row.names(designMatrix) <- meta$Sample

    return(list(dm = designMatrix, tests = pairwise.tests, meta = meta))
}



#' Sort Differential Expression Table
#'
#' Sort a table, as returned by \code{\link{calculate_DEGs}}, descending by
#' differentially expressed genes (DEGs), by regulation in general (up and down)
#' or particularly by upregulation.
#'
#' The differential expression table returned by \code{\link{calculate_DEGs}}
#' usually contains multiple columns for multiple tools, which tested all given
#' genes for differential expression (TRUE: gene is differentially expressed;
#' FALSE: gene is not differentially expressed). By default, the table is sorted
#' by positive hits (TRUEs) in these columns. The tool with the highest number
#' of DEGs is sorted first, then the tool with the second most hits, and so on.
#'
#' @export
#' @param degTable A data.frame. Usually a single list element from a list
#'   returned by \code{\link{calculate_DEGs}}.
#' @param tools Vector with names/tools of performed differential expression
#'   tests.
#' @param remInf  Logical. Remove values of infinity?
#' @param orderByLog Logical. Order by absolute logarithmic fold change?
#' @param orderUpReg Logical. Order by upregulation?
#' @return Sorted degTable.
#' @examples
#' \dontrun{
#' tools <- c("DESeq", "DESeq2")
#' deg_res <- calculate_DEGs("some data", tools)
#' sortedDEG <- sortDiffExp(deg_res$DEGs[1], tools)
#' }
sortDiffExp <- function
(
    degTable,
    tools,
    remInf = FALSE,
    orderByLog = FALSE,
    orderUpReg = FALSE
) {
    # make sure to compensate for different caps in names
    found <- tolower(colnames(degTable)) %in% tolower(tools)
    tools <- colnames(degTable)[found]

    degTable <- degTable[order(rowSums(as.matrix(degTable[, tools])), decreasing = TRUE),]

    # additional sorting and operations
    if (orderByLog)
        degTable <- degTable[order(abs(degTable[, "log2FC"]), decreasing = TRUE),]
    if (orderUpReg)
        degTable <- degTable[order(degTable[, "log2FC"], decreasing = TRUE),]
    if (remInf)
        degTable <- degTable[!(is.infinite(degTable[, "log2FC"])),]

    return(degTable)
}


# Run DESeq
#
# Find significantly, differentially expressed genes using DESeq.
# Assumes negative binomial distribution of counts.
#
# @rdname run_DESeq
# @param counts Integer matrix of counting data. Must be RAW.\cr
#   Genes in rows, samples in columns.
# @param conds Character vector of at most two different
#   conditions for each column of 'counts'.
# @param method defaults to per-condition. See DESeq or DESeq manual:
#   'estimateDispersions'.
# @param sharingMode defaults to maximum. See DESeq or DESeq manual:
#   'estimateDispersions'.
# @param fitType Defaults to parametric. See DESeq or DESeq manual:
#   'estimateDispersions'.
# @seealso DESeq documentation. Available online at:
#   \url{https://bioconductor.org/packages/release/bioc/manuals/DESeq/man/DESeq.pdf}.
# @return A data.frame. See \code{\link{DESeq::nbinomTest}}.
.run_DESeq <- function
(
    counts,
    conds,
    method = "per-condition",
    sharingMode = "maximum",
    fitType = "parametric"
) {
    writeLines("DESeq ...")
    requireNamespace("locfit", quietly = TRUE)
    cond1 <- unique(conds)[1]
    cond2 <- unique(conds)[2]
    cds <- DESeq::newCountDataSet(countData = counts, conditions = conds)
    cds <- DESeq::estimateSizeFactors(cds)

    # estimateDispersions() could throw error:
    # ------------------------------------------------
    # Error in parametricDispersionFit(means, disps) :
    # Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')
    # Additional: Warning message:
    #    glm.fit: algorithm did not converge
    # ------------------------------------------------
    cds <- tryCatch(
        DESeq::estimateDispersions(
            cds,
            method = method,
            sharingMode = sharingMode,
            fitType = fitType
        ),
        # method <- c( "pooled", "per-condition", "blind" )
        error = function(error_message) {
            warning(paste("DEseq:", error_message))
            return(NA)
        }
    )
    if (!isS4(cds)) {
        res <- list("padj" = NA)
        return(res)
    }

    res <- DESeq::nbinomTest(cds, cond1, cond2)
    return(res)
}



# Run DESeq2
#
# Find significantly, differentially expressed genes using DESeq2.
# Assumes negative binomial distribution of counts.
#
# @rdname run_DESeq2
# @param counts Integer matrix of counting data. Must be RAW. Genes in rows,
#   samples in columns.
# @param conds Character vector of at most two different conditions for each
#   column of 'counts'.
# @param cpus Defaults to 1. Maximum number of cores (threads) to use.
# @seealso DESeq2 documentation. Available online at:
#   \url{https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf}.
# @return A DESeqResults object. See \code{\link{DESeq2::results}}.
.run_DESeq2 <- function
(
    counts,
    conds,
    cpus = 1
) {
    writeLines("DESeq2 ...")

    # construct data object
    cond1 <- unique(conds)[1]
    cond2 <- unique(conds)[2]
    colData <- data.frame(condition = factor(conds))
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData, formula( ~ condition))

    # collapse technical replicates
    # ...
    # subset of relevant columns of count data
    dds <- dds[, dds$condition == cond1 | dds$condition == cond2]
    # remove "unused" time-elements
    dds$condition <- droplevels(dds$condition)
    # "control" is the FIRST level in the treatment factor --> default log2FC are calculated as TREATMENT over CONTROL
    # NOTE: cond2 is usually "control"
    dds$condition <- relevel(dds$condition, cond2)

    # run differential expression pipeline
    if (cpus > 1) {
        if (cpus >= 4) cpus = 4
        res <- DESeq2::results(DESeq2::DESeq(dds,
                                             parallel = T,
                                             BPPARAM = BiocParallel::MulticoreParam(cpus)))
    } else {
        res <- DESeq2::results(DESeq2::DESeq(dds, parallel = F))
    }

    res$padj[is.na(res$padj)] <- 1
    return(res)
}


# contained in 'samr' package
# Run SAMseq
#
# Find significantly, differentially expressed genes using SAMseq.
#
# @rdname run_SAMseq
# @inheritParams DESeq2
# @param censoring.status Default to NULL. See \code{\link{SAMseq}}.
# @param resp.type Defaults to "Two class unpaired". See \code{\link{SAMseq}}.
# @param nperms Defaults to 100. Number of permutations used to estimate false
#   discovery rates. See \code{\link{SAMseq}}.
# @param nresamp Defaults to 20. Number of resamples used to construct test
#   statistic. See \code{\link{SAMseq}}.
# @param fdr.output Defaults to 1. FDR cut-off for output of significant genes.
#   See \code{\link{SAMseq}}.
# @seealso SAMseq manual. Available online at:
#   \url{https://statweb.stanford.edu/~tibs/SAM/sam.pdf}.
# @return A list.
.run_SAMseq <- function
(
    counts,
    conds,
    censoring.status = NULL,
    resp.type = "Two class unpaired",
    nperms = 100,
    nresamp = 20,
    fdr.output = 1
) {
    writeLines("SAMseq ...")
    # convert character vector of conds to a 1,2 vector
    class <- rep(0, length(conds))
    class[conds == unique(conds)[1]] <- 1
    class[conds == unique(conds)[2]] <- 2
    class <- as.numeric(class)

    SAMseq <- NULL
    SAMseq$test <- samr::SAMseq(
        x = counts,
        y = class,
        censoring.status = censoring.status,
        resp.type = resp.type,
        geneid = rownames(counts),
        genenames = rownames(counts),
        nperms = nperms,
        nresamp = nresamp,
        fdr.output = fdr.output
    )

    SAMseq$result.table <- rbind(
        SAMseq$test$siggenes.table$genes.up,
        SAMseq$test$siggenes.table$genes.lo
    )
    geneInds <- match(SAMseq$result.table[, 1], rownames(counts))
    SAMseq$score <- rep(0, nrow(counts))
    SAMseq$score[geneInds] <- as.numeric(SAMseq$result.table[, 3])
    SAMseq$FDR <- rep(1.0, nrow(counts))
    SAMseq$FDR[geneInds] <- as.numeric(SAMseq$result.table[, 5]) / 100

    return(SAMseq)
}


# Run PoissonSeq
#
# Find significantly, differentially expressed genes using PoissonSeq.
#
# @rdname run_PoissonSeq
# @inheritParams DESeq2
# @param SigFDR Signficance threshold.
# @seealso PoissonSeq documentation. Available online at: \url{https://cran.r-project.org/web/packages/PoissonSeq/PoissonSeq.pdf}.
# @return A list.
.run_PoissonSeq <- function
(
    counts,
    conds,
    SigFDR
) {
    writeLines("PoissonSeq ...")
    requireNamespace("PoissonSeq", quietly = TRUE)
    cond1 = unique(conds)[1]
    cond2 = unique(conds)[2]
    y = rep(0, length(conds))
    y[conds == cond1] = 1
    y[conds == cond2] = 2
    y = as.numeric(y)

    ## if there is only one sample in either condition than no logFC will be reported
    if (sum(conds == cond1) == 1 | sum(conds == cond2) == 1)
        condition.type <- 'quant'
    else
        condition.type <- 'twoclass'

    dat <- list(
        n = counts,
        y = y,
        type  = condition.type,
        pair  = FALSE,
        gname = as.character(rownames(counts))
    )
    #there is a number of addition arguments ?PS.Main
    para <- list(ct.sum = 5, npermu = 500)

    tmpres <- PoissonSeq::PS.Main(dat = dat, para = para)
    rownames(tmpres) <- tmpres$gname
    res <- {}
    res$deg <- tmpres$gname[tmpres$fdr < SigFDR]
    res$padj <- tmpres[rownames(counts), "fdr"]
    res$padj[is.na(res$padj)] <- 1 #sets p.value of untested genes to 1
    return(res)
}


# Run edgeR
#
# Find significantly, differentially expressed genes using edgeR.
#
# @rdname run_edgeR
# @inheritParams DESeq2
# @param remove.zeros Defaults to FALSE. If TRUE, remove all rows with counts
#    of zero.
# @seealso edgeR manual. Available online at
#   \url{https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf}.
# @return A list.
.run_edgeR <- function
(
    counts,
    conds,
    remove.zeros = FALSE
) {
    writeLines("edgeR...")
    dge <- edgeR::DGEList(
        counts = counts,
        group = conds,
        remove.zeros = remove.zeros
    )
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    dge <- edgeR::estimateCommonDisp(dge) #Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags.
    dge <- edgeR::estimateTagwiseDisp(dge, trend = "movingave") #Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood.
    de.common <- edgeR::exactTest(dge)
    res <- {}
    res$dge <- dge
    res$de.common <- de.common
    res$padj <- p.adjust(de.common$table$PValue, method = "BH") # Benjamini & Hochberg adjustment
    return(res)
}


# run baySeq
#
# Find significantly, differentially expressed genes using baySeq.
#
# @rdname run_baySeq
# @inherit DESeq2
# @param libsizes Default to colSums(counts). Number of reads per sample.
# @param groups Defaults to NULL. See \code{\link{baySeq}}.
# @param seglens See \code{\link{baySeq}}.
# @param samplesize Defaults to 5000. See \code{\link{baySeq}}.
# @param cpus Defaults to 1. Maximum number of cores (threads) to use.
# @seealso baySeq manual. Available online at \url{https://www.bioconductor.org/packages/devel/bioc/manuals/baySeq/man/baySeq.pdf}.
# @return A 'countData' object.
.run_baySeq <- function
(
    counts,
    conds,
    libsizes = colSums(counts),
    groups = NULL,
    seglens,
    samplesize = 5000,
    cpus = 1
) {
    writeLines("baySeq ...")

    if (cpus == 1) {
        cl <- NULL
    } else {
        cl <- parallel::makeCluster(cpus, "SOCK")
    }

    # convert character vector to integer vector
    cond1 <- unique(conds)[1]
    cond2 <- unique(conds)[2]
    replicates <- rep(0, length(conds))
    replicates[conds == cond1] <- 1
    replicates[conds == cond2] <- 2
    replicates <- as.numeric(replicates)

    # make default groups
    if (is.null(groups)) {
        groups <- list(
            NDE = rep(1, length(replicates)),
            DE  = c(rep(1, sum(conds == cond1)) , rep(2, sum(conds == cond2)))
        )
    }

    simCount <- as.matrix(counts)
    CD <- new(
        baySeq::.__C__countData,
        data = simCount,
        replicates = replicates,
        libsizes = as.integer(libsizes),
        groups = groups,
        seglens = seglens
    )

    CD@annotation <- data.frame(name = rownames(counts))
    # nba:
    CDP.NBML <- baySeq::getPriors.NB(
        CD,
        samplesize = samplesize,
        estimation = "QL",
        cl = cl
    )

    if (compareVersion(as.character(packageVersion("baySeq")), "2.0.0") >= 0)
        CDPost.NBML <- baySeq::getLikelihoods(CDP.NBML, pET = "BIC", cl = cl)
    else
        CDPost.NBML <- baySeq::getLikelihoods.NB(CDP.NBML, pET = "BIC", cl = cl)
    CDPost.NBML@estProps
    # tC <- topCounts(CDPost.NBML, group = 2, number = dim(CDPost.NBML)[1])
    # return(tC)

    if (cpus != 1) {
        parallel::stopCluster(cl = cl)
    }

    return(CDPost.NBML)
 }


# Run NOISeq
#
# Find significantly, differentially expressed genes using NOISeq.\cr
# Data-adaptive and nonparametric. Only biological replicates are supported at
# the moment.
#
# @rdname run_NOISeq
# @inherit DESeq2
# @param genelength Vector of gene lengths.
# @param k Default to 0.5. Counts equal to 0 are replaced by k.
# @param norm Defaults to "tmm". Normalization method. Choose 'n' if 'counts'
#   is already normalized. See ?\code{\link{oiseq}}.
# @param odds Defaults to 0.99. Proposed threshold for biological replicates.
#   Here, the probability of differential expression would be equivalent to 1 - FDR.
# @param lc Defaults to 1. Length correction is done by dividing expression by
#   length^lc. See \code{\link{noiseq}}.
# @seealso NOISeq manual. Available online at:
#   \url{https://www.bioconductor.org/packages/devel/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf}.
# @return A list.
.run_NOISeq <- function
(
    counts,
    conds,
    genelength,
    k = 0.5,
    norm = "tmm",
    odds = 0.99,
    lc = 1
) {
    writeLines("NOISeq ...")

    # HACK: genelength vector must be named
    if (is.null(names(genelength)))
        names(genelength) <- rownames(counts)
    factors <- data.frame(factor = factor(conds))
    noidata <- NOISeq::readData(
        data = counts,
        length = genelength,
        factors = factors
    )

    # NOTE: noires is a S4 class. Access using '@' operator instead of '$'.
    noires <- NOISeq::noiseqbio(
        noidata,
        k = k,
        norm = norm,
        factor = "factor",
        lc = lc
    ) # control_mean, treatment_mean, theta, prob, log2FC, length
    ## This code could be used for non-biological replicates, but the output has
    ## little informative value noires <- noiseq(noidata, k = k, norm = norm,
    ## replicates <- rep.type, factor = "factor", lc = lc)
    # control_mean, treatment_mean, M, D, prob, ranking, length selecting
    # differentially expressed genes
    res <- {}
    res$table <- noires@results[[1]]
    res$deg <- NOISeq::degenes(noires, q = odds)

    res$table$prob[is.na(res$table$prob)] <- 0.0

    return(res)
    ### A counting table.
}


# Run limma voom
#
# Find significantly, differentially expressed genes using limma voom.
#
# @rdname run_limmavoom
# @inherit DESeq2
# @note Voom is theoretically more powerful than limma-trend for RNA-seq data.
# @seealso limma manual. Available online at: \url{https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf}.
# @return A list.
.run_limmavoom <- function
(
    counts,
    conds
) {
    writeLines("limma voom ...")

    dge <- edgeR::DGEList(counts = counts, group = conds)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    voom.data <- limma::voom(dge, design = model.matrix( ~ factor(conds)))

    voom.data$genes <- rownames(counts)
    voom.fitlimma <- limma::lmFit(voom.data, design = model.matrix( ~ factor(conds))) # fit linear model
    voom.fitbayes <- limma::eBayes(voom.fitlimma) # compute t-statistics using empirical Bayes
    voom.pvalues <- voom.fitbayes$p.value[, 2]
    voom.adjpvalues <- p.adjust(voom.pvalues, method = "BH") # Benjamini & Hochberg

    result <- NULL
    result$data <- voom.data
    result$fitlimma <- voom.fitlimma
    result$fitbayes <- voom.fitbayes
    result$pvalues <- voom.pvalues
    result$adjpvalues <- voom.adjpvalues
    return(result)
}


# Pairwise Test Differential Expression
#
# Test for differential expression between a pair of conditions using one or
# more statistical tools.
#
# @rdname pairwise_testdiffexpr
# @param counts Matrix. Genes as rows, samples as columns.
# @param colinds Column indices.
# @param testconds Character vector.
# @param libsizes Library sizes. length(libsizes) == ncol(counts) must be TRUE.
# @param genelength Vector of gene lengths.
#   length(genelength) == nrow(counts) must be True.
# @param SigFDR Number. Defaults to 0.01. Adjusted p-value cut-off. Also known as
#   cut-off for significant false-discovery rate.
# @param logfcCut Numeric(1) > 0 or NA. Defaults to NA. Log2 fold changed
#   cut-off. If set, both pValCut and logfcCut must be satisfied to accept a
#   gene as hit. Additionally, boundaries will be marked in volcano plots.
#   This should usually be set to between 0.25 and 1.0.
# @param tools Character(n). By default, c("DESeq", "DESeq2", "edgeR", "limma").
#   The following tools are available (not case sensitive): \cr
#   DESeq2, limma, edgeR, SAMseq, NOISeq, baySeq, PoissonSeq
# @param odds For NOISeq only. Must be in [0..1]. Only odds above
#   this value are considered significant. Defaults to 1-FDR.
# @inheritParams parallelizationParameters
# @return A data.frame.
.pairwise_testdiffexpr <- function
(
    counts,
    colinds,
    testconds,
    libsizes,
    genelength,
    SigFDR = 0.01,
    logfcCut = NA,
    logfcnorm = "mrn",
    tools = "DESeq2",
    odds = 1 - SigFDR,
    cpus = 1
) {
    ## sanity checks
    if (!is.na(logfcCut) && (!is.numeric(logfcCut) || logfcCut < 0))
        stop("Log2 cut-off must be a non-negative number!")
    if (!is.numeric(SigFDR) || SigFDR <= 0)
        stop("SigFDR must be a number greater zero!")
    if (!(logfcnorm %in% c("mrn", "tpm", "rpkm")))
        stop("logfcnorm must be either 'mrn' or 'tpm'!")

    ## initiate data
    tools <- tolower(tools)
    nrreps <- length(colinds) / 2
    if (nrreps < 2)
        warning(
            "You have no replicates. These tools are therefore not available: limma, edgeR, NOISeq, SAMseq",
            immediate. = TRUE
        )

    stopifnot(length(genelength) == nrow(counts))

    # reorder columns
    testcounts <- counts[, colinds]
    testsizes  <- libsizes[colinds]
    testlength <- genelength
    # remove genes with zero counted reads
    testlength <- testlength[rowSums(testcounts) > 0]
    testcounts <- testcounts[rowSums(testcounts) > 0,]

    tpms <- get_tpm(testcounts, testlength, testsizes)
    colnames(tpms) <- colnames(testcounts)
    rpkms <- get_rpkm(testcounts, testlength, testsizes)
    colnames(tpms) <- colnames(testcounts)
    mrns <- get_mrn(testcounts)
    colnames(mrns) <- colnames(testcounts)

    # check, of there are multiple files available for the first and second condition
    rleConds <- rle(testconds)
    if (length(rleConds$values) > 2)
        stop("More than two conditions!")
    conds <- rleConds$values

    # if 2 or more samples per condition, determine row means
    if (rleConds$lengths[1] > 1) {
        meanAtpm <- rowMeans(tpms[, testconds == conds[1]]+1)
        meanArpkm <- rowMeans(rpkms[, testconds == conds[1]]+1)
        meanAmrn <- rowMeans(mrns[, testconds == conds[1]]+1)
    } else {
        meanAtpm <- tpms[, testconds == conds[1]]+1
        meanArpkm <- rpkms[, testconds == conds[1]]+1
        meanAmrn <- mrns[, testconds == conds[1]]+1
    }

    if (rleConds$lengths[2] > 1) {
        meanBtpm <- rowMeans(tpms[, testconds == conds[2]]+1)
        meanBrpkm <- rowMeans(rpkms[, testconds == conds[2]]+1)
        meanBmrn <- rowMeans(mrns[, testconds == conds[2]]+1)
    } else {
        meanBtpm <- tpms[, testconds == conds[2]]+1
        meanBrpkm <- rpkms[, testconds == conds[2]]+1
        meanBmrn <- mrns[, testconds == conds[2]]+1
    }

    log2_fc_tpm <- log2(meanAtpm / meanBtpm)
    log2_fc_tpm[is.nan(log2_fc_tpm)] <- 0
    fc_rpkm <- meanArpkm / meanBrpkm
    log2_fc_rpkm <- log2(meanArpkm / meanBrpkm)
    log2_fc_rpkm[is.nan(log2_fc_rpkm)] <- 0
    log2_fc_mrn <- log2(meanAmrn / meanBmrn)
    log2_fc_mrn[is.nan(log2_fc_mrn)] <- 0

    # 0/0   FC = 0      log2FC = 0
    # 0/x   FC = 0      log2FC = -Inf
    # x/0   FC = Inf    log2FC = Inf

    # table with count from different tools
    res_table <- data.frame(
        id = rownames(testcounts),
        mean_A_mrn = meanAmrn,
        mean_B_mrn = meanBmrn,
        log2_fc_mrn = log2_fc_mrn,
        mean_A_tpm = meanAtpm,
        mean_B_tpm = meanBtpm,
        log2_fc_tpm = log2_fc_tpm,
        mean_A_rpkm = meanArpkm,
        mean_B_rpkm = meanBrpkm,
        fc_rpkm = fc_rpkm,
        log2_fc_rpkm = log2_fc_rpkm
    )

    # results returned from tools

    if ("deseq" %in% tools) {
        if (nrreps > 1) {
            resdeseg <- .run_DESeq(
                testcounts,
                conds = testconds,
                method = "per-condition"
            ) # conditional pooling of replicates
        } else {
            resdeseg <- .run_DESeq(
                testcounts,
                conds = testconds,
                method = "blind",
                sharingMode = "fit-only"
            ) # Treat all samples as if they all have their own condition.
        }
        # mark significant genes
        sign_deseq <- resdeseg$id[!(is.na(resdeseg$padj)) & resdeseg$padj < SigFDR]
        sign_deseq <- is.element(res_table$id, sign_deseq)
        res_table <- cbind(res_table, DESeq = sign_deseq)
        res_table <- cbind(res_table, DESeq_adj_pval = resdeseg$padj)
    }

    if ("deseq2" %in% tools) {
        resdeseg2 <- .run_DESeq2(counts = testcounts, conds = testconds, cpus = cpus)

        # mark significant genes
        sign_deseq2 <- rownames(resdeseg2)[!(is.na(resdeseg2$padj)) &
                                               resdeseg2$padj < SigFDR]
        sign_deseq2 <- is.element(res_table$id, sign_deseq2)
        res_table <- cbind(res_table, DESeq2 = sign_deseq2) # named column
        res_table <- cbind(res_table, DESeq2_adj_pval = resdeseg2$padj)
    }

    if ("noiseq" %in% tools && nrreps >= 2) {
        res_noiseq <- .run_NOISeq(
            counts = testcounts,
            conds = testconds,
            genelength = testlength,
            k = 0.5,
            norm = "tmm",
            lc = 1,
            odds = odds
        )

        log2_noiseq <- res_noiseq$deg$log2FC
        signNOISeq <- is.element(res_table$id, rownames(res_noiseq$deg))
        res_table <- cbind(res_table, NOISeq = signNOISeq)
        res_table <- cbind(
            res_table,
            "NOISeq_1-odds" = (1 - res_noiseq$table$prob)
        ) # 1 - odds is better for sorting
    }

    # NOTE: FDRs differ for each run
    if ("samseq" %in% tools && nrreps >= 2) {
        res_samseq <- .run_SAMseq(testcounts, conds = testconds)

        SAMseq <- rownames(res_samseq$test$samr.obj$x)[res_samseq$FDR < SigFDR]
        signSAMseq <- is.element(res_table$id, SAMseq)
        res_table <- cbind(res_table, SAMseq = signSAMseq)
        res_table <- cbind(res_table, SAMseq_FDR = res_samseq$FDR)
    }

    if ("limma" %in% tools && nrreps >= 2) {
        res_limma <- .run_limmavoom(testcounts, conds = testconds)
        signLimma <- names(res_limma$adjpvalues)[res_limma$adjpvalues < SigFDR]
        signLimma <- is.element(res_table$id, signLimma)
        res_table <- cbind(res_table, Limma = signLimma)
        res_table <- cbind(res_table, Limma_adj_pval = res_limma$adjpvalues)
    }

    if ("edger" %in% tools && nrreps >= 2) {
        res_edgeR <- .run_edgeR(testcounts, conds = testconds, remove.zeros = FALSE)
        signEdgeR <- rownames(res_edgeR$de.common$table[res_edgeR$padj < SigFDR, ])
        signEdgeR <- is.element(res_table$id, signEdgeR)
        res_table <- cbind(res_table, EdgeR = signEdgeR)
        res_table <- cbind(res_table, EdgeR_adj_pval = res_edgeR$padj)
    }

    if ("bayseq" %in% tools) {
        res_bayseq <- .run_baySeq(
            counts = testcounts,
            conds  = testconds,
            libsizes = testsizes,
            seglens = as.double(testlength),
            samplesize = 5000,
            cpus = cpus
        )

        bayseq_select <- baySeq::topCounts(
            res_bayseq,
            group = 2,
            number = dim(res_bayseq)[1]
        )
        res_bayseq <- (as.character(
            bayseq_select$annotation[which(bayseq_select$FDR.DE < SigFDR)]
        ))
        signbaySeq <- is.element(res_table$id, res_bayseq)
        res_table <- cbind(res_table, baySeq = signbaySeq)
        res_table <- cbind(res_table, baySeq_FDR = bayseq_select$FDR.DE)
    }

    if ("poissonseq" %in% tools) {
        res_poissonseq <- .run_PoissonSeq(
            counts = testcounts,
            conds = testconds,
            SigFDR = SigFDR
        )
        padj_poissoinseq <- res_poissonseq$padj
        signPoissonSeq <- is.element(res_table$id, res_poissonseq$deg)
        res_table <- cbind(res_table, PoissonSeq = signPoissonSeq)
        res_table <- cbind(res_table, PoissonSeq_adj_pval = padj_poissoinseq)
    }

    if (!is.na(logfcCut)) {
        # apply log2 fold change threshold
        if (logfcnorm == "mrn")
            fold_hits <- abs(res_table$log2_fc_mrn) >= logfcCut
        else # must be "tpm" now
            fold_hits <- abs(res_table$log2_fc_tpm) >= logfcCut
        tCol <- tolower(colnames(res_table)) %in% tolower(tools)
        tCol <- which(tCol == TRUE)
        res_table[, tCol] <- res_table[, tCol] & fold_hits
    }

    return(res_table)
}



#' Calculate Differentially Expressed Genes
#'
#' Calculate differentially expressed genes (DEGs) by comparing two conditions.
#' Create a DEG table and various plots to visualize the results.
#'
#' Each column of the design matrix represents an independent test
#' of two conditions (control vs. treatment). All other conditions are
#' ignored for this test. This also means that all samples belonging to "control"
#' will be grouped together as biological replicates. If you use technical
#' replicates, you should merge the corresponding rows by summing them up.
#' \cr\cr
#' RPKM, TPM, fold changes, etc. are calculated based on the 'counts' matrix.
#' \cr\cr
#' Each DEG test tool calculates its own fold changes and uses its own
#' normalization methods. Therefore, their results are hard (or impossible) to
#' compare.
#' \cr\cr
#' For each tool defined by 'tools', 2 columns will be created in the output DEG
#' table. The first column is always a boolean vector indicating whether a
#' gene's expression is significantly (based on the specific tool) regulated.
#' The second column are FDR, p-value or other statistical parameters, which
#' were used for this decision.
#'
#' @export
#' @param counts Matrix. Genes in rows, samples in columns.
#' @param geneLengths Integer vector. Lengths of genes. Make sure that the
#'   order is correct, i.e. rownames(geneLengths) == rownames(counts).
#' @param libSizes Library sizes. length(libsizes) must be equal to ncol(counts).
#' @param designMatrix Matrix describing all comparisons, i.e. tests for
#'   differential expression, to perform. One column is one comparison, which
#'   will result in one DEG table. For each column, two conditions are expected:
#'   control and treatment. The remaining rows, i.e. samples, are ignored for
#'   this comparison.
#' @param pValCut Number. Defaults to 0.01. Maximum FDR or adjusted p-value at
#'    which to reject the null hypothesis (i.e. significance threshold).
#' @param logfcCut Numeric(1) > 0 or NA. Defaults to NA. Log2 fold change
#'   cut-off. If set, both pValCut and logfcCut must be satisfied to accept a
#'   gene as DEG. Additionally, boundaries will be marked in volcano plots.
#'   This should usually be set to between 0.25 and 1.0.
#' @param logfcnorm Either "mrn", "tpm" or "rpkm". MRN stands for median ratio
#'   normalization as implemented in the EBSeq package. It is very similar to the
#'   VST transformation implemented in DESeq2. MRN normalization is more rebust
#'   against library size effects.
#' @param tools Character(n). By default, c("DESeq", "DESeq2", "edgeR", "limma").
#'   The following tools are available additionally (not case sensitive):
#'   SAMseq, NOISeq, baySeq, PoissonSeq.
#' @param outDir Path to output directory. Defaults to "stats".
#' @param prefix Prefix string to add to all output files. Defaults to the empty
#'   string "".
#' @param no.xls If TRUE, do not write results as xls files.
#' @param no.csv If TRUE, do not write results as csv files.
#' @param anno Optional. Annotation data.frame to merge into resulting
#'   differential gene expression data.frame(s). Merging is done based on the
#'   first columns.
#' @param all.x Only has an effect if 'anno' is defined. Logical indicating
#'   whether all rows of the annotation data.frame should be added, even if they
#'   could not be matched against the DEG data.frame. See ?\code{merge} for more
#'   information. Defaults to FALSE.
#' @param stop.on.error Logical indicating whether the whole execution should
#'   stop if an error is encountered for any comparison. By default (FALSE), the
#'   resulting comparison is returned as an error object in the result list, but
#'   the calculation of subsequent comparisons continues.
#' @inheritParams parallelizationParameters
#' @return
#' \describe{
#'   \item{DEG}{
#'     A list of data.frames. Each list entry corresponds to one column of the
#'     design matrix. More precisely, it corresponds to one specific test of two
#'     conditions among all 'tools'. For each tool using FDR or p-values, 'pValCut'
#'     will be used to decide which genes are significantly, differentially
#'     expressed (as indicated by the boolean vector of each tool).
#'   }
#'   \item{tool}{The name of the tool.}
#' }
#' @examples
#' calculate_DEGs
calculate_DEGs <- function
(
    counts,
    geneLengths,
    libSizes,
    designMatrix,
    pValCut = 0.01,
    logfcCut = NA,
    logfcnorm = "mrn",
    tools = c("DESeq", "DESeq2", "edgeR", "limma"),
    outDir = "stats",
    prefix = "",
    no.xls = FALSE,
    no.csv = FALSE,
    anno = "",
    all.x = FALSE,
    stop.on.error = FALSE,
    cpus = 1,
    workers = NA
) {
    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
    # correct tool names
    available <- c(
        "DESeq",
        "DESeq2",
        "Limma",
        "EdgeR",
        "baySeq",
        "NOISeq",
        "SAMseq",
        "PoissonSeq"
    )
    found <- tolower(available) %in% tolower(tools)
    tools <- available[found]

    if (!(logfcnorm %in% c("mrn", "tpm", "rpkm")))
        stop("Invalid logfcnorm. Choose either 'mrn', 'tpm' or 'rpkm'")

    # sanity check of design matrix
    .check_design_matrix(designMatrix, counts)

    # check for empty and low-count samples
    countSizes <- colSums(counts)
    if (sum(countSizes == 0) > 0) {
        zero_lib <- which(countSizes == 0)
        warning(
            "calculate_DEGs: One or more libraries have zero counts and are excluded: ",
            paste(zero_lib, collapse = ",")
        )
        if (length(zero_lib) == ncol(counts))
            stop("Add libraries has zero counts! Abort.")
        counts <- counts[, which(countSizes > 0)]
        designMatrix <- designMatrix[which(countSizes > 0), ]
    }
    if (sum(countSizes < 1000) > 0) {
        warning(
            "calculate_DEGs: One or more libraries have less than 1000 counts: ",
            paste(which(countSizes < 1000), collapse = ",")
        )
    }

    calcDEGs <- function(i) {
        writeLines("")
        writeLines(paste(
            "Working on comparison",
            colnames(designMatrix)[i],
            paste0("(", i, "/", ncol(designMatrix), ")")
        ))

        # find genes which are top ranked for differentially expression
        # correct with adj. correction
        DEGs <- .pairwise_testdiffexpr(
            counts    = counts,
            colinds   = reorder_idx(colnames(counts), c(
                rownames(designMatrix)[designMatrix[, i] == "treatment"],
                rownames(designMatrix)[designMatrix[, i] == "control"]
            )),
            testconds = c(
                rep("treatment", sum(designMatrix[, i] == "treatment")),
                rep("control" , sum(designMatrix[, i] == "control"))
            ),
            libsizes   = libSizes,
            genelength = geneLengths,
            SigFDR     = pValCut,
            logfcCut   = logfcCut,
            logfcnorm  = logfcnorm,
            tools      = tools,
            cpus       = allocCPUS$per_call
        )

        return(DEGs)
    }

    writeOut <- function(i){
        #     logfc = DEGs$restable$log2FC
        #     names(logfc) = DEGs$DeSeq$id
        DEGs <- allDEGs[[i]]
        used_tools = intersect(tools, colnames(DEGs))
        tools_with_hits = colSums(DEGs[tools]) > 0
        if (sum(tools_with_hits) < length(used_tools)) {
            if (length(used_tools[!tools_with_hits]) == 1) {
                message(used_tools[!tools_with_hits], " has no hits and is excluded from graphs.")
            } else {
                message(
                    paste(used_tools[!tools_with_hits], sep = ", "),
                    " have no hits and are excluded from graphs."
                )
            }
        }
        tools_with_hits <- used_tools[tools_with_hits]

        comparison <- colnames(designMatrix)[i]

        # NOTE: MA plot does not seem to work independent of stable normalization.
        #       But only DESeq2 supplies an appropriate function ~ rlog or vst.
        norm <- toupper(logfcnorm)
        patA <- paste0("mean_A_", logfcnorm)
        patB <- paste0("mean_B_", logfcnorm)
        mA <- grep(patA, colnames(DEGs))
        mB <- grep(patB, colnames(DEGs))
        stopifnot(length(mA) == 1 && length(mB) == 1)
        make_MA_plot(
            file.path(outDir, paste0(prefix, comparison, "_MA")),
            counts = as.matrix(DEGs[,c(mA, mB)]),
            title = paste("MA Plot", "-", comparison, "-", norm)
        )

        if (length(tools_with_hits) > 0) {
            if (length(tools_with_hits) >= 1 && length(tools_with_hits) <= 4) {
                venn <- make_Venn_diagramms(
                    file = file.path(outDir, paste0(prefix, comparison, "_venn")),
                    degTable = DEGs,
                    tools = tools_with_hits
                )
            }

            # the function will throw a warning for less than 2 tools
            if (length(tools_with_hits) > 1) {
                inter <- make_Intersection_Barplots(
                    file = file.path(outDir, paste0(prefix, comparison, "_bars")),
                    degTable = DEGs,
                    tools = tools_with_hits
                )
            }

            fCol = grep(paste0("log2_fc_", logfcnorm), colnames(DEGs))
            stopifnot(length(fCol) == 1)
            make_volcano_plot(
                file = file.path(outDir, paste0(prefix, comparison, "_volc")),
                degTable = DEGs,
                tools = tools_with_hits,
                title = paste("Vulcano Plot -", comparison),
                FDR = pValCut,
                fCol = fCol,
                log2fold = logfcCut
            )

            DEGs <- sortDiffExp(DEGs, tools = tools_with_hits)
        } else
            warning(paste("No hit by any tool for comparison: ", comparison))

        # add annotation information, if requested
        if (anno != "" && !is.null(anno) && !is.na(anno)) {
            DEGs <- merge(
                DEGs,
                anno,
                by.x = 1,
                by.y = 1,
                all.x = all.x,
                all.y = FALSE
            )
        }

        if (!no.csv) {
            csv_out <- file.path(outDir, paste0(prefix, comparison, ".csv"))
            write.table(
                DEGs,
                file = csv_out,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE
            )
            writeLines(paste("Results at:", csv_out))
        }

        if (!no.xls) {
            xls_out <- file.path(outDir, paste0(prefix, comparison, ".xls"))
            WriteXLS::WriteXLS(
                DEGs,
                ExcelFileName = xls_out,
                SheetNames = comparison,
                AdjWidth = TRUE,
                BoldHeaderRow = TRUE,
                FreezeRow = 1,
                FreezeCol = 1,
                row.names = FALSE
            )
            writeLines(paste("Results at:", xls_out))
        }

        return(TRUE)
    }

    allocCPUS <- allocate_cpus(cpus, ncol(designMatrix), workers)
    print(allocCPUS)
    bio_par <- BiocParallel::MulticoreParam(
        workers = allocCPUS$in_parallel,
        tasks = ncol(designMatrix),
        progressbar = TRUE,
        stop.on.error = stop.on.error
    )
    allDEGs <- tryCatch({
        BiocParallel::bplapply(1:ncol(designMatrix), calcDEGs, BPPARAM = bio_par)
    }, error = identity)
    .stop_at_error(allDEGs, from = match.call()[[1]], warn = !stop.on.error)
    names(allDEGs) <- colnames(designMatrix)

    writeLines("Start plotting and writing out ...")
    bio_par <- BiocParallel::MulticoreParam(
        workers = 1,
        tasks = 1,
        stop.on.error = stop.on.error
    )
    allRes <- tryCatch({
        BiocParallel::bplapply(1:length(allDEGs), writeOut, BPPARAM = bio_par)
    }, error = identity)
    .stop_at_error(allRes, from = match.call()[[1]], warn = !stop.on.error)

    if (is.na(logfcCut)) {
        desc <- paste0("adjusted p-value cut-off = ", pValCut)
    } else {
        desc <- paste0("adjusted p-value cut-off = ", pValCut, "; ", "log2 fold-change cut-off = ", logfcCut)
    }

    return(list(
        DEGs  = allDEGs,
        desc  = desc,
        tools = tools,
        tool  = "DEG"
    ))
}
