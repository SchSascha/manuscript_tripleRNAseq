####-
##-====================================================
##-========================  VISUALIZATION ========================
##-====================================================
####-


#' Make Hierarchical Clustering Plot
#'
#' Given a matrix of counts (rows) per sample (columns), calculate the
#' column-wise distances and perform hierarchical clustering. The plot
#' is saved as PDF file.
#'
#' @export
#' @param file Prefix path. Appropriate file extensions are added automatically.
#' @param counts Numeric matrix. Samples in columns, genes in rows.
#' @param main Main title for the plot. Defaults to the empty string "".
#' @param conds Character(n). If supplied, the plot will be colored based on
#'   the conditions of each sample.\cr
#'   Assumption: conds[i] for colnames(counts)[i]
#' @param wScale Width scale of image. width = wScale * |samples|. Increase for
#'   a bigger plot. Defaults to 0.1875.
#' @param bottom Bottom margin of PDF file. Defaults to 7.
#' @param method Clustering method. UPGMA by default. See ?\code{hclust} for
#'    details.
#' @inheritDotParams graphics::plot
#' @return An object of type 'dist'. Basically the distance matrix calculated
#'   from 'counts' and used for clustering.
#' @examples
#' # get some data to cluster
#' counts <- matrix(replicate(20, rnorm(20)), ncol = 5)
#' # column names must be unique!
#' colnames(counts) <- c("S1.1", "S1.2", "S3.1", "S3.2", "S3.3")
#' conds <- c("A", "A", "B", "B", "B")
#' # make plot
#' # make_hclust_plot("test", counts, "Test title", conds = conds)
make_hclust_plot <- function
(
    file,
    counts,
    main   = "",
    conds  = NA,
    wScale = 0.1875,
    bottom = 7,
    method = "ave",
    ...
){
    samples <- colnames(counts)
    dists = dist(t(counts))
    hc <- hclust(dists, method)
    dend <- as.dendrogram(hc)
    width <- max(7, as.integer(ncol(counts) * wScale))
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

    # ensure valid file ending
    if (!grepl("\\.pdf$", file))
        file <- paste0(file, ".pdf")

    pdf(file = file, width = width)
    par(mar = c(bottom, 5, 5, 3))

    if (!is.na(conds) && "dendextend" %in% installed.packages()) {
        if (is.list(conds) || !is.na(conds)) {
            uniq_conds <- unique(conds)
            cols <- rainbow(length(uniq_conds))
            idx <- lapply(1:length(uniq_conds), function(i) which(uniq_conds[i] == conds))
            col_for_label <- vector(mode = "character", ncol(counts))
            for (i in 1:length(idx))
                col_for_label[idx[[i]]] <- cols[i]
            # order of leaves change for dendrogram
            order <- reorder_idx(colnames(counts), labels(dend))
            dendextend::labels_colors(dend) <- col_for_label[order]
            dendextend::color_branches(dend)
            par(cex=0.6)
            plot(dend, ...)
            par(cex=1)
            title(main = main)
            mtext("Samples")
            legend("topright", legend = uniq_conds, fill = cols)
            dev.off()
            return(dists)
        }
    } else if (mget("__warn_hclust__", ifnotfound = TRUE, envir = .GlobalEnv) == TRUE){
        # only show this warnings once
        assign("__warn_hclust__", FALSE, envir = environment())
        warning("Colorless plot. For more colors in hclust, please install package 'dendextend'.")
    }
}





#' Make Heat Hierarchical Clustering
#'
#' Given a matrix of counts (rows) per sample (columns), calculate the
#' column-wise distances and perform hierarchical clustering. Additionally,
#' pairwise distances are presented as heat map. The plot is safed as PDF file.
#'
#' @export
#' @param overwrite Logical indicating whether to overwrite existing output
#'   files. Default to FALSE.
#' @inheritParams make_hclust_plot
#' @inheritParams .make_DESeq2_transform
#' @inheritDotParams pheatmap::pheatmap
#' @return An object of type 'dist'. Basically the distance matrix calculated
#'   from 'counts' and used for clustering.
#' @examples
#' # get some data to cluster
#' dat <- as.integer(abs(replicate(20, rnorm(20))) * 1000)
#' counts <- matrix(dat, ncol = 5)
#' # column names must be unique!
#' colnames(counts) <- c("S1.1", "S1.2", "S3.1", "S3.2", "S3.3")
#' conds <- c("A", "A", "B", "B", "B")
#' # make plot
#' make_heat_clustering_plot("heat", counts, main = "Test title", conds = conds)
make_heat_clustering_plot <- function
(
    file,
    counts,
    designMatrix = NA,
    conds = NA,
    norm = "auto",
    geneLen = NA,
    overwrite = FALSE,
    main = "",
    ...
){
    # ensure valid file ending
    if (!grepl("\\.pdf$", file))
        file <- paste0(file, ".pdf")

    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    if (!overwrite && file.exists(file))
        stop("Cannot make heat clustering: output file already exists! ", shQuote(file))

    dists <- .distances_by_DESeq2(
        counts = counts,
        designMatrix = designMatrix,
        conds = conds,
        norm = norm,
        geneLen = geneLen,
        err_msg = "Cannot make heat clustering:"
    )

    # convert and name distance matrix
    distMat <- as.matrix(dists)
    rownames(distMat) <- attr(dists, "Labels")
    colnames(distMat) <- NULL
    colors <- grDevices::colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)

    # write out plot
    pdf(file, onefile = FALSE)
    pheatmap::pheatmap(
        mat = distMat,
        clustering_distance_rows = dists,
        clustering_distance_cols = dists,
        col = colors,
        main = main,
        ...
    )
    dev.off()
    return(dists)
}



#' Make MA Plot
#'
#' An MA plot is a scatter plot created by plotting
#' intensity ratio (of A vs. B) by their average intensity.
#'
#' Either 'counts' or 'degTable' must be set to a valid expression object. Log2
#' transformation is applied to count data. Plots are mostly created by
#' functions of the affy package.
#'
#' @export
#' @param counts A two column matrix of expression values, 'A vs B'. If supplied,
#'   'degTable' is ignored.
#' @param degTable Single list element of result of \code{\link{calculate_DEGs}}.
#'   A data.frame or named matrix.
#' @param title Main title for the plot.
#' @param dotSize Size of dots in the plot. Defaults to 0.25.
#' @inheritParams make_hclust_plot
#' @inheritDotParams affy::ma.plot
#' @references Gautier L, Cope L, Bolstad BM and Irizarry RA (2004).
#'   "affy-analysis of Affymetrix GeneChip data at the probe level."
#'   Bioinformatics, 20(3), pp. 307-315. ISSN 1367-4803, doi:
#'   10.1093/bioinformatics/btg405.
#' @return Path to output PDF file.
#' @examples
#' # get some data to cluster
#' dat <- as.integer(abs(replicate(10, rnorm(10))) * 1000)
#' counts <- matrix(dat, ncol = 2)
#' colnames(counts) <- c("S1", "S2")
#' # make plot
#' make_MA_plot("MA", counts = counts)
#'
#' \dontrun{
#' deg_res <- calculate_DEGs("with/valid/input")
#' make_MA_plot("MA", degTable = deg_res$DEGs[[1]])
#' }
make_MA_plot <- function
(
    file,
    counts   = NA,
    degTable = NA,
    title    = "",
    dotSize  = 1,
    ...
){
    if (is.null(dim(counts))) {
        if (is.null(dim(degTable)))
            stop("Cannot make MA plots: arguments 'counts' and 'degTable' are both NA! See ?make_MA_plot")
        if (!is.numeric(degTable[,2]) || !is.numeric(degTable[,3]))
            stop("Cannot make MA plots with degTable: count data must be numeric!")
        R <- log2(degTable[,2])
        G <- log2(degTable[,3])
    } else {
        if (ncol(counts) != 2)
            stop(paste("Cannot make MA plots: counts must have exactly two columns! But it has:", ncol(counts)))
        if (!is.numeric(counts))
            stop("Cannot make MA plots: counts must be numeric! See ?make_MA_plot")
        R <- log2(counts[,1])
        G <- log2(counts[,2])
    }
    M <-        R - G
    A <- 0.5 * (R + G)

    # find NA
    m_na <- is.na(M)
    a_na <- is.na(A)
    # find Inf
    m_inf <- is.infinite(M)
    a_inf <- is.infinite(A)
    # remove Inf and NA rows
    keep <- !(m_na | a_na | m_inf | a_inf)
    M <- M[keep]
    A <- A[keep]

    # start plotting
    dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
    file <- paste0(file, ".pdf")

    pdf(file)
    affy::ma.plot(
        A    = A,
        M    = M,
        cex  = dotSize,
        # xlab = "average intensity (A)", # microarray
        # ylab = "intensity ratio (M)",
        xlab = "Average Log2 Count (A)",       # RNA-seq
        ylab = "Log2 Count Ratio (M)",
        # main = if (title != "") title else "Intensity Ratio by Average Intensities (MA)", # microarray
        # main = if (title != "") title else "Count Ratio by Average Counts (MA)",          # RNA-seq
        main = if (title != "") title else "MA Plot",
        plot.method = "smoothScatter",
        ...
    )
#     mtext(sub)
    dev.off()
    return(file)
}


#' Make Volcano Plot
#'
#' Volcano plots are used to detect interesting genes. They also provide a brief
#' overview for the (differentially expressed) genes of two conditions. Adjusted
#' p-values are plotted against log2 fold changes. The plot should look similar
#' to a volcanic eruption.
#'
#' Already a single gene with a very low p-value, compared to the other genes in
#' 'degTable', will result in a volcano plot, were most of the data points are
#' concentrated next to the x-axis and almost indistinguishable. Scaling the
#' plot via the 'scale' argument may prevent this mostly unwanted behavior.
#'
#' @export
#' @param tools Vector with names of statistical tools. Usually the ones
#'   used for \code{\link{calculate_DEGs}}. If NA, this will be set to all the
#'   tools found in 'degTable'. Tools are found based on indices defined by
#'   'pCols'. \strong{If set, this will overwrite pCols}.
#' @param gCol Column index (1) for gene names. Defaults to 1.
#' @param pCols Column indices for p-values. If NA, the function will
#'   use any column containing the substring "adj_pval". E.g., a column named
#'   "DESeq2_adj_pval" would be used.
#' @param fCol Column index for log2 fold changes. If NA, the function will
#'   use a column named "log2_fc". For other data.frames, supply the index to
#'   the corresponding column instead.
#' @param interactive If TRUE, a HTML file with an interactive volcano plot
#'   will be created additionally. Defaults to TRUE.
#' @param infMode Character(1). Either 'finite' or 'ignore'. If 'finite', all
#'   values of Inf or -Inf will be set to max(log-fold)+1 and min(log-fold)-1,
#'   respectively. If 'ignore', Inf and -Inf values will be out of bounds and
#'   not shown. Defaults to 'finite'.
#' @param scale Percentage of lowest p-values to scale down. Between 0 and 1.
#'   If NA, no scaling is performed. Otherwise, p-values are sorted and
#'   'scale'-percent of the top values are set to a fixed value. This value is
#'   determined by using the lowest p-value of the remaining ones plus some
#'   offset. The resulting points will be located above a horizontal, colored
#'   line. Defaults to NA. See details.
#' @param FDR FDR cut-off. Values below this threshold are considered
#'   significant. Defaults to 0.01.
#' @param log2fold Numeric(1) >= 0. Log2 fold change threshold. Increasing this
#'   value removes statistical hits but improves the FDR for remaining genes.
#'   Removed hits are shown in red. If 0, this value is ignored. Defaults to 0.
#' @inheritParams make_MA_plot
#' @inheritDotParams graphics::plot
#' @return Path to output file(s). If interactive is TRUE, the second list
#' element will be the path to the resulting interactive HTML file.
#' @examples
#' \dontrun{
#' deg_res <- calculate_DEGs("with/valid/input")
#' make_volcano_plot("volc", deg_res$DEGs[[1]])
#' }
#' make_volcano_plot
make_volcano_plot <- function
(
    file,
    degTable,
    tools = NA,
    title = "",
    gCol = 1,
    pCols = NA,
    fCol = NA,
    interactive = TRUE,
    infMode = "finite",
    scale = NA,
    dotSize = 0.25,
    FDR = 0.01,
    log2fold = 0,
    ...
){
    # sanity check
    if (!(is.data.frame(degTable) || is.matrix(degTable)) || is.na(colnames(degTable)))
        stop("volcano plot: 'degTable' must be a matrix or data.frame and contain column names!")
    if (is.na(fCol) || is.null(fCol))
        fCol <- which(colnames(degTable) == "log2_fc")

    if ((is.list(tools) ||is.vector(tools)) && !is.na(tools) && !is.null(tools)) {
        inds <- tolower(colnames(degTable)) %in% tolower(tools)
        inds <- sapply(1:length(inds), function(i) if (inds[i] == TRUE) i else 0)
        inds <- inds[which(inds > 0)]
        tools <- colnames(degTable)[inds]
        pCols <- inds + 1 # HACK: p_values are always in the next column
    } else {
        if (is.na(pCols) || is.null(pCols))
            pCols <- grep("adj_pval|1-odds", names(degTable))
        tools <- gsub("adj_pval|1-odds|_", "", names(degTable)[pCols])
    }

    # subtitle string
    sub_title <- paste(tools, collapse = ", ")
    p_adj_sep <- degTable[pCols]
    p_adj_sep[which(is.na(p_adj_sep), arr.ind=TRUE)] <- 1
    # get row-wise maximums. This is about 100 times fast than apply(..., max())
    p_adj_comb <- do.call(function(...) pmax(..., na.rm=TRUE), p_adj_sep)
    sign_df <- as.data.frame(p_adj_sep < FDR)
    num_sign_df <- rowSums(sign_df)
    if (max(num_sign_df) < length(tools))
        message("Intersection of all tools is zero.")
    labels <- c(
        "Insignificant",
        "All Tools No Fold",
        "Only Fold Change",
        "Not all Tools",
        "All Tools & FC"
    )
    deg_df <- data.frame(degTable[gCol], degTable[fCol], p_adj_comb)
    names(deg_df) <- c("Genes", "Fold Change", "FDR")
    deg_df["group"] <- labels[1]
    abs_folds <- abs(unlist(deg_df[["Fold Change"]], use.names = FALSE))

    if (is.na(log2fold))
        log2fold <- 0
    else if (log2fold < 0)
        stop("volcano: Cannot make plot with negative, absolute log2 fold change value!")

    if (log2fold > 0) {
        # significant, but log2 fold change too low - suspicious
        sigAlone <- p_adj_comb < FDR & abs_folds < log2fold
        deg_df[sigAlone, "group"] <- labels[2]

        # large enough log2 fold change but p too high
        foldAlone <- p_adj_comb > FDR & abs_folds > log2fold
        deg_df[foldAlone, "group"] <- labels[3]
    }

    # large enough log2 fold change but not all tools hit
    sigOne <- num_sign_df > 0 & num_sign_df < ncol(sign_df) & abs_folds > log2fold
    deg_df[sigOne, "group"] <- labels[4]

    # sufficient log2 fold change and p-value
    sigAll <- p_adj_comb < FDR & abs_folds > log2fold
    deg_df[sigAll, "group"] <- labels[5]

    # completely uninteresting
    if (log2fold)
        sigNone <- !(sigAlone | foldAlone | sigOne | sigAll)
    else
        sigNone <- !(sigOne | sigAll)

    # make static plot
    # adjust limits and infinity outlier
    neg_inf <- deg_df[["Fold Change"]] == -Inf
    pos_inf <- deg_df[["Fold Change"]] == Inf
    mi_fold <- min(deg_df[["Fold Change"]][!neg_inf], Inf, na.rm = TRUE)
    ma_fold <- max(deg_df[["Fold Change"]][!pos_inf], -Inf, na.rm = TRUE)
    # extreme case: only infinity fold changes
    if (is.infinite(abs(mi_fold))) {
        warning("Minimum fold change of all genes is infinite.")
        mi_fold <- -1
    }
    if (is.infinite(abs(ma_fold))) {
        warning("Minimum fold change of all genes is infinite.")
        mi_fold <- 1
    }

    if (infMode == "finite") {
        deg_df[["Fold Change"]][neg_inf] <- mi_fold - 1
        deg_df[["Fold Change"]][pos_inf] <- ma_fold + 1
        xlim <- c(mi_fold - 1, ma_fold + 1)
    } else if (infMode == "ignore") {
        xlim <- c(mi_fold, ma_fold)
    } else
        stop("Unknown mode for infMode: ", infMode)

    colors <- c("black", "red", "darkorange", "blue", "green")

    # prevent infinity log(p-values) & rescale plot
    if (is.na(scale)) {
        p_thresh <- NA
        p_adj_comb[p_adj_comb < 1e-300] <- 1e-300
    } else {
        if (scale <= 0 || scale >= 1) {
            warning("Scale argument for volcano plots must be greater zero and lower one. Setting it to 0.01")
            scale <- 1e-2
        }
        p_thresh <- sort(p_adj_comb)[max(1, as.integer(scale * length(p_adj_comb)))]
        p_thresh <- max(p_thresh, 1e-290)
        off_P <- 0.98 # percentage of Y-axis used for P offset
        off_set <- p_thresh ** (1/off_P) - p_thresh
        p_adj_comb[p_adj_comb < p_thresh] <- p_thresh + off_set
    }
    p_adj_comb <- -log10(p_adj_comb)

    dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
    pdf_file <- paste0(file, "_", paste(tools, collapse = "_"), ".pdf")
    pdf(pdf_file)
    plot(
        x = deg_df[["Fold Change"]],
        y = p_adj_comb,
        xaxt = "n",
        type = "n", # consider values, but do not plot points (yet)
        pch = 20,
        cex = dotSize,
        main = if (title != "") title else "Volcano Plot",
        xlim = xlim,
        ylim = c(0, max(p_adj_comb) * 1.01),
        xlab = "Log2 Fold Change",
        ylab = "-Log10( Adj. P-Value )",
        col = colors[1],
        ...
    )
    mtext(sub_title)

    # make better axis labels and colors
    x_ticks <- seq(from = mi_fold, to = ma_fold)
    x_ticks <- as.integer(x_ticks)
    axis(1, at = x_ticks, labels = as.character(x_ticks), las=2)
    if (infMode == "finite") {
        axis(
            side = 1,
            at = c(mi_fold - 1, ma_fold + 1),
            labels = c("-Inf", "+Inf"),
            col.ticks = "orange",
            col.axis = "orange",
            las=2
        )
    }

    # axis for extreme p-values
    if (!is.na(p_thresh))
        axis(
            side = 4, at = -log10(p_thresh + off_set),
            labels = paste0("> ", round(-log10(p_thresh + 0.5 * off_set), digits = 2)),
            col = "magenta",
            col.axis = "magenta"
        )

    # Add infinity borders
    if (infMode == "finite")
        abline(v = c(mi_fold-0.5, ma_fold+0.5), col = "orange", lty = 1)
    if (!is.na(p_thresh))
        abline(h = -log10(p_thresh + 0.5 * off_set), col = "magenta", lty = 1)

    # add dots
    points(
        deg_df[["Fold Change"]][sigNone],
        p_adj_comb[sigNone],
        pch = 20,
        cex = dotSize,
        col = colors[1]
    )
    if (log2fold > 0) {
        points(
            deg_df[["Fold Change"]][sigAlone],
            p_adj_comb[sigAlone],
            pch = 20,
            cex = dotSize,
            col = colors[2]
        )
        points(
            deg_df[["Fold Change"]][foldAlone],
            p_adj_comb[foldAlone],
            pch = 20,
            cex = dotSize,
            col = colors[3]
        )
    }
    points(
        deg_df[["Fold Change"]][sigOne],
        p_adj_comb[sigOne],
        pch = 20,
        cex = dotSize,
        col = colors[4]
    )
    points(
        deg_df[["Fold Change"]][sigAll],
        p_adj_comb[sigAll],
        pch = 20,
        cex = dotSize,
        col = colors[5]
    )

    # Add informative lines
    abline(v = c(log2fold, -log2fold), col = "black", lty = 2)
    abline(h = -log10(FDR), col = "blue", lty = 2)

    legend(x = "topright", legend = labels, col = colors, pch = 20)
    dev.off()

    html_file <- NULL
    if (interactive) {
        # The result will be a functional HTML
        if ("htmlwidgets" %in% installed.packages() && "plotly" %in% installed.packages()){
            if (log2fold > 0) {
                colors <- c("red", "green", "black", "darkorange", "blue")
            } else {
                colors <- c("green", "black", "blue")
            }

            p <- plotly::plot_ly()
            p <- plotly::add_trace(
                p,
                data = deg_df,
                x = deg_df[["Fold Change"]],
                y = p_adj_comb,
                text = deg_df$Gene,
                mode = "markers",
                type = "scatter",
                colors = colors,
                color = deg_df$group
            )
            p <- plotly::layout(
                p,
                title = if (title != "") paste(title, sub_title, sep="\n") else paste0("volcano Plot\n", sub_title),
                xaxis = list(title = "Log2-FC", showgrid = TRUE),
                yaxis = list(title = "-log10(Adj. P-Value)")
            )
            if (infMode == "finite") {
                dx <- c(mi_fold-0.5, ma_fold+0.5)
                dy <- c(0, max(p_adj_comb))
                p <- plotly::add_trace(
                    p,
                    x = dx[1],
                    y = dy,
                    type = "scatter",
                    mode = "lines",
                    name = "infinity"
                )
                p <- plotly::add_trace(
                    p,
                    x = dx[2],
                    y = dy,
                    type = "scatter",
                    mode = "lines",
                     name = "infinity"
                 )
            }
            #abline(v = c(mi_fold-0.5, ma_fold+0.5), col = "orange", lty = 1)
            #abline(v = c(log2fold, -log2fold), col = "black", lty = 2)
            #abline(h = -log10(FDR), col = "blue", lty = 2)
            outDir <- dirname(file)
            html_file <- paste0(
                basename(file),
                "_",
                paste(tools, collapse = "_"),
                ".html"
            )
            htmlwidgets::saveWidget(
                plotly::as_widget(p),
                html_file,
                selfcontained = FALSE,
                background = "gray"
            )
            if (outDir != ".") {
                file.rename(html_file, file.path(outDir, html_file))
                htmlDir <- paste0(tools::file_path_sans_ext(html_file), "_files")
                file.copy(htmlDir, outDir, recursive = TRUE)
                unlink(htmlDir, recursive = TRUE, force = TRUE)
            }
        } else if(mget("__warn_volcano__", ifnotfound = TRUE, envir = .GlobalEnv) == TRUE) {
            # only show this warnings once in a session
            assign("__warn_volcano__", FALSE, envir = environment())
            warning("For interactive volcano plots, please install packages: 'plotly' & 'htmlwidgets'")
        }
    }

    return(c(pdf_file, html_file))
}


#' Make Intersection Bar Plots
#'
#' Create an intersection bar plot based on the intersections of gene lists, e.g.
#' differentially expressed genes. Especially for more than 4 lists (sets) of
#' genes, they are a much clearer alternative to Venn diagrams.
#'
#' @export
#' @inheritParams make_volcano_plot
#' @param degTable Single list element of result of \code{\link{calculate_DEGs}}.
#'   If NA, a valid 'DEGlist' object must be given. Defaults to NA.
#' @param tools Vector with names of statistical tools. Usually the ones
#'   used for \code{\link{calculate_DEGs}}. If NA, this will be set to all the
#'   tools detected in 'degTable' or 'DEGlist'. Defaults to NA.
#' @param DEGlist A names list. A vector of genes per tool. Can be created
#'   using \code{\link{make_DEG_list}}. If NA, 'degTable' must be given.
#'   Defaults to NA.
#' @param showEmpty Logical(1). If TRUE, empty sets will be shown.
#'   Defaults to FALSE.
#' @param nintersections Integer(1). Only the top 'nsets' sets will be shown.
#'   If NA, all sets will be shown. Defaults to 20.
#' @param colorless If FALSE, the plot will only be black/white. Defaults to
#'   FALSE.
#' @param format Output file format. PDF, SVG and PNG are supported. Defaults to
#'   PDF.
#' @param as.log Logical (1) indicating whether to use log2 transformation.
#'   Defaults to FALSE.
#' @param x.label Custom label for x-axis. Defaults to "Set Size".
#' @param y.label Custom label for y-axis. Defaults to "Number of Genes".
#' @param height Image height. Only for PNG format. Defaults to 2000.
#' @param width Image width. Only for PNG format. Defaults to 3000.
#' @param dpi Image resolution. Only for PNG. Defaults to 300.
#' @return A named list with keywords:\cr
#' \describe{
#'   \item{DEGlist}{Named list of DEG genes per tool.}
#'   \item{file   }{Path to output file.}
#'   \item{UpSet  }{A data frame. The result of the conversion of degTable
#'                  using \code{UpSetR::fromList.}}
#' }
#' @references Jake R Conway, Alexander Lex, Nils Gehlenborg UpSetR: An R
#'   Package for the Visualization of Intersecting Sets and their Properties
#'   doi: https://doi.org/10.1093/bioinformatics/btx364
#' @examples
#' deg_list <- list(tool1 = c("A", "B", "C", "D"),
#'                  tool2 = c("A", "B", "C", "D", "F"),
#'                  tool3 = c("A", "B", "C", "E"))
#' make_Intersection_Barplots("inter1", DEGlist = deg_list)
#'
#' \dontrun{
#' deg_res <- calculate_DEGs("with/valid/input")
#' make_Intersection_Barplots("inter2", deg_res$DEGs[[1]])
#' }
make_Intersection_Barplots <- function
(
    file,
    degTable = NA,
    tools = NA,
    DEGlist = NA,
    showEmpty = FALSE,
    nintersections = 20,
    #   TODO PACKAGE: good default value?
    colorless = FALSE,
    format = "pdf",
    as.log = FALSE,
    x.label = "Set Size",
    y.label = "Number of Genes",
    height = 2000,
    width = 3000,
    dpi = 300
){
    DEGlist <- .get_DEGlist(degTable, tools, DEGlist)
    if (length(names(DEGlist)) < 2) {
        warning("Intersection Bar Plots for less than 2 tools are not supported.")
        return()
    }
    file <- paste0(file, "_", paste(names(DEGlist), collapse = "_"))

    upSetFrame <- UpSetR::fromList(DEGlist)

    if (format == "pdf") {
        file <- paste0(file, ".pdf")
        pdf(file, onefile = FALSE, width = 11, height = 7)
        # onefile = FALSE is a workaround for an empty first page in PDF format
        # https://github.com/hms-dbmi/UpSetR/issues/90
    } else if (format == "svg") {
        file <- paste0(file, ".svg")
        svg(file, width = 11, height = 7)
    } else if (format == "png"){
        file <- paste0(file, ".png")
        png(file, width = width, height = height, res = dpi)
        # TODO png goes not work in R320 on server Uranos
        # make_Intersection_Barplots(..., format = "png") gives:
        # ...
        # Error in grid.Call(L_stringMetric, as.graphicsAnnot(x$label)) :
        #    X11 Schrift -adobe-helvetica-%s-%s-*-*-%d-*-*-*-*-*-*-*, Typ 1 in Groesse 12 konnte nicht geladen werden
    } else
        stop(paste("unknown output file format:", format))

    UpSetR::upset(
        upSetFrame,
        order.by="freq",
        nsets = length(names(DEGlist)),
        nintersects = nintersections,
        scale.sets = if (as.log) "log2" else "identity",
        text.scale = c(2, 1.5, 2, 1.5, 1.25, 2),
        number.angles = 30,
        main.bar.color = if (colorless) "gray" else "blue"
    )
    dev.off()

    return(list(DEGlist = DEGlist, file = file, upset = upSetFrame))
}



#' Make Venn Diagrams For Differentially Expressed Genes
#'
#' Create a Venn diagram based on the intersections of gene lists, e.g.
#' differentially expressed genes (DEGs).
#'
#' @export
#' @inheritParams make_Intersection_Barplots
#' @param height Image height. Defaults to 3000.
#' @param width Image width. Defaults to 3000.
#' @param dpi Image resolution. Defaults to 300.
#' @return A list with keywords:\cr
#' \describe{
#'   \item{venn   }{Venn diagram object.}
#'   \item{DEGlist}{A named list with gene names per tools.}
#' }
#' @examples
#' deg_list <- list(tool1 = c("A", "B", "C", "D"),
#'                  tool2 = c("A", "B", "C", "D", "F"),
#'                  tool3 = c("A", "B", "C", "E"))
#' make_Venn_diagramms("venn1", DEGlist = deg_list)
#'
#' \dontrun{
#' deg_res <- calculate_DEGs("with/valid/input")
#' make_Venn_diagramms("inter2", deg_res$DEGs[[1]])
#' }
make_Venn_diagramms <- function
(
    file,
    degTable = NA,
    tools = NA,
    DEGlist = NA,
    height = 3000,
    width = 3000,
    dpi = 300
){
    # fix any issue regarding capital letters
    DEGlist <- .get_DEGlist(degTable, tools, DEGlist)
    tools <- names(DEGlist)

    # numGenes <- sapply(1:length(tools), function(i)
    # length(degTable$id[degTable[, tools[i]]]))
    numGenes <- sapply(1:length(DEGlist), function(i) length(DEGlist[[i]]))
    if (0 %in% numGenes)
        warning(
            paste(tools[numGenes %in% 0], "has no hits and skipped for plotting!\n"),
            immediate. = TRUE
        )

    tools <- tools[numGenes > 0]
    if (length(tools) > 1) {
        if (length(tools) > 4) {
            # NOTE: For more than 4 sets, split them into groups of 4.
            for (i in 1:(length(tools) - 3)) {
                .make_venn_diffexpr(
                    degTable,
                    tools[i:(i + 3)],
                    file = paste0(file, "_", paste(tools[i:(i + 3)], collapse = "_"), ".pdf"),
                    height = height,
                    width = width,
                    dpi = dpi
                )
            }
        } else {
            .make_venn_diffexpr(
                degTable = degTable,
                tools = tools,
                DEGlist = DEGlist,
                file = paste0(file, "_", paste(tools, collapse = "_"), ".pdf"),
                height = height,
                width = width,
                dpi = dpi
            )
        }
    } # else: just one set
}


# Helper Function for make_venn_diagrams
#
# @rdname make_venn_diffexpr
# @inherit make_Venn_diagramms
.make_venn_diffexpr <- function
(
    file,
    degTable = NA,
    tools = NA,
    DEGlist = NA,
    height = 3000,
    width = 3000,
    dpi = 500
){
    # For each tool, look at their top hits (boolean vector).
    # And then grab the corresponding gene names.
    DEGlist <- .get_DEGlist(degTable, tools, DEGlist)

    # make colorful venn.
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") # no log files
    vp <- VennDiagram::venn.diagram(
        x = DEGlist[tools],
        filename = NULL,
        main = "Venn Diagram",
        main.fontface = "bold",
        sub = paste0("total number of DEGs: ", length(unique(unlist(DEGlist)))),
        fill = rainbow(length(DEGlist[tools])),
        # alpha = 0.3,
        height = height,
        width = width,
        resolution = dpi,
        print.mode = c("raw", "percent")
    )
    pdf(file)
    grid::grid.draw(vp)
    dev.off()

    return(list(venn = vp, DEGlist = DEGlist[tools]))
}


#' Make Overlap Matrix
#'
#' Calculate a binary matrix representing overlapping sets. The input is a list
#' element of the DEG table returned by the function \code{\link{calculate_DEGs}}.
#'
#' @export
#' @inheritParams make_Intersection_Barplots
#' @return #' A matrix with genes as rows and one column per tool. If a value
#'   of 1 for column j ind row i, the gene i was a significant hit for tool j.
#' @examples
#' deg_list <- list(tool1 = c("A", "B", "C", "D"),
#'                  tool2 = c("A", "B", "C", "D", "F"),
#'                  tool3 = c("A", "B", "C", "E"))
#' make_overlap_matrix(DEGlist = deg_list)
#'
#' \dontrun{
#' deg_res <- calculate_DEGs("with/valid/input")
#' make_overlap_matrix(deg_res$DEGs[[1]])
#' }
make_overlap_matrix <- function
(
    degTable = NA,
    tools = NA,
    DEGlist = NA
){
    DEGlist <- .get_DEGlist(degTable, tools, DEGlist)
    UpSetR::fromList(DEGlist)
}


#' Make Correlation Plots
#'
#' Calculate and visualize the pairwise euclidean distances between all samples
#' in a count matrix.
#'
#' If a design matrix is supplied, two plots are created per comparison, i.e.
#' only pairwise distances within the same condition are calculated and plotted.
#' Instead of individual points, the density is drawn, i.e. the overall location
#' and concentration of points.
#'
#' @export
#' @param dat Matrix or data.frame of count data. Genes in rows, samples in
#' columns.
#' @param outDir Path to directory folder. Defaults to the working directory.
#' @param prefix prefix string that is attached to every output file.
#' @param designMatrix Optional. Design matrix with samples in rows and
#'   comparisons in columns. Defaults to NA.
#' @param overwrite Logical indicating whether to overwrite existing output
#'   files. Defaults to FALSE.
#' @return Character(1). Path to output file.
#' @examples
#' # get some data to cluster
#' dat <- as.integer(abs(replicate(20, rnorm(20))) * 1000)
#' counts <- matrix(dat, ncol = 5)
#' # column names must be unique!
#' colnames(counts) <- c("S1.1", "S1.2", "S3.1", "S3.2", "S3.3")
#' # make plot
#' make_correlation_plots(counts, prefix = "cor")
make_correlation_plots <- function
(
    dat,
    outDir = getwd(),
    prefix = "",
    designMatrix = NA,
    overwrite = FALSE
){
    # log scale is better - and fix values of infinite by adding 1
    l_dat <- as.data.frame(log10(dat+1))
    # create each pairwise plot within one condition
    # TODO: sometimes, special signs are ignored which leads to naming errors. Handle such errors here!
    remove_signs  <- "[][~!$%(){}+\"?,']"
    replace_signs <- "[~@#$%^&*+:<>,/;=-]"
    colnames(l_dat) <- gsub(remove_signs , "" , colnames(l_dat))
    colnames(l_dat) <- gsub(replace_signs, "_", colnames(l_dat))
    labels <- colnames(l_dat)

    # select columns to correlate pairwise
    if (!is.matrix(designMatrix) || is.na(designMatrix)) {
        cols <- 1:ncol(l_dat)
    } else {
        message("Creating correlation plots based on Design Matrix...")
        res <- list()
        # START OF RECURSION: for each comparison and treatment and control...
        for (c in 1:ncol(designMatrix)) {
            message("Processing: ", colnames(designMatrix)[c])
            # select all columns used for the current test
            treatCols <- reorder_idx(
                colnames(l_dat),
                rownames(designMatrix)[designMatrix[, c] == "treatment"]
            )
            contCols  <- reorder_idx(
                colnames(l_dat),
                rownames(designMatrix)[designMatrix[, c] == "control"]
            )
            res[[2*c-1]] <- make_correlation_plots(
                dat = dat[,treatCols],
                outDir = outDir,
                prefix = paste0(prefix, colnames(designMatrix)[c], "_treatment"),
                overwrite = overwrite
            )
            res[[2*c]] <- make_correlation_plots(
                dat = dat[,contCols],
                outDir = outDir,
                prefix = paste0(prefix, colnames(designMatrix)[c], "_control"),
                overwrite = overwrite
            )
        }

        return(unlist(res))
    }

    # prepare and check output path
    outFile <- file.path(outDir, paste0(prefix, ".pdf"))
    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
    if (file.exists(outFile) && !overwrite)
        stop(paste("Output file already exists:", outFile))

    spl <- list()
    k <- 1
    font_size <- max(1, 18.0 - 2*length(cols))
    for (j in cols) {
        for (i in cols) {
            if (i == j) {
                pear <- 1.0
                sp <- ggplot2::ggplot(
                        l_dat,
                        ggplot2::aes_string(
                            x = paste0("`", labels[[i]], "`"),
                            y = paste0("`", labels[[j]], "`")
                        )
                    ) +
                    ggplot2::geom_blank() +
                    ggplot2::geom_abline(
                        intercept = 0,
                        slope = 1,
                        color = "red",
                        size = 0.25
                    ) +
                    ggplot2::theme_bw() +
                    ggplot2::scale_x_continuous(position = "top") +
                    ggplot2::labs(
                        title = paste("Pearson:", round(pear, 3)),
                        x = if (j==1) labels[[i]] else "",
                        y = if (i==1) labels[[j]] else ""
                    ) +
                    ggplot2::theme(
                        text         = ggplot2::element_text(size = font_size),
                        plot.title   = ggplot2::element_text(size = font_size),
                        axis.title.x = ggplot2::element_text(size = font_size),
                        axis.title.y = ggplot2::element_text(size = font_size)
                    )
            } else {
                pear <- cor(l_dat[,i], l_dat[,j], method="pearson")
                # make plot
                sp <- ggplot2::ggplot(
                        l_dat,
                        ggplot2::aes_string(
                            x = paste0("`", labels[[i]], "`"),
                            y = paste0("`", labels[[j]], "`")
                        )
                    ) +
                    ggplot2::labs(
                        title = paste("Pearson:", round(pear, 3)),
                        x = if (j==1) labels[[i]] else "",
                        y = if (i==1) labels[[j]] else ""
                    ) +
                    ggplot2::geom_abline(
                        intercept = 0,
                        slope = 1,
                        color = "red",
                        size = 0.25
                    ) +
                    # NOTE: print all dots
                    # geom_point(size=0.25) +
                    # NOTE: density plot saves a lot of time for opening the resulting PDF
                    # stat_density_2d(aes(fill = ..level..), geom="polygon") +
                    # NOTE: heat scale ~ to much for each plot
                    # scale_fill_gradient(low="blue", high="red") +
                    ggplot2::geom_density_2d() +
                    ggplot2::geom_smooth(color = "blue", size = 0.25, se=FALSE) +
                    ggplot2::theme_bw() +
                    # NOTE: x-axis on top of plots
                    ggplot2::scale_x_continuous(position = "top") +
                    ggplot2::theme(
                        text         = ggplot2::element_text(size = font_size),
                        plot.title   = ggplot2::element_text(size = font_size),
                        axis.title.x = ggplot2::element_text(size = font_size),
                        axis.title.y = ggplot2::element_text(size = font_size)
                    )
            }

            spl[[k]] <- sp
            k <- k + 1
        }
    }

    # write multiple scatter plots into a single plot
    message("Saving plot to: ", shQuote(outFile))
    ggplot2::ggsave(
        outFile,
        gridExtra::arrangeGrob(
            grobs  = spl,
            ncol   = length(cols),
            bottom = "Pearson >= 0.95 is good; Pearson >= 0.90 is still ok; Pearson < 0.90 is not ok"
        ),
        width  = min(30, 3+length(cols)*2),
        height = min(30, 3+length(cols)*2))
    return(outFile)
}



####-
##-====================================================
##-=============================  PCA ==============================
##-====================================================
####-


#' Make DESeq2 transform
#'
#' FOR INTERNAL USE ONLY.
#' Given a matrix of count data with samples in columns and genes as rows,
#' apply 'rlog' or 'vst' logarithm transformations using the DESeq2 package.
#'
#' @param err_msg Message to show on error, e.g. when this function is called
#'   from another function.
#'
#' @param counts Numeric matrix with m genes (rows) and n columns
#'   (samples). DESeq2 \code{\link{rlog}} or \code{\link{vst}} transformation is
#'   applied based on the number of samples. rlog for less than 50 samples, vst
#'   otherwise.
#' @param designMatrix Character matrix with n rows (samples) and k columns
#'   (comparisons). If NA, 'conds' must be supplied. Else,
#'   'conds' is ignored.
#' @param conds Character(n). This vector defines the group of each sample.
#'   Samples of the same group are assigned the same color. If NA,
#'   'designMatrix' must be supplied. For technical reasons, do not use
#'   "control" as condition label.
#' @param norm Normalization method. For 'auto', 'rlog' or 'vst' is chosen based
#'   on the number of samples. These can only be used on integer matrices. For
#'   'rpkm' or 'tpm', either RPKM or TPM normalization is applied. This requires
#'   'geneLen'.
#'   If norm = 'cpm', samples are normalized for library size and
#'   multiplied by 10^6.
#'   If norm = 'mrn', samples are scaled by their geometric mean (median). This
#'   method is favored when vst or rlog are not applicable. It is used for
#'   volcano plots by default. Internally, log2(mrn(x)+1) values are calculated
#'   to make results comparable with other log-based methods like vst or rlog.
#'   If norm = 'none', an unnormalized DESeqDataSet object
#'   is returned.
#' @param geneLen Vector(m) of gene lengths. Only required for TPM or RPKM
#'   normalization.
#'
#' @return A DESeqDataSet object with DESeq2::vst transform (ncol(counts) > 50)
#'   or DESeq2::rlog transform (otherwise).
#' @examples
#' \dontrun{
#' counts <- matrix(1, ncol = 2, nrow = 10)
#' colnames(counts) <- c("S1", "S2")
#' rownames(counts) <- paste0("gene_", 1:10)
#' conds <- c(rep("A", 1), rep("B", 1))
#' .make_DESeq2_transform(counts = counts, conds = conds)
#' }
.make_DESeq2_transform <- function
(
    counts,
    designMatrix = NA,
    conds = NA,
    norm = "auto",
    geneLen = NA,
    err_msg = ".make_DESeq2_transform:"
){
    # check data integrity
    if (is.na(counts) || is.null(counts) || !(is.matrix(counts) || is.data.frame(counts)))
        stop(
            err_msg,
            " Cannot plot PCA: 'counts' must be either, a matrix or data.frame!"
        )
    if (!is.matrix(designMatrix) && !is.data.frame(designMatrix) && (is.na(designMatrix) || is.null(designMatrix)) && (is.na(conds) || is.null(conds)))
        stop(
            err_msg,
            " 'designMatrix' and conds are both NA or NULL. One of both must be supplied!"
        )

    # define condition vector, i.e. the condition of each sample
    if (is.matrix(designMatrix) || is.data.frame(designMatrix))
        conds <- conditions_from_design(designMatrix)

    if ("control" %in% conds) {
        warning("'control' cannot be used as condition. It is renamed to 'contr' instead.")
        conds[conds == "control"] <- "contr"
    }

    colData <- S4Vectors::DataFrame(condition = factor(conds))
    rownames(colData) <- colnames(counts)
    # remove genes with zero counts in all samples
    inds <- rowSums(counts) > 0
    counts <- counts[inds,]

    # DESeq2 only accepts integer matrices!
    if (norm == "auto" || norm == "rlog" || norm == "vst")
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData, formula(~ condition))

    # RPKM and TPM require the geneLen attribute
    if ( !is.atomic(geneLen) || !is.na(geneLen))
        geneLen <- geneLen[inds]
    if ( (norm == "rpkm" || norm == "tpm") && is.atomic(geneLen) && (is.na(geneLen) || is.null(geneLen) ))
        stop("Cannot use RPKM or TPM normalization on counts: function attribute 'geneLen' is NA!")

    # NOTE: rlog and vst compensate gene counts of zero
    if (norm == "auto" && ncol(dds) < 50 || norm == "rlog") {
        # NOTE: 'blind=FALSE' means experimental design is not used
        message("Using rlog transformation.")
        return(DESeq2::rlog(dds, blind = FALSE))
    } else if (norm == "auto" || norm == "vst") {
        message("Using VST transformation.")
        # NOTE: vst() is a fast wrapper for varianceStabilizingTransformation(),
        # by subsetting to a smaller number of genes (nsub = 1000), but only
        # available since version 1.12.0
        if(packageVersion("DESeq2") >= 1.12) {
            return(DESeq2::vst(dds, blind = FALSE))
        } else {
            return(DESeq2::varianceStabilizingTransformation(dds, blind = FALSE))
        }
    } else if (norm == "rpkm") {
        rpkm <- get_rpkm(
            counts,
            geneLen,
            colSums(counts)
        )
        se <- SummarizedExperiment::SummarizedExperiment(rpkm, colData=colData)
        return(DESeq2::DESeqTransform(se))
    } else if (norm == "tpm") {
        tpm <- get_tpm(
            counts,
            geneLen,
            colSums(counts)
        )
        se <- SummarizedExperiment::SummarizedExperiment(tpm, colData=colData)
        return(DESeq2::DESeqTransform(se))
    } else if (norm == "cpm") {
        cpm <- t(t(counts) / colSums(counts))*10^6
        se <- SummarizedExperiment::SummarizedExperiment(cpm, colData=colData)
        return(DESeq2::DESeqTransform(se))
    } else if (norm == "mrn") {
        mrn <- log2(get_mrn(counts) + 1)
        se <- SummarizedExperiment::SummarizedExperiment(mrn, colData=colData)
        return(DESeq2::DESeqTransform(se))
    } else if (norm == "none") {
        se <- SummarizedExperiment::SummarizedExperiment(counts, colData=colData)
        return(DESeq2::DESeqTransform(se))
    } else
        stop("Wrong value for 'norm'!")
}


# Distances by DESeq2
# FOR INTERNAL USE ONLY.
#
# Given a matrix of count data with samples in columns and genes as rows,
# calculate the pairwise distances between samples.\cr
# 'rlog' or 'vst' logarithm transformations are performed using the DESeq2
# package.
#
# @rdname distances_by_DESeq2
# @inheritParams .make_DESeq2_transform
# @return A numeric matrix (pairwise sample distances).
.distances_by_DESeq2 <- function
(
    counts,
    designMatrix,
    conds,
    norm,
    geneLen,
    err_msg = ".distances_by_DESeq2"
){
    tdds <- .make_DESeq2_transform(
      counts = counts,
      designMatrix = designMatrix,
      conds = conds,
      norm = norm,
      geneLen = geneLen,
      err_msg = err_msg
    )
    return(dist(t(SummarizedExperiment::assay(tdds))))
}



#' Make Principle Component Analysis Plot
#'
#' Perform a Principle Component Analysis (PCA) based on count data.
#'
#' Conditions are extracted from a design matrix or can be supplied using the
#' character vector 'conds'. It must be of the same length as the number of used
#' samples in counts.
#' \cr\cr
#' Samples are grouped based on their conditions. Each group has a unique color.
#' A secondary grouping can be defined by using 'shapes'. Each "secondary" group
#' has a unique shape than.
#' \cr\cr
#' The PCA plot is build using DESeq2 DataSets and \code{ggplot2}.
#'
#' @export
#' @inheritParams .make_DESeq2_transform
#' @param file Character. Output file name.
#' @param shapes Character(n). Optional. Secondary group information. If
#'   supplied, data point shapes are changed based on this group information in
#'   addition to the colors.
#' @param main Main title of the plot.
#' @param colName The title for groups as defined by colors
#'   (shown in the legend).
#' @param shapeName The title for groups as defined by shapred
#'   (shown in legend).
#' @param dotSize Size of output dots (data points). Defaults to 3.
#' @param custom_colors Character(n). Optional. This vector overwrite the
#'   colors chosen by ggplot. If set, it must have the same size as the
#'   conditions.
#' @param add_eclipse Logical. If TRUE, for each condition, the center of
#'   replicates will be surrounded by an eclipse describing the variance.
#'   Defaults to TRUE.
#' @param overwrite Logical indicating whether to overwrite existing output
#'   files. Defaults to FALSE.
#' @return A list with keywords:\cr
#' \describe{
#'   \item{tdds  }{DESeqDataSet object with rlog or vst transformation as used for PCA.}
#'   \item{pca   }{PCA data object}
#'   \item{biplot}{PCA data object as used for plotting the biplot}
#'   \item{g     }{ggplot object}
#' }
#' @examples
#' require("airway")
#' data("airway")
#' counts <- assay(airway)
#' # trim data to reduce runtime
#' counts <- counts[1:500,]
#' shapes <- as.character(colData(airway)$cell)
#' conds  <- as.character(colData(airway)$dex)
#' pca_file <- "PCA.pdf"
#' pca_res <- make_PCA_plot(file = pca_file,
#'                          counts = counts,
#'                          conds = conds,
#'                          norm = "auto",
#'                          shapes = shapes,
#'                          overwrite = TRUE)
make_PCA_plot <- function
(
    file,
    counts,
    designMatrix = NA,
    conds = NA,
    shapes = NULL,
    norm = "auto",
    main = "PCA biplot",
    colName = "groups",
    shapeName = "shapes",
    geneLen = NA,
    dotSize = 3,
    custom_colors = NA,
    add_eclipse = TRUE,
    overwrite = FALSE
){
    # display plot or safe plot to file?
    printToFile <- TRUE
    if (is.na(file) || is.null(file)) printToFile <- FALSE

    if (printToFile) {
        dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
        # ensure valid file ending
        if (!grepl("\\.pdf$", file))
            file <- paste0(file, ".pdf")
        if (!overwrite && file.exists(file))
            stop(
                "Cannot create PCA plot: output file already exists! ",
                shQuote(file)
            )
    }

    # apply rlog or vst transformation
    tdds <- .make_DESeq2_transform(
        counts = counts,
        designMatrix = designMatrix,
        conds = conds,
        norm = norm,
        geneLen = geneLen,
        err_msg = "Cannot plot PCA:"
    )

    # calculate PCA data (adapted from DESeq2 source code)

    # calculate the variance for each gene
    rv <- genefilter::rowVars(SummarizedExperiment::assay(tdds))

    # select the 500 genes by variance
    ntop <- 500 # DESeq2 default value
    select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

    # perform a PCA on the data in assay(x) for the selected genes
    # no scaling because tdds is already transformed
    pca <- prcomp(t(SummarizedExperiment::assay(tdds)[select,]))

    # calculate eigenvalues and variances
    eigenvalues <- (pca$sdev) ^ 2
    variance <- eigenvalues * 100 / sum(eigenvalues)
    # variance_cumulative <- cumsum(variance)

    # create PCA data for biplot
    biplot <- DESeq2::plotPCA(tdds, intgroup = "condition", returnData = TRUE)
    # percentVar <- round(100 * attr(biplot, "percentVar"))

    # create biplot
    if (length(shapes) < 2 || is.na(shapes) || is.null(shapes)) {
        aes <- ggplot2::aes(
            x = biplot$PC1,
            y = biplot$PC2,
            color = biplot$group
        )
    } else {
        aes <- ggplot2::aes(
            x = biplot$PC1,
            y = biplot$PC2,
            color = biplot$group,
            shape = factor(shapes)
        )
    }
    g <- ggplot2::ggplot(biplot, aes) +
        ggplot2::labs(
            title = main,
            subtitle = paste(round(sum(variance[1:2]), digits = 1), "%", "cumulative variance"),
            x = paste0("PC1: ", round(variance[1], digits = 1), "% variance"),
            y = paste0("PC2: ", round(variance[2], digits = 1), "% variance")
        ) +
        ggplot2::geom_point(size = dotSize) +
        ggplot2::theme(aspect.ratio=1) +
        ggplot2::scale_colour_discrete(name = colName) +
        ggplot2::scale_shape_discrete(name = shapeName)

    # add variance eclipses, if demanded
    if (add_eclipse)
        g <- g + ggplot2::stat_ellipse(ggplot2::aes(group = group))

    # change colors to specified values
    if (length(custom_colors) > 1 || !is.na(custom_colors)) {
        g <- ggplot2::scale_color_manual(custom_colors)
    }

    # display plot to current device
    if (printToFile) pdf(file)

    # plot ggplot biplot
    plot(g)

    # plot PCA summary table
    ncomp <- ncol(pca$x)
    if(ncomp > 10) ncomp <- 10 # print not more than 10 components
    plot.new()
    gridExtra::grid.table(round(
        # omit first row, which is standard deviation
        # standard deviation is usally not as important as variance for PCA
        # but values can become quite large and use too much horizontal space
        summary(pca)$importance[2:3, 1:ncomp],
        # round to two digits, again to save horizontal space
        digits = 2
    ))

    # plot contribution to variance for each component
    barplot(
        variance,
        names.arg = 1:length(variance),
        main = paste(
            "Variances:",
            "A good dimension reduction is achieved when the first few PCs",
            "account for a large proportion of the variability (80-90%)",
            sep = "\n"
        ),
        xlab = "Principal Components",
        ylab = "Percentage of variances"
    )

    # plot eigenvalues for each component
    barplot(
        eigenvalues,
        names.arg = 1:length(eigenvalues),
        main = paste(
            "Eigenvalues:",
            "An eigenvalue > 1 indicates that PCs account for more variance",
            "than accounted by one of the original variables",
            sep = "\n"
        ),
        xlab = "Principal Components",
        ylab = "Eigenvalue"
    )

    if (printToFile) dev.off()

    return(list(tdds = tdds, pca = pca, biplot = biplot, g = g))
}
