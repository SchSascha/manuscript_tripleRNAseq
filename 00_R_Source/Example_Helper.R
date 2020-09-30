#' Make Example Count Matrix
#'
#' Return a simple integer matrix, which is similar to a matrix returned by
#'   \code{featureCounts}.
#'
#' @export
#' @return Integer(m x n).
#' @examples makeExampleCountsMatrix()
makeExampleCountsMatrix <- function() {
    utils::data("counts_min", package = "Geo2RNAseq", envir = environment())
    return(get("counts_min", envir = environment()))
}


#' Make Example Design Matrix
#'
#' Return a simple design matrix, as used GEO2RNAseq. Samples in rows, tests
#' (comparisons) in columns. Each column represents a test, which will be
#' performed in the differential gene expression (DEG) analysis.
#'
#' Entries in the design matrix are either "treatment", "control" or "none". If
#' possible, column names should have the form "<A>_VS_<B>", where <A> and
#' <B> are variable parts. <A> represents the treatment, <B> represents the
#' control. If defined consistently, the condition of each sample can be
#' extracted automatically.
#'
#' @export
#' @return Character(m x n).
#' @examples makeExampleDesignMatrix()
makeExampleDesignMatrix <- function() {
    design_matrix <- matrix("none", ncol = 3, nrow = 6)
    colnames(design_matrix) <- c("DEX_VS_noDEX", "DEX_VS_2DEX", "2DEX_VS_no_DEX")
    rownames(design_matrix) <- paste0(
        "Sample_", c("DEX.1", "DEX.2", "Con.1", "2DEX.1", "Con.2", "2DEX.2")
    )
    design_matrix[c(1,2), 1] <- "treatment"
    design_matrix[c(3,5), 1] <- "control"
    design_matrix[c(1,2), 2] <- "treatment"
    design_matrix[c(4,6), 2] <- "control"
    design_matrix[c(4,6), 3] <- "treatment"
    design_matrix[c(3,5), 3] <- "control"
    return(design_matrix)
}


#' Make Example DEG Table
#'
#' Return an example data.frame with differentially expressed genes (DEGs) as
#' returned by \code{\link{calculate_DEGs}}.
#'
#' @export
#' @return A data.frame with specific columns.
#' @examples makeExampleDEGTable()
makeExampleDEGTable <- function() {
    utils::data("deg_min", package = "Geo2RNAseq", envir = environment())
    return(get("deg_min", envir = environment()))
}


#' Get getGEOdata Example Return Object
#'
#' Returns the return value of \code{\link{getGEOdata}} for GSE52778. Used in
#' vignette.
#'
#' @export
#' @param extDir Directory with external package data. Usually the subdirectory
#'    "extdata" inside the installation directory of the package.
#' @return Named list.
# @examples getGeoDemoDat
#' @examples
#' pkgDir <- system.file("extdata", package = "Geo2RNAseq")
#' getGeoDemoDat(pkgDir)
getGeoDemoDat <- function(extDir) {
    utils::data("geo_dat", package = "Geo2RNAseq", envir = environment())
    gdd <- get("geo_dat", envir = environment())
    gdd$METAfile <- file.path(extDir, "GEO", "sapiens_NA_GSE52778.tsv")
    gdd$SDRFfile <- file.path(extDir, "GEO", "sapiens_NA_GSE52778_SDRF.tsv")
    gdd$IDFfile  <- file.path(extDir, "GEO", "sapiens_NA_GSE52778_IDF.tsv")
    return(gdd)
}
