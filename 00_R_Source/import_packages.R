#' @importFrom baySeq getLikelihoods getPriors.NB getLikelihoods.NB topCounts
#' @importFrom DESeq estimateSizeFactors estimateDispersions nbinomTest newCountDataSet
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq
#' @importFrom edgeR DGEList calcNormFactors estimateCommonDisp estimateTagwiseDisp exactTest
#' @importFrom futile.logger flog.threshold
# @importFrom GenomicRanges exonsBy
#' @importFrom gdata last read.xls
#' @importFrom genefilter rowVars
#' @importFrom graphics abline axis barplot legend mtext par plot plot.new points title
#' @importFrom grDevices colorRampPalette dev.copy2pdf dev.off pdf png rainbow svg
#' @importFrom grid grid.draw
#' @importFrom gridExtra grid.table
#' @importFrom gtools combinations
#' @importFrom limma voom lmFit eBayes
# @importFrom IRanges PartitioningByEnd
#' @importFrom methods new
#' @importFrom NOISeq degenes readData noiseqbio
#' @importFrom parallel makeCluster mclapply stopCluster
# @importFrom pasillaBamSubset dm3_chr4 untreated1_chr4 untreated3_chr4
#' @importFrom pheatmap pheatmap
#' @importFrom RCurl getURL url.exists
#' @importFrom RSQLite dbConnect dbDisconnect
#' @importFrom Rsamtools asSam asBam countBam scanBamFlag scanBamHeader ScanBamParam sortBam
#' @importFrom Rsubread featureCounts
#' @importFrom R.utils gunzip
#' @importFrom samr SAMseq
#' @importFrom ShortRead FastqStreamer
#' @importFrom stats as.dendrogram cor dist formula hclust model.matrix prcomp p.adjust relevel
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment assay
#' @importFrom tools file_path_sans_ext
#' @importFrom UpSetR fromList upset
#' @importFrom utils globalVariables capture.output combn read.csv read.table write.table zip compareVersion download.file flush.console installed.packages packageVersion tar unzip
#' @importFrom VennDiagram venn.diagram
#' @importFrom WriteXLS WriteXLS
#' @import Biostrings
#' @import locfit
#' @import RColorBrewer
# @import PoissonSeq
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "3.2.0") utils::globalVariables(c("."))


#' Minimal Counting Table
#'
#' This data set is used as minimum example of a reads-per-gene counting table.
#'
#' @name counts_min
#' @aliases counts_min
#' @docType data
#' @author Bastian Seelbinder \email{bastian.seelbinder@leibniz-hki.de}
#' @keywords data
NULL

#' Minimal Differential Expression Table
#'
#' This data set is used as minimum example of a gene significance table.
#'
#' @name deg_min
#' @aliases deg_min
#' @docType data
#' @author Bastian Seelbinder \email{bastian.seelbinder@leibniz-hki.de}
#' @keywords data
NULL

#' getGEOdata Return Value
#'
#' This data set is used to show the return value of \code{\link{getGEOdata}}
#' without having to call the function. This is used in the vignette.
#'
#' @name geo_dat
#' @aliases geo_dat
#' @docType data
#' @author Bastian Seelbinder \email{bastian.seelbinder@leibniz-hki.de}
#' @keywords data
NULL
