library("tidyverse")
source("scripts/lib_annotate.R")
source("scripts/lib_overlaps.R")


# @param matrix. samples in columns
# @param threshold sample-wise cutoff value
# @param as.lower If TRUE, a values must be equal or below the given threshold. Otherwise, values must be equal or greater.
get_gene_list_per_sample <- function(mat, threshold, as.lower = F) {
    n <- ncol(mat)
    l <- lapply(1:n, function(i) if (as.lower) rownames(mat)[mat[,i] <= threshold[i]] else rownames(mat)[mat[,i] >= threshold[i]])
    names(l) <- colnames(mat)
    return(l)
}


load("R_workspace.RData")

# re-annotate counts_sp1 genes to entrez
columns(org.Hs.eg.db)
entrez <- mapIds(org.Hs.eg.db, keys = rownames(counts_sp1), column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
entrez <- na.omit(entrez)
counts_sp1 <- counts_sp1[names(entrez), ]
rownames(counts_sp1) <- entrez

# load metabolically relavent genes
metabolic_entrez <- read.csv("../doc/recon301_genes.tsv", sep = "\t") %>% unlist() %>% sub("\\..+", "", .)

# filter
keep <- rownames(counts_sp1) %in% metabolic_entrez
counts_sp1_sel <- counts_sp1[keep,]


# filter zero coverage genes
counts_sp1_sel <- counts_sp1_sel[rowSums(counts_sp1_sel) > 0,]
counts_sp1_sel <- counts_sp1_sel[,grep("Afu_alone", colnames(counts_sp1_sel), invert = T)]
colnames(counts_sp1_sel) <- sub("Donor", "D", colnames(counts_sp1_sel))
colnames(counts_sp1_sel) <- gsub("plus", "+", colnames(counts_sp1_sel))
colnames(counts_sp1_sel) <- gsub(".\\+.", "+", colnames(counts_sp1_sel))
colnames(counts_sp1_sel) <- gsub("4h30min", "4.5h", colnames(counts_sp1_sel))
n <- ncol(counts_sp1_sel)

sapply(1:n, function(i) hist(log(counts_sp1_sel[,i]+1), breaks = 1000, freq = F))
quant_df <- sapply(1:n, function(i) quantile(counts_sp1_sel[,i], c(.20, .25, .30, .8, .75, .70)))
colnames(quant_df) <- colnames(counts_sp1_sel)
quant_df

write.table(quant_df, "../results/metabolic/quant_df.tsv", sep = "\t", quote = F)


# < 30% and > 70% seems to be most symmetric

l.low <- get_gene_list_per_sample(counts_sp1_sel, quant_df[3,], as.lower = T)
l.up <- get_gene_list_per_sample(counts_sp1_sel, quant_df[6,], as.lower = F)
names(l.low) <- paste("low.", names(l.low))
names(l.up) <- paste("up.", names(l.up))
l <- c(l.low, l.up)
WriteXLS::WriteXLS(lapply(l, data.frame), "../results/metabolic/genes.q30.xls")
make_upset_plot(l, file = "../results/metabolic/upset.q30.pdf", save = T, main = "q = 30%", keep.order = T)

l.low <- get_gene_list_per_sample(counts_sp1_sel, quant_df[2,], as.lower = T)
l.up <- get_gene_list_per_sample(counts_sp1_sel, quant_df[5,], as.lower = F)
names(l.low) <- paste("low.", names(l.low))
names(l.up) <- paste("up.", names(l.up))
l <- c(l.low, l.up)
WriteXLS::WriteXLS(lapply(l, data.frame), "../results/metabolic/genes.q25.xls")
make_upset_plot(l, file = "../results/metabolic/upset.q25.pdf", save = T, main = "q = 25%")

l.low <- get_gene_list_per_sample(counts_sp1_sel, quant_df[1,], as.lower = T)
l.up <- get_gene_list_per_sample(counts_sp1_sel, quant_df[4,], as.lower = F)
names(l.low) <- paste("low.", names(l.low))
names(l.up) <- paste("up.", names(l.up))
l <- c(l.low, l.up)
WriteXLS::WriteXLS(lapply(l, data.frame), "../results/metabolic/genes.q20.xls")
make_upset_plot(l, file = "../results/metabolic/upset.q20.pdf", save = T, main = "q = 20%")
