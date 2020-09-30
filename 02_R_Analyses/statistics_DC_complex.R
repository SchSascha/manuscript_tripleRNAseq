###############################
### DESeq2 complex test for CMV

library("DESeq2")
library("BiocParallel")
library("WriteXLS")
library("tidyverse")
library("gplots")

register(MulticoreParam(6))

get_items_per_intersection <- function(l) attr(venn(l, show.plot = F, ), "intersections")

load("R_workspace.RData")
outDir <- "../results/stats_DC_complex"
dir.create(outDir)

WriteXLS("meta", "../results/stats_DC_complex/meta.xls")

# test difference in singleinfections FOR TIME
# model: y ~ pathogen + time (ALL samples)
# Q1: compensating for pathogens, is there a systematic change between longer and shorter cultivation

meta_tmp <- meta %>% filter(type %in% c("alone", "si_AFU", "si_CMV"))
counts_sp1_tmp <- counts_sp1[,meta_tmp$samples]
rownames(meta_tmp) <- meta_tmp$samples
meta_tmp <- meta_tmp[,-1]
# rewrite time - we only want baseline & after
meta_tmp$growth[meta_tmp$growth == "4h30min"] <- "2h"
meta_tmp$type <- droplevels(meta_tmp$type)
meta_tmp$growth <- droplevels(meta_tmp$growth)
levels(meta_tmp$growth) <- c("0h", "xh")
dds <- DESeqDataSetFromMatrix(countData = counts_sp1_tmp,
                              colData = meta_tmp,
                              design = ~ donor + type + growth)
dds <- DESeq(dds, parallel = T)
res1 <- results(dds, contrast=c("growth", "xh", "0h"))
res1 <- res1[order(res1$padj),]
res1 <- as.data.frame(res1)
res1$sig <- res1$padj < 0.01
WriteXLS(
    "res1",
    ExcelFileName = file.path(outDir, "d+t+g.time.si.0h_vs_xh.xls"),
    row.names = T,
    col.names = T,
    BoldHeaderRow = T,
    FreezeRow = 1,
    FreezeCol = 1,
    AdjWidth = T
)

tmp <- res1 %>% mutate(gene = rownames(.)) %>% filter(padj < 0.01)
WriteXLS(
    "tmp",
    ExcelFileName = file.path(outDir, "d+t+g.time.si.0h_vs_xh.p001.xls"),
    row.names = T,
    col.names = T,
    BoldHeaderRow = T,
    FreezeRow = 1,
    FreezeCol = 1,
    AdjWidth = T
)



# test difference in singleinfections
# model: y ~ pathogen + time (CMV | AFU samples)
# Q1: compensating for time, is there a systematic change in pathogens
meta_tmp <- meta %>% filter(type %in% c("alone", "si_AFU"))
counts_sp1_tmp <- counts_sp1[,meta_tmp$samples]
rownames(meta_tmp) <- meta_tmp$samples
meta_tmp <- meta_tmp[,-1]
meta_tmp$type <- droplevels(meta_tmp$type)
meta_tmp$growth <- droplevels(meta_tmp$growth)
dds <- DESeqDataSetFromMatrix(countData = counts_sp1_tmp,
                              colData = meta_tmp,
                              design = ~ donor + type + growth)
dds <- DESeq(dds, parallel = T)
res1 <- results(dds, contrast=c("type", "si_AFU", "alone"))
res1 <- res1[order(res1$padj),]
res1 <- as.data.frame(res1)
res1$sig <- res1$padj < 0.01
WriteXLS(
  "res1",
  ExcelFileName = file.path(outDir, "d+t+g.type.si_afu_vs_alone.xls"),
  row.names = T,
  col.names = T,
  BoldHeaderRow = T,
  FreezeRow = 1,
  FreezeCol = 1,
  AdjWidth = T
)

res1 <- gdata::read.xls(file.path(outDir, "d+t+g.type.si_afu_vs_alone.xls"))
tmp <- res1 %>% mutate(gene = rownames(.)) %>% filter(padj < 0.01 & abs(log2FoldChange) >= 1.0)
WriteXLS(
    "tmp",
    ExcelFileName = file.path(outDir, "d+t+g.type.si_afu_vs_alone.p001.lfc1.xls"),
    row.names = T,
    col.names = T,
    BoldHeaderRow = T,
    FreezeRow = 1,
    FreezeCol = 1,
    AdjWidth = T
)



res1 %>% mutate(gene = rownames(.)) %>% filter(padj < 0.01) %>% dim
res1 %>% mutate(gene = rownames(.)) %>% filter(padj < 0.01 & abs(log2FoldChange) >= 1.0) %>% dim



# without Donors, to be sure
dds <- DESeqDataSetFromMatrix(countData = counts_sp1_tmp,
                              colData = meta_tmp,
                              design = ~ type + growth)
dds <- DESeq(dds, parallel = T)
res2 <- results(dds, contrast=c("type", "si_AFU", "alone"))
res2 <- res2[order(res2$padj),]
res2 <- as.data.frame(res2)
res2$sig <- res2$padj < 0.01
WriteXLS(
    "res2",
    ExcelFileName = file.path(outDir, "t+g.type.si_afu_vs_alone.xls"),
    row.names = T,
    col.names = T,
    BoldHeaderRow = T,
    FreezeRow = 1,
    FreezeCol = 1,
    AdjWidth = T
)


# check overlap
l <- get_items_per_intersection(list(with_donor = res1 %>% mutate(gene = rownames(.)) %>% filter(padj < 0.01 & abs(log2FoldChange) >= 1.0) %>% .$gene,
                                     not_donor  = res2 %>% mutate(gene = rownames(.)) %>% filter(padj < 0.01 & abs(log2FoldChange) >= 1.0) %>% .$gene))
l %>% map(length)




#######
## Model to test for differences in condition (without time)
meta_tmp <- meta %>% filter(condition != "none")
counts_sp1_tmp <- counts_sp1[,meta_tmp$sample]
rownames(meta_tmp) <- meta_tmp$samples
meta_tmp <- meta_tmp[,-1]
meta_tmp$type <- droplevels(meta_tmp$type)
dds <- DESeqDataSetFromMatrix(countData = counts_sp1_tmp,
                              colData = meta_tmp,
                              design = ~ donor + type)
dds <- DESeq(dds, parallel = T)
test_names1 <- grep("type_", resultsNames(dds), value = T)
# all tests are against alone by default
# we want co vs si as well
test_names2 <- list(c("type", "co_Both", "si_AFU"), c("type", "co_Both", "si_CMV"))



deg_res_list <- plyr::llply(test_names1, function(x) {
  res <- results(dds, name = x)
  res <- res[order(res$padj),]
  res <- as.data.frame(res)
  res$sig <- res$padj < 0.01
  return(res)
})
# names(deg_res_list) <- test_names

deg_res_list <- c(deg_res_list, plyr::llply(test_names2, function(x) {
    res <- results(dds, contrast = x)
    res <- res[order(res$padj),]
    res <- as.data.frame(res)
    res$sig <- res$padj < 0.01
    return(res)
}))

names(deg_res_list) <- c(test_names1, lapply(test_names2, paste, collapse = "_"))


deg_res_list_p <- plyr::llply(deg_res_list, function(df) return(df[ !is.na(df$padj) & df$padj < 0.01,]))

WriteXLS(
  "deg_res_list_p",
  ExcelFileName = file.path(outDir, "donor+type.type.all.p001.xls"),
  row.names = T,
  col.names = T,
  BoldHeaderRow = T,
  FreezeRow = 1,
  FreezeCol = 1,
  AdjWidth = T
)

deg_res_list_lfc1 <- plyr::llply(deg_res_list, function(df) return(df[!is.na(df$padj) & df$padj < 0.01 & abs(df$log2FoldChange) >= 1.0,]))
WriteXLS(
  "deg_res_list_lfc1",
  ExcelFileName = file.path(outDir, "donor+type.type.all.p001.lfc1.xls"),
  row.names = T,
  col.names = T,
  BoldHeaderRow = T,
  FreezeRow = 1,
  FreezeCol = 1,
  AdjWidth = T
)





#############
### CMV
####

meta_tmp <- meta %>% filter(type %in% c("alone", "si_CMV"))
counts_sp1_tmp <- counts_sp1[,meta_tmp$samples]
rownames(meta_tmp) <- meta_tmp$samples
meta_tmp <- meta_tmp[,-1]
meta_tmp$type <- droplevels(meta_tmp$type)
meta_tmp$growth <- droplevels(meta_tmp$growth)
dds <- DESeqDataSetFromMatrix(countData = counts_sp1_tmp,
                              colData = meta_tmp,
                              design = ~ donor + type + growth)
dds <- DESeq(dds, parallel = T)
res1 <- results(dds, contrast=c("type", "si_CMV", "alone"))
res1 <- res1[order(res1$padj),]
res1 <- as.data.frame(res1)
res1$sig <- res1$padj < 0.01
WriteXLS(
  "res1",
  ExcelFileName = file.path(outDir, "d+t+g.type.si_cmv_vs_alone.xls"),
  row.names = T,
  col.names = T,
  BoldHeaderRow = T,
  FreezeRow = 1,
  FreezeCol = 1,
  AdjWidth = T
)

res1 <- gdata::read.xls(file.path(outDir, "d+t+g.type.si_cmv_vs_alone.xls"))
tmp <- res1 %>% mutate(gene = rownames(.)) %>% filter(padj < 0.01 & abs(log2FoldChange) >= 1.0)
WriteXLS(
    "tmp",
    ExcelFileName = file.path(outDir, "d+t+g.type.si_cmv_vs_alone.p001.lfc1.xls"),
    row.names = T,
    col.names = T,
    BoldHeaderRow = T,
    FreezeRow = 1,
    FreezeCol = 1,
    AdjWidth = T
)


# without Donors, to be sure
dds <- DESeqDataSetFromMatrix(countData = counts_sp1_tmp,
                              colData = meta_tmp,
                              design = ~ type + growth)
dds <- DESeq(dds, parallel = T)
res2 <- results(dds, contrast=c("type", "si_CMV", "alone"))
res2 <- res2[order(res2$padj),]
res2 <- as.data.frame(res2)
res2$sig <- res2$padj < 0.01
WriteXLS(
    "res2",
    ExcelFileName = file.path(outDir, "t+g.type.si_cmv_vs_alone.xls"),
    row.names = T,
    col.names = T,
    BoldHeaderRow = T,
    FreezeRow = 1,
    FreezeCol = 1,
    AdjWidth = T
)


# check overlap
l <- get_items_per_intersection(list(with_donor = res1 %>% mutate(gene = rownames(.)) %>% filter(padj < 0.01 & abs(log2FoldChange) >= 1.0) %>% .$gene,
                                     not_donor  = res2 %>% mutate(gene = rownames(.)) %>% filter(padj < 0.01 & abs(log2FoldChange) >= 1.0) %>% .$gene))



# AFU | CMV vs alone overlap
resx1 <- gdata::read.xls(file.path(outDir, "d+t+g.type.si_cmv_vs_alone.xls"))
resx1 <- resx1 %>% mutate(gene = rownames(.)) %>% filter(padj < 0.01 & abs(log2FoldChange) >= 1.0)

resx2 <- gdata::read.xls(file.path(outDir, "d+t+g.type.si_afu_vs_alone.xls"))
resx2 <- resx2 %>% mutate(gene = rownames(.)) %>% filter(padj < 0.01 & abs(log2FoldChange) >= 1.0)

l <- get_items_per_intersection(list(cmv = resx1$X, afu = resx2$X))




######## Salmon tests for comparison
library("tximport")
library("readr")
library("data.table")
library("tidyverse")
library("DESeq2")
library("WriteXLS")

load("R_workspace.RData")

# create tx2gene from annotation
anno_file <- "../../../../sbidata/shared2/Combined_Genomes_forTripleRNAseqPreprocessing/Hsapiens_GRCh38_89_2017-05-07_Human_herpesvirus_5_strain_TB40E_clone_TB40-BAC4/anno/hs-GR38.89_afu-Af293_cmv-5.gtf"
anno <- data.table::fread(anno_file, sep = "\t", sep2 = ";")
anno <- anno$V9
anno <- gsub("\"", "", anno)
tx2gene <- anno %>% strsplit(split = ";") %>%
    sapply(., function(x) grep("(transcript_id)|(gene_id)", x, value = T))
tx2gene <- tx2gene %>% t
tx2gene <- tx2gene %>% sub(" transcript_id ", "", .) %>% sub("gene_id ", "", .)
tx2gene <- as.data.frame(tx2gene)
tx2gene <- tx2gene[,2:1]

res_dirs <- list.dirs(path = "../results/salmon_u/", full.names = TRUE, recursive = F)
files <- paste0(res_dirs, "/quant.sf")
stopifnot(file.exists(files))
sal_counts_sp1 <- tximport(files, type="salmon", tx2gene = tx2gene)

stopifnot(files %>% gsub(".trimo.quant/quant.sf", "", .) %>% gsub(".*(//)", "", .) ==
SDRF_info$table$Data.File %>% gsub(".*/", "", .) %>% gsub(".fq.gz", "", .))



# keep only DC genes
tx2gene <- tx2gene %>% filter(grepl("ENST", V2))

meta_tmp <- meta %>% filter(type %in% c("alone", "si_AFU", "si_CMV"))
counts_sp1_tmp <- sal_counts_sp1
keep <- rownames(counts_sp1_tmp$counts) %in% (tx2gene$V1 %>% unique)
counts_sp1_tmp$counts <- counts_sp1_tmp$counts[keep, meta$samples %in% meta_tmp$samples]
counts_sp1_tmp$abundance <- counts_sp1_tmp$abundance[keep, meta_tmp$samples]
counts_sp1_tmp$length <- counts_sp1_tmp$length[keep, meta_tmp$samples]
rownames(meta_tmp) <- meta_tmp$samples
meta_tmp <- meta_tmp[,-1]
# rewrite time - we only want baseline & after
meta_tmp$growth[meta_tmp$growth == "4h30min"] <- "2h"
meta_tmp$type <- droplevels(meta_tmp$type)
meta_tmp$growth <- droplevels(meta_tmp$growth)
levels(meta_tmp$growth) <- c("0h", "xh")
ddsTxi <- DESeqDataSetFromTximport(counts_sp1_tmp,
                                   colData = meta_tmp,
                                   design = ~ donor + type + growth)

dds <- DESeq(ddsTxi, parallel = T)
sal_res1 <- results(dds, contrast=c("growth", "xh", "0h"))
sal_res1 <- sal_res1[order(sal_res1$padj),]
sal_res1 <- as.data.frame(sal_res1)
sal_res1$sig <- sal_res1$padj < 0.01
WriteXLS(
    "sal_res1",
    ExcelFileName = file.path(outDir, "sal_d+t+g.time.si.0h_vs_xh.xls"),
    row.names = T,
    col.names = T,
    BoldHeaderRow = T,
    FreezeRow = 1,
    FreezeCol = 1,
    AdjWidth = T
)

# check overlap with featureCounts data
l1 <- get_items_per_intersection(list(feat = res1 %>% mutate(gene = rownames(.)) %>% filter(padj < 0.01 & abs(log2FoldChange) >= 1.0) %>% .$gene,
                                      salm  = sal_res1 %>% mutate(gene = rownames(.)) %>% filter(padj < 0.01 & abs(log2FoldChange) >= 1.0) %>% .$gene))
l2 <- get_items_per_intersection(list(feat = res1 %>% mutate(gene = rownames(.)) %>% filter(padj < 0.01) %>% .$gene,
                                      salm  = sal_res1 %>% mutate(gene = rownames(.)) %>% filter(padj < 0.01) %>% .$gene))

