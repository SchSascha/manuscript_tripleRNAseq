###############################
### DESeq2 complex test for CMV

library("DESeq2")
library("BiocParallel")
library("WriteXLS")
register(MulticoreParam(10))

outDir <- "../results/stats_CMV_complex"
load("R_workspace.RData")

meta <- data.frame(samples = rownames(design_matrix_sp3), condition = conds_sp3)
meta$type <- "none"
meta$type[grepl("SingleInfection_CMV", meta$condition)] <- "si"
meta$type[grepl("CoInfection_CMV_first", meta$condition)] <- "co_CMV"
meta$type[grepl("Coinfection_Afu_first", meta$condition)] <- "co_Afu"
meta$type[meta$condition == "CoInfection"] <- "co_Both"
meta$type <- factor(meta$type)
meta$type <- relevel(meta$type, "si")

meta$growth <- "none"
meta$growth[grepl("(0h)|(CMV_first)", meta$condition)] <- "0h"
meta$growth[meta$condition == "CoInfection"] <- "0h"
meta$growth[grepl("(2h)|(Afu_first)", meta$condition)] <- "2h"
meta$growth <- factor(meta$growth)
meta$growth <- relevel(meta$growth, "0h")

meta$donor <- gsub("_.*", "", meta$samples)
meta$donor <- factor(meta$donor)

meta$comb <- "none"
meta$comb[grepl("Single", meta$condition)] <- "si"
meta$comb[grepl("CoInfection", meta$condition, ignore.case = T)] <- "co"
meta$comb <- factor(meta$comb)
meta$comb <- relevel(meta$comb, "si")


# simpel test with all covariates
meta2 <- meta
rownames(meta2) <- meta$samples
meta2 <- meta2[,-1]

counts_sp3_2 <- counts_sp3[,meta2$condition != "none"]
meta2 <- meta2[meta2$condition != "none",]

dds <- DESeqDataSetFromMatrix(countData = counts_sp3_2,
                              colData = meta2,
                              design = ~ donor + type + growth)
dds <- DESeq(dds)
res1 <- results(dds, contrast=c("type", "co_Both", "si"))
res1 <- res1[order(res1$padj),]
res1 <- as.data.frame(res1)
res1$sig <- res1$padj < 0.05
WriteXLS(
    "res1",
    ExcelFileName = file.path(outDir, "type.d+t+g.co_Both_vs_si.xls"),
    SheetNames = "co_Both_vs_si",
    row.names = T,
    col.names = T,
    BoldHeaderRow = T,
    FreezeRow = 1,
    FreezeCol = 1,
    AdjWidth = T
)

res1_2 <- results(dds, list("growth_2h_vs_0h"))
res1_2 <- res1_2[order(res1_2$padj),]
res1_2 <- as.data.frame(res1_2)
res1_2$sig <- res1_2$padj < 0.05
WriteXLS(
    "res1_2",
    ExcelFileName = file.path(outDir, "growth.d+t+g.growth_2h_vs_0h.xls"),
    SheetNames = "co_Both_vs_si",
    row.names = T,
    col.names = T,
    BoldHeaderRow = T,
    FreezeRow = 1,
    FreezeCol = 1,
    AdjWidth = T
)



dds2 <- DESeqDataSetFromMatrix(countData = counts_sp3_2,
                               colData = meta2,
                               design = ~ type + growth)
dds2 <- DESeq(dds2)
res2 <- results(dds2, contrast=c("type", "co_Both", "si"))
res2 <- res2[order(res2$padj),]
res2 <- as.data.frame(res2)
res2$sig <- res2$padj < 0.05
WriteXLS(
    "res2",
    ExcelFileName = file.path(outDir, "type.t+g.co_Both_vs_si.xls"),
    SheetNames = "co_Both_vs_si",
    row.names = T,
    col.names = T,
    BoldHeaderRow = T,
    FreezeRow = 1,
    FreezeCol = 1,
    AdjWidth = T
)



dds3 <- DESeqDataSetFromMatrix(countData = counts_sp3_2,
                               colData = meta2,
                               design = ~ donor + comb + growth + comb:growth)
dds3 <- DESeq(dds3, parallel = T)
res3 <- results(dds3, list("combco.growth2h"))
res3 <- res3[order(res3$padj),]
res3 <- as.data.frame(res3)
res3$sig <- res3$padj < 0.05
WriteXLS(
    "res3",
    ExcelFileName = file.path(outDir, "type.d+c*g.combco.growth2h.xls"),
    SheetNames = "co_Both_vs_si",
    row.names = T,
    col.names = T,
    BoldHeaderRow = T,
    FreezeRow = 1,
    FreezeCol = 1,
    AdjWidth = T
)



###############
## make Box-plot
library(stringr)
library(tidyverse)
library(magrittr)

g <- ncounts_sp3_2 %>% as.data.frame %>% mutate(gene = rownames(.)) %>% gather("condition", "reads", -gene) %>%
    mutate(condition = ifelse(stringr::str_count(string = .$condition, pattern = "plus") == 1, "SI", "CO")) %>%
    ggplot(aes(x = gene, y = reads, color = condition, group = condition)) +
    geom_line() +
    ggtitle("# Normalized Reads per Condition") +
    labs(x = "Gene", y = "# norm. Reads") +
    theme_bw() +
    theme(axis.text.x=element_blank())
g
ggsave("../results/other_plots/CMV_normalized_read_distribution.grouped.pdf", dpi = 600, width = 10)

g <- ncounts_sp3_2 %>% as.data.frame %>%
    mutate(Gene = rownames(.)) %>%
    filter(Gene %in% na.omit(rownames(res1)[res1$sig])) %>%
    gather("condition", "Reads", -Gene) %>%
    mutate(condition = ifelse(stringr::str_count(string = .$condition, pattern = "plus") == 1, "SI", "CO")) %>%
    mutate(Gene = sub("cds-ABV", "", Gene)) %>%
    ggplot(aes(x = Gene, y = Reads, fill = condition)) +
    geom_boxplot(notch = T) +
    scale_y_log10() +
    ggtitle("Top 12 Most significant CMV Genes") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
g
ggsave("../results/other_plots/CMV_complex_top_box.log.pdf", g, dpi = 400)


# CMV - multifactor test with EdgeR
# library("edgeR")
# meta2 <- droplevels(meta2)
# y <- DGEList(counts_sp3_2)
# y$samples <- cbind(y$samples, meta2)
# y <- calcNormFactors(y)
# design <- model.matrix(~ donor + growth + type, meta2)
# y <- estimateDisp(y, design)
# fit <- glmQLFit(y, design)
#
# qlf <- glmQLFTest(fit) # co_both vs intercept
# topTags(qlf)
