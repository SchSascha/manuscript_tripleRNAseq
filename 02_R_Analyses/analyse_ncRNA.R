#!R
library("refGenome")
library("tidyverse")
library("broom")
library("readxl")
library("writexl")

get_sign_genes_geo2 <- function(deg_res, with_fold = F, sigP = 0.05) {
    l <- lapply(deg_res, function(x) {
        # select tools
        tools <- grep("_adj_pval", x %>% colnames, value = T, fixed = T)
        idx <- x[,tools] %>% {. < sigP} %>% {rowSums(.) == length(tools)} %>% which
        x <- x[idx,]
        return(x)
    })

    names(l) <- names(deg_res)
    return(l)
}

simp_names <- function(s) {
    s <- gsub("Co[Ii]nfection" , "CO", s)
    s <- gsub("SingleInfection", "SI", s)
    s <- gsub("SingleCulture"  , "SC", s)
    return(s)
}


# create ensemblGenome object for storing Ensembl genomic annotation data
ens <- ensemblGenome()

# read GTF file into ensemblGenome object
read.gtf(ens, "../doc/Homo_sapiens.GRCh38.89.ncRNA.gtf")
head(ens)

# counts all annotations on each seqname
tableSeqids(ens)

# summarise features in GTF file
tableFeatures(ens) # 40125 exons

# create table of genes
genes <- getGenePositions(ens)
head(genes)

# biotypes
dim(genes) # 10406 nc
table(genes$gene_biotype) %>% tidy

# 3prime_overlapping_ncRNA         31
# bidirectional_promoter_lncRNA     8
# lincRNA                        7514
# macro_lncRNA                      1
# snoRNA                          943
# snRNA                          1909

write_xlsx(genes, "../results/stats/Hsapiens_ncRNA/ncRNA_anno.xlsx")


#### Check which/how many ncRNA are differentially expressed
dc_counts <- read.csv("../results/stats/Hsapiens_counts.csv", sep = "\t")
sum(genes$gene_id %in% rownames(dc_counts[rowSums(dc_counts) >= 1,]))
sum(genes$gene_id %in% rownames(dc_counts[rowSums(dc_counts) >= 10,]))
sum(genes$gene_id %in% rownames(dc_counts[rowSums(dc_counts) >= 100,]))
sum(genes$gene_id %in% rownames(dc_counts[rowSums(dc_counts) >= 1000,]))


#### DEG
dc_deg <- readRDS("../results/DEG_overlaps/deg_res_sp1.rds")
names(dc_deg$DEGs) <- names(dc_deg$DEGs) %>% simp_names()
l <- get_sign_genes_geo2(dc_deg$DEGs)

nc_list <- lapply(l, function(df) {
    idx <- which(df$id %in% genes$gene_id)
    if (idx %>% length == 0)
        return(data.frame())
    else
        df[idx,]
})

linc_list <- lapply(l, function(df) {
    idx <- which(df$id %in% (genes %>% filter(gene_biotype == "lincRNA") %>% .$gene_id))
    if (idx %>% length == 0)
        return(data.frame())
    else
        df[idx,]
})

stats_df <- data.frame(ndeg = sapply(l, nrow), ndeg_nc = sapply(nc_list, nrow), ndeg_linc = sapply(linc_list, nrow))
stats_df <- stats_df %>% mutate("test" = rownames(stats_df)) %>% select("test", everything())

dir.create("../results/stats/Hsapiens_ncRNA")
write_xlsx(stats_df, "../results/stats/Hsapiens_ncRNA/stats.xlsx")
write_xlsx(nc_list, "../results/stats/Hsapiens_ncRNA/deg_ncRNA.xlsx")
write_xlsx(linc_list, "../results/stats/Hsapiens_ncRNA/deg_lincRNA.xlsx")
