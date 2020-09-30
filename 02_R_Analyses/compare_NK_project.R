#!/bin/env R
library("plyr")
library("data.table")
library("futile.logger")
library("tidyverse")
library("VennDiagram")
library("ggplot2")
library("UpSetR")
library("gplots") # to calculate venn overlap
library("AnnotationDbi")
library("org.Hs.eg.db")


get_sign_genes_geo2 <- function(x, with_fold = F, sigP = 0.01) {
    # select tools
    tools <- grep("_adj_pval", x %>% colnames, value = T, fixed = T)
    if (!with_fold) {
        deg_genes <- x[,tools] %>% {. < sigP} %>% {rowSums(.) == length(tools)} %>% which %>% names
    } else {
        inds <- x[,tools] %>% {. < sigP} %>% {rowSums(.) == length(tools)} %>% which
        up <- x$log2_fc_mrn[inds] >= 0
        deg_genes <- paste0(rownames(x)[inds], ifelse(up, "+", "-"))
    }

    return(deg_genes)
}

get_sign_genes_nk_list <- function(deg_res, with_fold = F, sigP = 0.01) {
    l <- lapply(deg_res, function(x) {
        # select tools
        tools <- grep("_adj_pval", x %>% colnames, value = T, fixed = T)
        if (!with_fold) {
            deg_genes <- x[,tools] %>% {. < sigP} %>% {rowSums(.) == length(tools)} %>% which %>% names
        } else {
            inds <- x[,tools] %>% {. < sigP} %>% {rowSums(.) == length(tools)} %>% which
            up <- x$log2_rpkm[inds] >= 0
            deg_genes <- paste0(rownames(x)[inds], ifelse(up, "+", "-"))
        }

        return(deg_genes)
    })

    names(l) <- names(deg_res)
    return(l)
}

map_gene_names <- function(genes, orgdb, from_type = "ENSEMBL", to_type = "SYMBOL") {
    mapped_gene <- mapIds(org.Hs.eg.db, keys = genes, keytype = from_type, column = to_type)
    # original names of 'genes' are stored in names(mapped_gene)
    use_original <- is.na(mapped_gene)
    mapped_gene[use_original] <- names(mapped_gene)[use_original]
    return(mapped_gene)
}


get_items_per_intersection <- function(l) attr(venn(l, show.plot = F, ), "intersections")

rem_sign <- function(s) sub("[-|+]$", "", s)

make_venn_plot <- function(l, file = "", save = F, main = "Venn Diagram") {
    if (save && file == "")
        stop("Please supply file name if you want to save the figure!")
    if (length(l) > 5)
        stop("Cannot make Venn diagram with more than 5 sets! Use UpSet instead!")

    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") # no log files
    vp <- VennDiagram::venn.diagram(
        x = l,
        filename = NULL,
        main = main,
        main.fontface = "bold",
        sub = paste0("total number of Items: ", length(unique(unlist(l)))),
        fill = rainbow(length(l)),
        height = 3000,
        width = 3000,
        resolution = 500,
        print.mode = if (length(l) < 5) c("raw", "percent") else "raw"
    )

    if (save) {
        pdf(file)
        grid::grid.draw(vp)
        dev.off()
    } else {
        grid::grid.newpage()
        grid::grid.draw(vp)
    }

    return(vp)
}



##############
# Part 1: DC
# basically, we want to know if there is a difference between treating DC cells with AFu compared to treating NK cells with AFu

############
#>#># MOL 4
#######
deg_res_sp1 <- readRDS("../results/DEG_overlaps/deg_res_sp1.rds")
xf_deg <- deg_res_sp1$DEGs$SingleInfection_Afu_0h_vs_SingleCulture_DC
xf_best <- get_sign_genes_geo2(xf_deg, with_fold = T)

nk_df <- data.table::fread("../other_data/NK92_4_VS_HS.csv") %>% as.data.frame
nk_best <- nk_df %>% filter(abs(log2FC) >= 1.0)
signs <- ifelse(nk_best$log2FC >= 0, "+", "-")
nk_best <- paste(nk_best$id, signs, sep = "")

# general overlap of annotation
l <- list(NK = nk_df$id, DC = xf_deg$id)
inter <- get_items_per_intersection(l)
sapply(inter, length)
make_venn_plot(l, file = "../results/NK_compare/gene_list_overlap.DC.pdf", main = "Overlap of Annotation", save = T)

# best of both using signed genes
l <- list(DC = xf_best, NK = nk_best)
inter <- get_items_per_intersection(l)
sapply(inter, length)
make_venn_plot(l, file = "../results/NK_compare/DC.NK_noRep_best_sign_overlap.pdf", main = "Overlap of top genes with fc sign", save = T)
write_lines(inter$`DC:NK` %>% sub("[+-]", "", .) %>% map_gene_names(org.Hs.eg.db),
            "../results/NK_compare/DC.NK_noRep_best_sign_overlap.genelist.txt")

# best of both using unsigned genes
xf_best <- get_sign_genes_geo2(xf_deg, with_fold = F)
nk_best <- nk_df %>% filter(abs(log2FC) >= 1.0) %>% .$id
l <- list(DC = xf_best, NK = nk_best)
inter <- get_items_per_intersection(l)
sapply(inter, length)
make_venn_plot(l, file = "../results/NK_compare/DC.NK_noRep_best_noSign_overlap.pdf", main = "Overlap of top genes without fc sign", save = T)
write_lines(inter$`DC:NK` %>% sub("[+-]", "", .) %>% map_gene_names(org.Hs.eg.db),
            "../results/NK_compare/DC.NK_noRep_best_noSign_overlap.genelist.txt")


#># For Testing: using data with lower incubation time
xf_deg <- deg_res_sp1$DEGs$SingleInfection_Afu_4h30min_vs_SingleCulture_DC
xf_best <- get_sign_genes_geo2(xf_deg, with_fold = T)
l <- list(DC = xf_best, NK = nk_best)
inter <- get_items_per_intersection(l)
sapply(inter, length)
make_venn_plot(l, file = "../results/NK_compare/DC.AFU4.5h.NK_noRep_best_sign_overlap.pdf", main = "Overlap of top genes with fc sign", save = T)
write_lines(inter$`DC:NK` %>% sub("[+-]", "", .) %>% map_gene_names(org.Hs.eg.db),
            "../results/NK_compare/DC.AFU4.5h.NK_noRep_best_sign_overlap.genelist.txt")



############
#>#># MOL 0.5
#######

deg_res_sp1 <- readRDS("../results/DEG_overlaps/deg_res_sp1.rds")
xf_deg <- deg_res_sp1$DEGs$SingleInfection_Afu_0h_vs_SingleCulture_DC
xf_best <- get_sign_genes_geo2(xf_deg, with_fold = T)

nk_df <- data.table::fread("../other_data/NK92_0.5_VS_HS.csv") %>% as.data.frame
nk_best <- nk_df %>% filter(abs(log2FC) >= 1.0)
signs <- ifelse(nk_best$log2FC >= 0, "+", "-")
nk_best <- paste(nk_best$id, signs, sep = "")

# best of both using signed genes
l <- list(DC = xf_best, NK = nk_best)
inter <- get_items_per_intersection(l)
sapply(inter, length)
make_venn_plot(l, file = "../results/NK_compare/mol_0.5.DC.NK_noRep_best_sign_overlap.pdf", main = "Overlap of top genes with fc sign", save = T)
write_lines(inter$`DC:NK` %>% sub("[+-]", "", .) %>% map_gene_names(org.Hs.eg.db),
            "../results/NK_compare/mol_0.5.DC.NK_noRep_best_sign_overlap.genelist.txt")

# best of both using unsigned genes
xf_best <- get_sign_genes_geo2(xf_deg, with_fold = F)
nk_best <- nk_df %>% filter(abs(log2FC) >= 1.0) %>% .$id
l <- list(DC = xf_best, NK = nk_best)
inter <- get_items_per_intersection(l)
sapply(inter, length)
make_venn_plot(l, file = "../results/NK_compare/mol_0.5.DC.NK_noRep_best_noSign_overlap.pdf", main = "Overlap of top genes without fc sign", save = T)
write_lines(inter$`DC:NK` %>% sub("[+-]", "", .) %>% map_gene_names(org.Hs.eg.db),
            "../results/NK_compare/mol_0.5.DC.NK_noRep_best_noSign_overlap.genelist.txt")




##############
# Part 2: AFu
# basically, we want to know if there is a difference between treating AFu with NK cells or with DC

deg_res_sp2 <- readRDS("../results/DEG_overlaps/deg_res_sp2.rds")
xf_deg <- deg_res_sp2$DEGs$SingleInfection_Afu_0h_vs_SingleCulture_Afu_0h
xf_best <- get_sign_genes_geo2(xf_deg, with_fold = T)


#########################
## No Replicate Version
#####


nk_df <- data.table::fread("../other_data/NK92_4_VS_AF.csv") %>% as.data.frame
nk_best <- nk_df %>% filter(abs(log2FC) >= 1.0)
signs <- ifelse(nk_best$log2FC >= 0, "+", "-")
nk_best <- paste(nk_best$id, signs, sep = "")

# general overlap of annotation
l <- list(NK = nk_df$id, DC = xf_deg$id)
inter <- get_items_per_intersection(l)
sapply(inter, length)
make_venn_plot(l, file = "../results/NK_compare/gene_list_overlap.AFU.pdf", main = "Overlap of Annotation", save = T)

# best of both using signed genes
l <- list(DC = xf_best, NK = nk_best)
inter <- get_items_per_intersection(l)
sapply(inter, length)
make_venn_plot(l, file = "../results/NK_compare/AFU.NK_noRep_best_sign_overlap.pdf", main = "Overlap of top genes with fc sign", save = T)
write_lines(inter$`DC:NK` %>% sub("[+-]", "", .),
            "../results/NK_compare/AFU.NK_noRep_best_sign_overlap.genelist.txt")

# best of both using unsigned genes
xf_best <- get_sign_genes_geo2(xf_deg, with_fold = F)
nk_best <- nk_df %>% filter(abs(log2FC) >= 1.0) %>% .$id
l <- list(DC = xf_best, NK = nk_best)
inter <- get_items_per_intersection(l)
sapply(inter, length)
make_venn_plot(l, file = "../results/NK_compare/AFU.NK_noRep_best_noSign_overlap.pdf", main = "Overlap of top genes without fc sign", save = T)
write_lines(inter$`DC:NK` %>% sub("[+-]", "", .),
            "../results/NK_compare/AFU.NK_noRep_best_noSign_overlap.genelist.txt")


############
#>#># MOL 4
#######

nk_df <- data.table::fread("../other_data/NK92_0.5_VS_AF.csv") %>% as.data.frame
nk_best <- nk_df %>% filter(abs(log2FC) >= 1.0)
signs <- ifelse(nk_best$log2FC >= 0, "+", "-")
nk_best <- paste(nk_best$id, signs, sep = "")

# best of both using signed genes
l <- list(DC = xf_best, NK = nk_best)
inter <- get_items_per_intersection(l)
sapply(inter, length)
make_venn_plot(l, file = "../results/NK_compare/mol_0.5.AFU.NK_noRep_best_sign_overlap.pdf", main = "Overlap of top genes with fc sign", save = T)
write_lines(inter$`DC:NK` %>% sub("[+-]", "", .),
            "../results/NK_compare/mol_0.5.AFU.NK_noRep_best_sign_overlap.genelist.txt")


# best of both using unsigned genes
xf_best <- get_sign_genes_geo2(xf_deg, with_fold = F)
nk_best <- nk_df %>% filter(abs(log2FC) >= 1.0) %>% .$id
l <- list(DC = xf_best, NK = nk_best)
inter <- get_items_per_intersection(l)
sapply(inter, length)
make_venn_plot(l, file = "../results/NK_compare/mol_0.5.AFU.NK_noRep_best_noSign_overlap.pdf", main = "Overlap of top genes without fc sign", save = T)
write_lines(inter$`DC:NK` %>% sub("[+-]", "", .),
            "../results/NK_compare/mol_0.5.AFU.NK_noRep_best_noSign_overlap.genelist.txt")



#########################
## Replicate Version
#####

deg_res_sp2 <- readRDS("../results/DEG_overlaps/deg_res_sp2.rds")
xf_deg <- deg_res_sp2$DEGs$SingleInfection_Afu_0h_vs_SingleCulture_Afu_0h
xf_best <- get_sign_genes_geo2(xf_deg, with_fold = T)

# ASP_US aspergillus unstimulated for control
# DC_asp1 for direct comparison
# NK_asp1 for NK comparison
# NKDC_Asp1 for ? comparison

sheets <- readxl::excel_sheets("../other_data/Statistics.xlsx")
nk_df <- lapply(sheets, function(s) readxl::read_xlsx("../other_data/Statistics.xlsx", sheet = s))
names(nk_df) <- sheets
sheets
nk_df <- nk_df[c("DC_Asp_VS_Asp_US", "NK_Asp_VS_Asp_US", "NKDC_Asp_VS_Asp_US")]
nk_df <- lapply(nk_df, function(x) x[,-1])
nk_df <- lapply(nk_df, function(x) x[,c("id", "log2FC", "deseqpadj", "deseqpadj2", "res_limma$adjpvalues", "p_edgeR")])
cn <- c("id", "log2_rpkm", "DESeq_adj_pval", "DESeq2_adj_pval", "limma_adj_pval", "edgeR_adj_pval") # NOTE: fc is RPKM, not MRN. But I do not want to rewrite the stuff above...
nk_df <- lapply(nk_df, function(x) {x <- as.data.frame(x); colnames(x) <- cn; rownames(x) <- x$id; x})

# with signed genes
xf_best <- get_sign_genes_geo2(xf_deg, with_fold = T)
nk_best <- get_sign_genes_nk_list(nk_df, with_fold = T)

for(test in names(nk_best)) {
    l <- list(DC = xf_best, NK = nk_best[[test]])
    inter <- get_items_per_intersection(l)
    sapply(inter, length)
    make_venn_plot(l,
                   save = T,
                   file = paste0("../results/NK_compare/AFU.NK_Rep_best_sign.", test, ".pdf"),
                   main = paste0("Overlap Top Rep NK | XF - ", test))
    write_lines(inter$`DC:NK` %>% sub("[+-]", "", .),
                paste0("../results/NK_compare/AFU.NK_Rep_best_sign.", test, ".genelist.txt"))
}


# with unsigned genes
xf_best <- get_sign_genes_geo2(xf_deg, with_fold = F)
nk_best <- get_sign_genes_nk_list(nk_df, with_fold = F)

for(test in names(nk_best)) {
    l <- list(DC = xf_best, NK = nk_best[[test]])
    inter <- get_items_per_intersection(l)
    sapply(inter, length)
    make_venn_plot(l,
                   save = T,
                   file = paste0("../results/NK_compare/AFU.NK_Rep_best_noSign.", test, ".pdf"),
                   main = paste0("Overlap Top Rep NK | XF - no sign - ", test))
    write_lines(inter$`DC:NK` %>% sub("[+-]", "", .),
                paste0("../results/NK_compare/AFU.NK_Rep_best_noSign.", test, ".genelist.txt"))
}


# no deseq version
xf_best <- get_sign_genes_geo2(xf_deg[,c(-12,-13)], with_fold = F)
nk_df <- lapply(nk_df, function(x) x[,-3])
nk_best <- get_sign_genes_nk_list(nk_df, with_fold = F)

for(test in names(nk_best)) {
    l <- list(DC = xf_best, NK = nk_best[[test]])
    inter <- get_items_per_intersection(l)
    sapply(inter, length)
    make_venn_plot(l,
                   save = T,
                   file = paste0("../results/NK_compare/AFU.NK_Rep_best_noSign.nodeseq.", test, ".pdf"),
                   main = paste0("Overlap Top Rep NK | XF - no DESeq - no sign - ", test))
    write_lines(inter$`DC:NK` %>% sub("[+-]", "", .),
                paste0("../results/NK_compare/AFU.NK_Rep_best_noSign.nodeseq.", test, ".genelist.txt"))
}


# NK load version - for proof checking
nk_df <- readxl::read_xlsx("../other_data/Statistics.xlsx", sheet = "NK_Asp_VS_Asp_US")
nk_df <- nk_df[,c("id", "log2FC", "deseqpadj", "deseqpadj2", "res_limma$adjpvalues", "p_edgeR")]
cn <- c("id", "log2_rpkm", "DESeq_adj_pval", "DESeq2_adj_pval", "limma_adj_pval", "edgeR_adj_pval") # NOTE: fc is RPKM, not MRN. But I do not want to rewrite the stuff above...
nk_df <- lapply(list(nk_df), function(x) {x <- as.data.frame(x); colnames(x) <- cn; rownames(x) <- x$id; x})

nk_best <- get_sign_genes_nk_list(nk_df, with_fold = T)[[1]]
xf_best <- get_sign_genes_geo2(xf_deg, with_fold = T)
l <- list(DC = xf_best, NK = nk_best)
make_venn_plot(l, file = "../results/NK_compare/AFU.NK_noRep_besasdfrlap.pdf", main = "Overlap of top genes without fc sign", save = F)



