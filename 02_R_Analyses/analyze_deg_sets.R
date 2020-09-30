#!/usr/bin/env R
library("plyr")
library("tidyverse")
library("ggplot2")

get_sign_genes_geo2 <- function(deg_res, with_fold = F, sigP = 0.01, lfc = 1) {
    l <- lapply(deg_res, function(x) {
        # select tools
        tools <- grep("_adj_pval", x %>% colnames, value = T, fixed = T)
        print(tools)
        keep1 <- x[,tools] %>% {. <= sigP} %>% {rowSums(.) == length(tools)}
        keep2 <- abs(x$log2_fc_mrn) >= lfc
        keep <- keep1 & keep2
        return(x[keep,])
    })

    names(l) <- names(deg_res)
    return(l)
}



deg_res_sp1 <- readRDS("../results/DEG_overlaps/deg_res_sp1.rds")
best <- get_sign_genes_geo2(deg_res_sp1$DEGs, lfc = 1)
lapply(best, dim)

best <- get_sign_genes_geo2(deg_res_sp1$DEGs, lfc = 0)
lapply(best, dim)


################################################################
### DEG plots


deg_res_sp1 <- readRDS("../results/DEG_overlaps/deg_res_sp1.rds")
deg_res_sp2 <- readRDS("../results/DEG_overlaps/deg_res_sp2.rds")
deg_res_sp3 <- readRDS("../results/DEG_overlaps/deg_res_sp3.rds")
best_sp1 <- get_sign_genes_geo2(deg_res_sp1$DEGs, lfc = 1) %>% sapply(., nrow) %>% tibble(Comparison = names(.), n = ., ) %>% mutate(Species = "Human")
best_sp2 <- get_sign_genes_geo2(deg_res_sp2$DEGs, lfc = 1) %>% sapply(., nrow) %>% tibble(Comparison = names(.), n = ., ) %>% mutate(Species = "A. fumigatus")
best_sp3 <- get_sign_genes_geo2(deg_res_sp3$DEGs, lfc = 1) %>% sapply(., nrow) %>% tibble(Comparison = names(.), n = ., ) %>% mutate(Species = "CMV")
best <- rbind(best_sp1, best_sp2, best_sp3)

ggplot(best, aes(x = Comparison, y = n, fill = Species, group = Species)) +
    #facet_grid( ~ Species) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1.0, vjust = 0.5))


best_sp1 <- get_sign_genes_geo2(deg_res_sp1$DEGs, lfc = 1, sigP = 0.01) %>% sapply(., nrow) %>% tibble(Comparison = names(.), n = ., ) %>% mutate(Species = "Homo sapiens")
best_sp2 <- get_sign_genes_geo2(deg_res_sp2$DEGs, lfc = 0, sigP = 0.05) %>% sapply(., nrow) %>% tibble(Comparison = names(.), n = ., ) %>% mutate(Species = "Aspergillus fumigatus")
best_sp3 <- get_sign_genes_geo2(deg_res_sp3$DEGs, lfc = 0, sigP = 0.05) %>% sapply(., nrow) %>% tibble(Comparison = names(.), n = ., ) %>% mutate(Species = "Cytomegalovirus")
best <- rbind(best_sp1, best_sp2, best_sp3)
best <- mutate(best, Species = factor(Species, levels = c("Homo sapiens", "Aspergillus fumigatus", "Cytomegalovirus")))

# relabel
best$comp_bak <- best$Comparison
best$Comparison <- sub("SingleCulture_DC", "DC", best$Comparison, fixed = T)
best$Comparison <- sub("SingleCulture_Afu_0h", "Afu 0h", best$Comparison, fixed = T)
best$Comparison <- sub("SingleCulture_Afu_4h30min", "Afu 4.5h", best$Comparison, fixed = T)
best$Comparison <- sub("Co[iI]nfection_Afu_first", "DC + CMV 2h + Afu 0h", best$Comparison, fixed = F)
best$Comparison <- sub("CoInfection_CMV_first", "DC + CMV 0h + Afu 4.5h", best$Comparison, fixed = T)
best$Comparison <- sub("SingleInfection_Afu_0h", "DC + Afu 0h", best$Comparison, fixed = T)
best$Comparison <- sub("SingleInfection_Afu_4h30min", "DC + Afu 4.5h", best$Comparison, fixed = T)
best$Comparison <- sub("SingleInfection_CMV_0h", "DC + CMV 0h", best$Comparison, fixed = T)
best$Comparison <- sub("SingleInfection_CMV_2h", "DC + CMV 2h", best$Comparison, fixed = T)

best$Comparison <- sub("CoInfection", "DC + CMV 0h + Afu 0h", best$Comparison, fixed = T)
best$Comparison <- sub("coinfection", "DC + CMV + Afu", best$Comparison, fixed = T)
best$Comparison <- sub("singleinfection", "DC + {CMV | Afu}", best$Comparison, fixed = T)

best$Comparison <- sub("_vs_", "  vs  ", best$Comparison, fixed = T)
best$Comparison <- sub("_VS_", "  vs  ", best$Comparison, fixed = T)

best$Comparison <- factor(best$Comparison,
                          levels = c(
                              "DC + Afu 0h  vs  DC",
                              "DC + Afu 4.5h  vs  DC",
                              "DC + CMV 0h  vs  DC",
                              "DC + CMV 2h  vs  DC",
                              "DC + CMV 0h + Afu 0h  vs  DC",
                              "DC + CMV 0h + Afu 4.5h  vs  DC",
                              "DC + CMV 2h + Afu 0h  vs  DC",

                              "DC + Afu 0h  vs  Afu 0h",
                              "DC + Afu 4.5h  vs  Afu 4.5h",
                              "DC + CMV 0h + Afu 0h  vs  Afu 0h",
                              "DC + CMV 2h + Afu 0h  vs  Afu 0h",
                              "DC + CMV 0h + Afu 4.5h  vs  Afu 4.5h",

                              "DC + CMV 0h + Afu 0h  vs  DC + Afu 0h",
                              "DC + CMV 2h + Afu 0h  vs  DC + Afu 0h",
                              "DC + CMV 0h + Afu 4.5h  vs  DC + Afu 4.5h",
                              "DC + CMV 0h + Afu 0h  vs  DC + CMV 0h",
                              "DC + CMV 0h + Afu 4.5h  vs  DC + CMV 0h",
                              "DC + CMV 2h + Afu 0h  vs  DC + CMV 2h",

                              "DC + CMV 0h + Afu 4.5h  vs  DC + CMV 0h + Afu 0h",
                              "DC + CMV 2h + Afu 0h  vs  DC + CMV 0h + Afu 0h",
                              "DC + CMV 2h + Afu 0h  vs  DC + CMV 0h + Afu 4.5h",
                              "DC + CMV + Afu  vs  DC + {CMV | Afu}"
                          ),
                          labels = c(
                              "moDC + Afu 0h  vs  moDC",
                              "moDC + Afu 4.5h  vs  moDC",
                              "moDC + CMV 0h  vs  moDC",
                              "moDC + CMV 2h  vs  moDC",
                              "moDC + CMV 0h + Afu 0h  vs  moDC",
                              "moDC + CMV 0h + Afu 4.5h  vs  moDC",
                              "moDC + CMV 2h + Afu 0h  vs  moDC",

                              "moDC + Afu 0h  vs  Afu 0h",
                              "moDC + Afu 4.5h  vs  Afu 4.5h",
                              "moDC + CMV 0h + Afu 0h  vs  Afu 0h",
                              "moDC + CMV 2h + Afu 0h  vs  Afu 0h",
                              "moDC + CMV 0h + Afu 4.5h  vs  Afu 4.5h",

                              "moDC + CMV 0h + Afu 0h  vs  moDC + Afu 0h",
                              "moDC + CMV 2h + Afu 0h  vs  moDC + Afu 0h",
                              "moDC + CMV 0h + Afu 4.5h  vs  moDC + Afu 4.5h",
                              "moDC + CMV 0h + Afu 0h  vs  moDC + CMV 0h",
                              "moDC + CMV 0h + Afu 4.5h  vs  moDC + CMV 0h",
                              "moDC + CMV 2h + Afu 0h  vs  moDC + CMV 2h",

                              "moDC + CMV 0h + Afu 4.5h  vs  moDC + CMV 0h + Afu 0h",
                              "moDC + CMV 2h + Afu 0h  vs  moDC + CMV 0h + Afu 0h",
                              "moDC + CMV 2h + Afu 0h  vs  moDC + CMV 0h + Afu 4.5h",
                              "moDC + CMV + Afu  vs  moDC + {CMV | Afu}"
                          ))


g <- ggplot(best, aes(x = Comparison, y = n, group = Species)) +
    facet_wrap( ~ Species, scales = "free") +
    geom_bar(stat = "identity", fill = "gray", color = "black") +
    theme_bw() +
    scale_fill_discrete(guide = F) +
    theme(text = element_text(family = "sans"),
          axis.text.x = element_text(size = 8, angle = 90, hjust = 1.0, vjust = 0.5))
g
ggsave("../results/other_plots/num_deg_overview.pdf", dpi = 300, width = 9)

ggplot(best, aes(x = Comparison, y = n, fill = Species, group = Species)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1.0, vjust = 0.5))



################################################################
### Co-Infection specific genes in H. sapiens
## NOTES
###
# The analysis is composed of two parts
# Part 1 is here - a simple interesect of all DEG considereing both CoI vs SI tests
# Part 2 is an additional statistical test between all CoI vs all SI samples

source("scripts/lib_annotate.R")

sig_deg_sp1 <- get_sign_genes_geo2(deg_res_sp1$DEGs, lfc = 1, sigP = 0.01)
names(sig_deg_sp1)

of_interest <- sig_deg_sp1[c("CoInfection_vs_SingleInfection_Afu_0h", "CoInfection_vs_SingleInfection_CMV_0h")]
of_interest_genes <- llply(of_interest, function(x) as.character(x$id))
coi_genes <- intersect(of_interest_genes[[1]], of_interest_genes[[2]])

# convert to symbol for better interpretation
coi_genes_sym <- to_symbol(coi_genes)

# get original rows to compare
# Basically: A gene is coinfection specific if the fold change is positive (or, at least, has the same direction)
coi_genes_deg <- ldply(of_interest, .id = "comparison", function(df) {
    df[names(coi_genes_sym),] %>%
        mutate(symbol = coi_genes_sym) %>%
        dplyr::select(id, symbol, everything())
}) %>%
    arrange(symbol, comparison) %>%
    ddply("symbol", function(df) {
        df$CoI_specific <- (df$log2_fc_mrn[1] * df$log2_fc_mrn[2]) > 0
        df %>% dplyr::select(comparison, id, symbol, CoI_specific, everything())
    })


dir.create("../results/DEG_overlaps/CoI_specific_genes", showWarnings = F)
writexl::write_xlsx(coi_genes_deg, "../results/DEG_overlaps/CoI_specific_genes/CoI_specific_genes.xlsx")



################################################################
### check out old data

cmv_dat <- data.table::fread("../results/stats/Hsapiens/Hsapiens_CoInfection_vs_SingleInfection_Afu_0h.csv", sep = "\t", sep2 = ";")
cmv_dat <- as.data.frame(cmv_dat)
best <- get_sign_genes_geo2(list(cmv_dat), lfc = 0)
lapply(best, dim)

cmv_dat <- data.table::fread("../results/stats/Hsapiens/Hsapiens_CoInfection_vs_SingleCulture_DC.csv", sep = "\t", sep2 = ";")
cmv_dat <- as.data.frame(cmv_dat)
best <- get_sign_genes_geo2(list(cmv_dat), lfc = 0)
lapply(best, dim)

#sum((cmv_dat$DESeq2_adj_pval < 0.01) & (cmv_dat$DESeq_adj_pval < 0.01) & (cmv_dat$Limma_adj_pval < 0.01) & (cmv_dat$EdgeR_adj_pval < 0.01))

cmv_dat <- data.table::fread("../results/stats/Hsapiens/Hsapiens_CoInfection_vs_SingleInfection_CMV_0h.csv", sep = "\t", sep2 = ";")
cmv_dat <- as.data.frame(cmv_dat)
best <- get_sign_genes_geo2(list(cmv_dat), lfc = 0)
lapply(best, dim)

cmv_dat <- data.table::fread("../results/stats/Hsapiens/Hsapiens_CoInfection_vs_SingleInfection_Afu_0h.csv", sep = "\t", sep2 = ";")
cmv_dat <- as.data.frame(cmv_dat)
best <- get_sign_genes_geo2(list(cmv_dat), lfc = 0)
lapply(best, dim)

cmv_dat <- data.table::fread("../results/stats/Hsapiens/Hsapiens_CoInfection_vs_SingleInfection_CMV_0h.anno.csv", sep = "\t", sep2 = ";")
cmv_dat <- as.data.frame(cmv_dat)
best <- get_sign_genes_geo2(list(cmv_dat), lfc = 0)
lapply(best, dim)

