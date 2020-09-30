#!/usr/bin/env R

library("readxl")
library("tidyverse")

hs <- read_excel("../results/stats/Hsapiens_gene_biotype_stats.xlsx") %>%
    group_by(group) %>%
    summarise_if(is.numeric, sum) %>%
    mutate(Specie = "hs") %>%
    rename(type = group) %>%
    select(Specie, type, everything())
af <- read_excel("../results/stats/Afumigatus_genetype_stats.abs.xlsx") %>%
    mutate(Specie = "af")

hc <- read_excel("../results/stats/CMV_counts.xls") %>%
    mutate(type = "mRNA") %>%
    group_by(type) %>%
    summarise_if(is.numeric, sum) %>%
    mutate(Specie = "cmv")

# combine
dist.abs <- rbind(hs, af, hc)
writexl::write_xlsx(dist.abs, "../results/stats/read_distribution.abs.comb.xlsx")

dist.perc <- t(t(dist.abs[,-c(1:2)]) / colSums(dist.abs[,-c(1:2)])) * 100
dist.perc <- cbind(dist.abs[,1:2], dist.perc)
writexl::write_xlsx(dist.abs, "../results/stats/read_distribution.perc.comb.xlsx")

# group per specie
# HS accounts for 99% of reads in some infection samples
dist.perc %>%
    group_by(Specie) %>%
    summarise_if(is.numeric, sum) %>% View()
