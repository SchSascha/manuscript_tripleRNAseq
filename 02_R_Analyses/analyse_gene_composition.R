#!/usr/bin/env R

library("biomaRt")
library("tidyverse")
library("ggplot2")
library("scales")
library("ggforce") # for 'geom_arc_bar'
library("magrittr")


annotate <- function(df) {
    if (nrow(df) == 0)
        return(df)

    df$id <- rownames(df)
    ids <- df$id
    res <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'wikigene_description', 'gene_biotype'),
                 filters = 'ensembl_gene_id',
                 values = ids,
                 mart = ensembl)
    colnames(res)[1] <- "id"
    res_nodup <- res[!duplicated(res$id),]

    anno_df <- merge(df, res_nodup, by = "id", all.x = T)
    return(anno_df)
}





#######################
### Homo Sapiens
####

counts_sp1 <- read.csv("../results/stats/Hsapiens_counts.csv", sep = "\t")

hs_anno <- rtracklayer::readGFF("../other_data/Homo_sapiens.GRCh38.89.gtf")

# add gene function
listMarts()
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
attributes[1:5,]

unique(counts_sp1_anno$gene_biotype)

counts_sp1_sample_prop <- counts_sp1_anno %>% select(-id) %>% group_by(gene_biotype) %>% summarise_if(is.numeric, sum)

# They are too diverse & hardly visible in a pie chart
# Therefore, group them together
# order matters!
# see definitions on embl gene biotypes:
# use: https://www.gencodegenes.org/pages/biotypes.html
# either: http://vega.archive.ensembl.org/info/about/gene_and_transcript_types.html
remap_map_df <- counts_sp1_sample_prop$gene_biotype %>% data.frame("gene_biotype" = .)
remap_map_df$group[grepl("^rRNA", remap_map_df$gene_biotype)] <- "rRNA"
remap_map_df$group[remap_map_df$gene_biotype %in% c(
    "macro_lncRNA",
    "3prime_overlapping_ncRNA",
    "bidirectional_promoter_lncRNA",
    "sense_intronic",
    "sense_overlapping",
    "antisense",
    "non_coding",
    "lincRNA"
)] <- "lncRNA"

remap_map_df$group[remap_map_df$gene_biotype %in% c("miRNA")] <- "miRNA"
remap_map_df$group[remap_map_df$gene_biotype %in% c("Mt_rRNA", "Mt_tRNA")] <- "mitoRNA"
remap_map_df$group[remap_map_df$gene_biotype %in% c("protein_coding")] <- "mRNA"
remap_map_df$group[remap_map_df$gene_biotype %in% c("snoRNA")] <- "snoRNA"
remap_map_df$group[remap_map_df$gene_biotype %in% c("snRNA")] <- "snRNA"
remap_map_df$group[remap_map_df$gene_biotype %in% c("misc_RNA")] <- "miscRNA"
remap_map_df$group[grepl("pseudogene", remap_map_df$gene_biotype)] <- "pseudogene"
remap_map_df$group[is.na(remap_map_df$group)] <- "other"

# add this info to counting table
counts_sp1_sample_prop_grp <- left_join(counts_sp1_sample_prop, remap_map_df) %>% select(gene_biotype, group, everything())
writexl::write_xlsx(counts_sp1_sample_prop_grp, "../results/stats/Hsapiens_gene_biotype_stats.xlsx")

# calculate percentages
counts_sp1_sample_perc <- counts_sp1_sample_prop_grp %>% mutate_if(is.numeric, .funs = function(x) x/sum(x))
writexl::write_xlsx(counts_sp1_sample_perc, "../results/stats/Hsapiens_gene_biotype_stats.perc.xlsx")
counts_sp1_perc <- counts_sp1_sample_perc %>% select(-gene_biotype) %>% group_by(group) %>% summarise_if(is.numeric, sum)
writexl::write_xlsx(counts_sp1_perc, "../results/stats/Hsapiens_gene_biotype_grouped_stats.perc.xlsx")

# Convert to long for plotting
counts_sp1_perc_long <- gather(counts_sp1_perc, key = "Sample", value = "perc", -group)
# the mean of rows of tss columns results in the same invariance
counts_sp1_perc_long <- counts_sp1_perc_long %>% group_by(group) %>% summarise(value = mean(perc))
sum(counts_sp1_perc_long$value)
counts_sp1_perc_long <- mutate(counts_sp1_perc_long, Label = paste0(round(value * 100, digits = 2), "%"))
counts_sp1_perc_long$group <- factor(counts_sp1_perc_long$group, levels = c("mRNA",
                                                                            "pseudogene",
                                                                            "lncRNA",
                                                                            "miRNA",
                                                                            "miscRNA",
                                                                            "mitoRNA",
                                                                            "snoRNA",
                                                                            "snRNA",
                                                                            "rRNA",
                                                                            "other"))

blank_theme <- theme_minimal() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold")
    )

# design colors palette
custom_scale <- rep("#000000", length(counts_sp1_perc_long$group))

custom_scale[which(levels(counts_sp1_perc_long$group) == "other")] <- "#999999" # light gray
custom_scale[which(levels(counts_sp1_perc_long$group) == "rRNA")] <- "#666666" # light gray
custom_scale[which(!(levels(counts_sp1_perc_long$group) %in% c("rRNA", "other")))] <- RColorBrewer::brewer.pal(8, "Paired")



cols <- function(a) image(1:length(a), 1, as.matrix(1:length(a)), col = a, axes=FALSE , xlab="", ylab="")
cols(custom_scale)
cols(RColorBrewer::brewer.pal(12, "Paired"))




g <- ggplot(counts_sp1_perc_long, aes(x = "", y = value, fill = group)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    coord_polar(theta = "y") +
    theme(axis.text.x=element_blank()) +
    blank_theme +
    scale_y_continuous(breaks = cumsum(counts_sp1_perc_long$value) - counts_sp1_perc_long$value / 2, labels = counts_sp1_perc_long$Label)
g



df <- counts_sp1_perc_long %>%
    mutate(end = 2 * pi * cumsum(value)/sum(value),
           start = lag(end, default = 0),
           middle = 0.5 * (start + end),
           hjust = ifelse(middle > pi, 1, 0),
           vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))


g <- ggplot(df) + blank_theme +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                     start = start, end = end, fill = group)) +
    geom_text(aes(x = 1.05 * sin(middle), y = 1.05 * cos(middle), label = Label,
                  hjust = hjust, vjust = vjust)) +
    coord_fixed() +
    scale_x_continuous(limits = c(-1.5, 1.4),  # Adjust so labels are not cut off
                       name = "", breaks = NULL, labels = NULL) +
    scale_y_continuous(limits = c(-1.1, 1.1),      # Adjust so labels are not cut off
                       name = "", breaks = NULL, labels = NULL) +
    scale_fill_manual(values = custom_scale) +
    ggtitle("Human")
g
ggsave(g, filename = "../results/other_plots/Hsapiens_gene_stats.pie.pdf", dpi = 300)
writexl::write_xlsx(counts_sp1_perc_long, "../results/stats/Hsapiens_gene_biotype_grouped_long.xlsx")







#######################
### A. fumigatus
####
library("rtracklayer")

counts_sp2 <- read.csv("../results/stats/Afumigatus_novo_counts.csv", sep = "\t")

sp2_anno <- readGFF("../other_data/A_fumigatus_Af293_current_features.gff")
sp2_anno <- as.tibble(sp2_anno)
sp2_anno <- filter(sp2_anno, type != "chromosome")
sp2_anno <- filter(sp2_anno, type != "gene" & !is.na(Name))

counts_sp2 <- counts_sp2 %>% mutate(id = rownames(.)) %>% select(id, everything()) %>% as.tibble()
stopifnot(counts_sp2$id %in% sp2_anno$Name)

counts_sp2_anno <- merge(counts_sp2, sp2_anno %>% select(Name, type), all.x = T, by.x = "id", by.y = "Name")
counts_sp2_typed <- counts_sp2_anno %>% select(-id) %>% group_by(type) %>% summarise_if(is.numeric, sum)
writexl::write_xlsx(counts_sp2_typed, "../results/stats/Afumigatus_genetype_stats.abs.xlsx")

# calculate percentages
counts_sp2_sample_perc <- counts_sp2_typed %>% mutate_if(is.numeric, .funs = function(x) x/sum(x))
writexl::write_xlsx(counts_sp2_sample_perc, "../results/stats/Afumigatus_genetype_stats.perc.xlsx")

counts_sp2_sample_perc_long <- gather(counts_sp2_sample_perc, "sample", "value", -type)
counts_sp2_sample_perc_long <- counts_sp2_sample_perc_long %>% group_by(type) %>% summarise(value = mean(value))
counts_sp2_sample_perc_long <- mutate(counts_sp2_sample_perc_long, Label = paste0(round(value * 100, digits = 2), "%"))

df <- counts_sp2_sample_perc_long %>%
    mutate(end = 2 * pi * cumsum(value)/sum(value),
           start = lag(end, default = 0),
           middle = 0.5 * (start + end),
           hjust = ifelse(middle > pi, 1, 0),
           vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))


g <- ggplot(df) + blank_theme +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                     start = start, end = end, fill = type)) +
    geom_text(aes(x = 1.05 * sin(middle), y = 1.05 * cos(middle), label = Label,
                  hjust = hjust, vjust = vjust)) +
    coord_fixed() +
    scale_x_continuous(limits = c(-1.5, 1.4),  # Adjust so labels are not cut off
                       name = "", breaks = NULL, labels = NULL) +
    scale_y_continuous(limits = c(-1.1, 1.1),      # Adjust so labels are not cut off
                       name = "", breaks = NULL, labels = NULL) +
    scale_fill_manual(values = custom_scale) +
    ggtitle("A. fumigatus")
g
ggsave(g, filename = "../results/other_plots/Afumigatus_gene_stats.pie.all.pdf", dpi = 300)


# this time: only samples with AFU in it
counts_sp2_sample_perc_long <- gather(counts_sp2_sample_perc, "sample", "value", -type)
counts_sp2_sample_perc_long <- filter(counts_sp2_sample_perc_long, grepl("Afu", sample))
counts_sp2_sample_perc_long <- counts_sp2_sample_perc_long %>% group_by(type) %>% summarise(value = mean(value))
counts_sp2_sample_perc_long <- mutate(counts_sp2_sample_perc_long, Label = paste0(round(value * 100, digits = 2), "%"))

df <- counts_sp2_sample_perc_long %>%
    mutate(end = 2 * pi * cumsum(value)/sum(value),
           start = lag(end, default = 0),
           middle = 0.5 * (start + end),
           hjust = ifelse(middle > pi, 1, 0),
           vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))


g <- ggplot(df) + blank_theme +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                     start = start, end = end, fill = type)) +
    geom_text(aes(x = 1.05 * sin(middle), y = 1.05 * cos(middle), label = Label,
                  hjust = hjust, vjust = vjust)) +
    coord_fixed() +
    scale_x_continuous(limits = c(-1.5, 1.4),  # Adjust so labels are not cut off
                       name = "", breaks = NULL, labels = NULL) +
    scale_y_continuous(limits = c(-1.1, 1.1),      # Adjust so labels are not cut off
                       name = "", breaks = NULL, labels = NULL) +
    scale_fill_manual(values = custom_scale) +
    ggtitle("A. fumigatus")
g
ggsave(g, filename = "../results/other_plots/Afumigatus_gene_stats.pie.only_afu.pdf", dpi = 300)


# based on all
counts_sp2_sample_perc_long <- plyr::ddply(counts_sp2_typed, "type", function(x) sum(x[-1]))
colnames(counts_sp2_sample_perc_long) <- c("type", "value")
counts_sp2_sample_perc_long <- mutate(counts_sp2_sample_perc_long, value = value / sum(value))
counts_sp2_sample_perc_long <- mutate(counts_sp2_sample_perc_long, Label = paste0(round(value * 100, digits = 2), "%"))

df <- counts_sp2_sample_perc_long %>%
    mutate(end = 2 * pi * cumsum(value)/sum(value),
           start = lag(end, default = 0),
           middle = 0.5 * (start + end),
           hjust = ifelse(middle > pi, 1, 0),
           vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))


g <- ggplot(df) + blank_theme +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                     start = start, end = end, fill = type)) +
    geom_text(aes(x = 1.05 * sin(middle), y = 1.05 * cos(middle), label = Label,
                  hjust = hjust, vjust = vjust)) +
    coord_fixed() +
    scale_x_continuous(limits = c(-1.5, 1.4),  # Adjust so labels are not cut off
                       name = "", breaks = NULL, labels = NULL) +
    scale_y_continuous(limits = c(-1.1, 1.1),      # Adjust so labels are not cut off
                       name = "", breaks = NULL, labels = NULL) +
    scale_fill_manual(values = custom_scale) +
    ggtitle("A. fumigatus")
g
ggsave(g, filename = "../results/other_plots/Afumigatus_gene_stats.pie.accum.pdf", dpi = 300)

