# Description: This script loads necessary data and plots RNA types
#

library("tidyverse")
library("reshape2")
library("ggpubr")
library("ggsci")
library("scales")
library("xlsx")
library("RColorBrewer")
library("here")

setwd(here())

##### get file
df <- read_tsv( "mapping_stats.csv", col_types = c("fcnnnnnnnnnnnnnnnnnn")
  )

##### get file
df <- read_tsv(  "mapping_stats.csv", col_types = c("fcnnnnnnnnnnnnnnnnnn")
)

##### subtract relevant information for ggplot
### get percentage mapped reads
df_pct <- df[,c(1, 6:8)]
condition <- str_replace( df_pct$X1, "Donor[1-4]_[0-9]{2}_", "" )
str(df_pct)
condition <- factor(condition) %>% fct_inorder()
df_pct$X1 <- condition
names(df_pct)[1] <- "condition"
df_pct2 <- melt(df_pct)
df_pct2$variable <- factor(str_replace( df_pct2$variable, "percentage_reads_mapped_", "" ), 
                           levels = c( "Hsapiens", "Afumigatus", "CMV" ),
                           labels = c( "Homo sapiens", "Aspergillus fumigatus", "CMV" ))

#### reads mapped in exons ####
df_exons <- cbind(df_pct[1], df[,c(15:17)])
# scale sum of exon mapped reads to 100%
tmpSum <- rowSums(df_exons[,2:4])
df_exons[,2:4] <- df_exons[,2:4] / tmpSum * 100
rownames(df_exons) <- NULL
df_exons <- melt(df_exons)
df_exons$variable <- factor(str_replace( df_exons$variable, "percentage_reads_mapping_in_exons_", "" ), 
                            levels = c( "Hsapiens", "Afumigatus", "CMV" ) )

#### include type of RNA
# first, get all the mean values of exon mappings based on df_exons
mean_exon_coverage <- rbind(
  aggregate(value ~ condition, filter(df_exons,variable=="Hsapiens"), mean),
  aggregate(value ~ condition, filter(df_exons,variable=="Afumigatus"), mean),
  aggregate(value ~ condition, filter(df_exons,variable=="CMV"), mean)
)
stopifnot(df_exons %>% filter(variable=="Hsapiens" & condition == "DC_alone") %>% pull(value) %>% mean ==
            mean_exon_coverage$value[1])
stopifnot(df_exons %>% filter(condition == "DC_plus_CMV_2h") %>% pull(value) %>% sum == 400)

mean_exon_coverage <- cbind( c(replicate(10,"Hsapiens"), replicate(10,"Afumigatus"), replicate(10,"CMV")), 
                             mean_exon_coverage )
names(mean_exon_coverage) <- c("variable", "condition", "value")
stopifnot(mean_exon_coverage %>% filter(condition == "DC_alone") %>% pull(value) %>% sum == 100)
stopifnot(mean_exon_coverage %>% filter(condition == "DC_plus_CMV_2h") %>% pull(value) %>% sum == 100)

#### get percentages for type ####
df_RNAtypes_hsa.raw <- read.xlsx("Hsapiens_gene_biotype_grouped_stats.perc.xlsx", 1)
df_RNAtypes_afu.raw <- read.xlsx("Afumigatus_genetype_stats.perc.xlsx", 1)

df_RNAtypes_hsa <- df_RNAtypes_hsa.raw
colnames(df_RNAtypes_hsa) <- c("type", as.character(condition))
df_RNAtypes_hsa <- melt(df_RNAtypes_hsa)
df_RNAtypes_hsa <- cbind(organism = replicate(400, "HSapiens"), df_RNAtypes_hsa)
# sanity check
df_RNAtypes_hsa %>% filter(variable == "DC_plus_CMV_0h") %>% pull(value) %>% sum # should be 4


df_RNAtypes_afu <- df_RNAtypes_afu.raw
df_RNAtypes_afu$type <- fct_recode(df_RNAtypes_afu$type, other = "tRNA") 
colnames(df_RNAtypes_afu) <- c("type", as.character(condition))
df_RNAtypes_afu <- melt(df_RNAtypes_afu)
df_RNAtypes_afu <- cbind(organism = replicate(280, "Afumigatus"), df_RNAtypes_afu)
# sanity check
df_RNAtypes_afu %>% filter(variable == "DC_alone") %>% pull(value) %>% sum
df_RNAtypes_afu %>% filter(variable == "Afu_alone_0h") %>% pull(value) %>% sum

df_RNAtypes_cmv <- df_exons %>% filter(variable == "CMV")
df_RNAtypes_cmv <- cbind(organism = replicate(40, "CMV"), df_RNAtypes_cmv)
df_RNAtypes_cmv[3] <- df_RNAtypes_cmv[,2]
df_RNAtypes_cmv[2] <- factor(replicate(40, "mRNA"))
names(df_RNAtypes_cmv) <- names(df_RNAtypes_afu)
df_RNAtypes_cmv$value = 1
# sanity check -> should also be sum = 4 for 4 replicates
df_RNAtypes_cmv %>% filter(variable == "DC_alone") %>% pull(value) %>% sum
df_RNAtypes_cmv %>% filter(variable == "Afu_alone_0h") %>% pull(value) %>% sum


# modify relative coverage based on mean_exon_coverage
for (r in 1:4) {
  for (i in 1:10) {
    df_RNAtypes_hsa$value[(r-1)*100+(i-1)*10+c(1:10)] <-
      df_RNAtypes_hsa$value[(r-1)*100+(i-1)*10+c(1:10)]*mean_exon_coverage$value[i]
    
    df_RNAtypes_afu$value[(r-1)*70+(i-1)*7+c(1:7)] <- 
      df_RNAtypes_afu$value[(r-1)*70+(i-1)*7+c(1:7)]*mean_exon_coverage$value[i+10]
    
    df_RNAtypes_cmv$value[(r-1)*10+(i-1)+1] <-
      df_RNAtypes_cmv$value[(r-1)*10+(i-1)+1]*mean_exon_coverage$value[i+20]
  }
}

df_exons_all_types <- rbind(df_RNAtypes_hsa, df_RNAtypes_afu, df_RNAtypes_cmv)
names(df_exons_all_types) <- c("variable", "type", "condition", "value")
df_exons_all_types$type
df_exons_all_types$type <- fct_relevel(df_exons_all_types$type, c("pseudogene","miscRNA", "other"), after = 11)
df_exons_all_types$type <- fct_relevel(df_exons_all_types$type, c("mitoRNA","miRNA", "lncRNA"), after = 2)

df_exons_all_types$condition <- df_exons_all_types$condition %>%
  fct_recode(
    "moDC" = "DC_alone",
    "moDC + Afu 0h" = "DC_plus_Afu_0h",
    "moDC + Afu 4.5h" = "DC_plus_Afu_4h30min",
    "moDC + CMV 0h" = "DC_plus_CMV_0h",
    "moDC + CMV 2h" = "DC_plus_CMV_2h",
    "moDC + CMV 0h + Afu 0h" = "DC_plus_CMV_0h_plus_Afu_0h",
    "moDC + CMV 0h + Afu 4.5h" = "DC_plus_CMV_0h_plus_Afu_4h30min",
    "moDC + CMV 2h + Afu 0h" = "DC_plus_Afu_0h_plus_CMV_2h",
    "Afu 0h" = "Afu_alone_0h",
    "Afu 4.5h" = "Afu_alone_4h30min"
  ) %>%
  fct_relevel("moDC + CMV 2h", after = 4)


#### plot ####
p_RNAtypes <- ggbarplot(df_exons_all_types, x = "condition", 
                        y = "value",
                        facet.by = "variable",
                        color = "type", fill = "type"
                        , palette = "Paired"
                        , panel.labs = list(variable = 
                                              c("Homo sapiens", "Aspergillus fumigatus", "Cytomegalovirus"))
                        ,legend = "right",
                        # combine = T, 
                        # merge = "flip",
                        add = c("mean"),
                        xlab = F
) + 
  ylab("Reads mapped in exons (in %)") +
  rotate_x_text(45) +
  grids(axis = "y", linetype = "solid") +
  theme(
    # strip.text.x = element_text(
    # size = 12, 
    # color = "red", 
    # family = "Arial",
    # face = c("italic")
    # ) 
  )
p_RNAtypes
p_RNAtypes <- p_RNAtypes +
  # rremove("xlab") +
  rremove("legend.title") %>%
  ggpar(legend.title = "", font.family = "Arial") +
  grids(axis = "y", linetype = "solid")
# pdf
p_RNAtypes %>%
  ggexport(filename = paste0(PLOT_PATH, "Fig2_vi.pdf"), width = 9, height = 6) # width = 12, height = 7.5)
# png
p_RNAtypes %>%
  ggexport(filename = paste0(PLOT_PATH, "Fig2_iii.png"), width = 2400, height = 1500, res = 250)

# sanity check
sum(df_exons_all_types %>% filter(condition=="DC + Afu 0h") %>% pull(value))/4
sum(df_exons_all_types %>% filter(condition=="DC + Afu 4.5h") %>% pull(value))/4
sum(df_exons_all_types %>% filter(condition=="Afu 0h") %>% pull(value))/4
# <- all fine
sum(df_exons_all_types %>% filter(condition=="DC + CMV 0h + Afu 0h") %>% pull(value))/4
sum(df_exons_all_types %>% filter(condition=="DC + CMV 0h") %>% pull(value))/4

# only CMV
p_exon_CMV <- filter(df_exons_all_types, variable == "CMV") %>% 
  ggbarplot(x = "condition", 
            y = "value",
            facet.by = "variable",
            # color = "variable",
            panel.labs = list(variable = c("Homo sapiens", "Aspergillus fumigatus", "Cytomegalovirus")),
            # combine = T, 
            # merge = "flip",
            fill = '#a6cee3',
            color = '#a6cee3',
            add = c("mean"),
            lab.size = 16
  ) + 
  ylab("") +
  xlab("") +
  grids(axis = "y", linetype = "solid")
p_exon_CMV +
  rotate_x_text(45)
# pdf
(p_exon_CMV + rremove("x.text") +
    # rremove("x.ticks") +
    grids(axis = "y", linetype = "solid")) %>%
  ggexport(filename = "Fig2_CMV.pdf", width = 5, height = 4)
# png
(p_exon_CMV + rremove("x.text")) %>%
  ggexport(filename = "Fig2_CMV.png", width = 1000, height = 800, res = 300)
