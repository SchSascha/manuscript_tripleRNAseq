### Content: generate heatmap for selected genes
#

library("tidyverse")
library("pheatmap")
library("RColorBrewer")
library("here")

setwd(here())

# load data and gene symbols
dat <- read_tsv("./dat/heatmap_genes_expr_mrn.csv")
dat_meta <- read_tsv("./dat/heatmap_genes_symbols.csv")

# wrangle data sets
dat_long <- melt(dat)
dat_long$variable <- str_replace( dat_long$variable, "Donor[1-4]_[0-9]{2}_", "" )
dat_long$variable <- factor(dat_long$variable)
dat_long$X1 <- factor(dat_long$X1)

dat2 <- dat_long %>% filter(variable == "DC_alone" | 
                      variable == "DC_plus_Afu_0h" |
                      variable == "DC_plus_CMV_0h" |
                      variable == "DC_plus_CMV_0h_plus_Afu_0h"
                    )
dat2$variable <- droplevels(dat2$variable)

str(dat2)
levels(dat2$variable)
dat2$variable <- factor(dat2$variable,
                        levels = levels(dat2$variable),
                        labels = c("DC", "DC + Afu 0h", "DC + CMV 0h", "DC + (CMV, Afu) 0h"))

# heatmap of mean values
hm_dat1 <- dat2 %>% filter(variable == "DC")
hm_dat1.mean <- tapply(hm_dat1$value, hm_dat1$X1,mean)
hm_dat2 <- dat2 %>% filter(variable == "DC + Afu 0h")
hm_dat2.mean <- tapply(hm_dat2$value, hm_dat2$X1,mean)
hm_dat3 <- dat2 %>% filter(variable == "DC + CMV 0h")
hm_dat3.mean <- tapply(hm_dat3$value, hm_dat3$X1,mean)
hm_dat4 <- dat2 %>% filter(variable == "DC + (CMV, Afu) 0h")
hm_dat4.mean <- tapply(hm_dat4$value, hm_dat4$X1,mean)

hm_dat <- as.matrix(cbind(hm_dat1.mean, hm_dat2.mean, hm_dat3.mean, hm_dat4.mean))
row.names(hm_dat) <- dat_meta$`figure symbol`[match(row.names(hm_dat), dat_meta$Ensembl)]

# heatmap
pdf("heatmap_selectedGenes.pdf", width = 4, height = 9)
pheatmap(hm_dat
         , scale = "row"
         , cluster_cols = F
         , border_color = "white"
         , cellwidth = 15
         , cellheight = 15
         , labels_col = c("DC alone", "DC + Afu 0h", "DC + CMV 0h", "DC + (Afu, CMV) 0h")
)
dev.off()

