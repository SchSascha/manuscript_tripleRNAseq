# loading and plotting qPCR data vs respective RNA-seq data

library("tidyverse")
library("ggpubr")
library("ggsci")
library("reshape2")
library("cowplot")
library("here")

setwd(here())

DAT_FILE <- "dat/expr_qpcr.tsv"

#### functions #####################################################
#
# IN
#   dat <- data frame with all genes
#
get_stats <- function(dat) {
  res_df <- tibble(gene_symbol = NA, method = NA, condition1 = NA, condition2 = NA, value = NA)
  for ( g in levels(dat$gene_symbol) ) {
    tmp_df <- dat %>% filter(gene_symbol == g)
    tmp_df <- tmp_df[!is.na(tmp_df$value),]
    sampleCount <- table(tmp_df$condition)
    lvl2drop <- names(sampleCount[sampleCount <= 1])
    tmp_df <- tmp_df[!tmp_df$condition %in% lvl2drop,]
    tmp_res <- pairwise.t.test(x = tmp_df$maxNormVal, g = tmp_df$condition, p.adjust.method = "BH")
    
    tmp_res <- melt(tmp_res$p.value)
    
    tmp_res <- cbind( gene_symbol = g, method = dat$method[1], tmp_res )
    names(tmp_res) <- c("gene_symbol", "method", "condition1", "condition2", "value")
    
    res_df <- rbind( res_df, tmp_res )
  }
  
  res_df <- res_df[!is.na(res_df$value),]
  
  return(res_df)
}

#
## get_relativeNorm_vals
# PRE:
#   ALAS is reference measurement
# IN:
#   dat <- Donor specific dataset
#   ref_condition <- condition to be referenced to for plot normalisation
# OUT:
#   dat including further column "normVal" with normalised values to given ref_condition
get_relativeNorm_vals <- function(dat, ref_condition) {
  dat$normVal <- NA
  ALAS_ref <- dat %>% filter(gene_symbol == "ALAS", condition == ref_condition) %>% pull(value)
  for (i in 1:dim(dat)[1]) {
    ALAS_val <- dat %>% filter(gene_symbol == "ALAS", condition == dat$condition[i]) %>% pull(value)
    tmpGene_ref <- dat %>% 
      filter(gene_symbol == dat$gene_symbol[i], condition == ref_condition) %>% 
      pull(value)
    dat$normVal[i] <- 2^(tmpGene_ref - dat$value[i]) / 2^(ALAS_ref - ALAS_val)
  }
  return(dat)
}

# normalize relative to maximum value
# IN 
#   df
# OUT
#   df including additional column "maxNormVal"
get_max_norm <- function(dat, valueColName = "normVal") {
  dat$maxNormVal <- NA
  
  for (i in 1:dim(dat)[1]) {
    maxVal <- dat %>% filter(gene_symbol == dat$gene_symbol[i]) %>% pull(valueColName) %>% max(na.rm = T)
    if (maxVal == 0)
      dat$maxNormVal[i] <- dat %>% select(valueColName) %>% slice(i) %>% pull()
    else
      dat$maxNormVal[i] <- dat %>% select(valueColName) %>% slice(i) %>% pull() / maxVal
  }
  return(dat)
}

# ================================================================ #


#### get file ####
tbl <- read_tsv( DAT_FILE, col_types = c("fffffn") )
str(tbl)
tbl$method <- fct_relevel(tbl$method, "qpcr", after = 1)
summary(tbl)
table(tbl$gene_symbol)

tbl$gene_symbol <- tbl$gene_symbol %>% 
  fct_recode("IL-15" = "IL15",
             "IFN-\u03B3" = "IFNG")
levels(tbl$gene_symbol)
str(tbl)

# ================================================================ #


#### normalise data ####
# get mean for all expr data
# normalise to ALAS control for qPCR

# qpcr
tbl_qpcr <- tbl %>% filter(method == "qpcr")
summary(tbl_qpcr)

tbl_qpcr2_by <- by( data = tbl_qpcr, 
                    INDICES = tbl_qpcr$Donor,
                    FUN = get_relativeNorm_vals, "DC + (CMV, Afu) 0h"
)

tbl_qpcr2 <- rbind(tbl_qpcr2_by[[1]],tbl_qpcr2_by[[2]],tbl_qpcr2_by[[3]],tbl_qpcr2_by[[4]])
str(tbl_qpcr2)
summary(tbl_qpcr2)
tbl_qpcr_mean <- aggregate(normVal ~ gene_symbol + condition + organism + method, 
                           data = tbl_qpcr2, 
                           FUN = mean)
str(tbl_qpcr_mean)
tbl_qpcr_mean <- get_max_norm(tbl_qpcr_mean)


# expr
tbl_expr <- tbl %>% filter(method == "expr")

tbl_expr_mean <- aggregate(value ~ gene_symbol + condition + organism + method, 
                           data = tbl_expr, 
                           FUN = mean)
names(tbl_expr_mean) <- c(names(tbl_expr_mean)[1:4],"normVal") # TODO; here: fast and dirty, consider revising for cleaner code
tbl_expr_mean <- get_max_norm(tbl_expr_mean)
# sanity check - manually check
stopifnot(tbl_expr %>% filter(gene_symbol == "TLR3", condition == "DC + CMV 0h") %>% pull(value) %>% mean ==
            tbl_expr_mean %>% filter(gene_symbol == "TLR3", condition == "DC + CMV 0h") %>% pull(normVal)
)

# ================================================================ #


#### plotting ####

## prepare data
tbl_qpcr3_by <- by( data = tbl_qpcr2, 
                    INDICES = tbl_qpcr$Donor,
                    FUN = get_max_norm
)
tbl_qpcr3 <- rbind(tbl_qpcr3_by[[1]],tbl_qpcr3_by[[2]],tbl_qpcr3_by[[3]],tbl_qpcr3_by[[4]])

tbl_expr2_by <- by( data = tbl_expr, 
                    INDICES = tbl_expr$Donor,
                    FUN = get_max_norm, valueColName = "value"
)
tbl_expr2 <- rbind(tbl_expr2_by[[1]],tbl_expr2_by[[2]],tbl_expr2_by[[3]],tbl_expr2_by[[4]])


tbl_4_plot_2 <- rbind(tbl_qpcr3 %>% select(gene_symbol, organism, Donor, condition, method, value, maxNormVal),
                      tbl_expr2 %>% select(gene_symbol, organism, Donor, condition, method, value, maxNormVal)) %>%
  filter(organism == "af" | organism == "hs" | organism == "cmv")

tbl_4_plot_2$method <- fct_recode(tbl_4_plot_2$method, "RNA-seq" = "expr", qPCR = "qpcr")
str(tbl_4_plot_2)
summary(tbl_4_plot_2)
levels(tbl_4_plot_2$gene_symbol)
# DC -> moDC
levels(tbl_4_plot_2$condition)
tbl_4_plot_2$condition <- fct_recode(tbl_4_plot_2$condition,
                                     "moDC" = "DC",
                                     "moDC + CMV 0h" = "DC + CMV 0h",
                                     "moDC + Afu 0h" = "DC + Afu 0h",
                                     "moDC + (CMV, Afu) 0h" = "DC + (CMV, Afu) 0h"
)
tbl_4_plot_2$gene_symbol <- fct_drop(tbl_4_plot_2$gene_symbol)
summary(tbl_4_plot_2)
str(tbl_4_plot_2)
levels(tbl_4_plot_2$gene_symbol)
table(tbl_4_plot_2$gene_symbol)

## add statistics
# multiple correction per gene separately for RNA-seq and qPCR data
qpcr_stats <- get_stats( tbl_4_plot_2 %>% filter(method == "qPCR") )
rnaseq_stats <- get_stats( tbl_4_plot_2 %>% filter(method == "RNA-seq") )

expr_stats <- rbind(qpcr_stats, rnaseq_stats)

## get plots
p_stats <- list()
for (g in 1:length(levels(tbl_4_plot_2$gene_symbol))) {
  # get individual gene
  curr_gene <- levels(tbl_4_plot_2$gene_symbol)[g]
  pl_dat_tmp <- tbl_4_plot_2 %>% filter(gene_symbol == curr_gene)
  
  # get stats both for qPCR and RNA-seq
  curr_stat <- expr_stats %>% filter(gene_symbol == curr_gene & method == "qPCR")
  numTests <- dim(curr_stat)[1]
  names(curr_stat)
  names(curr_stat)[3:4] <- c("group1", "group2")
  curr_stat <- curr_stat %>% mutate(y.position = seq(1.1, 1.1+(numTests-1)*0.125,0.125))
  curr_stat$pStars <- NA
  curr_stat$pStars[curr_stat$value <= 0.001] <- "\u2217\u2217\u2217"
  curr_stat$pStars[curr_stat$value <= 0.01 & curr_stat$value > 0.001] <- "\u2217\u2217"
  curr_stat$pStars[curr_stat$value <= 0.05 & curr_stat$value > 0.01] <- "\u2217"
  curr_stat$pStars[curr_stat$value > 0.05] <- "n.s."
  
  curr_stat2 <- expr_stats %>% filter(gene_symbol == curr_gene & method == "RNA-seq")
  numTests <- dim(curr_stat2)[1]
  names(curr_stat2)
  names(curr_stat2)[3:4] <- c("group1", "group2")
  curr_stat2 <- curr_stat2 %>% 
    mutate(y.position = seq(max(curr_stat$y.position)+0.1,
                            max(curr_stat$y.position)+0.1+(numTests-1)*0.125,0.125)) 
  curr_stat2$pStars <- NA
  curr_stat2$pStars[curr_stat2$value <= 0.001] <- "\u2217\u2217\u2217"
  curr_stat2$pStars[curr_stat2$value <= 0.01 & curr_stat2$value > 0.001] <- "\u2217\u2217"
  curr_stat2$pStars[curr_stat2$value <= 0.05 & curr_stat2$value > 0.01] <- "\u2217"
  curr_stat2$pStars[curr_stat2$value > 0.05] <- "n.s." 
  
  curr_stat <- curr_stat[curr_stat$pStars != "n.s.",]
  curr_stat2 <- curr_stat2[curr_stat2$pStars != "n.s.",]
  
    p_stats[[g]] <- ggbarplot(pl_dat_tmp, 
                              x = "condition", 
                              y = "maxNormVal",
                              color = "method",
                              fill = "method",
                              facet.by = "gene_symbol",
                              position = position_dodge(0.8),
                              add = "mean_se"
    ) %>%
      ggpar(font.family = "Arial", ylim = c(0, 2.5)) +
      xlab("") +
      ylab("") +
      # ylab("relative RNA-seq and qPCR expression values") + # uncomment for y axis label
      theme_pubr() +
      rremove("legend") +
      rremove("x.text") + # comment out to get x-axis labels 
      # rotate_x_text(90) +
      grids(axis = "y", linetype = "solid") 

    if (dim(curr_stat)[1] > 0) {
      p_stats[[g]] <- p_stats[[g]] + 
      stat_pvalue_manual(curr_stat, label = "pStars",
                         position = position_nudge(0.2,0), 
                         tip.length = 0, label.size = 3,
                         color = "#00bfc4")
    }
    
    if (dim(curr_stat2)[1] > 0 ) {
      p_stats[[g]] <- p_stats[[g]] + 
      stat_pvalue_manual(curr_stat2, label = "pStars", 
                         position = position_nudge(-0.2,0), 
                         tip.length = 0, label.size = 3,
                         color = "#f8766d")
    }
  
    if (g <= 13) {
      p_stats[[g]] <- p_stats[[g]] +
        rremove("x.text") +
        rremove("x.ticks") +
        grids(axis = "y", linetype = "solid") 
    } else
      p_stats[[g]] <- p_stats[[g]] +
      rremove("x.text") + # comment out to get x-axis labels 
      # rotate_x_text(90) +
      grids(axis = "y", linetype = "solid")
}

# put plot together...
p_stats_final <- plot_grid(plotlist = p_stats, align = "hv", ncol = 5)
p_stats_final %>%
  ggsave( filename = "pl_expression.pdf", width = 10, height = 10,
            device = cairo_pdf)
