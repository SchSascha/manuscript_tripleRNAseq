# Plot cytokine data

library("ggpubr")
library("readxl")
library("RColorBrewer")
library("cowplot")
library("tidyverse")
library("here")

setwd(here())

source("tripleSeq_funs.R")

#### functions #####################################################

#
# get_relative_normed to a reference abundance, e.g. COI
# IN:z
#   dat <- Donor specific dataset
#   ref <- condition to be referenced to for plot normalisation
# testing:
# DAT <- elisa
# DAT <- elisa_suppl
# REF = "DC + CMV 0h + Afu 0h"
get_relativeVals <- function(DAT, REF = "DC + CMV 0h + Afu 0h") {
  DAT$relVal <- NA
  vars  <- levels(DAT$variable)
  conds <- levels(DAT$condition)
  for ( v in 1:length( vars ) ) {
    ref <- DAT %>% filter( variable == vars[v] & condition == REF )
    # alternatively create the average reference already here; 
    # PRO: if ref is missing for one repl it's not critical
    ref$value <- ref %>% pull(value) %>% mean(na.rm = T)
    for (c in 1:length( conds ) ) {
      # ids <- DAT %>% filter( variable == vars[v] & condition == conds[c] )
      ids <- which( DAT$variable == vars[v] & DAT$condition == conds[c] )
      for (i in 1:4) {
        DAT$relVal[ids[i]] <- DAT$value[ids[i]] / ref$value[i]
      }
    }
    
  }
  return(DAT)
}

#
# IN
#   df - stats df with conditions still with numbers
# OUT
#   df
reorder_stats <- function(df) {
  res_df <- tibble()
  for (g in levels(df$variable)) {
    tmp_df <- df %>% filter(variable == g)
    testOrder <- tmp_df$condB %>% str_extract("^[0-9]") %>% as.numeric()
    testOrder[testOrder==8] <- 50
    testOrder[testOrder==5] <- 60
    testOrder[testOrder==6] <- 70
    testOrder[testOrder==7] <- 80
    testOrder[testOrder<10] <- testOrder[testOrder<10]*10
    testOrder2 <- tmp_df$condA %>% str_extract("^[0-9]") %>% as.numeric()
    testOrder2[testOrder2==8] <- 50
    testOrder2[testOrder2==5] <- 60
    testOrder2[testOrder2==6] <- 70
    testOrder2[testOrder2==7] <- 80
    testOrder2[testOrder2>10] <- testOrder2[testOrder2>10]/10
    testOrder <- testOrder + testOrder2
    tmp_df <- tmp_df[order(testOrder),]
    
    res_df <- rbind(res_df, tmp_df)
  }
  return(res_df)
  
}

# ================================================================ #


# ================================================================ #


#### get data ####
elisa <- read_tsv("plex.dat_extConditionname.csv", col_types = c("fffnf"))
str(elisa)

elisa.filtered <- elisa %>%
  filter(
    variable == "I-TAC (46)" | 
      variable == "IL-1beta (18)" |
      variable == "IFN-beta (30)" |
      variable == "IP-10 (22)" |
      variable == "IL-8 (27)" |
      variable == "IFN-gamma (43)" |
      variable == "TNF-alpha (45)" 
    )

elisa.filtered$variable <- droplevels(elisa.filtered$variable)
elisa.filtered$variable <- fct_recode(elisa.filtered$variable,
                                      "IFN-\u03b2" = "IFN-beta (30)",
                                      "IFN-\u03b3" = "IFN-gamma (43)",
                                      "IL-1\u03b2" = "IL-1beta (18)",
                                      "IL-8" = "IL-8 (27)",
                                      "CXCL10" = "IP-10 (22)",
                                      "CXCL11" = "I-TAC (46)",
                                      "TNF-\u03b1" = "TNF-alpha (45)" 
                                      )
elisa.filtered$variable <- fct_relevel(elisa.filtered$variable,
                                       "IFN-\u03b2",
                                       "IFN-\u03b3",
                                       "IL-1\u03b2",
                                       "IL-8",
                                       "CXCL10",
                                       "CXCL11",
                                       "TNF-\u03b1")
levels(elisa.filtered$variable)

elisa.filtered <- elisa.filtered %>% 
  filter(condition != "Afu 0h" & condition != "Afu 4.5h" )
elisa.filtered$condition <- droplevels(elisa.filtered$condition)
elisa.filtered$condition <- fct_relevel(elisa.filtered$condition, 
                                        "DC + CMV 2h", after = 4)
levels(elisa.filtered$condition)
str(elisa.filtered)

# ================================================================ #


#### prepare data ##################################################

elisa.filtered.rel <- get_relativeVals(elisa.filtered)

levels(elisa$condition2)
elisa_suppl <- elisa %>% filter(condition2 != "SC Afu")
summary(elisa)
summary(elisa_suppl)
elisa_suppl$condition2 <- droplevels(elisa_suppl$condition2)
levels(elisa_suppl$condition)
elisa_suppl$condition  <- droplevels(elisa_suppl$condition)
levels(elisa_suppl$variable)
elisa_suppl$variable <- fct_recode(elisa_suppl$variable,
                             "CXCL11" = "I-TAC (46)",
                             "IFN-\u03b1" = "IFN-alpha (48)",
                             "IFN-\u03b2" = "IFN-beta (30)",
                             "IFN-\u03b3" = "IFN-gamma (43)",
                             "IL-10" = "IL-10 (28)",
                             "IL-12p70" = "IL-12p70 (34)",
                             "IL-17A" = "IL-17A (36)",
                             "IL-1\u03b1" = "IL-1alpha (62)",
                             "IL-1\u03b2" = "IL-1beta (18)",
                             "IL-2" = "IL-2 (19)",
                             "IL-23" = "IL-23 (63)",
                             "IL-6" = "IL-6 (25)",
                             "IL-8" = "IL-8 (27)",
                             "CXCL10" = "IP-10 (22)",
                             "CCL5" = "RANTES (42)",
                             "TNF-\u03b1" = "TNF-alpha (45)" 
)
elisa_suppl$variable <- fct_relevel(elisa_suppl$variable,
                              "IL-1\u03b1",
                              "IL-1\u03b2",
                              "IL-2",
                              "IL-6",
                              "IL-8",
                              "IL-10",
                              "IL-12p70",
                              "IL-17A",
                              "IL-23",
                              "IFN-\u03b1",
                              "IFN-\u03b2",
                              "IFN-\u03b3",
                              "CXCL10",
                              "CXCL11",
                              "CCL5",
                              "TNF-\u03b1"
)
elisa_suppl$condition <- fct_relevel(elisa_suppl$condition, 
                                        "DC + CMV 2h", after = 4)
elisa_suppl$condition2 <- fct_recode(elisa_suppl$condition2, SC = "SC moDC")

# DC -> moDC
levels(elisa_suppl$condition)
elisa_suppl$condition <- fct_recode(elisa_suppl$condition,
                                    "moDC" = "DC",
                                    "moDC + Afu 0h" = "DC + Afu 0h",
                                    "moDC + Afu 4.5h" = "DC + Afu 4.5h",
                                    "moDC + CMV 0h" = "DC + CMV 0h",
                                    "moDC + CMV 2h" = "DC + CMV 2h",
                                    "moDC + CMV 0h + Afu 0h" = "DC + CMV 0h + Afu 0h",
                                    "moDC + CMV 0h + Afu 4.5h" = "DC + CMV 0h + Afu 4.5h",
                                    "moDC + CMV 2h + Afu 0h" = "DC + CMV 2h + Afu 0h"
)

## do relative plot for all cytokines
levels(elisa_suppl$condition)
elisa_suppl_rel <- get_relativeVals(DAT = elisa_suppl, REF = "moDC + CMV 0h + Afu 0h")

# ================================================================ #


#### load stats ####
stats_fn <- "stats.ttest.single.tsv"
stats_dat <- read_tsv(stats_fn, col_types = c("fffnnnn"))
str(stats_dat)

stats_dat$variable <- stats_dat$variable %>% fct_recode(
                                   "CXCL11" = "I-TAC (46)",
                                   "IFN-\u03b1" = "IFN-alpha (48)",
                                   "IFN-\u03b2" = "IFN-beta (30)",
                                   "IFN-\u03b3" = "IFN-gamma (43)",
                                   "IL-10" = "IL-10 (28)",
                                   "IL-12p70" = "IL-12p70 (34)",
                                   "IL-17A" = "IL-17A (36)",
                                   "IL-1\u03b1" = "IL-1alpha (62)",
                                   "IL-1\u03b2" = "IL-1beta (18)",
                                   "IL-2" = "IL-2 (19)",
                                   "IL-23" = "IL-23 (63)",
                                   "IL-6" = "IL-6 (25)",
                                   "IL-8" = "IL-8 (27)",
                                   "CXCL10" = "IP-10 (22)",
                                   "CCL5" = "RANTES (42)",
                                   "TNF-\u03b1" = "TNF-alpha (45)"
)

# reorder stats for visual order of stats bars
stats_dat <- reorder_stats(stats_dat)

levels(stats_dat$condA)
levels(stats_dat$condB)
levels(elisa_suppl_rel$condition)
stats_dat$condA <- stats_dat$condA %>% fct_recode(
  "moDC + Afu 0h" = "2. DC + Afu",
  "moDC + Afu 4.5h" = "3. DC + Afu",
  "moDC + CMV 0h" = "4. DC + CMV",
  "moDC + CMV 0h + Afu 0h" = "5. DC + Afu + CMV",
  "moDC + CMV 0h + Afu 4.5h" = "6. DC + Afu + CMV",
  "moDC + CMV 2h + Afu 0h" = "7. DC + Afu + CMV",
  "moDC + CMV 2h" = "8. DC + CMV"
)
stats_dat$condB <- stats_dat$condB %>% fct_recode(
  "moDC" = "1. DC alone",
  "moDC + Afu 0h" = "2. DC + Afu",
  "moDC + Afu 4.5h" = "3. DC + Afu",
  "moDC + CMV 0h" = "4. DC + CMV",
  "moDC + CMV 2h" = "8. DC + CMV"
)
names(stats_dat)[2:3] <- c("group1", "group2")

# ================================================================ #


#### plot #######################################

# plot for Fig.5 selection or for supplement (all analytes)?
# TO BE MODIFIED
ALL_ANALYTES <- T # F: only those that are in main manuscript, T: all measured analytes
WO_REDUNDANCY <- F # F: show all, T: show only analytes that are not shown in main manuscript 

# Which data set?
elisa_dat <- elisa_suppl_rel

# Figure 5 selection:
if (!ALL_ANALYTES) {
  elisa_dat <-
    elisa_suppl_rel %>% filter(
      elisa_suppl_rel$variable %in% c(
        "IFN-\u03b2",
        "IFN-\u03b3",
        "IL-1\u03b2",
        "IL-8",
        "CXCL10",
        "CXCL11",
        "TNF-\u03b1"
      )
    )
  elisa_dat$variable <- elisa_dat$variable %>% fct_drop()
}

# kick out redundancy if requested:
if (ALL_ANALYTES & WO_REDUNDANCY) {
  elisa_dat <-
    elisa_suppl_rel %>% filter(
      elisa_suppl_rel$variable %in% c(
        "IL-1α",
        "IL-2",
        "IL-6",
        "IL-10",
        "IL-12p70",
        "IL-17A",
        "IL-23",
        "IFN-α",
        "CCL5"
      )
    )
  elisa_dat$variable <- elisa_dat$variable %>% fct_drop()
}

p_stats <- list()
p_stats.woLabel <- list()
for (g in 1:length(levels(elisa_dat$variable))) {
  # get individual gene
  curr_gene <- levels(elisa_dat$variable)[g]
  pl_dat_tmp <- elisa_dat %>% filter(variable == curr_gene)
 
  p_stats[[g]] <- ggbarplot(pl_dat_tmp, x = "condition", 
                            y = "relVal",
                            facet.by = "variable",
                            color = "condition2",
                            fill = "condition2",
                            position = position_dodge(0.9),
                            add = c("mean_se")
  ) %>%
    ggpar(font.family = "Arial", legend.title = "") +
    xlab("") +
    ylab("") +
    theme_pubr() +
    rremove("legend")
  
  p_stats.woLabel[[g]] <- p_stats[[g]]
  
  if (ALL_ANALYTES) {
    if (WO_REDUNDANCY) {
      g_bottom_threshold <- 6
    } else {
      g_bottom_threshold <- 12
    }
  } else {
    g_bottom_threshold <- 4
  }
  
  # if (g <= g_bottom_threshold) {
  if (g <= g_bottom_threshold) { # for filtered cytos for Fig 5
    p_stats[[g]] <- p_stats[[g]] +
      rremove("x.text") +
      rremove("x.ticks") +
      grids(axis = "y", linetype = "solid") 
  } else
    p_stats[[g]] <- p_stats[[g]] +
    rremove("x.text") +
    rotate_x_text(90) +
    grids(axis = "y", linetype = "solid") +
    ylab("relative cytokine abundance")
  
  p_stats.woLabel[[g]] <- p_stats[[g]] +
    rremove("x.text") +
    rremove("x.ticks") +
    grids(axis = "y", linetype = "solid") 
  
  # get limit info
  tmp_pl_obj <- ggplot_build(p_stats[[g]])
  y_lim <- tmp_pl_obj$layout$panel_scales_y[[1]]$range$range[2]
  y_lim <- ceiling(y_lim*10)/10 
  # get stats
  curr_stat <- stats_dat %>% filter(variable == curr_gene)
  numTests <- dim(curr_stat)[1]
  curr_stat <- curr_stat %>% 
    mutate(y.position = seq(y_lim, y_lim+(numTests-1)*(ceiling(y_lim)/10),ceiling(y_lim)/10))
  curr_stat$pStars <- NA
  curr_stat$pStars[curr_stat$padj <= 0.001] <- "\u2217\u2217\u2217"
  curr_stat$pStars[curr_stat$padj <= 0.01 & curr_stat$padj > 0.001] <- "\u2217\u2217"
  curr_stat$pStars[curr_stat$padj <= 0.05 & curr_stat$padj > 0.01] <- "\u2217"
  curr_stat$pStars[curr_stat$padj <= 0.1 & curr_stat$padj > 0.05] <- "(\u2217)"
  curr_stat$pStars[curr_stat$padj > 0.1] <- "n.s."
  
  curr_stat <- curr_stat[curr_stat$pStars != "n.s.",]
  
  if (dim(curr_stat)[1] > 0) {
    p_stats[[g]] <- p_stats[[g]] + 
      stat_pvalue_manual(curr_stat, label = "pStars",
                         tip.length = 0, label.size = 2)
    p_stats.woLabel[[g]] <- p_stats.woLabel[[g]] + 
      stat_pvalue_manual(curr_stat, label = "pStars",
                         tip.length = 0, label.size = 2)
  }
  
}

# aspect ratio depending on number of data to plot
if (ALL_ANALYTES) {
  if (WO_REDUNDANCY) {
    p.width <- 11
    p.height <- 7.5
  } else {
    p.width <- 10
    p.height <- 10
  }
} else {
  p.width <- 11
  p.height <- 7.5
}

# labeling makes plots small, hence two versions w/ and w/o labels 
p_stats_final <- plot_grid(plotlist = p_stats, align = "hv", ncol = 3) 
p_stats_fina.woLabel <- plot_grid(plotlist = p_stats.woLabel, align = "hv", ncol = 4)

