#!/usr/bin/env R
# Geo2RNAseq-0.100.1
#
#library("Geo2RNAseq")
source("R/Plotting.R")
library("ggplot2")

load("R_workspace.RData")

# new labels
conds_sp1_bak <- conds_sp1
conds_sp1 <- conds_sp1_bak
conds_sp1 <- sub("SingleCulture_DC", "DC", conds_sp1, fixed = T)
conds_sp1 <- sub("SingleInfection_Afu_0h", "DC + Afu 0h", conds_sp1, fixed = T)
conds_sp1 <- sub("SingleInfection_Afu_4h30min", "DC + Afu 4.5h", conds_sp1, fixed = T)
conds_sp1 <- sub("SingleInfection_CMV_0h", "DC + CMV 0h", conds_sp1, fixed = T)
conds_sp1 <- sub("SingleInfection_CMV_2h", "DC + CMV 2h", conds_sp1, fixed = T)
conds_sp1 <- sub("CoInfection_Afu_first", "DC + CMV 2h + Afu 0h", conds_sp1, fixed = T)
conds_sp1 <- sub("CoInfection_CMV_first", "DC + CMV 0h + Afu 4.5h", conds_sp1, fixed = T)
conds_sp1 <- sub("CoInfection", "DC + CMV 0h + Afu 0h", conds_sp1, fixed = T)
conds_sp1 <-
    factor(
        conds_sp1,
        levels = c(
            "DC",
            "DC + Afu 0h",
            "DC + Afu 4.5h",
            "DC + CMV 0h",
            "DC + CMV 2h",
            "DC + CMV 0h + Afu 0h",
            "DC + CMV 0h + Afu 4.5h",
            "DC + CMV 2h + Afu 0h",
            "none"
        )
    )


conds_sp2_bak <- conds_sp2
conds_sp2 <- conds_sp2_bak
conds_sp2 <- sub("SingleCulture_Afu_0h", "Afu 0h", conds_sp2, fixed = T)
conds_sp2 <- sub("SingleCulture_Afu_4h30min", "Afu 4.5h", conds_sp2, fixed = T)
conds_sp2 <- sub("SingleInfection_Afu_0h", "DC + Afu 0h", conds_sp2, fixed = T)
conds_sp2 <- sub("SingleInfection_Afu_4h30min", "DC + Afu 4.5h", conds_sp2, fixed = T)
conds_sp2 <- sub("CoInfection_Afu_first", "DC + CMV 2h + Afu 0h", conds_sp2, fixed = T)
conds_sp2 <- sub("CoInfection_CMV_first", "DC + CMV 0h + Afu 4.5h", conds_sp2, fixed = T)
conds_sp2 <- sub("CoInfection", "DC + CMV 0h + Afu 0h", conds_sp2, fixed = T)
conds_sp2 <-
    factor(
        conds_sp2,
        levels = c(
            "Afu 0h",
            "Afu 4.5h",
            "DC + Afu 0h",
            "DC + Afu 4.5h",
            "DC + CMV 0h + Afu 0h",
            "DC + CMV 0h + Afu 4.5h",
            "DC + CMV 2h + Afu 0h",
            "none"
        )
    )


conds_sp3_bak <- conds_sp3
conds_sp3 <- conds_sp3_bak
conds_sp3 <- sub("SingleInfection_CMV_0h", "DC + CMV 0h", conds_sp3, fixed = T)
conds_sp3 <- sub("SingleInfection_CMV_2h", "DC + CMV 2h", conds_sp3, fixed = T)
conds_sp3 <- sub("Coinfection_Afu_first", "DC + CMV 2h + Afu 0h", conds_sp3, fixed = T)
conds_sp3 <- sub("CoInfection_CMV_first", "DC + CMV 0h + Afu 4.5h", conds_sp3, fixed = T)
conds_sp3 <- sub("CoInfection", "DC + CMV 0h + Afu 0h", conds_sp3, fixed = T)

## best treat related

pca_res_sp1_tmp <- make_PCA_plot(
    file         = NULL, #file.path("../results/PCA", name_sp1, paste0(name_sp1, "_PCA_mrn")),
    counts       = mrn_sp1[, conds_sp1 != "none"],
    designMatrix = NA,
    conds        = conds_sp1[conds_sp1 != "none"],
    shapes       = shapes_sp1[conds_sp1 != "none"],
    norm         = "none",
    main         = paste0("PCA - ", name_sp1),
    colName      = "",
    shapeName    = "",
    overwrite    = FORCE_OVERWRITE
)
g <- pca_res_sp1_tmp$g + theme_bw() + ggplot2::theme(aspect.ratio=1)
ggsave(filename = file.path("../results/PCA", paste0(name_sp1, ".nolog.mrn.pdf")),
       plot = g,
       dpi = 300)

pca_res_sp1_tmp <- make_PCA_plot(
    file         = NULL, #file.path("../results/PCA", name_sp1, paste0(name_sp1, "_PCA_mrn")),
    counts       = counts_sp1[, conds_sp1 != "none"] + 1,
    designMatrix = NA,
    conds        = conds_sp1[conds_sp1 != "none"],
    shapes       = shapes_sp1[conds_sp1 != "none"],
    norm         = "none",
    main         = paste0("PCA - ", name_sp1),
    colName      = "",
    shapeName    = "",
    overwrite    = FORCE_OVERWRITE
)
g <- pca_res_sp1_tmp$g + theme_bw() + ggplot2::theme(aspect.ratio=1)
ggsave(filename = file.path("../results/PCA", paste0(name_sp1, ".nolog.pdf")),
       plot = g,
       dpi = 300)

### best treat related end


pca_res_sp1 <- make_PCA_plot(
    file         = NULL, #file.path("../results/PCA", name_sp1, paste0(name_sp1, "_PCA_mrn")),
    counts       = log2(mrn_sp1[, conds_sp1 != "none"]+1),
    designMatrix = NA,
    conds        = conds_sp1[conds_sp1 != "none"],
    shapes       = shapes_sp1[conds_sp1 != "none"],
    norm         = "none",
    main         = paste0("PCA - ", name_sp1),
    colName      = "",
    shapeName    = "",
    overwrite    = FORCE_OVERWRITE
)
g <- pca_res_sp1$g + theme_bw() + ggplot2::theme(aspect.ratio=1)
ggsave(filename = file.path("../results/PCA", paste0(name_sp1, ".mrn.pdf")),
       plot = g,
       dpi = 300)

pca_res_sp2 <- make_PCA_plot(
    file         = NULL,
    counts       = log2(mrn_sp2[, conds_sp2 != "none"]+1),
    designMatrix = NA,
    conds        = conds_sp2[conds_sp2 != "none"],
    shapes       = shapes_sp2[conds_sp2 != "none"],
    norm         = "none",
    main         = paste0("PCA - ", name_sp2),
    colName      = "",
    shapeName    = "",
    overwrite    = FORCE_OVERWRITE
)
g <- pca_res_sp2$g + theme_bw() + ggplot2::theme(aspect.ratio=1)
ggsave(filename = file.path("../results/PCA", paste0(name_sp2, ".mrn.pdf")),
       plot = g,
       dpi = 300)

pca_res_sp3 <- make_PCA_plot(
    file         = NULL,
    counts       = log2(mrn_sp3[, conds_sp3 != "none"]+1),
    designMatrix = NA,
    conds        = conds_sp3[conds_sp3 != "none"],
    shapes       = shapes_sp3[conds_sp3 != "none"],
    norm         = "none",
    main         = paste0("PCA - ", name_sp3),
    colName      = "",
    shapeName    = "",
    overwrite    = FORCE_OVERWRITE
)
g <- pca_res_sp3$g + theme_bw() + ggplot2::theme(aspect.ratio=1)
ggsave(filename = file.path("../results/PCA", paste0(name_sp3, ".mrn.pdf")),
       plot = g,
       dpi = 300)



# custom PCA plots
# DC - DC vs single vs coinfection
nconds <- rep("none", length(conds_sp1_bak))
nconds[grep("DC", conds_sp1_bak)] <- "DC"
nconds[grep("SingleInfection_Afu", conds_sp1_bak)] <- "DC + Afu"
nconds[grep("SingleInfection_CMV", conds_sp1_bak)] <- "DC + CMV"
nconds[grep("CoInfection", conds_sp1_bak)] <- "DC + Afu + CMV"
nconds <- factor(nconds, levels = c("none", "DC", "DC + Afu", "DC + CMV", "DC + Afu + CMV"))

pca_res_sp4 <- make_PCA_plot(
    file         = NULL,# file.path("../results/other_plots/", paste0(name_sp1, "_PCA_GroupedInfection_mrn")),
    counts       = log2(mrn_sp1[, conds_sp1_bak != "none"]+1),
    designMatrix = NA,
    conds        = nconds[nconds != "none"],
    shapes       = shapes_sp1[conds_sp1_bak != "none"],
    norm         = "none",
    main         = paste0("PCA - ", name_sp1),
    colName      = "",
    shapeName    = "",
    overwrite    = FORCE_OVERWRITE
)
g <- pca_res_sp4$g + theme_bw() + ggplot2::theme(aspect.ratio=1)
ggsave(filename = file.path("../results/PCA", paste0(name_sp1, ".grouped.mrn.pdf")),
       plot = g,
       dpi = 300)


# AF - group infection
nconds <- conds_sp2_bak
nconds[grep("Infection", conds_sp2_bak)] <- "DC + Afu +/- CMV"
nconds[grep("SingleCulture_Afu_0h", conds_sp2_bak)] <- "Afu 0h"
nconds[grep("SingleCulture_Afu_4h30min", conds_sp2_bak)] <- "Afu 4.5h"
nconds <- factor(nconds, levels = c("none", "Afu 0h", "Afu 4.5h", "DC + Afu +/- CMV"))

pca_res_sp5 <- make_PCA_plot(
    file         = NULL, #file.path("../results/other_plots/", paste0(name_sp2, "_PCA_GroupedInfection_mrn")),
    counts       = log2(mrn_sp2[, conds_sp2 != "none"]+1),
    designMatrix = NA,
    conds        = nconds[nconds != "none"],
    shapes       = shapes_sp2[nconds != "none"],
    norm         = "none",
    main         = paste0("PCA - ", name_sp2),
    colName      = "",
    shapeName    = "",
    overwrite    = FORCE_OVERWRITE
)
g <- pca_res_sp5$g + theme_bw() + ggplot2::theme(aspect.ratio=1)
ggsave(filename = file.path("../results/PCA", paste0(name_sp2, ".grouped.mrn.pdf")),
       plot = g,
       dpi = 300)

# CMV - group single vs coinfection
nconds <- conds_sp3_bak
nconds[grep("Single", conds_sp3_bak)] <- "DC + CMV"
nconds[grep("Coinfection", conds_sp3_bak, ignore.case = T)] <- "DC + Afu + CMV"
nconds <- factor(nconds, levels = c("none", "DC + CMV", "DC + Afu + CMV"))

pca_res_sp6 <- make_PCA_plot(
    file         = NULL, #file.path("../results/other_plots/", paste0(name_sp3, "_PCA_GroupedInfection_mrn")),
    counts       = log2(mrn_sp3[, conds_sp3 != "none"]+1),
    designMatrix = NA,
    conds        = nconds[nconds != "none"],
    shapes       = shapes_sp3[nconds != "none"],
    norm         = "none",
    main         = paste0("PCA - log2 MRN + 1 - ", name_sp3),
    colName      = "",
    shapeName    = "",
    overwrite    = FORCE_OVERWRITE
)
g <- pca_res_sp6$g + theme_bw() + ggplot2::theme(aspect.ratio=1)
ggsave(filename = file.path("../results/PCA", paste0(name_sp3, ".grouped.mrn.pdf")),
       plot = g,
       dpi = 300)

pca_res_sp8 <- make_PCA_plot(
    file         = NULL, #file.path(statsDir, name_sp3, paste0(name_sp3, "_PCA_none")),
    counts       = counts_sp3[, conds_sp3_bak != "none"],
    designMatrix = NA,
    conds        = conds_sp3_bak[conds_sp3_bak != "none"],
    shapes       = shapes_sp3[conds_sp3_bak != "none"],
    norm         = "none",
    main         = paste0("PCA - none - ", name_sp3),
    colName      = "",
    shapeName    = "",
    overwrite    = FORCE_OVERWRITE
)

