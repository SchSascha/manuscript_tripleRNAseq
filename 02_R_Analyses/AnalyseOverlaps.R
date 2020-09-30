#!R
library("plyr")
library("tidyverse")
library("magrittr")
library("VennDiagram")
library("ggplot2")
library("UpSetR")
library("gplots") # to calculate venn overlap
library("futile.logger")

rds_dir <- "../results/DEG_overlaps/"
deg_res_sp1 <- readRDS(file.path(rds_dir, "deg_res_sp1.rds"))
deg_res_sp2 <- readRDS(file.path(rds_dir, "deg_res_sp2.rds"))
deg_res_sp2_nod <- readRDS(file.path(rds_dir, "deg_res_sp2_nodeseq.rds"))
deg_res_sp2_nod <- deg_res_sp2_nod[[1]]
deg_res_sp3 <- readRDS(file.path(rds_dir, "deg_res_sp3.rds"))
deg_res_sp1 <- deg_res_sp1$DEGs
deg_res_sp2 <- deg_res_sp2$DEGs
deg_res_sp3 <- deg_res_sp3$DEGs


get_sign_genes_geo2 <- function(deg_res, with_fold = F, sigP = 0.01) {
    l <- lapply(deg_res, function(x) {
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
    })

    names(l) <- names(deg_res)
    return(l)
}


#' Make UpSet aka Intersection Bar Plots
#'
#' @param main Title of figure
#' @param l names list. Each named entry shall be a set. The intersections between these sets are calculated within this function.
#' @param file File to the save the plot to.
#' @param save If T, the plot will be saved to 'file'
#' @param nintersections The maximum number of intersections (bars) to display
#' @param colorless If F, the bars will be colored gray
make_upset_plot <- function(l, file = "", save = F, main = "UpSet", nintersections = 20, font.size = 2, number.angles = -30, colorless = F) {
    upSetFrame <- UpSetR::fromList(l)

    if (save) {
        if (file == "")
            stop("Please supply file name if you want to save the figure!")
        if (!grepl(".pdf$", file)) file <- paste0(file, ".pdf")
        pdf(file, onefile = FALSE, width = 11, height = 7)
    }

    UpSetR::upset(
        upSetFrame,
        order.by="freq",
        nsets = length(names(l)),
        nintersects = nintersections,
        scale.sets = "identity",
        text.scale = c(2, 1.5, 2, 1.5, 1.25, font.size),
        number.angles = number.angles,
        main.bar.color = if (colorless) "gray" else "blue"
    )
    grid.text(main,x = 0.65, y=0.95, gp=gpar(fontsize=20))

    if (save)
        dev.off()
}


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


get_items_per_intersection <- function(l) attr(venn(l, show.plot = F, ), "intersections")

rem_sign <- function(s) sub("[-|+]$", "", s)


#### Part 1 - Intersections of DEG genes (only names)


with_fold_change_direction <- F
if (with_fold_change_direction) {
    # consider fold change
    degs_sp1 <- get_sign_genes_geo2(deg_res_sp1, with_fold = T)
    degs_sp2 <- get_sign_genes_geo2(deg_res_sp2, with_fold = T, sigP = 0.05)
    degs_sp2_nod <- get_sign_genes_geo2(deg_res_sp2_nod, with_fold = T, sigP = 0.05)
    degs_sp3 <- get_sign_genes_geo2(deg_res_sp3, with_fold = T, sigP = 0.05)
} else {
    # ignore fold change
    degs_sp1 <- get_sign_genes_geo2(deg_res_sp1)
    degs_sp2 <- get_sign_genes_geo2(deg_res_sp2, sigP = 0.05)
    degs_sp2_nod <- get_sign_genes_geo2(deg_res_sp2_nod, sigP = 0.05)
    degs_sp3 <- get_sign_genes_geo2(deg_res_sp3, sigP = 0.05)
}



# sanity checks - these numbers are verified
stopifnot( sapply(degs_sp1, length) == c(241, 523, 500, 423, 556, 208, 1, 0, 2, 1269, 827, 1363, 953, 1705, 2156, 1860) )
stopifnot( sapply(degs_sp2, length) == c(1, 0, 1, 1, 0, 0, 280, 233, 215, 166, 186) )
stopifnot( sapply(degs_sp3, length) == c(1, 2, 0, 0, 0, 0, 11) )

save = F

# DC cells
## all against 'DC alone' - 'SingleCulture_DC'
### single infection
inds <- grepl('SingleCulture_DC', names(degs_sp1)) & grepl('SingleInfection', names(degs_sp1), ignore.case = T)
l <- degs_sp1[inds]
names(l) <- l %>% names %>% sub("_vs_SingleCulture_DC", "", .)
make_venn(l, file = "DC-SingleCulture_DC-SingleInfection.pdf", save = save, main = "DC alone vs SingleInfection")
saveRDS(get_items_per_intersection(l), "DC-SingleCulture_DC-SingleInfection.rds")

### coinfection
inds <- grepl('SingleCulture_DC', names(degs_sp1)) & grepl('CoInfection', names(degs_sp1), ignore.case = T)
l <- degs_sp1[inds]
names(l) <- l %>% names %>% sub("_vs_SingleCulture_DC", "", .)
make_venn(l, file = "DC-SingleCulture_DC-Coinfection.pdf", save = save, main = "DC alone vs CoInfection")
saveRDS(get_items_per_intersection(l), "DC-SingleCulture_DC-Coinfection.rds")

### all
inds <- grepl('SingleCulture_DC', names(degs_sp1))
l <- degs_sp1[inds]
names(l) <- l %>% names %>% sub("_vs_SingleCulture_DC", "", .)
#make_venn(l, file = "DC-SingleCulture_DC-all.pdf", save = save, main = "DC alone vs All")
make_upset_plot(l, "DC-SingleCulture_DC-all.pdf", save = save, nintersections = 30, number.angles = -50, font.size = 1.5)
saveRDS(get_items_per_intersection(l), "DC-SingleCulture_DC-all.rds")


### Use intersections to answer particular questions of interest

# All sets with AFu involved
items <- get_items_per_intersection(l)
s1 <- grep("SingleInfection_CMV", names(items), invert = T, ignore.case = T)
saveRDS(items[s1], "DC-All_AFu.rds")

# All sets with CMV involved
s2 <- grep("SingleInfection_AFu", names(items), invert = T, ignore.case = T)
saveRDS(items[s2], "DC-All_CMV.rds")

# Which are present regardless of pathogen and growth
i3 <- items[["SingleInfection_CMV_0h:SingleInfection_Afu_4h30min:SingleInfection_Afu_0h:SingleInfection_CMV_2h:CoInfection_CMV_first:CoInfection:CoInfection_Afu_first"]]
saveRDS(items[s3], "DC-All_intersect.rds")

# single AFU  vs CMV vs Co ---> unique per group
a <- which("SingleInfection_Afu_4h30min:SingleInfection_Afu_0h" == names(items))
stopifnot(length(a) == 1)
b <- which(names(items) == "CoInfection_CMV_first:CoInfection:CoInfection_Afu_first")
stopifnot(length(b) == 1)
c <- which(names(items) == "SingleInfection_CMV_0h:SingleInfection_CMV_2h")
stopifnot(length(c) == 1)
saveRDS(items[c(a, b, c)], "DC-AFu_vs_CMV_vs_Both.rds")

# unique genes for SI_CMV | SI_AFU | CO
s1 <- grep("SingleInfection_CMV", names(l), ignore.case = T)
write.table(l[s1] %>% unlist %>% rem_sign %>% unique, "DC-vs_SingleCulture-Si_CMV-union.csv", quote = F, row.names = F, col.names = F)
s2 <- grep("SingleInfection_AFu", names(l), ignore.case = T)
write.table(l[s2] %>% unlist %>% rem_sign %>% unique, "DC-vs_SingleCulture-Si_AFu-union.csv", quote = F, row.names = F, col.names = F)
s3 <- grep("CoInfection", names(l), ignore.case = T)
write.table(l[s3] %>% unlist %>% rem_sign %>% unique, "DC-vs_SingleCulture-Co-union.csv", quote = F, row.names = F, col.names = F)
write.table(l[s3] %>% lapply(rem_sign) %>% Reduce(intersect, .), "DC-vs_SingleCulture-Co-inter.csv", quote = F, row.names = F, col.names = F)


# AFu cells
## Time 0h
inds <- grepl('SingleCulture_Afu_0h', names(degs_sp2))
l <- degs_sp2[inds]
names(l) <- l %>% names %>% sub("_vs_SingleCulture_Afu_0h", "", .)
make_venn_plot(l, file = "AFU-SingleCulture_Afu_0h.pdf", save = save, main = "Afu 0h alone")
saveRDS(get_items_per_intersection(l), "AFU-SingleCulture_Afu_0h.rds")

## Time 4.5h
inds <- grepl('SingleCulture_Afu_4h30min', names(degs_sp2))
l <- degs_sp2[inds]
names(l) <- l %>% names %>% sub("_vs_SingleCulture_Afu_4h30min", "", .)
make_venn_plot(l, file = "AFU-SingleCulture_Afu_4h30min.pdf", save = save, main = "Afu 4.5h alone")
saveRDS(get_items_per_intersection(l), "AFU-SingleCulture_Afu_4.5h.rds")

## together
inds <- grepl('SingleCulture_Afu', names(degs_sp2))
l <- degs_sp2[inds]
names(l) <- l %>% names %>% sub("_vs_SingleCulture_Afu_0h", "", .)
names(l) <- l %>% names %>% sub("_vs_SingleCulture_Afu_4h30min", "", .)
make_upset_plot(l, "AFU-SingleCulture_AFu-all.upset.pdf", save = save, nintersections = 30, number.angles = -50, font.size = 1.5)
make_venn_plot(l, file = "AFU-SingleCulture_AFu-all.venn.pdf", save = save, main = "Afu alone")
saveRDS(get_items_per_intersection(l), "AFU-SingleCulture_Afu-all.rds")



# AFu cells - NO DESEQ
## Time 0h
inds <- grepl('SingleCulture_Afu_0h', names(degs_sp2_nod))
l <- degs_sp2_nod[inds]
names(l) <- l %>% names %>% sub("_vs_SingleCulture_Afu_0h", "", .)
make_venn_plot(l, file = "AFU-NOD-SingleCulture_Afu_0h.pdf", save = save, main = "Afu 0h alone")
saveRDS(get_items_per_intersection(l), "AFU-NOD-SingleCulture_Afu_0h.rds")

## Time 4.5h
inds <- grepl('SingleCulture_Afu_4h30min', names(degs_sp2_nod))
l <- degs_sp2_nod[inds]
names(l) <- l %>% names %>% sub("_vs_SingleCulture_Afu_4h30min", "", .)
make_venn_plot(l, file = "AFU-NOD-SingleCulture_Afu_4h30min.pdf", save = save, main = "Afu 4.5h alone")
saveRDS(get_items_per_intersection(l), "AFU-NOD-SingleCulture_Afu_4.5h.rds")

## together
inds <- grepl('SingleCulture_Afu', names(degs_sp2_nod))
l <- degs_sp2_nod[inds]
names(l) <- l %>% names %>% sub("_vs_SingleCulture_Afu_0h", "", .)
names(l) <- l %>% names %>% sub("_vs_SingleCulture_Afu_4h30min", "", .)
make_upset_plot(l, "AFU-NOD-SingleCulture_AFu-all.upset.005.pdf", save = save, nintersections = 30, number.angles = -50, font.size = 1.5)
make_venn_plot(l, file = "AFU-NOD-SingleCulture_AFu-all.venn.005.pdf", save = save, main = "Afu alone")
saveRDS(get_items_per_intersection(l), "AFU-NOD-SingleCulture_Afu-all.rds")

## just co vs si
inds <- grepl('SingleCulture_Afu', names(degs_sp2_nod))
l <- degs_sp2_nod[inds]
names(l) <- l %>% names %>% sub("_vs_SingleCulture_Afu_0h", "", .)
names(l) <- l %>% names %>% sub("_vs_SingleCulture_Afu_4h30min", "", .)
l %>%unlist %>% unique %>% length() # 3044 might be to much for STRING ~> change from 5% fdr to 1% fdr
#write.table(l %>%unlist %>% unique %>% rem_sign, "AFU-NOD-vs_SingleCulture-combined-p001.csv", quote = F, row.names = F, col.names = F)


# convert genes names to other known IDs for STRING DB
select_best_alias <- function(id, map) {
    if (length(id) > 1)
        return(sapply(id, function(x) select_best_alias(x, map)))

    row <- which(map[,1] == id)
    aliase <- map[row,] %>% unlist
    if (T %in% grepl("AFUB", aliase)) {
        c <- grep("AFUB", aliase, value = T)[1]
        return(c)
    } else {
        return(aliase[2])
    }

    # else if (T %in% grepl("CADAFU", aliase)) {
    #     c <- grep("CADAFU", aliase, value = T)[1]
    #     return(c)
    #
}

alias_mapping <- read.table("afu_alias.tsv", sep = "\t", as.is = T)
rownames(alias_mapping) <- alias_mapping[[1]]

write.table(l %>%unlist %>% unique %>% rem_sign %>% as.character %>% {select_best_alias(., map = alias_mapping)},
            "AFU-NOD-vs_SingleCulture-union-p001.alias.csv", quote = F, row.names = F, col.names = F)
write.table(get_items_per_intersection(l)[[1]] %>% rem_sign %>% as.character %>% {select_best_alias(., map = alias_mapping)},
            "AFU-NOD-vs_SingleCulture-inter-p001.alias.csv", quote = F, row.names = F, col.names = F)
write.table(l %>%unlist %>% unique %>% rem_sign %>% as.character,
            "AFU-NOD-vs_SingleCulture-union-p001.csv", quote = F, row.names = F, col.names = F)
write.table(get_items_per_intersection(l)[[1]] %>% rem_sign %>% as.character,
            "AFU-NOD-vs_SingleCulture-inter-p001.csv", quote = F, row.names = F, col.names = F)
# Differences between single AFu and CoInfection
write.table(l[grep("CoInfection", l %>% names)] %>% unlist %>% rem_sign %>% unique,
            "AFU-NOD-vs_SingleCulture-CoInfection-p001.csv", quote = F, row.names = F, col.names = F)
write.table(l[grep("SingleInfection", l %>% names)] %>% unlist %>% rem_sign %>% unique,
            "AFU-NOD-vs_SingleCulture-SiInfection_AFu-p001.csv", quote = F, row.names = F, col.names = F)


# CMV cells
save = T

inds <- grepl('single', names(degs_sp3), ignore.case = T) & grepl('single', names(degs_sp3), ignore.case = T)
l <- degs_sp3[inds]
#names(l) <- l %>% names %>% sub("_vs_SingleInfection_CMV_0h", "", .)
make_venn_plot(l, file = "CMV-vs_single_infection_0h.pdf", save = save, main = "CMV single infection vs co-infections")

# add complex test results
cmv_comp_df <- readxl::read_excel("../results/stats_CMV_complex/type.d+t+g.co_Both_vs_si.xls")
l <- c(l, "CMV d+t+g si vs co" = (cmv_comp_df %>% filter(sig == T) %>% .[1]))
names(l)[5] <- "CMV d+t+g si vs co"
make_upset_plot(l, "CMV.upset.pdf", save = save, main = "CMV DEG overlaps")

