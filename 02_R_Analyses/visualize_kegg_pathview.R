library("pathview")
library("data.table")
library("tidyverse")
library("org.Hs.eg.db")
library("SetRank")

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

out <- "../results/pathview/"
dir.create(out)
deg_res_sp1 <- readRDS("../results/DEG_overlaps/deg_res_sp1.rds")

ens2sym <- SetRank::createIDConverter("org.Hs.eg.db", from = "ENSEMBL", to = "SYMBOL")
sym2ens <- SetRank::createIDConverter("org.Hs.eg.db", from = "SYMBOL", to = "ENSEMBL")
sig_deg_res_sp1 <- get_sign_genes_geo2(deg_res_sp1$DEGs, lfc = 0)

# for combined pathview maps, make sure all of them have the same gene set (which they dont!)
shared_gene_ids <- lapply(deg_res_sp1$DEGs, function(x) x$id) %>% purrr::reduce(intersect)
deg_res_sp1_shared <- lapply(deg_res_sp1$DEGs, function(x) x %>% filter(id %in% shared_gene_ids) %>% arrange(id))

# for combined pathview maps, I select multiple DEG tables
comb_cols <- c("CoInfection_vs_SingleCulture_DC", "CoInfection_vs_SingleInfection_CMV_0h", "CoInfection_vs_SingleInfection_Afu_0h")
names(comb_cols) <- comb_cols

# FYI for checking with ensembl id refers to which gene of interest
#!  df <- data.frame(ens = dat %>% names)
#!  df$sym <- apply(df, 1, ens2sym)



#######################
## Fold Changes Combined plots for publication

dat <- lapply(comb_cols, function(col) deg_res_sp1_shared[[col]]$log2_fc_mrn)
dat <- do.call(cbind, dat)
colnames(dat) <- comb_cols
rownames(dat) <- deg_res_sp1_shared$CoInfection_vs_SingleCulture_DC$id

pathview(gene.data = dat, pathway.id = "04620", species = "hsa", gene.idtype = "ensembl", same.layer = T,
         kegg.dir = out, out.suffix = "combined", limit = list(cpd = 1.5), bins = c(cpd = 20), multi.state = T)

pathview(gene.data = dat, pathway.id = "04623", species = "hsa", gene.idtype = "ensembl",
         kegg.dir = out, out.suffix = "combined", limit = list(cpd = 1.5), bins = c(cpd = 20))

pathview(gene.data = dat, pathway.id = "04625", species = "hsa", gene.idtype = "ensembl",
         kegg.dir = out, out.suffix = "combined", limit = list(cpd = 1.5), bins = c(cpd = 20))




#######################
## REL-A | REL | REL-B

dat <- sig_deg_res_sp1$SingleInfection_Afu_0h_vs_SingleCulture_DC$log2_fc_mrn
names(dat) <- sig_deg_res_sp1$SingleInfection_Afu_0h_vs_SingleCulture_DC$id
pathview(gene.data = dat, pathway.id = "04064", species = "hsa", gene.idtype = "ensembl",
         kegg.dir = out, out.suffix = "si_afu_vs_alone.p001")


dat <- sig_deg_res_sp1$SingleInfection_CMV_0h_vs_SingleCulture_DC$log2_fc_mrn
names(dat) <- sig_deg_res_sp1$SingleInfection_CMV_0h_vs_SingleCulture_DC$id
pathview(gene.data = dat, pathway.id = "04064", species = "hsa", gene.idtype = "ensembl",
         kegg.dir = out, out.suffix = "si_cmv_vs_alone.p001")



# will all genes
dat <- deg_res_sp1$DEGs$SingleInfection_Afu_0h_vs_SingleCulture_DC$log2_fc_mrn
names(dat) <- deg_res_sp1$DEGs$SingleInfection_Afu_0h_vs_SingleCulture_DC$id
pathview(gene.data = dat, pathway.id = "04064", species = "hsa", gene.idtype = "ensembl",
         kegg.dir = out, out.suffix = "si_afu_vs_alone", limit = list(cpd = 1.5))

dat <- deg_res_sp1$DEGs$SingleInfection_CMV_0h_vs_SingleCulture_DC$log2_fc_mrn
names(dat) <- deg_res_sp1$DEGs$SingleInfection_CMV_0h_vs_SingleCulture_DC$id
pathview(gene.data = dat, pathway.id = "04064", species = "hsa", gene.idtype = "ensembl",
         kegg.dir = out, out.suffix = "si_cmv_vs_alone", limit = list(cpd = 1.5))

dat <- deg_res_sp1$DEGs$CoInfection_vs_SingleCulture_DC$log2_fc_mrn
names(dat) <- deg_res_sp1$DEGs$CoInfection_vs_SingleCulture_DC$id
pathview(gene.data = dat, pathway.id = "04064", species = "hsa", gene.idtype = "ensembl",
         kegg.dir = out, out.suffix = "co_both_vs_alone", limit = list(cpd = 1.5))

dat <- deg_res_sp1$DEGs$CoInfection_Afu_first_vs_SingleCulture_DC$log2_fc_mrn
names(dat) <- deg_res_sp1$DEGs$CoInfection_Afu_first_vs_SingleCulture_DC$id
pathview(gene.data = dat, pathway.id = "04064", species = "hsa", gene.idtype = "ensembl",
         kegg.dir = out, out.suffix = "co_afu_vs_alone", limit = list(cpd = 1.5))

dat <- deg_res_sp1$DEGs$CoInfection_CMV_first_vs_SingleCulture_DC$log2_fc_mrn
names(dat) <- deg_res_sp1$DEGs$CoInfection_CMV_first_vs_SingleCulture_DC$id
pathview(gene.data = dat, pathway.id = "04064", species = "hsa", gene.idtype = "ensembl",
         kegg.dir = out, out.suffix = "co_cmv_vs_alone", limit = list(cpd = 1.5))




######################
## TLR3


# val := dat["ENSG00000164342"]

# val: -1.4
for (lim in c(1.5)) {
    dat <- deg_res_sp1$DEGs$SingleInfection_Afu_0h_vs_SingleCulture_DC$log2_fc_mrn
    names(dat) <- deg_res_sp1$DEGs$SingleInfection_Afu_0h_vs_SingleCulture_DC$id
    pathview(gene.data = dat, pathway.id = "04620", species = "hsa", gene.idtype = "ensembl",
             kegg.dir = out, out.suffix = "si_afu_vs_alone", limit = list(cpd = lim), bins = c(cpd = 20))

    # val: 3.7
    dat <- deg_res_sp1$DEGs$SingleInfection_CMV_0h_vs_SingleCulture_DC$log2_fc_mrn
    names(dat) <- deg_res_sp1$DEGs$SingleInfection_CMV_0h_vs_SingleCulture_DC$id
    pathview(gene.data = dat, pathway.id = "04620", species = "hsa", gene.idtype = "ensembl",
             kegg.dir = out, out.suffix = "si_cmv_vs_alone", limit = list(cpd = lim), bins = c(cpd = 20))

    # val: 1.9
    dat <- deg_res_sp1$DEGs$CoInfection_vs_SingleCulture_DC$log2_fc_mrn
    names(dat) <- deg_res_sp1$DEGs$CoInfection_vs_SingleCulture_DC$id
    pathview(gene.data = dat, pathway.id = "04620", species = "hsa", gene.idtype = "ensembl",
             kegg.dir = out, out.suffix = "co_both_vs_alone", limit = list(cpd = lim), bins = c(cpd = 20))

    # val: 1.46
    dat <- deg_res_sp1$DEGs$CoInfection_Afu_first_vs_SingleCulture_DC$log2_fc_mrn
    names(dat) <- deg_res_sp1$DEGs$CoInfection_Afu_first_vs_SingleCulture_DC$id
    pathview(gene.data = dat, pathway.id = "04620", species = "hsa", gene.idtype = "ensembl",
             kegg.dir = out, out.suffix = "co_afu_vs_alone", limit = list(cpd = lim), , bins = c(cpd = 20))

    # val: 2.39
    dat <- deg_res_sp1$DEGs$CoInfection_CMV_first_vs_SingleCulture_DC$log2_fc_mrn
    names(dat) <- deg_res_sp1$DEGs$CoInfection_CMV_first_vs_SingleCulture_DC$id
    pathview(gene.data = dat, pathway.id = "04620", species = "hsa", gene.idtype = "ensembl",
             kegg.dir = out, out.suffix = "co_cmv_vs_alone", limit = list(cpd = lim), bins = c(cpd = 20))

    # val: 3.3
    dat <- deg_res_sp1$DEGs$CoInfection_vs_SingleInfection_Afu_0h$log2_fc_mrn
    names(dat) <- deg_res_sp1$DEGs$CoInfection_vs_SingleInfection_Afu_0h$id
    pathview(gene.data = dat, pathway.id = "04620", species = "hsa", gene.idtype = "ensembl",
             kegg.dir = out, out.suffix = "co_both_vs_si_afu", limit = list(cpd = lim), bins = c(cpd = 20))

    # val: -1.78
    dat <- deg_res_sp1$DEGs$CoInfection_vs_SingleInfection_CMV_0h$log2_fc_mrn
    names(dat) <- deg_res_sp1$DEGs$CoInfection_vs_SingleInfection_CMV_0h$id
    pathview(gene.data = dat, pathway.id = "04620", species = "hsa", gene.idtype = "ensembl",
             kegg.dir = out, out.suffix = "co_both_vs_si_cmv", limit = list(cpd = lim), bins = c(cpd = 20))
}



######################
## Dectin-1
## synonym: CLEC7A
# suppressed in AFU
# requires a bit of time to be fully depressed (co_cmv_first shows a positive signal)

# val := dat["ENSG00000172243"]

# val: -0.74
dat <- deg_res_sp1$DEGs$SingleInfection_Afu_0h_vs_SingleCulture_DC$log2_fc_mrn
names(dat) <- deg_res_sp1$DEGs$SingleInfection_Afu_0h_vs_SingleCulture_DC$id
pathview(gene.data = dat, pathway.id = "04625", species = "hsa", gene.idtype = "ensembl",
         kegg.dir = out, out.suffix = "si_afu_vs_alone", limit = list(cpd = 1.5), bins = c(cpd = 20))

# val: 0.016
dat <- deg_res_sp1$DEGs$SingleInfection_CMV_0h_vs_SingleCulture_DC$log2_fc_mrn
names(dat) <- deg_res_sp1$DEGs$SingleInfection_CMV_0h_vs_SingleCulture_DC$id
pathview(gene.data = dat, pathway.id = "04625", species = "hsa", gene.idtype = "ensembl",
         kegg.dir = out, out.suffix = "si_cmv_vs_alone", limit = list(cpd = 1.5), bins = c(cpd = 20))

# val: -0.62
dat <- deg_res_sp1$DEGs$CoInfection_vs_SingleCulture_DC$log2_fc_mrn
names(dat) <- deg_res_sp1$DEGs$CoInfection_vs_SingleCulture_DC$id
pathview(gene.data = dat, pathway.id = "04625", species = "hsa", gene.idtype = "ensembl",
         kegg.dir = out, out.suffix = "co_both_vs_alone", limit = list(cpd = 1.5), bins = c(cpd = 20))

# val: -0.59
dat <- deg_res_sp1$DEGs$CoInfection_Afu_first_vs_SingleCulture_DC$log2_fc_mrn
names(dat) <- deg_res_sp1$DEGs$CoInfection_Afu_first_vs_SingleCulture_DC$id
pathview(gene.data = dat, pathway.id = "04625", species = "hsa", gene.idtype = "ensembl",
         kegg.dir = out, out.suffix = "co_afu_vs_alone", limit = list(cpd = 1.5), , bins = c(cpd = 20))

# val: 0.28
dat <- deg_res_sp1$DEGs$CoInfection_CMV_first_vs_SingleCulture_DC$log2_fc_mrn
names(dat) <- deg_res_sp1$DEGs$CoInfection_CMV_first_vs_SingleCulture_DC$id
pathview(gene.data = dat, pathway.id = "04625", species = "hsa", gene.idtype = "ensembl",
         kegg.dir = out, out.suffix = "co_cmv_vs_alone", limit = list(cpd = 1.5), bins = c(cpd = 20))

# val: 0.14
dat <- deg_res_sp1$DEGs$CoInfection_vs_SingleInfection_Afu_0h$log2_fc_mrn
names(dat) <- deg_res_sp1$DEGs$CoInfection_vs_SingleInfection_Afu_0h$id
pathview(gene.data = dat, pathway.id = "04625", species = "hsa", gene.idtype = "ensembl",
         kegg.dir = out, out.suffix = "co_both_vs_si_afu", limit = list(cpd = 1.5), bins = c(cpd = 20))

# val: -0.62
dat <- deg_res_sp1$DEGs$CoInfection_vs_SingleInfection_CMV_0h$log2_fc_mrn
names(dat) <- deg_res_sp1$DEGs$CoInfection_vs_SingleInfection_CMV_0h$id
pathview(gene.data = dat, pathway.id = "04625", species = "hsa", gene.idtype = "ensembl",
         kegg.dir = out, out.suffix = "co_both_vs_si_cmv", limit = list(cpd = 1.5), bins = c(cpd = 20))




######################
## Just a lot of other pathways
###


pathway_names <- c("05163", "04668", "04630")
for (pn in pathway_names) {
    for (comp in names(deg_res_sp1$DEGs)) {
        dat <- deg_res_sp1$DEGs[[comp]]$log2_fc_mrn
        names(dat) <- deg_res_sp1$DEGs[[comp]]$id
        pathview(gene.data = dat, pathway.id = pn, species = "hsa", gene.idtype = "ensembl",
                 kegg.dir = out, out.suffix = comp, limit = list(cpd = 1.5), bins = c(cpd = 20))
    }
}


pathway_names <- c("04623", "04621")
for (pn in pathway_names) {
  for (comp in names(deg_res_sp1$DEGs)) {
    dat <- deg_res_sp1$DEGs[[comp]]$log2_fc_mrn
    names(dat) <- deg_res_sp1$DEGs[[comp]]$id
    pathview(gene.data = dat, pathway.id = pn, species = "hsa", gene.idtype = "ensembl",
             kegg.dir = out, out.suffix = comp, limit = list(cpd = 1.5), bins = c(cpd = 20),
             kegg.native = F)
  }
}

