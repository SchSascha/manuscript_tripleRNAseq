library("biomaRt")
library("tidyverse")
library("AnnotationDbi")
library("org.Hs.eg.db")


sortDiffExp <- function(degTable, tools) {
    # make sure to compensate for different caps in names
    tools_idx <- which(tolower(colnames(degTable)) %in% tolower(tools))
    # by bools first, by adj. p-val second
    # bool index must be inverse for it to work!
    degTable <- degTable[order((!degTable[, tools_idx]) %>% as.matrix %>% rowSums,
                               degTable[, tools_idx+1] %>% as.matrix %>% rowSums,
                               decreasing = FALSE), ]
    return(degTable)
}


annotate <- function(df) {
    if (nrow(df) == 0)
        return(df)

    ids <- df$id
    res <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'wikigene_description', 'gene_biotype'),
                 filters = 'ensembl_gene_id',
                 values = ids,
                 mart = ensembl)
    res <- filter(res, !duplicated(ensembl_gene_id)) # it can still happen...
    colnames(res)[1] <- "id"

    return(merge(df, res, all.x = T))
}


annotate_deg_file <- function(file, do_sort = T) {
    if (grepl("\\.csv$", file)) {
        pages <- "unnamed"
        df_list <- list(read.csv(file, sep = "\t"))
    } else if (grepl("\\.xls$", file)) {
        pages <- readxl::excel_sheets(file)
        df_list <- lapply(pages, function(l) readxl::read_xls(file, sheet = l))
    } else if (grepl("\\.xlsx$", file)) {
        pages <- readxl::excel_sheets(file)
        df_list <- lapply(pages, function(l) readxl::read_xlsx(file, sheet = l))
    }

    anno_df_list <- lapply(df_list, function(df){
        merged_df <- annotate(df)

        if (do_sort) {
            tools <- grep("_adj_pval", colnames(merged_df), value = T)
            tools <- gsub("_adj_pval", "", tools)
            merged_df <- sortDiffExp(merged_df, tools)
        }
    })
    names(anno_df_list) <- pages

    return(anno_df_list)
}


listMarts()
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
attributes[1:5,]
