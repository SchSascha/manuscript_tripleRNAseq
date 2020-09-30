# R341
# Geo2RNAseq version  0.99.12
library("Geo2RNAseq")

# get info #####################################################################

# Fix SDRF -> add files
SDRF_file <- "../doc/FungiNet_template_RNA-seq_V2.1_DC.csv"
SDRF_info <- parse_SDRF(SDRF_file, skip = 8, col.name.row=7)

# determine order of samples in SDRF file
samp_names <- substr(SDRF_info$samples,
                     start  = nchar("Donor1_01_1"),
                     stop = nchar("Donor4_06_DC_plus_CMV_0h_plus_Afu_4h30min"))
fq <- dir("../raw/") # sample names
fq_ren <- gsub("\\.fq\\.gz", "", fq)
fq_ren <- sub("ID-.*-", "", fq_ren)
fq_ren <- gsub("\\+", "plus", fq_ren)
fq_ren <- sub("4\\.5h", "4h30min", fq_ren)

# Reorder samples for match the meta-data order
library("stringr")
to_reverse <- str_count(fq_ren,"_plus_") == 2 & !grepl("DC_plus_Afu_0h_plus_CMV_2h", fq_ren, fixed = T)
tmp <- strsplit(fq_ren[to_reverse], "_plus_", fixed = T)
fq_fixed <- sapply(tmp, function(t) paste(t[1], t[3], t[2], sep = "_plus_"))
fq_ren[to_reverse] <- fq_fixed
# the resulting order is identical
stopifnot(fq_ren == samp_names)

# add file info
SDRF_info$table$Data.File <- paste0("../raw/", fq)
# add single-end specification
SDRF_info$table$Library.Strategy <- paste0(SDRF_info$table$Library.Strategy, "; single-end")
# write meta-data out
write_meta_table(SDRF_info$table, as.xls=T, outDir="../doc/")


# load updated file
new_SDRF_file <- "../doc/meta_data.xls"
SDRF_info <- parse_meta_xls(new_SDRF_file)

# continue as usual
SDRF           <- SDRF_info$table
number_samples <- SDRF_info$nrSamples
samples        <- SDRF_info$samples
design_matrix  <- SDRF_info$designMatrix

# global settings ##############################################################
MAX_CPUS <- 75 # allocate no more than <int> CPUs <-- fill in!
# If TRUE, existing files will be overwritten without questioning.
# If FALSE, most methods will skip step for existing output files.
# In these cases, the method will return 'call[x] = "NOT USED"'
FORCE_OVERWRITE <- FALSE # <-- fill in!

# gene ids start with ...
geneID_sp1 <- "ENSG"
geneID_sp2 <- "Afu"
geneID_sp3 <- "cds-ABV"

name_sp1 <- "Hsapiens"
name_sp2 <- "Afumigatus"
name_sp3 <- "CMV"


# paired or unpaired reads?
paired <- F

# files ########################################################################

# sp1 = species 1
# sp2 = species 2
# comb = combined = sp1+sp2

# caveat: TopHat2/HISAT2 are looking for the genome file in the index directory
# move it there or create a link if the genome is located in its own directory
d <- "/sbidata/shared2/Combined_Genomes_forTripleRNAseqPreprocessing/Hsapiens_GRCh38_89_2017-05-07_Human_herpesvirus_5_strain_TB40E_clone_TB40-BAC4/fasta/"
genome_sp1  <- file.path(d, "Homo_sapiens.GRCh38.dna.primary_assembly.fa")
genome_sp2  <- file.path(d, "A_fumigatus_Af293_current_chromosomes.fasta")
genome_sp3  <- file.path(d, "EF999921.1.fasta")
genome_comb <- file.path(d, "../index_hisat2/hs-GR38.89_afu-Af293_cmv-5.fa")

# <gene_id "something";> must be first attribute
anno_sp1  <- file.path(d, "../anno/Homo_sapiens.GRCh38.89.mod.gtf")
anno_sp2  <- file.path(d, "../anno/A_fumigatus_Af293_current_features.mod.gtf")
anno_sp3  <- file.path(d, "../anno/EF999921.1.mod.gtf")
anno_comb <- "/sbidata/shared2/Combined_Genomes_forTripleRNAseqPreprocessing/Hsapiens_GRCh38_89_2017-05-07_Human_herpesvirus_5_strain_TB40E_clone_TB40-BAC4/anno/hs-GR38.89_afu-Af293_cmv-5.gtf"

# build index with make_HiSat2_index() or make_Tophat_index()
index_comb <- "/sbidata/shared2/Combined_Genomes_forTripleRNAseqPreprocessing/Hsapiens_GRCh38_89_2017-05-07_Human_herpesvirus_5_strain_TB40E_clone_TB40-BAC4/index_hisat2/hs-GR38.89_afu-Af293_cmv-5.fa"



SDRF$Genome.Annotation.Source    <- paste("Ensembl", "Aspergillus Genome Database", "GenBank", sep = "; ")
SDRF$Genome.Annotation.Version   <- paste("GRCh38_89", "s03-m05-r09", "EF999921.1", sep = "; ")
SDRF$Genome.Annotation.Date      <- paste("2017-05-07", "2018-05-06", "2016-07-26", sep = "; ")
SDRF$Genome.Annotation.Data.File <- genome_comb

SDRF$Structural.Annotation.Source    <- paste("Ensembl", "Aspergillus Genome Database", "GenBank", sep = "; ")
SDRF$Structural.Annotation.Version   <- paste("GRCh38_89", "s03-m05-r09", "EF999921.1", sep = "; ")
SDRF$Structural.Annotation.Date      <- paste("2017-05-07", "2018-05-06", "2016-07-26", sep = "; ")
SDRF$Structural.Annotation.Data.File <- anno_comb


# raw FASTQ files ##############################################################

raw_fastq_files <- SDRF_info$seFiles

fastqDir <- "../results/fastq"
stopifnot(file.exists(raw_fastq_files))

# quality control 1 (before trimming) ##########################################

qualityDir <- "../results/quality"
qualityRawDir <- file.path(qualityDir, "raw")
if (length(dir(qualityRawDir)) > 0)
    warning(paste0("Directory \"", qualityRawDir, "\" not empty! Overwrite? ", FORCE_OVERWRITE), immediate. = TRUE)

if (FORCE_OVERWRITE || length(dir(qualityRawDir)) == 0)
{
    writeLines("FastQC - raw ...")
    fq_res <- run_FastQC(raw_fastq_files, outDir = qualityRawDir, cpus = MAX_CPUS, extend = T)
    SDRF$Data.Processing.Software <- "FastQC" # according to FungiNet Lookup-Table
    SDRF$Description <- fq_res$version
    SDRF$Parameters.And.Values <- fq_res$calls
}

# trimming #####################################################################

trimmedDir <- fastqDir
windowsizetrimming <- 15
qualcuttrimming <- 25
phred <- "-phred33"
leading <- 3
trailing <- 3
minlen <- 30
adapters <- NA
#trimoJar <- "/opt/ngs/Trimmomatic-0.36/trimmomatic-0.36.jar"

if (length(list.files(trimmedDir, pattern = "\\.trimo(\\.pe)*\\.fastq$")) > 0)
    warning(paste0("Directory \"", trimmedDir, "\" not empty! Overwrite? ", FORCE_OVERWRITE), immediate. = TRUE)

trimming_res <- run_Trimmomatic(
    files      = raw_fastq_files,
    is.paired  = paired,
    outDir     = trimmedDir,
    cpus       = MAX_CPUS,
    windowsize = windowsizetrimming,
    qualcut    = qualcuttrimming,
    phred      = phred,
    leading    = leading,
    trailing   = trailing,
    minlen     = minlen,
    adapters   = adapters,
    overwrite  = FORCE_OVERWRITE
)

input <- trimming_res$input
surviving <- trimming_res$surviving

# update SDRF/metadata
software_nr <- 1
rows <- which(trimming_res$calls != "NOT USED")
col <- paste0("Data.Processing.Software.", software_nr)
SDRF[rows, col] <- "trimmomatic" # according to FungiNet Lookup-Table
col <- paste0("Description.", software_nr)
SDRF[rows, col] <- basename(sub(".*/(.+)\\.jar.*", "\\1", trimming_res$calls[1])) # name and version of tool
col <- paste0("Parameters.And.Values.", software_nr)
SDRF[rows, col] <- trimming_res$calls # complete call

trimmed_fastq_files <- trimming_res$files

# sanity check. Can be ignored during experimental sessions
stopifnot(length(trimmed_fastq_files) == number_samples * (paired+1))

fastq_files <- trimmed_fastq_files

## number of reads after trimming based on Trimmomatic output --> grep
number_raw <- input
number_trimmed <- surviving

names(number_raw) <- basename(raw_fastq_files)
names(number_trimmed) <- basename(trimmed_fastq_files)

# percentage of trimmed reads
perc_trimmed <- sort(1 - number_trimmed / number_raw) * 100
fastq_files <- trimmed_fastq_files

# quality control 2 (after trimming) ###########################################

qualityDir <- "../results/quality"
qualitySubDir <- file.path(qualityDir, "trimmed")
if (length(dir(qualitySubDir)) > 0)
    warning(paste0("Directory \"" , qualitySubDir, "\" not empty! overwrite? ", FORCE_OVERWRITE))

if (FORCE_OVERWRITE || length(dir(qualitySubDir)) == 0)
{
    writeLines("FastQC - trimmed ...")
    fq_res_trim <- run_FastQC(fastq_files, outDir = qualitySubDir, cpus = MAX_CPUS, extend = F)
}

# mapping ######################################################################

index_map <- index_comb
anno_map  <- anno_comb

# STRAND SPECIFIC STRAND SPECIFIC STRAND SPECIFIC STRAND SPECIFIC STRAND SPECIFIC
# STRAND SPECIFIC STRAND SPECIFIC STRAND SPECIFIC STRAND SPECIFIC STRAND SPECIFIC
# STRAND SPECIFIC STRAND SPECIFIC STRAND SPECIFIC STRAND SPECIFIC STRAND SPECIFIC
# STRAND SPECIFIC STRAND SPECIFIC STRAND SPECIFIC STRAND SPECIFIC STRAND SPECIFIC

mapper <- "hisat2"
mapDir <- "../results/mapping"
bamDir <- file.path(mapDir, "bamfiles")
convert_sam <- TRUE # should SAM files be converted to BAM files if BAM files are not found?

samDir <- file.path(mapDir, "samfiles")
bamDir <- file.path(mapDir, "bamfiles")

if (FORCE_OVERWRITE || length(dir(samDir)) == 0)
{
    mapping_res <- run_Hisat2(
        files     = fastq_files,
        index     = index_map,
        outDir    = mapDir,
        is.paired = paired,
        addArgs   = "--rna-strandness F", # only forward reads
        splice    = T,
        cpus      = MAX_CPUS,
        overwrite = FORCE_OVERWRITE
    )

    ## NOTE: use that only if the run_Hisat2() parameter 'as.bam' is FALSE.
    # bam_files <- sam_to_bam(mapping_res$mapFiles, sort=TRUE, bamDir = bamDir)
    bam_files <- mapping_res$files

    software_nr <- software_nr + 1
    rows <- which(mapping_res$calls != "NOT USED")
    col <- paste0("Data.Processing.Software.", software_nr)
    SDRF[rows, col] <- "HiSat2" # according to FungiNet Lookup-Table
    col <- paste0("Description.", software_nr)
    SDRF[rows, col] <- mapping_res$version
    col <- paste0("Parameters.And.Values.", software_nr)
    SDRF[rows, col] <- mapping_res$calls

} else {
    if (convert_sam) {
        sam_files <- list.files(samDir, pattern = "\\.sam$", full.names = TRUE)
        if (length(sam_files) != number_samples)
            stop(paste("Invalid number of SAM files. Expected", number_samples, "but got", length(sam_files)))
        bam_files <- sam_to_bam(mapping_res$mapFiles, sort=TRUE, bamDir = bamDir)
    } else
        bam_files <- list.files(bamDir, pattern = "\\.bam$", full.names = TRUE)

    if (length(bam_files) != number_samples)
        warning(paste("Found", length(bam_files), "BAM files, but expected", number_samples, ". Maybe SAMtools was interrupted. Try to convert again? ", convert_sam), immediate. = TRUE)
}

save.image("mapping.rda")


# counting #####################################################################
anno_count <- anno_comb


# featureCounts out is on stderr --> have to use 2>&1
# in addition, "featureCounts -v" is surrounded by empty lines

countDir <- "../results/counting"
if (length(dir(countDir)) > 0)
    warning(paste("Directory", countDir, "not empty! Overwrite?", FORCE_OVERWRITE), immediate. = TRUE)
if (FORCE_OVERWRITE || length(dir(countDir)) == 0)
{
    res_counting <- run_featureCounts(
        files             = bam_files,
        annotation        = anno_count,
        isGTF             = TRUE,
        IDtype            = "gene_id",
        featureType       = "exon",
        outDir            = countDir,
        isPairedEnd       = paired,
        allowMultiOverlap = FALSE,
        cpus              = MAX_CPUS
    )

    software_nr <- software_nr + 1
    rows <- which(res_counting$calls != "NOT USED")
    col <- paste0("Data.Processing.Software.", software_nr)
    SDRF[rows, col] <- "featureCounts" #  according to FungiNet Lookup-Table
    col <- paste0("Description.", software_nr)
    SDRF[rows, col] <- res_counting$version
    col <- paste0("Parameters.And.Values.", software_nr)
    SDRF[rows, col] <- res_counting$calls

} else {
    res_counting <- read.featureCounts.files(countDir)
}


# stranded counting
countDir2 <- "../results/counting_stranded"
res_counting2 <- run_featureCounts(
    files             = bam_files,
    annotation        = anno_count,
    isGTF             = TRUE,
    IDtype            = "gene_id",
    featureType       = "exon",
    outDir            = countDir2,
    isPairedEnd       = paired,
    strandSpecific    = 1,
    allowMultiOverlap = FALSE,
    cpus              = MAX_CPUS
)


save.image("counting.rda")

# for dual, the counting variables need additional processing

# generate count file for each FASTQ/BAM/SAM file
# (needed for MultiQCs featureCounts report)
for (i in 1:nrow(res_counting$summary)) {
    # add "Status" line to beginning of table
    table <- c(rownames(res_counting$summary)[i], res_counting$summary[i, ])
    names(table) <- c("Status", names(res_counting$summary[i, ]))

    # write table to file
    write.table(
        table,
        file = file.path(countDir, paste0(table[1], ".counts.summary")),
        sep = "\t",
        quote = FALSE,
        col.names = FALSE
    )
}

counts <- res_counting$counts
countfile <- res_counting$countFile
sumfile <- res_counting$sumFile

counts_sp1 <- counts[grep(geneID_sp1, rownames(counts)), ]
counts_sp2 <- counts[grep(geneID_sp2, rownames(counts)), ]
counts_sp3 <- counts[grep(geneID_sp3, rownames(counts)), ]

colnames(counts_sp1) <- samples
colnames(counts_sp2) <- samples
colnames(counts_sp3) <- samples

# get gene lengths from featureCounts result
gene_length_sp1 <- res_counting$anno$Length[grep(geneID_sp1, res_counting$anno$GeneID)]
gene_length_sp2 <- res_counting$anno$Length[grep(geneID_sp2, res_counting$anno$GeneID)]
gene_length_sp3 <- res_counting$anno$Length[grep(geneID_sp3, res_counting$anno$GeneID)]

names(gene_length_sp1) <- res_counting$anno$GeneID[grep(geneID_sp1, res_counting$anno$GeneID)]
names(gene_length_sp2) <- res_counting$anno$GeneID[grep(geneID_sp2, res_counting$anno$GeneID)]
names(gene_length_sp3) <- res_counting$anno$GeneID[grep(geneID_sp3, res_counting$anno$GeneID)]

statsDir <- "../results/stats"

# write counts to CSV file
write_count_table(
    file = file.path(statsDir, paste0(name_sp1, "_counts")),
    counts = counts_sp1
)
write_count_table(
    file = file.path(statsDir, paste0(name_sp2, "_counts")),
    counts = counts_sp2
)
write_count_table(
    file = file.path(statsDir, paste0(name_sp3, "_counts")),
    counts = counts_sp3
)

# write counts to XLS file
write_count_table(
    file       = file.path(statsDir, paste0(name_sp1, "_counts")),
    counts     = counts_sp1,
    as.xls     = TRUE,
    sheetNames = basename(getwd())
)
write_count_table(
    file       = file.path(statsDir, paste0(name_sp2, "_counts")),
    counts     = counts_sp2,
    as.xls     = TRUE,
    sheetNames = basename(getwd())
)
write_count_table(
    file       = file.path(statsDir, paste0(name_sp3, "_counts")),
    counts     = counts_sp3,
    as.xls     = TRUE,
    sheetNames = basename(getwd())
)



# check for rRNA contamination #################################################

detect_high_coverage(counts_sp1)
detect_high_coverage(counts_sp2)
detect_high_coverage(counts_sp3)


# MultiQC ######################################################################
# TODO run some modules for sp1 and sp1 separately?

multiqc_call <- run_MultiQC(dir="../results/", outDir="../results", config = get_MultiQC_config())

# SAMtools #####################################################################

flagstatDir <- "mapping/stats"

bam_files <- sortBAMs(bam_files, overwrite = FORCE_OVERWRITE, delete.orig = F)
calls <- indexBAMs(bam_files, cpus = MAX_CPUS)


# mapping stats ################################################################

mapping_stats_df <- calc_triple_mapping_stats(
  number_raw,
  fastq_files,
  name_sp1,
  name_sp2,
  name_sp3,
  genome_sp1,
  genome_sp2,
  genome_sp3,
  anno_sp1,
  anno_sp2,
  anno_sp3,
  counts_sp1,
  counts_sp2,
  counts_sp3
)



# sort mapping stats by rownames in a more natural way
library("gtools")

# write mapping stats to CSV file
write_count_table(
    file   = file.path(statsDir, "mapping_stats"),
    counts = mapping_stats_df[mixedsort(rownames(mapping_stats_df)),],
    rnames = TRUE
)

# write mapping stats to XLS file
write_count_table(
    file       = file.path(statsDir, "mapping_stats"),
    counts     = mapping_stats_df[mixedsort(rownames(mapping_stats_df)),],
    as.xls     = TRUE,
    sheetNames = basename(getwd()),
    rnames     = TRUE
)
write_count_table(
    file       = file.path(statsDir, "mapping_stats"),
    counts     = mapping_stats_df,
    as.xls     = F,
    sheetNames = basename(getwd()),
    rnames     = TRUE
)

# RPKM and TPM values ##########################################################

rpkm_sp1 <- get_rpkm(counts_sp1, gene_length_sp1, colSums(counts_sp1))
rpkm_sp2 <- get_rpkm(counts_sp2, gene_length_sp2, colSums(counts_sp2))
rpkm_sp3 <- get_rpkm(counts_sp3, gene_length_sp3, colSums(counts_sp3))
tpm_sp1  <- get_tpm( counts_sp1, gene_length_sp1, colSums(counts_sp1))
tpm_sp2  <- get_tpm( counts_sp2, gene_length_sp2, colSums(counts_sp2))
tpm_sp3  <- get_tpm( counts_sp3, gene_length_sp3, colSums(counts_sp3))
mrn_sp1  <- get_mrn(counts_sp1)
mrn_sp2  <- get_mrn(counts_sp2)
mrn_sp3  <- get_mrn(counts_sp3+1)

# write RPKM and TPM values to CSV file
write_count_table(
    file   = file.path(statsDir, name_sp1, paste0(name_sp1, "_rpkm")),
    counts = rpkm_sp1
)
write_count_table(
    file   = file.path(statsDir, name_sp2, paste0(name_sp2, "_rpkm")),
    counts = rpkm_sp2
)
write_count_table(
    file   = file.path(statsDir, name_sp3, paste0(name_sp3, "_rpkm")),
    counts = rpkm_sp3
)
write_count_table(
    file   = file.path(statsDir, name_sp1, paste0(name_sp1, "_tpm")),
    counts = tpm_sp1
)
write_count_table(
    file   = file.path(statsDir, name_sp2, paste0(name_sp1, "_tpm")),
    counts = tpm_sp2
)
write_count_table(
    file   = file.path(statsDir, name_sp3, paste0(name_sp3, "_tpm")),
    counts = tpm_sp3
)
write_count_table(
    file   = file.path(statsDir, paste0(name_sp1, "_mrn")),
    counts = mrn_sp1
)
write_count_table(
    file   = file.path(statsDir, paste0(name_sp2, "_mrn")),
    counts = mrn_sp2
)
write_count_table(
    file   = file.path(statsDir, paste0(name_sp3, "_mrn")),
    counts = mrn_sp3
)



# clustering ###################################################################


# NOTE: experimental function. Check resulting vector!
conds_sp1 <- conditions_from_design(design_matrix_sp1)
conds_sp2 <- conditions_from_design(design_matrix_sp2)
conds_sp3 <- conditions_from_design(design_matrix_sp3)
bot <- 12

# hclust heat map version
make_heat_clustering_plot(
    file.path(statsDir, name_sp1, paste0(name_sp1, "_heat_hierarchical_clustering_counts.mrn")),
    main         = "Hierarchical Clustering Counts",
    counts       = counts_sp1,
    designMatrix = design_matrix_sp1,
    conds        = conds_sp1,
    norm         = "mrn",
    overwrite    = FORCE_OVERWRITE
)
# only DC samples
dc_ind <- grep("DC", colnames(counts_sp1))
make_heat_clustering_plot(
    file.path(statsDir, name_sp1, paste0(name_sp1, "_heat_hierarchical_clustering_counts.DC.mrn")),
    main         = "Hierarchical Clustering Counts",
    counts       = counts_sp1[,dc_ind],
    conds        = conds_sp1[dc_ind],
    norm         = "mrn",
    overwrite    = FORCE_OVERWRITE
)

# only AFU samples
make_heat_clustering_plot(
    file.path(statsDir, name_sp2, paste0(name_sp2, "_heat_hierarchical_clustering_counts.mrn")),
    main         = "Hierarchical Clustering Counts",
    counts       = counts_sp2,
    conds        = conds_sp2,
    norm         = "mrn",
    overwrite    = FORCE_OVERWRITE
)

af_ind <- grep("Afu", colnames(counts_sp2))
make_heat_clustering_plot(
    file.path(statsDir, name_sp2, paste0(name_sp2, "_heat_hierarchical_clustering_counts.AF.mrn")),
    main         = "Hierarchical Clustering Counts",
    counts       = counts_sp2[,af_ind],
    conds        = conds_sp2[af_ind],
    norm         = "mrn",
    overwrite    = FORCE_OVERWRITE
)

cmv_ind <- grep("CMV", colnames(counts_sp3))
make_heat_clustering_plot(
    file.path(statsDir, name_sp3, paste0(name_sp3, "_heat_hierarchical_clustering_counts.mrn")),
    main         = "Hierarchical Clustering Counts",
    counts       = counts_sp3[,cmv_ind],
    conds        = conds_sp3[cmv_ind],
    norm         = "mrn",
    overwrite    = FORCE_OVERWRITE
)



# PEARSON CORRELATION ##########################################################

# NOTE: Regression is performed. This takes a considerable amount of time
#       for large gene sets and/or samples. If in doubt, use only the first.
#make_correlation_plots(counts, "correlation", "RAW_",  design_matrix, overwrite=T)
#make_correlation_plots(tpm,    "correlation", "TPM_",  design_matrix, overwrite=T)
#make_correlation_plots(rpkm,   "correlation", "RPKM_", design_matrix, overwrite=T)
sel <- which(meta$condition == "SingleInfection_CMV_0h")
make_correlation_plots(mrn_sp3[,sel], "../results/correlation", "MRN_si_CMV_0h", overwrite=T)
sel <- which(meta$condition == "SingleInfection_CMV_2h")
make_correlation_plots(mrn_sp3[,sel], "../results/correlation", "MRN_si_CMV_2h", overwrite=T)

# PCA plots ####################################################################

# NOTE: if designMatrix is supplied, 'conds' is ignored. In that case, use 'designMatrix = NA'
conds_sp1 <- conditions_from_design(design_matrix)
SDRF_file_sp2 <- "../doc/FungiNet_template_RNA-seq_V2.1_AFu.csv"
SDRF_info_sp2 <- parse_SDRF(SDRF_file_sp2, skip = 8, col.name.row=7)
conds_sp2 <- conditions_from_design(SDRF_info_sp2$designMatrix)
SDRF_file_sp3 <- "../doc/FungiNet_template_RNA-seq_V2.1_CMV.csv"
SDRF_info_sp3 <- parse_SDRF(SDRF_file_sp3, skip = 8, col.name.row=7)
conds_sp3 <- conditions_from_design(SDRF_info_sp3$designMatrix)


colName <- "conditions"   # <-- fill in!
## in this vector, time series information can be supplied.
shapes <- gsub("_.*", "", SDRF_info$samples)
shapes_sp1 <- shapes        # <-- optional. fill in!
shapes_sp2 <- shapes        # <-- optional. fill in!
shapes_sp3 <- shapes        # <-- optional. fill in!
shapeName <- "Donor" # <-- optional. fill in!
# Fix conds
#conds_sp1 <- gsub("_.+$", "", conds_sp1)
#conds_sp2 <- gsub("_.+$", "", conds_sp2)

pca_res_sp1 <- make_PCA_plot(
    file         = file.path(statsDir, name_sp1, paste0(name_sp1, "_PCA_vst")),
    counts       = counts_sp1[, conds_sp1 != "none"],
    designMatrix = NA,
    conds        = conds_sp1[conds_sp1 != "none"],
    shapes       = shapes_sp1[conds_sp1 != "none"],
    norm         = "vst",
    main         = paste0("PCA - vst - ", name_sp1),
    colName      = colName,
    shapeName    = shapeName,
    overwrite    = FORCE_OVERWRITE
)



pca_res_sp2 <- make_PCA_plot(
    file         = file.path(statsDir, name_sp2, paste0(name_sp2, "_PCA_vst")),
    counts       = counts_sp2[, conds_sp2 != "none"],
    designMatrix = NA,
    conds        = conds_sp2[conds_sp2 != "none"],
    shapes       = shapes_sp2[conds_sp2 != "none"],
    norm         = "vst",
    main         = paste0("PCA - vst - ", name_sp2),
    colName      = colName,
    shapeName    = shapeName,
    overwrite    = FORCE_OVERWRITE
)
pca_res_sp3 <- make_PCA_plot(
    file         = file.path(statsDir, name_sp3, paste0(name_sp3, "_PCA_rlog")),
    counts       = counts_sp3[, conds_sp3 != "none"],
    designMatrix = NA,
    conds        = conds_sp3[conds_sp3 != "none"],
    shapes       = shapes_sp3[conds_sp3 != "none"],
    norm         = "rlog",
    main         = paste0("PCA - rlog - ", name_sp3),
    colName      = colName,
    shapeName    = shapeName,
    overwrite    = FORCE_OVERWRITE
)

#source("R/Helper.R") # for MRN method

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
pca_res_sp1$g + theme_bw() + coord_fixed()

pca_res_sp2 <- make_PCA_plot(
    file         = file.path("../results/PCA", name_sp2, paste0(name_sp2, "_PCA_mrn")),
    counts       = log2(mrn_sp2[, conds_sp2 != "none"]+1),
    designMatrix = NA,
    conds        = conds_sp2[conds_sp2 != "none"],
    shapes       = shapes_sp2[conds_sp2 != "none"],
    norm         = "none",
    main         = paste0("PCA - ", name_sp2),
    colName      = colName,
    shapeName    = shapeName,
    overwrite    = FORCE_OVERWRITE
)

pca_res_sp3 <- make_PCA_plot(
    file         = file.path("../results/PCA", name_sp3, paste0(name_sp3, "_PCA_mrn")),
    counts       = log2(mrn_sp3[, conds_sp3 != "none"]+1),
    designMatrix = NA,
    conds        = conds_sp3[conds_sp3 != "none"],
    shapes       = shapes_sp3[conds_sp3 != "none"],
    norm         = "none",
    main         = paste0("PCA - ", name_sp3),
    colName      = colName,
    shapeName    = shapeName,
    overwrite    = FORCE_OVERWRITE
)



# custom PCA plots
# DC - DC vs single vs coinfection
nconds <- rep("none", length(conds_sp1))
nconds[grep("DC", conds_sp1)] <- "DC"
nconds[grep("SingleInfection_Afu", conds_sp1)] <- "SingleInfection_Afu"
nconds[grep("SingleInfection_CMV", conds_sp1)] <- "SingleInfection_CMV"
nconds[grep("CoInfection", conds_sp1)] <- "Coinfection"
pca_res_sp4 <- make_PCA_plot(
    file         = file.path(statsDir, name_sp1, paste0(name_sp1, "_PCA_GroupedInfection_vst")),
    counts       = counts_sp1[, nconds != "none"],
    designMatrix = NA,
    conds        = nconds[nconds != "none"],
    shapes       = shapes_sp1[nconds != "none"],
    norm         = "vst",
    main         = paste0("PCA - vst - ", name_sp1),
    colName      = colName,
    shapeName    = shapeName,
    overwrite    = FORCE_OVERWRITE
)

pca_res_sp4 <- make_PCA_plot(
    file         = file.path("../results/other_plots/", paste0(name_sp1, "_PCA_GroupedInfection_mrn")),
    counts       = log2(mrn_sp1[, conds_sp1 != "none"]+1),
    designMatrix = NA,
    conds        = nconds[nconds != "none"],
    shapes       = shapes_sp1[conds_sp1 != "none"],
    norm         = "none",
    main         = paste0("PCA - log2 MRN + 1 - ", name_sp1),
    colName      = colName,
    shapeName    = shapeName,
    overwrite    = FORCE_OVERWRITE
)


# AF - group infection
nconds <- conds_sp2
nconds[grep("Infection", conds_sp2)] <- "Infection"

pca_res_sp5 <- make_PCA_plot(
    file         = file.path(statsDir, name_sp2, paste0(name_sp2, "_PCA_GroupedInfection_vst")),
    counts       = counts_sp2[, nconds != "none"],
    designMatrix = NA,
    conds        = nconds[nconds != "none"],
    shapes       = shapes_sp2[nconds != "none"],
    norm         = "vst",
    main         = paste0("PCA - vst - ", name_sp2),
    colName      = colName,
    shapeName    = shapeName,
    overwrite    = FORCE_OVERWRITE
)

pca_res_sp5 <- make_PCA_plot(
    file         = file.path("../results/other_plots/", paste0(name_sp2, "_PCA_GroupedInfection_mrn")),
    counts       = log2(mrn_sp2[, conds_sp2 != "none"]+1),
    designMatrix = NA,
    conds        = nconds[nconds != "none"],
    shapes       = shapes_sp2[nconds != "none"],
    norm         = "none",
    main         = paste0("PCA - log2 MRN + 1 - ", name_sp2),
    colName      = colName,
    shapeName    = shapeName,
    overwrite    = FORCE_OVERWRITE
)

# CMV - group single vs coinfection
nconds <- conds_sp3
nconds[grep("Single", conds_sp3)] <- "SingleInfection"
nconds[grep("CoInfection", conds_sp3, ignore.case = T)] <- "CoInfection"

pca_res_sp6 <- make_PCA_plot(
    file         = file.path(statsDir, name_sp3, paste0(name_sp3, "_PCA_GroupedInfection_rlog")),
    counts       = counts_sp3[, nconds != "none"],
    designMatrix = NA,
    conds        = nconds[nconds != "none"],
    shapes       = shapes_sp3[nconds != "none"],
    norm         = "rlog",
    main         = paste0("PCA - rlog - ", name_sp3),
    colName      = colName,
    shapeName    = shapeName,
    overwrite    = FORCE_OVERWRITE
)

pca_res_sp6 <- make_PCA_plot(
    file         = file.path("../results/other_plots/", paste0(name_sp3, "_PCA_GroupedInfection_mrn")),
    counts       = log2(mrn_sp3[, conds_sp3 != "none"]+1),
    designMatrix = NA,
    conds        = nconds[nconds != "none"],
    shapes       = shapes_sp3[nconds != "none"],
    norm         = "none",
    main         = paste0("PCA - log2 MRN + 1 - ", name_sp3),
    colName      = colName,
    shapeName    = shapeName,
    overwrite    = FORCE_OVERWRITE
)


pca_res_sp7 <- make_PCA_plot(
    file         = file.path(statsDir, name_sp3, paste0(name_sp3, "_PCA_grouped_none")),
    counts       = counts_sp3[, nconds != "none"],
    designMatrix = NA,
    conds        = nconds[nconds != "none"],
    shapes       = shapes_sp3[nconds != "none"],
    norm         = "none",
    main         = paste0("PCA - rlog - ", name_sp3),
    colName      = colName,
    shapeName    = shapeName,
    overwrite    = FORCE_OVERWRITE
)

pca_res_sp8 <- make_PCA_plot(
    file         = file.path(statsDir, name_sp3, paste0(name_sp3, "_PCA_none")),
    counts       = counts_sp3[, conds_sp3 != "none"],
    designMatrix = NA,
    conds        = conds_sp3[conds_sp3 != "none"],
    shapes       = shapes_sp3[conds_sp3 != "none"],
    norm         = "none",
    main         = paste0("PCA - none - ", name_sp3),
    colName      = colName,
    shapeName    = shapeName,
    overwrite    = FORCE_OVERWRITE
)

# differential expressed genes #################################################
pvalcut <- 0.01
log2cut <- NA
tools <- c("DESeq", "DESeq2", "edgeR", "limma")
colnames(counts) <- samples
deg_anno <- "" # <-- optional: fill in if you need additional annotation in DEG tables


# dual RNA-seq experiments (should) have two different design matrices
design_matrix_sp1 <- design_matrix
design_matrix_sp2 <- parse_SDRF("../doc/FungiNet_template_RNA-seq_V2.1_AFu.csv", skip = 8, col.name.row=7)$designMatrix
design_matrix_sp3 <- parse_SDRF("../doc/FungiNet_template_RNA-seq_V2.1_CMV.csv", skip = 8, col.name.row=7)$designMatrix


################
### Simple tests
deg_res_sp1 <- calculate_DEGs(
    counts       = counts_sp1,
    geneLengths  = gene_length_sp1,
    libSizes     = colSums(counts_sp1),
    designMatrix = design_matrix_sp1,
    pValCut      = pvalcut,
    logfcCut      = log2cut,
    tools        = tools,
    logfcnorm    = "mrn",
    outDir       = file.path(statsDir, name_sp1),
    prefix       = paste0(name_sp1, "_"),
    anno         = deg_anno,
    stop.on.error = FALSE,
    cpus          = MAX_CPUS,
    worker        = 10
)

deg_res_sp2 <- calculate_DEGs(
    counts       = counts_sp2,
    geneLengths  = gene_length_sp2,
    libSizes     = colSums(counts_sp2),
    designMatrix = design_matrix_sp2,
    pValCut      = 0.05,
    logfcCut      = log2cut,
    logfcnorm    = "mrn",
    tools        = tools,
    outDir       = file.path(statsDir, name_sp2),
    prefix       = paste0(name_sp2, "_"),
    anno         = deg_anno,
    stop.on.error = FALSE,
    cpus          = MAX_CPUS,
    worker        = 10
)


deg_res_sp2_nodeseq <- calculate_DEGs(
    counts       = counts_sp2,
    geneLengths  = gene_length_sp2,
    libSizes     = colSums(counts_sp2),
    designMatrix = design_matrix_sp2,
    pValCut      = 0.05,
    logfcCut      = log2cut,
    logfcnorm    = "mrn",
    tools        = c("DESeq2", "edgeR", "limma"),
    outDir       = file.path(statsDir, name_sp2, "nodeseq"),
    prefix       = paste0(name_sp2, "_"),
    anno         = deg_anno,
    stop.on.error = FALSE,
    cpus          = MAX_CPUS,
    worker        = 10
)


#############################
### additional tests for CMV - grouping
conds_sp3 <- conditions_from_design(design_matrix_sp3)
meta <- data.frame(samples = rownames(design_matrix_sp3), condition = conds_sp3)
design_matrix_sp3_2 <- createDesignMatrix(meta = meta,
                                        condCol = "condition",
                                        groupConds = list(
                                            coinfection = grep("coinfection", meta$condition, value = T, ignore.case = T),
                                            singleinfection = grep("singleinfection", meta$condition, value = T, ignore.case = T)
                                        ))$dm
# additional tests for CMV
design_matrix_sp3_3 <- design_matrix_sp3_2
colnames(design_matrix_sp3_3) <- "SingleInfection_CMV_0h_vs_SingleInfection_CMV_2h"
design_matrix_sp3_3[conds_sp3 == "SingleInfection_CMV_0h", 1] <- "treatment"
design_matrix_sp3_3[conds_sp3 == "SingleInfection_CMV_2h", 1] <- "control"

rownames(design_matrix_sp3_2) <- rownames(design_matrix_sp3)
design_matrix_sp3_2 <- design_matrix_sp3_2[,2, drop = F]

deg_res_sp3 <- calculate_DEGs(
    counts       = counts_sp3,
    geneLengths  = gene_length_sp3,
    libSizes     = colSums(counts_sp3),
    designMatrix = cbind(design_matrix_sp3, design_matrix_sp3_2, design_matrix_sp3_3),
    pValCut      = 0.05,
    logfcCut      = log2cut,
    logfcnorm    = "mrn",
    tools        = tools,
    outDir       = file.path(statsDir, name_sp3),
    prefix       = paste0(name_sp3, "_"),
    anno         = deg_anno,
    stop.on.error = FALSE,
    cpus          = MAX_CPUS,
    worker        = 10
)

saveRDS(deg_res_sp1, file = "../results/DEG_overlaps/deg_res_sp1.rds")
saveRDS(deg_res_sp2, file = "../results/DEG_overlaps/deg_res_sp2.rds")
saveRDS(deg_res_sp2_nodeseq, file = "../results/DEG_overlaps/deg_res_sp2_nodeseq")
saveRDS(deg_res_sp3, file = "../results/DEG_overlaps/deg_res_sp3.rds")



#############################
### additional tests for Hsapienss
library("tidyverse")

meta <- data.frame(samples = rownames(design_matrix_sp1), condition = conds_sp1)
meta$type <- "none"
meta$type[grepl("SingleCulture_DC", meta$condition, ignore.case = T)] <- "alone"
meta$type[grepl("SingleInfection_CMV", meta$condition, ignore.case = T)] <- "si_CMV"
meta$type[grepl("SingleInfection_AFU", meta$condition, ignore.case = T)] <- "si_AFU"
meta$type[grepl("CoInfection_CMV_first", meta$condition, ignore.case = T)] <- "co_CMV"
meta$type[grepl("Coinfection_Afu_first", meta$condition, ignore.case = T)] <- "co_AFU"
meta$type[meta$condition == "CoInfection"] <- "co_Both"
meta$type <- factor(meta$type)
meta$type <- relevel(meta$type, "alone")

meta$growth <- "none"
meta$growth[grepl("SingleCulture_DC", meta$condition)] <- "0h"
meta$growth[grepl("(0h)", meta$condition)] <- "0h"
meta$growth[grepl("(2h)", meta$condition)] <- "2h"
meta$growth[grepl("(4h30min)", meta$condition)] <- "4h30min"
meta$growth[meta$type %in% c("co_Both", "co_AFU", "co_CMV")] <- "0h"
meta$growth <- factor(meta$growth)
meta$growth <- relevel(meta$growth, "0h")

meta$donor <- gsub("_.*", "", meta$samples)
meta$donor <- factor(meta$donor)

design_matrix_sp1_2 <- createDesignMatrix(meta = meta,
                                          condCol = "type",
                                          condVal = "si_CMV",
                                          timeCol = "growth",
                                          timeVal = "0h",
                                          mode = "CombFixCond")$dm
design_matrix_sp1_2 <- design_matrix_sp1_2[,grepl("(si-CMV)|(si-AFU)", colnames(design_matrix_sp1_2))]
design_matrix_sp1_2 <- design_matrix_sp1_2[,grep("(none)|(si-AFU_2h)|(si-CMV_4h30min)", colnames(design_matrix_sp1_2), invert = T)]
rownames(design_matrix_sp1_2) <- rownames(design_matrix_sp1)

deg_res_sp1_2 <- calculate_DEGs(
    counts       = counts_sp1,
    geneLengths  = gene_length_sp1,
    libSizes     = colSums(counts_sp1),
    designMatrix = design_matrix_sp1_2,
    pValCut      = pvalcut,
    logfcCut      = log2cut,
    tools        = tools,
    logfcnorm    = "mrn",
    outDir       = file.path(statsDir, name_sp1),
    prefix       = paste0(name_sp1, "_"),
    anno         = deg_anno,
    stop.on.error = FALSE,
    cpus          = MAX_CPUS,
    worker        = 10
)


# CoInfection vs both single infections using all growth times
design_matrix_sp1_3 <- createDesignMatrix(meta = meta,
                                          condCol = "type",
                                          condVal = "CoI",
                                          groupConds = list("CoI" = c("co_CMV", "co_AFU", "co_Both"),
                                                            "SI" = c("si_CMV", "si_AFU")))
rownames(design_matrix_sp1_3$dm) <- design_matrix_sp1_3$meta$samples
design_matrix_sp1_3 <- design_matrix_sp1_3$dm[,2, drop=F]

deg_res_sp1_3 <- calculate_DEGs(
    counts       = counts_sp1,
    geneLengths  = gene_length_sp1,
    libSizes     = colSums(counts_sp1),
    designMatrix = design_matrix_sp1_3,
    pValCut      = pvalcut,
    logfcCut      = log2cut,
    tools        = tools,
    logfcnorm    = "mrn",
    outDir       = file.path(statsDir, name_sp1),
    prefix       = paste0(name_sp1, "_"),
    anno         = deg_anno,
    stop.on.error = FALSE,
    cpus          = 4,
    worker        = 10
)

saveRDS(deg_res_sp1_3, "../results/DEG_overlaps/deg_res_sp1_3.rds")





# CoInfection vs both single infections, but only at 0h growth
design_matrix_sp1_4 <- createDesignMatrix(meta = meta,
                                          condCol = "type",
                                          condVal = "CoI",
                                          timeCol = "growth",
                                          timeVal = "0h",
                                          mode = "CombFixTime",
                                          groupConds = list("CoI" = c("co_CMV", "co_AFU", "co_Both"),
                                                            "SI" = c("si_CMV", "si_AFU")))
rownames(design_matrix_sp1_4$dm) <- design_matrix_sp1_4$meta$samples
design_matrix_sp1_4 <- design_matrix_sp1_4$dm[, "CoI_0h_VS_SI_0h", drop = F]
# fix the issues that CoI 4.5h is considered 0h
design_matrix_sp1_4[grepl("4h30min", rownames(design_matrix_sp1_4)), 1] <- "none"
colnames(design_matrix_sp1_4)

deg_res_sp1_4 <- calculate_DEGs(
    counts       = counts_sp1,
    geneLengths  = gene_length_sp1,
    libSizes     = colSums(counts_sp1),
    designMatrix = design_matrix_sp1_4,
    pValCut      = pvalcut,
    logfcCut      = log2cut,
    tools        = tools,
    logfcnorm    = "mrn",
    outDir       = file.path(statsDir, name_sp1),
    prefix       = paste0(name_sp1, "_"),
    anno         = deg_anno,
    stop.on.error = FALSE,
    cpus          = 4,
    worker        = 10
)
saveRDS(deg_res_sp1_4, "../results/DEG_overlaps/deg_res_sp1_4.rds")




#####################
### Update Meta-Data

# it's not a setting for the DEG tools directly,
# hence writing pvalcut and log2cut to Description and not Parameters.And.Values
software_nr <- software_nr + 1
col <- paste0("Description.", software_nr)
SDRF[, col] <- deg_res_sp1$desc
col <- paste0("Data.Processing.Software.", software_nr)
SDRF[, col] <- paste(deg_res_sp1$tools, collapse = ", ")
col <- paste0("Storage.Date.", 3)
SDRF[, col] <- Sys.Date()

# fill in! <-- Archive directory of the final results
# SDRF$Storage.Location.3 <- "stored at HKI Backupserver" # old SDRF version
SDRF$Storage.Location...Server.1 <- "stored at HKI Jena: /sbiarchiv" # according to FungiNet Lookup-Table
SDRF$Storage.Location...Folder.1 <- "/sbiarchiv/xxx"
SDRF$Derived.Data.File <- "results.zip"

# save results #################################################################

# save all R objects
save.image("R_workspace.RData")

# TODO dual: save to two SDRF files???
# --> if yes: fill second SDRF right from the beginning?

# save updated (filled) SDRF table
write.table(
    SDRF,
    file      = paste0(tools::file_path_sans_ext(SDRF_file), "_final.csv"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
)

# zip important result files, e.g. for upload to FungiNetDB
archive(zipFile = "results.zip")

# optional: cleanup
system(paste("rm", "-r", mapDir, countDir, fastqDir, qualityDir))

# TODO standardize names of directory variables?
# mapDir, countDir, fastqDir, qualityDir?

# optional: move remaining files to archive directory
system("mv * xxx") # <-- fill in!
