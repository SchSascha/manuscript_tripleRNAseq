#!/bin/bash
CHROM_SIZES=/sbidata/shared2/Combined_Genomes_forTripleRNAseqPreprocessing/Hsapiens_GRCh38_89_2017-05-07_Human_herpesvirus_5_strain_TB40E_clone_TB40-BAC4/fasta/hs-GR38.89_afu-Af293_cmv-5.tab
parallel --plus -j 2 "bedtools genomecov -bg -ibam {} -g $CHROM_SIZES > ../results/coverage/{/..}.bedgraph" ::: ../results/mapping/bamfiles/*.sort.bam

# covert to bigwig format
# bedgraph must be case sensitive sorted first
export LC_COLLATE=C
parallel "sort -k1,1 -k2,2n {} > {.}.sorted.bedgraph" ::: ../results/coverage/*.bedgraph
parallel --plus -j 2 "bedGraphToBigWig {} $CHROM_SIZES {..}.bw" ::: ../results/coverage/*.sorted.bedgraph


