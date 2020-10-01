# manuscript_tripleRNAseq

Here, code to derive the main conclusions described in Seelbinder et al. (2020),"Triple RNA-seq reveals synergy in a human virus-fungus co-infection model", Cell reports (*accepted*) is presented.

# Contents
* **innatedb_curated_genes_2019-04-11.xls :** A list of genes (HGNC symbol) present in InnateDB at the time of processing (2019-04-11)
* **00_R_Source/ :** Source code used bei GEO2RNAseq to set pre-processing beginning with quality control and ending with statistical analysis of differentially expressed genes
* **01_R_PreProcessing/ :** Contains the main script to execute pre-processing
* **02_R_Analyses/ :** Contains scripts to analyse pre-processing output. It includes code for:
	* Gene correlation network analysis, including figures
	* Gene variance analysis, including heatmap figures
	* Read library composition analysis, including figures
	* Comparing statistical results across different experimental conditions
* **03_advanced_plots/ :** Contains the code for more complex plots. It includes code for:
	* visualization of RNA types
	* visualization of RNA-seq and qPCR data
	* visualization of cytokine measurements
	* KEGG pathway visualization
	* Principle Component Analysis across different experimental conditions
