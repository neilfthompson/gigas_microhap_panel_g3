# gigas_microhap_panel

#### Requirements 
- [amplitools](https://github.com/bensutherland/amplitools)       
- [simple_pop_stats](https://github.com/bensutherland/simple_pop_stats/tree/master)     


Microhaplotype analysis of marker panel developed by Sutherland et al. using new populations from USDA and existing samples from initial panel development effort.

Details of the initial panel can be found in Sutherland et al. (2023)


### Genotyping with standard hotspot (single SNP)
To analyze the VariantCaller output considering only the hotspot SNP (i.e., the target variant), follow the instructions on [amplitools](https://github.com/bensutherland/amplitools), as outlined in brief below.      
Instructions are also present in the following `gigas_microhap_panel/01_scripts/single_variant_hotspot_runall.R`       

Put the VariantCaller output file (i.e., `R_*USDA_OYSTER_*.xls`) in `02_input data`, then use the initiator to convert the tab-delimited file to a genepop file, including the flag `hotspot_only` as TRUE.      

Once the prepped matrix has been prepared, using the script `01_scripts/format_genepop.sh 02_input_data/prepped_matrices/<filename>` to convert the prepped matrix to a genepop.      

Copy the output genepop to the input folder of `simple_pop_stats`.    

Interactively run the script `gigas_microhap_panel/01_scripts/single_variant_hotspot_analysis_sps.R` for instructions on how to use amplitools and simple pop stats to analyze the data. This script will:      
- load the genepop
- assign population variable to each individual
- check for missing data and filter
- drop monomorphic or low MAF variants
- generate per locus statistics (i.e., mean FST and HOBS)
- plot samples in a PCA based on individual genotypes 
- convert to Rubias format
- prepare a CKMR object
- assess cutoffs for CKMR for sibship (offspring, parents) or parent-offspring relationships
- output putative parent-offspring or sibship pairs

