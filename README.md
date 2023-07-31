# gigas_microhap_panel

#### Requirements 
- [amplitools](https://github.com/bensutherland/amplitools)

Microhaplotype analysis of marker panel developed by Sutherland et al. using new populations from USDA and existing samples from initial panel development effort.

Details of the initial panel can be found in Sutherland et al. (2023)


### Genotyping with standard hotspot (single SNP)
To analyze the VariantCaller output considering only the hotspot SNP (i.e., the target variant), follow the instructions on [amplitools](https://github.com/bensutherland/amplitools), as outlined in brief below.      

Put the VariantCaller output file (i.e., `R_*USDA_OYSTER_*.xls`) in `02_input data`, then use the initiator to convert the tab-delimited file to a genepop file, including the flag `hotspot_only` as TRUE.      
