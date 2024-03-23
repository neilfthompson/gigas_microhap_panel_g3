# Compare hotspot SNP to top HOBS SNP, part 2
# note: requires that part 1 was run and '03_results/hotspot_data.txt' is present in amplitools
# Ben Sutherland
# initialized 2024-03-18

# Clear space
# rm(list=ls())

# Source amplitools via the initiator

# Load libraries
library(vcfR)

# Set user variables
vcf_filtered.FN            <- "14_extract_mhap/mpileup_calls_biallele_snp_q20_dp10_under_0.1_missing_w_AF_0.01.vcf" # filtered VCF
vcf_unfiltered.FN          <- "14_extract_mhap/mpileup_calls_biallele_snp.vcf" # unfiltered VCF
hotspot_data.FN <- "03_results/hotspot_data.txt" # created in part 1, has empirical names for expected hotspots

#### 01. Identify hotspot variants within the empirical VCF ####
## Load the hotspot data
hotspot_data.df <- read.delim2(file = hotspot_data.FN, header = T, sep = "\t")

## Load data unfiltered VCF to see how many hotspots were identified overall
my_vcf <- read.vcfR(file = vcf_unfiltered.FN)

# Convert vcf to genind
my_data.genind <- vcfR2genind(x = my_vcf)

# How many of the hotspot file sites were observed in the empirical data? 
empirical_sites.vec <- locNames(my_data.genind)
length(intersect(x = hotspot_data.df$expect.name, y = empirical_sites.vec)) # 437 in both

## Load data with filtered VCF to proceed (will match mhaps)
my_vcf <- read.vcfR(file = vcf_filtered.FN)

# Convert vcf to genind
my_data.genind <- vcfR2genind(x = my_vcf)

# How many of the hotspot file sites were observed in the empirical data? 
empirical_sites.vec <- locNames(my_data.genind)
length(intersect(x = hotspot_data.df$expect.name, y = empirical_sites.vec)) # 373 in both

# Which variants that are in the hotspot file were also empirically identified? 
comparison_loci <- intersect(x = hotspot_data.df$expect.name, y = empirical_sites.vec)
comparison_loci.df <- as.data.frame(comparison_loci)

write.table(x = comparison_loci.df, file = "03_results/comparison_loci.txt", col.names = F, row.names = F
            , sep = "\t"
            )

# NOTE: here should make sure the logic from step 1 didn't cause some hotspots to be missed (TODO)
#### NOTE: HERE COULD ALSO WHITTLE DOWN TO AMPLICON NAME, then this could be used for mhap too ####

# Next: go to hotspot_vs_top_hobs_03_.R
