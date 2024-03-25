# Exploring mhaps; requires that amplitargets is cloned and at same level as amplitools
# note: requires that part 3 was run and '03_results/prepared_filtered.RData' is present in simple_pop_stats
# Ben Sutherland
# 2024-03-20

# Clear space
# rm(list=ls())
# Source simple_pop_stats

# Prepared data from previous step
load(file = "03_results/prepared_filtered.RData")

head(comparison_loci.df) # amplicons in which the hotspot was empirically genotyped (comparables)
obj  # all loci

# Locus names?
head(as.data.frame(locNames(obj))) # obj locus names must match the comparison, but after that should be adjusted to not interfere with downstream processes

# Individual names? 
head(indNames(obj)) 


#### 01. Create object subsets for hotspot or top hobs ####
# Create hotspot object
obj_hotspot <- obj[, loc=comparison_loci.df$hotspot_locus]
obj_hotspot

# Remove atypical characters to not disrupt downstream analysis
locNames(obj_hotspot) <- gsub(pattern = "_", replacement = "", x = locNames(obj_hotspot))
locNames(obj_hotspot) <- gsub(pattern = ":", replacement = "", x = locNames(obj_hotspot))
locNames(obj_hotspot) <- gsub(pattern = "-", replacement = "", x = locNames(obj_hotspot))
head(locNames(obj_hotspot))


# Create top hobs object
obj_top_hobs <- obj[,loc=comparison_loci.df$max_hobs_id]
obj_top_hobs
head(locNames(obj_top_hobs))

# Remove atypical characters to not disrupt downstream analysis
locNames(obj_top_hobs) <- gsub(pattern = "_", replacement = "", x = locNames(obj_top_hobs))
locNames(obj_top_hobs) <- gsub(pattern = ":", replacement = "", x = locNames(obj_top_hobs))
locNames(obj_top_hobs) <- gsub(pattern = "-", replacement = "", x = locNames(obj_top_hobs))
head(locNames(obj_top_hobs))


#### 02. To make a rubias file, will need to create a stock code file ####
# Need to create a tab-delim stock code file in format of e.g., 
## row 1: collection	repunit
## row 2: boundary_bay	lower_mainland

# Here we will just create a df based on existing populations where collection = repunit
stock_code.df <- as.data.frame(unique(pop(obj)))
colnames(stock_code.df) <- "collection"
stock_code.df$repunit <- stock_code.df$collection
stock_code.df

# Write it out
write_delim(x = stock_code.df, file = "00_archive/stock_code.txt", delim = "\t", col_names = T)
micro_stock_code.FN <- "00_archive/stock_code.txt"
datatype <- "SNP" # required for genepop_to_rubias_SNP
# this is for annotate_rubias(), for an unknown reason it requires the name micro_stock_code.FN

# Also need 'indiv \t pop' popmap
annot_popmap_original.df <- read.delim2(file = "00_archive/my_data_ind-to-pop_annot.txt")
head(annot_popmap_original.df)
annot_popmap_original.df <- annot_popmap_original.df[, c("alt.ID", "pop")]
colnames(annot_popmap_original.df) <- c("indiv", "pop")
head(annot_popmap_original.df)
write.table(x = annot_popmap_original.df, file = "00_archive/sample_name_popmap.txt", quote = F, sep = "\t", row.names = F, col.names = T)

# Warning: I think because the dataset was read in as a VCF, the suffix is .1 not .01 for the second allele! (mhaps may be different)
data <- obj_top_hobs$tab
data[1:5,1:7]


##### Make rubias files #####
# Note: will need to source an updated genepop_to_rubias_SNP.r because it is using .1 instead of .01 (### TODO FIX ###)
popmap.FN <- "00_archive/sample_name_popmap.txt"

## hotspot
genepop_to_rubias_SNP(data = obj_hotspot, sample_type = "reference"
                      , custom_format = TRUE, micro_stock_code.FN = micro_stock_code.FN
                      , pop_map.FN = popmap.FN
)

# Obtain some variables to create filename
date <- format(Sys.time(), "%Y-%m-%d")
indiv.n <- nInd(obj_hotspot)
loci.n  <- nLoc(obj_hotspot)
rubias_custom.FN <- paste0("03_results/rubias_hotspot", "_", indiv.n, "_ind_", loci.n, "_loc_", date, ".txt")

# Copy to rename rubias output file
file.copy(from = "03_results/rubias_output_SNP.txt", to = rubias_custom.FN, overwrite = T)
print(paste0("The output is saved as ", rubias_custom.FN))

## top HOBS
genepop_to_rubias_SNP(data = obj_top_hobs, sample_type = "reference"
                      , custom_format = TRUE, micro_stock_code.FN = micro_stock_code.FN
                      , pop_map.FN = popmap.FN
)

# Obtain some variables to create filename
date <- format(Sys.time(), "%Y-%m-%d")
indiv.n <- nInd(obj_hotspot)
loci.n  <- nLoc(obj_hotspot)
rubias_custom.FN <- paste0("03_results/rubias_top_hobs", "_", indiv.n, "_ind_", loci.n, "_loc_", date, ".txt")

# Copy to rename rubias output file
file.copy(from = "03_results/rubias_output_SNP.txt", to = rubias_custom.FN, overwrite = T)
print(paste0("The output is saved as ", rubias_custom.FN))

# Save image
save.image("03_results/post-filters_prepared_for_parentage_rubias_built.RData")
print("Here you need to copy the above rubias file to amplitools results folder.")

file.copy(from = "03_results/rubias_hotspot_212_ind_372_loc_2024-03-23.txt", to = "../amplitools/03_results/")
file.copy(from = "03_results/rubias_top_hobs_212_ind_372_loc_2024-03-23.txt", to = "../amplitools/03_results/")

### Next, go to hotspot_vs_top_hobs_05

