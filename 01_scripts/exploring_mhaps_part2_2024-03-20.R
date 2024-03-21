# Exploring mhaps; requires that amplitargets is cloned and at same level as amplitools
# Ben Sutherland
# 2024-03-20

# Clear space
# rm(list=ls())
# Source simple_pop_stats

# Prepared data from exploring_mhaps part 1
load(file = "prepared_filtered.RData")

head(comparison.df)
obj

# Locus names?
locNames(obj) # obj locus names must match the comparison, but after that should be adjusted to not interfere with downstream processes

# Individual names? 
indNames(obj) # keeping these as the indiv name would be beneficial for downstream interpretation, but will need to rebuild popmap
simplify_names(df = obj, format = "amplitools")
obj <- obj_simplified
head(indNames(obj))


#### 01. Create object subsets for hotspot or top hobs ####
# Create hotspot object
obj_hotspot <- obj[, loc=comparison.df$hotspot_locus]
obj_hotspot
# now adjust locnames
locNames(obj_hotspot) <- gsub(pattern = "_", replacement = "", x = locNames(obj_hotspot))
locNames(obj_hotspot) <- gsub(pattern = ":", replacement = "", x = locNames(obj_hotspot))
locNames(obj_hotspot) <- gsub(pattern = "-", replacement = "", x = locNames(obj_hotspot))
head(locNames(obj_hotspot))

# Create top hobs object
obj_top_hobs <- obj[,loc=comparison.df$max_hobs_id]
obj_top_hobs
head(locNames(obj_top_hobs))

# now adjust locnames
locNames(obj_top_hobs) <- gsub(pattern = "_", replacement = "", x = locNames(obj_top_hobs))
locNames(obj_top_hobs) <- gsub(pattern = ":", replacement = "", x = locNames(obj_top_hobs))
locNames(obj_top_hobs) <- gsub(pattern = "-", replacement = "", x = locNames(obj_top_hobs))


#### 02. Create stock code file ####
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


#### 03. Recreate popmap, as we have adjusted to the true individual identifier for simplicity
# Note: the original popmap will not work, as we have switched to alt.id
head(indNames(obj_top_hobs))

popmap.df <- cbind(indNames(obj), as.character(pop(obj)))
popmap.df <- as.data.frame(popmap.df)
colnames(popmap.df) <- c("indiv", "pop")
head(popmap.df)
popmap.FN <- "02_input_data/popmap_for_rubias_conversion.txt"
write_delim(x = popmap.df, file = popmap.FN, delim = "\t")


# # Convert back to original form
# obj_top_hobs_inds.df <- as.data.frame(indNames(obj_top_hobs))
# colnames(obj_top_hobs_inds.df) <- "amplitools.ID"
# head(obj_top_hobs_inds.df)
# obj_top_hobs_inds.df <- separate(data = obj_top_hobs_inds.df, col = "amplitools.ID", into = c("run", "barcode", "indiv"), sep = "__", remove = F)
# head(obj_top_hobs_inds.df)
# 
# inds <- paste0(obj_top_hobs_inds.df$run, "_", obj_top_hobs_inds.df$barcode, ".fastq")
# head(indNames(obj_top_hobs))
# tail(indNames(obj_top_hobs))
# indNames(obj_top_hobs) <- inds
# head(indNames(obj_top_hobs))
# tail(indNames(obj_top_hobs))

#temp <- obj_top_hobs

# now it will match the popmap again
#locNames(temp) <- paste0("A", seq(1:nLoc(temp)))

#### ISSUE ####
# I think because the dataset was read in as a VCF, the suffix is .1 not .01 for the second allele
data <- obj_top_hobs$tab
data[1:5,1:7]

#obj <- obj_hotspot; prefix <- "hotspot"
obj <- obj_top_hobs; prefix <- "top_hobs"

# Note: will need to source an updated genepop_to_rubias_SNP.r because it is using .1 instead of .01
popmap.FN
genepop_to_rubias_SNP(data = obj, sample_type = "reference"
                      , custom_format = TRUE, micro_stock_code.FN = micro_stock_code.FN
                      , pop_map.FN = popmap.FN
)
print("Your output is available as '03_results/rubias_output_SNP.txt")

# Obtain some variables to create filename
date <- format(Sys.time(), "%Y-%m-%d")
indiv.n <- nInd(obj)
loci.n  <- nLoc(obj)
rubias_custom.FN <- paste0("03_results/rubias_", prefix, "_", indiv.n, "_ind_", loci.n, "_loc_", date, ".txt")

# Copy to rename rubias output file
file.copy(from = "03_results/rubias_output_SNP.txt", to = rubias_custom.FN, overwrite = T)
print(paste0("The output is saved as ", rubias_custom.FN))

# Go back and do the alternate one too

# Save image
save.image("03_results/post-filters_prepared_for_parentage_rubias_built.RData")

print("Here you need to copy the above rubias file to amplitools results folder.")



#### 02. CKMR-sim with above
# Clear space, and launch amplitools initiator (i.e., 01_scripts/00_initiator.R)

# Set user variables
#input_rubias.FN <- "03_results/rubias_hotspot_86_ind_337_loc_2024-03-20.txt"
input_rubias.FN <- "03_results/rubias_top_hobs_86_ind_337_loc_2024-03-20.txt"

parent_pop <- "VIU_F1"
offspring_pop <- "VIU_F2"
cutoff <- 5

# Prepare run folder
date <- format(Sys.time(), "%Y-%m-%d")
input_rubias_short.FN <- gsub(pattern = "03_results/", replacement = "", x = input_rubias.FN)
input_rubias_short.FN <- gsub(pattern = ".txt", replacement = "", x = input_rubias_short.FN)

run_folder.FN <- paste0("03_results/ckmr_input_", input_rubias_short.FN, "_"
                        , offspring_pop, "_vs_", parent_pop,"_", date
)
print("Making new result folder...")
print(run_folder.FN)
dir.create(run_folder.FN)


# Run ckmr on the input file
ckmr_from_rubias(input.FN = input_rubias.FN
                 , parent_pop = parent_pop
                 , offspring_pop = offspring_pop
                 , cutoff = 5
                 , output.dir = run_folder.FN
)

# # Filtered loci and pilot study null allele removed
# ckmr_from_rubias(input.FN = "03_results/rubias_output_SNP_filtered_and_null_pilot_drop.txt", parent_pop = "F0"
#                  , offspring_pop = "F1", cutoff = 5
# )

# Plot the output results
# graph_relatives(input.FN = "03_results/ckmr_input_rubias_135_ind_343_loc_2024-02-26_F1_vs_F0_2024-02-26/po_F0_vs_F1_pw_logl_5.txt", logl_cutoff = 5
#                 , drop_string = "", directed = F, plot_width = 8, plot_height = 8
# )

# # Before going to the next script, set some variables
# report.FN <- "03_results/parentage_with_all_filtered_loci/po_F0_vs_F1_pw_logl_5_report.txt"

save.image(paste0(run_folder.FN, "/ckmr_completed.RData"))

# Now go to 01_scripts/exploring_families.R






### OTHER TODO
# - multimappers number snp per locus
# - will need to indicate the size of the amplicons by average
#### POSSIBLE OTHER TODO
# - try to solve the issue of so few of the hotspots coming through (reduce filters?)


#### FINALLY NEXT: 
# genotype microhaps and then compare against the best SNP
