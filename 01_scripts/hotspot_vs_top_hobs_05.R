# hotspot vs. top hobs, parentage
# note: requires that part 4 was run and rubias files are present in amplitools
# Ben Sutherland
# 2024-03-23

# Clear space
# rm(list=ls())
# Source amplitools

# Set user variables
hotspot_rubias.FN <-"03_results/rubias_hotspot_212_ind_372_loc_2024-03-23.txt"
top_hobs_rubias.FN <-"03_results/rubias_top_hobs_212_ind_372_loc_2024-03-23.txt"
parent_pop <- "VIU_F1"
offspring_pop <- "VIU_F2"
cutoff <- 5

## Parentage with hotspot comparables
# Prepare run folder
date <- format(Sys.time(), "%Y-%m-%d")
input_rubias_short.FN <- gsub(pattern = "03_results/", replacement = "", x = hotspot_rubias.FN)
input_rubias_short.FN <- gsub(pattern = ".txt", replacement = "", x = input_rubias_short.FN)

run_folder.FN <- paste0("03_results/ckmr_input_", input_rubias_short.FN, "_"
                        , offspring_pop, "_vs_", parent_pop,"_", date
)
print("Making new result folder...")
print(run_folder.FN)
dir.create(run_folder.FN)

# Run ckmr on the input file
ckmr_from_rubias(input.FN = hotspot_rubias.FN
                 , parent_pop = parent_pop
                 , offspring_pop = offspring_pop
                 , cutoff = 5
                 , output.dir = run_folder.FN
)

## Parentage with top hobs comparables
# Prepare run folder
date <- format(Sys.time(), "%Y-%m-%d")
input_rubias_short.FN <- gsub(pattern = "03_results/", replacement = "", x = top_hobs_rubias.FN)
input_rubias_short.FN <- gsub(pattern = ".txt", replacement = "", x = input_rubias_short.FN)

run_folder.FN <- paste0("03_results/ckmr_input_", input_rubias_short.FN, "_"
                        , offspring_pop, "_vs_", parent_pop,"_", date
)
print("Making new result folder...")
print(run_folder.FN)
dir.create(run_folder.FN)

# Run ckmr on the input file
ckmr_from_rubias(input.FN = top_hobs_rubias.FN
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
