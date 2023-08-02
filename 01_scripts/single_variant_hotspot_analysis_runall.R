# Analyze the variantCaller output using amplitools
#  this script is a runall containing instructions on how to analyze the data using amplitools
#  please see amplitools instructions for full details (github.com/bensutherland/amplitools)
#  initialized 2023-07-31 (B. Sutherland)

# Put input file in 02_input_data

# Source the 00_initiator.R 

# Convert input to a prepped matrix
proton_to_genepop(hotspot_only = TRUE, neg_control = "BLANK")

# Go to terminal and finalize the genepop using
# 01_scripts/format_genepop.sh 02_input_data/prepped_matrices/R_2023_07_03_12_36_08_user_S5XL-00533-1229-USDA_OYSTER_20230702_gen_data.txt

# Copy the output 
# 02_input_data/prepped_genepops/R_2023_07_03_12_36_08_user_S5XL-00533-1229-USDA_OYSTER_20230702_gen_data.gen
# to simple_pop_stats input folder 02_input_data

# Move to the next script, gigas_microhap_panel/01_scripts/single_variant_hotspot_analysis_sps.R