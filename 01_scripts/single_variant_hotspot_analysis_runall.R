# Analyze the variantCaller output using amplitools
#  this script is a runall containing instructions on how to analyze the data using amplitools
#  please see amplitools instructions for full details (github.com/bensutherland/amplitools)
#  initialized 2023-07-31 (B. Sutherland)

# Put input file in 02_input_data

# Source the 00_initiator.R 

# Convert input to a prepped matrix
proton_to_genepop(hotspot_only = TRUE, neg_control = "BLANK")

# Go to terminal and finalize the genepop using
# 01_scripts/format_genepop.sh

# Copy the output to simple_pop_stats input folder
# once complete, you will have a file 'amplitools/03_results/cgig_all_rubias.txt'

# source 00_initiator.R again
