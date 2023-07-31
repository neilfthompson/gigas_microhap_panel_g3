# Analyze the variantCaller output using amplitools
#  this script is a runall containing instructions on how to analyze the data using amplitools
#  please see amplitools instructions for full details (github.com/bensutherland/amplitools)
#  initialized 2023-07-31 (B. Sutherland)

# Source the 00_initiator.R 

proton_to_genepop(hotspot_only = TRUE, neg_control = "BLANK")
