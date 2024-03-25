# Take output of mtype2 and prepare it for CKMR-sim
# Ben Sutherland (SBIO)
# 2024-03-21

# Clear space
# rm(list=ls())

# Source amplitools initiator

# Load libraries
#install.packages("tidyverse")
devtools::install_github("delomast/EFGLmh")
library(tidyverse)
library(EFGLmh)

# Set options
options(tibble.max_extra_cols = 10)


# Set variables
#annot.FN <- "~/Documents/00_sbio/GBMF_UBC_Pacific_oyster/amplicon_panel_mhap/simple_pop_stats/00_archive/my_data_ind-to-pop_annot.txt"
barcode_interp.FN <- "~/Documents/00_sbio/GBMF_UBC_Pacific_oyster/amplicon_panel_mhap/00_archive/GBMF_pilot_sample_ID_barcode_run_from_variantCaller_output_2023-07-27.txt"
sample_interp.FN <- "~/Documents/00_sbio/GBMF_UBC_Pacific_oyster/amplicon_panel_mhap/00_archive/GBMF_pilot_my_data_ind-to-pop_annot_2023-05-01.txt"


comp_loci.df <- read.table(file = "03_results/comparison_loci.txt", header = F)
head(comp_loci.df)
comp_loci.df <-separate(data = comp_loci.df, col = "V1", into = c("contig", "middle", "SNP.pos")
                    , sep = "_", remove = F)
head(comp_loci.df)
comp_loci.df$amplicon <- paste0(comp_loci.df$contig, "_", comp_loci.df$middle)
head(comp_loci.df)
comp_loci.df$amplicon <- gsub(pattern = "_|:|-", replacement = "", x = comp_loci.df$amplicon)

#### 01. Prepare mtype2 output for EFGLmh input ####
#mhap.FN <- "14_extract_mhap/genos.txt" # full mhaps
mhap.FN <- "14_extract_mhap/genos_hotspots.txt" # hotspot SNP only
mhap.df <- read.delim(file = mhap.FN, header = T)
dim(mhap.df)
head(mhap.df)
mhap.df$Indiv <- gsub(pattern = ".*/", replacement = "", x = mhap.df$Indiv)
mhap.df$Indiv <- gsub(pattern = ".sorted.bam", replacement = "", x = mhap.df$Indiv)
head(mhap.df)

output_mtype2.example <- mhap.df
head(output_mtype2.example)


#### 02. Annotate population and update individual names 
barcode_interp.df <- read.delim2(file = barcode_interp.FN)
dim(barcode_interp.df) # 575 lines
barcode_interp.df$matcher <- paste0(barcode_interp.df$Run.Name, "_", barcode_interp.df$Barcode)
head(barcode_interp.df)
barcode_interp.df <- barcode_interp.df[,c("Sample.Name", "matcher")]

sample_interp.df <- read.delim2(file = sample_interp.FN)
dim(sample_interp.df)  # 370 lines
head(sample_interp.df)


mhap.df <- merge(x = mhap.df, y = barcode_interp.df, by.x = "Indiv", by.y = "matcher", all.x = T)
head(mhap.df)
mhap.df$Indiv <- mhap.df$Sample.Name
unique(mhap.df$Indiv)
head(mhap.df)

mhap.df <- merge(x = mhap.df, y = sample_interp.df, by.x = "Indiv", by.y = "indiv", all.x = T)
head(mhap.df)

head(output_mtype2.example)



mhap.df <- mhap.df[,c("Indiv", "Locus", "Allele1", "Allele2", "Allele1_count", "Allele2_count", "p")]
head(mhap.df)


# Convert from long form to wide form
mhap.tbl <- as_tibble(mhap.df)
wideFormat <- mtype2wide(x = mhap.tbl)  
dim(wideFormat)
as.data.frame(wideFormat[1:5,1:5])
head(wideFormat$Indiv)

wideFormat.df <- as.data.frame(wideFormat)


# Bring the population ID in
head(sample_interp.df)
#wideformat.df <- merge(x = wideFormat.df, y = sample_interp.df, by.x = "Indiv", by.y = "indiv", all.x = T)
wideformat.df <- merge(x = sample_interp.df, y = wideFormat.df, by.x = "indiv", by.y = "Indiv", all.y = T)
wideformat.df[1:10, 1:10]
wideformat.df$Pop <- wideformat.df$pop
wideformat.df <- wideformat.df[-which(colnames(wideformat.df)=="pop")]
wideformat.df[1:10, 1:10]



# And finally...
mh.dat <- readInData(input = wideformat.df, pedigreeColumn = "Pop"
                    , nameColumn = "indiv") # into EFGLmh format
mh.dat

allele_rich <- aRich(mh.dat)
head(allele_rich)
hist(table(allele_rich$aRich))

### IDENTIFY MULTIMAPPERS HERE ####
str(mh.dat)

# Remove bad loci (too much missing?)
# Remove bad inds (too much missing?)

head(getLoci(x = mh.dat))
head(comp_loci.df)
# Which loci are in the mh.dat but not in the comparable loci? 
remove <- setdiff(x = getLoci(x = mh.dat), y = comp_loci.df$amplicon) 
mh.dat <- removeLoci(mh.dat, lociRemove = remove)

nloc <- length(getLoci(mh.dat))

# , then convert to rubias
parentage_baseline <- exportRubias_baseline(x = mh.dat, pops = c("VIU_F1", "VIU_F2"), repunit = "Pop", collection = "Pop")
parentage_baseline[1:5,1:5]



#write_delim(x = parentage_baseline, file = paste0("03_results/rubias_mhaps_", nloc, "_2024-03-23.txt"), delim = "\t")
write_delim(x = parentage_baseline, file = paste0("03_results/rubias_mhaps_hotspot", nloc, "_2024-03-23.txt"), delim = "\t")

# Set user variables
hotspot_rubias.FN <-"03_results/rubias_mhaps_hotspot373_2024-03-23.txt"
parent_pop <- "VIU_F1"
offspring_pop <- "VIU_F2"
cutoff <- 5

## Parentage with hotspot comparables
# Prepare run folder
date <- format(Sys.time(), "%Y-%m-%d")
input_rubias_short.FN <- gsub(pattern = "03_results/", replacement = "", x = hotspot_rubias.FN)
input_rubias_short.FN <- gsub(pattern = ".txt", replacement = "", x = input_rubias_short.FN)

run_folder.FN <- paste0("03_results/ckmr_input_", input_rubias_short.FN, "_"
                        , offspring_pop, "_vs_", parent_pop,"_", nloc, "_", date
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


# , then amplitools ckmr-sim wrapper
# May come across issues with the naming of .A1 and .A2, but this should be easy to fix, e.g., to .01 and .02


