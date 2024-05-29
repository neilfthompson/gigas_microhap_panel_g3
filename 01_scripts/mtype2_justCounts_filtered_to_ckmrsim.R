# Take output of mtype2 and prepare it for CKMR-sim
# Ben Sutherland (SBIO, VIU)
# 2024-05-23

# Clear space
# rm(list=ls())

# Source amplitools initiator

#### 00. Front matter ####
# Load libraries
#install.packages("tidyverse")
#devtools::install_github("delomast/EFGLmh")
library(tidyverse)
library(EFGLmh)

# Set options
options(tibble.max_extra_cols = 10)

# Set variables
barcode_interp.FN <- "~/Documents/00_sbio/GBMF_UBC_Pacific_oyster/amplicon_panel_mhap/00_archive/GBMF_pilot_sample_ID_barcode_run_from_variantCaller_output_2023-07-27.txt"
sample_interp.FN <- "~/Documents/00_sbio/GBMF_UBC_Pacific_oyster/amplicon_panel_mhap/00_archive/GBMF_pilot_my_data_ind-to-pop_annot_2023-05-01.txt"
mhap.FN <- "14_extract_mhap/filtered_genos.txt"

# Global variables
date <- format(Sys.time(), "%Y-%m-%d")


#### 01. Prepare mtype2 output for EFGLmh input ####
mhap.df <- read.delim(file = mhap.FN, header = T)
dim(mhap.df)
head(mhap.df)


#### 02. Annotate population and update individual names 
barcode_interp.df <- read.delim2(file = barcode_interp.FN)
dim(barcode_interp.df) # 575 lines
head(barcode_interp.df)

# Rename the run name as previously done for individual names
barcode_interp.df[grep(pattern = "OYR", x = barcode_interp.df$Run.Name), "Run.Name"] <- "OYR"
barcode_interp.df[grep(pattern = "Oyster", x = barcode_interp.df$Run.Name), "Run.Name"] <- "Oyster"
head(barcode_interp.df)

# This is what the genotypes are in the format of
head(mhap.df) # OYR_IonCode_0501.bam

barcode_interp.df$matcher <- paste0(barcode_interp.df$Run.Name, "_", barcode_interp.df$Barcode, ".bam")
head(barcode_interp.df)
barcode_interp.df <- barcode_interp.df[,c("Sample.Name", "matcher")]

# Are all the individual names present in the barcode interpretation?
table(mhap.df$Indiv %in% barcode_interp.df$matcher) # should all be TRUE

# Bring in the sample interp that has the indiv name and the population
sample_interp.df <- read.delim2(file = sample_interp.FN)
dim(sample_interp.df)  # 370 lines
head(sample_interp.df)

# Combine the genotypes with the true sample identifiers
mhap.df <- merge(x = mhap.df, y = barcode_interp.df, by.x = "Indiv", by.y = "matcher", all.x = T)
head(mhap.df)
unique(mhap.df$Sample.Name)

mhap.df$Indiv <- mhap.df$Sample.Name
unique(mhap.df$Indiv)
head(mhap.df)

# Combine the samples with the population info
mhap.df <- merge(x = mhap.df, y = sample_interp.df, by.x = "Indiv", by.y = "indiv", all.x = T)
head(mhap.df)

# Keep only the needed cols
mhap.df <- mhap.df[,c("Indiv", "Locus", "Allele1", "Allele2", "Allele1_count", "Allele2_count", "p")]
head(mhap.df)
dim(mhap.df)


#### Prepare wide form ####
# Convert from long form to wide form
mhap.tbl <- as_tibble(mhap.df)
wideFormat <- mtype2wide(x = mhap.tbl)  
dim(wideFormat)
as.data.frame(wideFormat[1:5,1:5])
head(wideFormat$Indiv)

wideFormat.df <- as.data.frame(wideFormat)

# Bring in the population IDs
head(sample_interp.df)

wideformat.df <- merge(x = sample_interp.df, y = wideFormat.df, by.x = "indiv", by.y = "Indiv", all.y = T)
wideformat.df[1:10, 1:10]
wideformat.df$Pop <- wideformat.df$pop
wideformat.df <- wideformat.df[-which(colnames(wideformat.df)=="pop")]
wideformat.df[1:10, 1:10]


#### Update identifiers ####
correct_ids <- TRUE

# Based on an earlier analysis of parentage, several individuals were determined to be 
#  mislabeled, or were from a mixed family. Therefore we will correct these individuals here
#  based on the variable above 'correct_ids=TRUE'
if(correct_ids==TRUE){
  
  print("Correcting IDs...")
  
  # Correct the known erroneous samples as follows
  wideformat.df$indiv <- gsub(pattern = "F10-01", replacement = "F9-101", x = wideformat.df$indiv)
  wideformat.df$indiv <- gsub(pattern = "F14-30", replacement = "F17-130", wideformat.df$indiv)
  wideformat.df$indiv <- gsub(pattern = "F2-11", replacement = "F3-111", wideformat.df$indiv)
  wideformat.df$indiv <- gsub(pattern = "F15-43", replacement = "F14-143", wideformat.df$indiv)
  wideformat.df$indiv <- gsub(pattern = "F15-10", replacement = "F14-110", wideformat.df$indiv)
  wideformat.df$indiv <- gsub(pattern = "F15-06", replacement = "F14-106", wideformat.df$indiv)
  
}else{
  
  print("Not changing IDs. Running analysis with known erroneous sample labels.")
  
}


#### Convert to EFGLmh format ####
mh.dat <- readInData(input = wideformat.df, pedigreeColumn = "Pop"
                    , nameColumn = "indiv") # into EFGLmh format
mh.dat
length(getInds(mh.dat)) # 367 inds

# Number of alleles per locus, before filtering
allele_rich <- aRich(mh.dat)
head(allele_rich)

#pdf(file = "03_results/frequency_of_alleles_per_amplicon_hist.pdf", width = 9, height = 4)
hist(allele_rich$aRich, breaks = 100, main = "", las = 1
     , xlab = "Number alleles per amplicon"
     )

text(x = 25, y = 45, labels = paste0("Total amplicons: ", nrow(allele_rich)))
#dev.off()

table(allele_rich$aRich==1) # number of monomorphs
# hist(table(allele_rich$aRich), breaks = 30, las = 1)

mh.dat

head(getLoci(x = mh.dat))
nloc <- length(getLoci(mh.dat))
nloc


#### Filtering ####
# Remove individuals with more than 30% missing genotypes
geno_success.df <- genoSuccess(x = mh.dat)
geno_success.df <- as.data.frame(geno_success.df)
head(geno_success.df)
hist(geno_success.df$success, breaks = 20, main = "", xlab = "Genotyping Rate (%)", las = 1)
abline(v = 0.7, lty = 3)
numInds(x = mh.dat)

inds_to_rem <- geno_success.df[geno_success.df$success < 0.7, "Ind"]
mh_filt.dat <- removeInds(mh.dat, inds = inds_to_rem)
numInds(x = mh_filt.dat)

# Remove loci with more than 30% missing genotypes
loci_success.df <- lociSuccess(x = mh_filt.dat)
loci_success.df <- as.data.frame(loci_success.df)
head(loci_success.df)

hist(loci_success.df$success, breaks = 20, main = "", xlab = "Genotyping Rate (%)", las = 1)
abline(v = 0.7, lty = 3)

loci_to_rem <- loci_success.df[loci_success.df$success < 0.7, "locus"]
length(loci_to_rem)
mh_filt.dat <- removeLoci(x = mh_filt.dat, lociRemove = loci_to_rem)
mh_filt.dat

# Keep only the target pops
drop_pops <- setdiff(x = unique(mh_filt.dat$genotypes$Pop), y = c("VIU_F1", "VIU_F2"))
mh_filt_parentage.dat <- removePops(x = mh_filt.dat, pops = drop_pops)

# Allelic richness in parentage pop only
allele_rich <- aRich(mh_filt_parentage.dat)
head(allele_rich)

#pdf(file = "03_results/frequency_of_alleles_per_amplicon_hist.pdf", width = 9, height = 4)
hist(allele_rich$aRich, breaks = 100, main = "", las = 1
     , xlab = "Number alleles per amplicon"
)

text(x = 10, y = 45, labels = paste0("Total amplicons: ", nrow(allele_rich)))
#dev.off()

table(allele_rich$aRich==1) # number of monomorphs
allele_rich.df <- as.data.frame(allele_rich)

monomorphs <- allele_rich.df[allele_rich.df$aRich==1, "locus"]
mh_filt_parentage_filt.dat <- removeLoci(x = mh_filt_parentage.dat, lociRemove = monomorphs)


# Convert to rubias and write out
nloc <- length(getLoci(mh_filt_parentage_filt.dat))
output.FN <- paste0("03_results/rubias_mhaps_", nloc, "_", date, ".txt")
parentage_baseline <- exportRubias_baseline(x = mh_filt_parentage_filt.dat, pops = c("VIU_F1", "VIU_F2"), repunit = "Pop", collection = "Pop")
parentage_baseline[1:5,1:10]
write_delim(x = parentage_baseline, file = output.FN, delim = "\t")


#### Parentage ####
# Set user variables
rubias.FN <-output.FN
parent_pop <- "VIU_F1"
offspring_pop <- "VIU_F2"
cutoff <- 5

## Parentage
# Prepare run folder
input_rubias_short.FN <- gsub(pattern = "03_results/", replacement = "", x = rubias.FN)
input_rubias_short.FN <- gsub(pattern = ".txt", replacement = "", x = input_rubias_short.FN)

run_folder.FN <- paste0("03_results/ckmr_input_", input_rubias_short.FN, "_"
                        , offspring_pop, "_vs_", parent_pop,"_", nloc, "_", date
)
print("Making new result folder...")
print(run_folder.FN)
dir.create(run_folder.FN)

# Run ckmr on the input file
ckmr_from_rubias(input.FN = rubias.FN
                 , parent_pop = parent_pop
                 , offspring_pop = offspring_pop
                 , cutoff = 5
                 , output.dir = run_folder.FN
)


# END #
