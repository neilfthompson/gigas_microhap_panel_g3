# Exploring mhaps; requires that amplitargets is cloned and at same level as amplitools
# Ben Sutherland
# 2024-03-18

# Clear space
# rm(list=ls())

# Set working directory
setwd("~/Documents/00_sbio/GBMF_UBC_Pacific_oyster/amplicon_panel_mhap/amplitools")

# Load libraries
library(vcfR)
library(adegenet)
library(tidyr)

# Set variables
vcf.FN <- "14_extract_mhap/mpileup_calls_SNP_only_biallelic_q20_dp10_Fmiss_0.1_w_AF_maf0.01.vcf"
hotspot.FN <- "../amplitargets/cgig/current/WGAG22008_BJS_OYRv01_Hotspot_A.bed"
regions.FN <- "../amplitargets/cgig/current/WGAG22008_BJS_OYRv01_Region_A.bed"


#### 01. Load amplicon panel files ####
# Load hotspot file
hotspot.df <- read.delim2(file = hotspot.FN, skip = 1, header = F)
head(hotspot.df)
colnames(hotspot.df) <- c("contig", "target.pos.start", "target.pos.end", "mname", "extra", "strand", "alleles", "type")
head(hotspot.df)

# Load regions file
regions.df <- read.delim2(file = regions.FN, skip = 1, header = F)
head(regions.df)
colnames(regions.df) <- c("contig", "amplicon.start", "amplicon.end", "SP.val", "extra", "assorted")
head(regions.df)


#### 02. Create hotspot identifier corresponding to VCF file identifier ####
# note: VCF file identifier: e.g., JH816222_1:110-279_7

## Identifier build requires hotspot and regions files
head(hotspot.df)
head(regions.df)

# Prepare regions file fields for merging
regions.df <- separate(data = regions.df, col = "assorted", into = c("GENE_ID", "Pool"), sep = ";", remove = T)
regions.df <- separate(data = regions.df, col = "GENE_ID", into = c("GENE_ID_prefix", "GENE_ID"), sep = "=", remove = T)
head(regions.df) # GENE_ID allows you to match to the hotspot file

# Prepare hotspot file fields for merging
hotspot.df$mname <- as.character(hotspot.df$mname)
str(hotspot.df)

# Confirm all prepared
nrow(hotspot.df) # 592
length(intersect(x = regions.df$GENE_ID, y = hotspot.df$mname)) # 590
setdiff(x = regions.df$GENE_ID, y = hotspot.df$mname) # something odd here, but only impacts two loci
setdiff(y = regions.df$GENE_ID, x = hotspot.df$mname) # something odd here, but only impacts two loci

# Merge
hotspot_data.df <- merge(x = hotspot.df, y = regions.df, by.x = "mname", by.y = "GENE_ID")
head(hotspot_data.df)
str(hotspot_data.df)
 

## Create identifier for hotspots to match VCF file naming (see above)
hotspot_data.df$expect.name <- paste0(gsub(pattern = "\\.", replacement = "_", x = hotspot_data.df$contig.x)
                                      , ":", hotspot_data.df$amplicon.start, "-", hotspot_data.df$amplicon.end
                                        , "_", hotspot_data.df$target.pos.end - hotspot_data.df$amplicon.start)
head(hotspot_data.df)


#### 03. Load empirical data ####
my_vcf <- read.vcfR(file = vcf.FN)

# Convert vcf to genind
my_data.genind <- vcfR2genind(x = my_vcf)
head(locNames(my_data.genind))
head(indNames(my_data.genind))

# How many of the hotspot file sites were observed in the empirical data? 
empirical_sites.vec <- locNames(my_data.genind)
length(intersect(x = hotspot_data.df$expect.name, y = empirical_sites.vec)) # 339 only in both
comparison_loci <- intersect(x = hotspot_data.df$expect.name, y = empirical_sites.vec)
head(comparison_loci)
# note: since the empirical approach is not identifying all of the hotspots, use the default approach
#  with the empirical chip to conduct standard hotspot workflow


#### 04. Annotate data, remove duplicates ####
# source simple_pop_stats without clearing space

# Rename data
obj <-  my_data.genind
obj

# Annotate samples
generate_popmap(df = obj)
# annotate it manually

# Create answer key
barcode_interp.FN  <- "../00_archive/GBMF_pilot_sample_ID_barcode_run_from_variantCaller_output_2023-07-27.txt"
sample_interp.FN   <- "../00_archive/GBMF_pilot_my_data_ind-to-pop_annot_2023-05-01.txt"

barcode_interp.df <- read.delim2(file = barcode_interp.FN)
head(barcode_interp.df)

sample_interp.df <- read.delim2(file = sample_interp.FN)
head(sample_interp.df)

# Only keep the relevant Run.Name (#TODO: THIS WILL CHANGE DEPENDING ON LOADED SAMPLES)
barcode_interp.df <- barcode_interp.df[grep(pattern = "315-Oyster", x = barcode_interp.df$Run.Name), ]

# Create popmap filler (build in amplitools format (run__barcode__indiv) due to chip and replicate info needed)
ind.present <- indNames(obj)
ind.present <- gsub(pattern = ".*IonCode", replacement = "IonCode", x = ind.present)
ind.present <- gsub(pattern = "\\.fastq", replacement = "", x = ind.present)
head(ind.present)

# Convert to df, retain order
popmap_filler.df <- as.data.frame(ind.present)
popmap_filler.df$order <- seq(1:nrow(popmap_filler.df))
head(popmap_filler.df)
dim(popmap_filler.df)

popmap_filler.df <-          merge(popmap_filler.df, barcode_interp.df, all.x = T, by.x = "ind.present", by.y = "Barcode"
                                     , sort = F)
head(popmap_filler.df)
dim(popmap_filler.df) # 192

popmap_filler.df   <-  merge(x = popmap_filler.df, y = sample_interp.df, by.x = "Sample.Name", by.y = "indiv"
                              , all.x = T, sort = F)
head(popmap_filler.df)
dim(popmap_filler.df) # 192

# Create individual identifier in amplitools format
popmap_filler.df$full.alt.ID <- paste0(popmap_filler.df$Run.Name, "__", popmap_filler.df$ind.present, "__", popmap_filler.df$Sample.Name)
head(popmap_filler.df)

# Reorder by popmap order
popmap_filler.df <- popmap_filler.df[order(popmap_filler.df$order), ]
head(popmap_filler.df)
write_delim(x = popmap_filler.df, file = "00_archive/popmap_filler.txt", delim = "\t")

head(x = cbind(ind.present, popmap_filler.df), n = 20)
tail(x = cbind(ind.present, popmap_filler.df), n = 5)

## Use the above to manually annotate the pop map

# Read in the annotated popmap
annotate_from_popmap(df = obj, popmap.FN = "00_archive/my_data_ind-to-pop_annot.txt"
                     , convert_to_alt_ID = TRUE
                     )

indNames(obj_annot)
pop(obj_annot)

# Drop blank samples (pop = NA)
keep <- indNames(obj_annot)[!is.na(pop(obj_annot))]
obj_annot <- obj_annot[keep]
table(pop(obj_annot))


## Get rid of replicate with lower genotyping rate
percent_missing_by_ind(df = obj_annot)
head(missing_data.df)

# Obtain the sample ID for sorting
missing_data.df <- separate(data = missing_data.df, col = "ind", into = c("run", "barcode", "indiv"), sep = "__", remove = F)
head(missing_data.df)

missing_data.df <- missing_data.df[with(missing_data.df, order(missing_data.df$indiv, missing_data.df$ind.per.missing, decreasing = F)), ]
head(missing_data.df)

# Remove the duplicated individual with more percent missing
missing_data.df <- missing_data.df[!duplicated(missing_data.df$indiv), ]
head(missing_data.df)
keep <- missing_data.df$ind # identifier
obj_annot <- obj_annot[keep]
obj_annot
table(pop(obj_annot))

# What is the missing data like? 
summary(missing_data.df$ind.per.missing)


#### 05. Filter ####
missing_data.df$ind.per.missing > 0.3 # no samples are missing more than 30% loci

# Filter loci based on missing data
obj.df <- genind2df(obj_annot)
obj.df[1:5,1:5]
obj.df <- t(obj.df)
obj.df[1:5,1:5]
obj.df <- obj.df[2:nrow(obj.df),] # remove pop row
obj.df[1:5,1:5]
dim(obj.df)
str(obj.df)

obj.df <- as.data.frame(obj.df)
dim(obj.df)
str(obj.df)
obj.df[1:5,1:5] # See top left of file
obj.df[(dim(obj.df)[1]-5):dim(obj.df)[1], (dim(obj.df)[2]-5):dim(obj.df)[2]] # See bottom right of file

# Add collector col
obj.df$marker.per.missing <- NA

for(i in 1:(nrow(obj.df))){
  
  # Per marker                      sum all NAs for the marker, divide by total number markers
  obj.df$marker.per.missing[i] <-  (sum(is.na(obj.df[i,]))-1) / (ncol(obj.df)-1) 
  
}

# Missing by marker
plot(100 * (1- obj.df$marker.per.missing), xlab = "Marker", ylab = "Genotyping rate (%)", las = 1
     , ylim = c(0,100)
     #, pch = plot_pch
     #, cex = plot_cex
)
abline(h = 70
       #, col = "grey60"
       , lty = 3)

# No markers to be removed

obj <- obj_annot

maf_filt(data = obj, maf = 0.01)
obj <- obj_maf_filt


#### 04. Pick highest MAF SNP per amplicon window ####
## Per locus statistics
per_locus_stats(data = obj)
head(per_loc_stats.df)
nrow(per_loc_stats.df)

head(per_loc_stats.df)
per_loc_hobs.df <- per_loc_stats.df[,c("mname", "Hobs")]
head(per_loc_hobs.df)

per_loc_hobs.df <- separate(data = per_loc_hobs.df, col = "mname", into = c("contig", "region_and_pos"), sep = ":", remove = F)
head(per_loc_hobs.df)

per_loc_hobs.df <- separate(data = per_loc_hobs.df, col = "mname", into = c("one", "two", "abs.pos"), sep = "_", remove = F)
head(per_loc_hobs.df)
per_loc_hobs.df$amplicon_name <- paste0(per_loc_hobs.df$one, "_", per_loc_hobs.df$two)
head(per_loc_hobs.df)

per_loc_hobs.df <- per_loc_hobs.df[,c("mname", "contig", "amplicon_name", "abs.pos", "Hobs")]
head(per_loc_hobs.df)
nrow(per_loc_hobs.df) 

# sort(table(per_loc_hobs.df$amplicon_name))
hist(table(per_loc_hobs.df$amplicon_name), breaks = 30, xlab = "SNPs per amplicon", main = "")

# remember, comparison_loci holds hotspot loci that were observed in the empirical dataset
head(comparison_loci)

comparison.df <- as.data.frame(comparison_loci)
head(comparison.df)
colnames(comparison.df) <- "hotspot_locus" 
head(comparison.df)

# Get the Hobs for these loci
comparison.df <- merge(x = comparison.df, y = per_loc_hobs.df, by.x = "hotspot_locus", by.y = "mname"
                    , sort = F
                    )
head(comparison.df)

comparison.df$max_hobs <- NA
comparison.df$max_hobs_id <- NA

head(comparison.df)

target_amplicon <- NULL ; slice <- NULL
for(i in 1:nrow(comparison.df)){
  
  target_amplicon <-  comparison.df$amplicon_name[i]
  
  slice <- per_loc_hobs.df[per_loc_hobs.df$amplicon_name==target_amplicon,]
  
  # Order the slice
  slice <- slice[with(slice, order(slice$Hobs, decreasing = T)), ]
  comparison.df$max_hobs[i] <- slice[1,"Hobs"]
  
  comparison.df$max_hobs_id[i] <- slice[1, "mname"]
  
}

head(comparison.df)

plot(x = comparison.df$Hobs, y = comparison.df$max_hobs
     , las = 1
     , ylab = "amplicon max HOBS"
     , xlab = "amplicon hotspot HOBS"
     )

save.image(file = "prepared_filtered.RData")


### OTHER TODO
# - multimappers number snp per locus
# - will need to indicate the size of the amplicons by average

#### NEXT: 
# - per amplicon, find the best HOBS, create obj and ckmr-sim with this
## BEST SNP

#### NEXT: 
# - create obj with amplicon max hobs and one with amplicon hotspot hobs
# - create rubias object with the above, convert to CKMR-sim, run analysis with either
# - try to solve the issue of so few of the hotspots coming through (reduce filters?)

#### FINALLY NEXT: 
# genotype microhaps and then compare against the best SNP


#### IF GO BACK
# should we use a combination of both chips, and sort by the largest seq library size per
# ioncode, then use that as the representative sample? 

