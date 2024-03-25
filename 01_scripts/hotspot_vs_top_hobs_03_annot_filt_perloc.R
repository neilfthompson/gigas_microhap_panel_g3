# Compare hotspot SNP to top HOBS SNP, part 3
# note: requires that part 2 was run and '03_results/comparison_loci.txt' is present in amplitools
# Ben Sutherland
# initialized 2024-03-18

# Clear space
# rm(list=ls())


# Source amplitools via the initiator

# Load libraries
library(vcfR)

# Set user variables
vcf_filtered.FN            <- "14_extract_mhap/mpileup_calls_biallele_snp_q20_dp10_under_0.1_missing_w_AF_0.01.vcf" # filtered VCF

# Load comparison loci (i.e., hotspot loci identified by de novo SNP identification)
comparison_loci.df <- read.delim2(file = "03_results/comparison_loci.txt", header = F, sep = "\t")
colnames(comparison_loci.df) <- "comparison.loci"
head(comparison_loci.df)
length(comparison_loci.df$comparison.loci)

## Load filtered VCF file
my_vcf <- read.vcfR(file = vcf_filtered.FN)

# Convert vcf to genind
my_data.genind <- vcfR2genind(x = my_vcf)


# Source simple_pop_stats, don't clear space #
#### Add population attribute, rename samples ####
# Rename data
obj <-  my_data.genind
obj

# Create blank popmap
generate_popmap(df = obj)

# Obtain filler information for the popmap
barcode_interp.FN  <- "../00_archive/GBMF_pilot_sample_ID_barcode_run_from_variantCaller_output_2023-07-27.txt"
sample_interp.FN   <- "../00_archive/GBMF_pilot_my_data_ind-to-pop_annot_2023-05-01.txt"

barcode_interp.df <- read.delim2(file = barcode_interp.FN)
head(barcode_interp.df) # sample names, barcodes, and run ID

sample_interp.df <- read.delim2(file = sample_interp.FN)
head(sample_interp.df) # sample names, population

# What individuals are present? 
ind.present <- indNames(obj)
head(ind.present) # format: runID_barcode.fastq
#ind.present <- gsub(pattern = ".*IonCode", replacement = "IonCode", x = ind.present)
#ind.present <- gsub(pattern = "\\.fastq", replacement = "", x = ind.present)
head(ind.present)

# Convert to df, retain order
popmap_filler.df <- as.data.frame(ind.present)
popmap_filler.df$order <- seq(1:nrow(popmap_filler.df))
head(popmap_filler.df)
dim(popmap_filler.df) # 212

# Merge with barcode interp to get sample names
head(popmap_filler.df)
head(barcode_interp.df)
# Create matching identifier on barcode interp file
barcode_interp.df$matcher <- paste0(barcode_interp.df$Run.Name, "_", barcode_interp.df$Barcode, ".fastq")
head(barcode_interp.df)
popmap_filler.df <- merge(popmap_filler.df, barcode_interp.df, all.x = T, by.x = "ind.present", by.y = "matcher"
                                   , sort = F)
head(popmap_filler.df)
dim(popmap_filler.df) # 212

# Merge with sample interp to get population attribute
head(popmap_filler.df)
head(sample_interp.df)
popmap_filler.df   <-  merge(x = popmap_filler.df, y = sample_interp.df, by.x = "Sample.Name", by.y = "indiv"
                             , all.x = T, sort = F)
head(popmap_filler.df)
dim(popmap_filler.df) # 212

# Create individual identifier in amplitools format
popmap_filler.df$full.alt.ID <- paste0(popmap_filler.df$Run.Name, "__", popmap_filler.df$Barcode, "__", popmap_filler.df$Sample.Name)
head(popmap_filler.df)

# Put back in original order
popmap_filler.df <- popmap_filler.df[order(popmap_filler.df$order), ]
write_delim(x = popmap_filler.df, file = "00_archive/popmap_filler.txt", delim = "\t")

#stop("Before proceeding, use the above text file to manually annotate the pop map")

# Read in the annotated popmap
annotate_from_popmap(df = obj, popmap.FN = "00_archive/my_data_ind-to-pop_annot.txt"
                     , convert_to_alt_ID = TRUE
)

# indNames(obj_annot)
# pop(obj_annot)

# Drop blank samples (pop = NA)
keep <- indNames(obj_annot)[!is.na(pop(obj_annot))]
obj_annot <- obj_annot[keep]
table(pop(obj_annot))


# Check per sample genotyping rate
percent_missing_by_ind(df = obj_annot)
head(missing_data.df)
summary(missing_data.df$ind.per.missing)
sd(missing_data.df$ind.per.missing)

#### 01. Filter ####
table(missing_data.df$ind.per.missing > 0.3) # no samples are missing more than 30% loci

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

summary(obj.df$marker.per.missing)
sd(obj.df$marker.per.missing)

# Plot geno rate per marker
pdf(file = "03_results/missing_data_per_locus.pdf", width = 7.5, height = 4.5)
plot(100 * (1- obj.df$marker.per.missing), xlab = "Locus (index)", ylab = "Genotyping rate (%)", las = 1
     , ylim = c(0,100)
     #, pch = plot_pch
     #, cex = plot_cex
)
abline(h = 70
       #, col = "grey60"
       , lty = 3)
dev.off()

# No markers to be removed


# Filter by MAF
maf_filt(data = obj_annot, maf = 0.01)
obj <- obj_maf_filt
obj



#### 02. Pick highest HOBS SNP per amplicon window ####
## Per locus statistics
per_locus_stats(data = obj)
head(per_loc_stats.df)
nrow(per_loc_stats.df)
per_loc_hobs.df <- per_loc_stats.df[,c("mname", "Hobs")] # keep only the necessary cols
head(per_loc_hobs.df)

# Split identifier mname to obtain name of amplicon
per_loc_hobs.df <- separate(data = per_loc_hobs.df, col = "mname", into = c("contig", "region_and_pos"), sep = ":", remove = F)
head(per_loc_hobs.df)

per_loc_hobs.df <- separate(data = per_loc_hobs.df, col = "mname", into = c("one", "two", "abs.pos"), sep = "_", remove = F)
head(per_loc_hobs.df)
per_loc_hobs.df$amplicon_name <- paste0(per_loc_hobs.df$one, "_", per_loc_hobs.df$two)
head(per_loc_hobs.df)

per_loc_hobs.df <- per_loc_hobs.df[,c("mname", "contig", "amplicon_name", "abs.pos", "Hobs")]
head(per_loc_hobs.df)
nrow(per_loc_hobs.df) 

# Observe the distribution of SNPs per amplicon
pdf(file = "03_results/SNP_per_amplicon.pdf", width = 7.5, height = 4.5)
hist(table(per_loc_hobs.df$amplicon_name), breaks = 30, xlab = "SNPs per amplicon", main = "", las = 1)
dev.off()

# remember, comparison_loci holds hotspot loci that were observed in the empirical dataset
head(comparison_loci.df)
colnames(comparison_loci.df) <- "hotspot_locus" 
head(comparison_loci.df)

# Get the Hobs for these loci
comparison_loci.df <- merge(x = comparison_loci.df, y = per_loc_hobs.df, by.x = "hotspot_locus", by.y = "mname"
                       , sort = F
)

### TODO: add an all.x flag to this? 

head(comparison_loci.df)

# Add empty vector to determine max hobs for each amplicon (any identified, filtered SNP)
comparison_loci.df$max_hobs <- NA
comparison_loci.df$max_hobs_id <- NA

head(comparison_loci.df)

target_amplicon <- NULL ; slice <- NULL
for(i in 1:nrow(comparison_loci.df)){
  
  target_amplicon <-  comparison_loci.df$amplicon_name[i]
  
  slice <- per_loc_hobs.df[per_loc_hobs.df$amplicon_name==target_amplicon,]
  
  # Order the slice
  slice <- slice[with(slice, order(slice$Hobs, decreasing = T)), ]
  comparison_loci.df$max_hobs[i] <- slice[1,"Hobs"]
  
  comparison_loci.df$max_hobs_id[i] <- slice[1, "mname"]
  
}

head(comparison_loci.df)

# For how many of the amplicons is the hotspot also the max HOBS?
table(comparison_loci.df$hotspot_locus==comparison_loci.df$max_hobs_id)
nrow(comparison_loci.df)

pdf(file = "03_results/max_hobs_vs_hotspot_hobs_per_amplicon.pdf", width = 7.5, height = 4.5)
plot(x = comparison_loci.df$Hobs, y = comparison_loci.df$max_hobs
     , las = 1
     , ylab = "amplicon max HOBS"
     , xlab = "amplicon hotspot HOBS"
)
dev.off()

save.image(file = "03_results/prepared_filtered.RData")

### Next, go to hotspot_vs_top_hobs_04



### OTHER TODO
# - multimappers number snp per locus
# - may need to indicate the size of the amplicons by average

