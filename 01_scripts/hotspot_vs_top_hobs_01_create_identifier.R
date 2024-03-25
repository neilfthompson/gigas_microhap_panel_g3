# Compare hotspot SNP to top HOBS SNP, part 1
# note: requires that amplitargets is cloned and at same level as amplitools
# Ben Sutherland
# initialized 2024-03-18

# Clear space
# rm(list=ls())

# Source amplitools via the initiator

# Load libraries

library(adegenet)
library(tidyr)

# Set user variables
hotspot.FN <- "../amplitargets/cgig/current/WGAG22008_BJS_OYRv01_Hotspot_A.bed" # Cgig_v.1.0 hotspot
regions.FN <- "../amplitargets/cgig/current/WGAG22008_BJS_OYRv01_Region_A.bed"  # Cgig_v.1.0 regions file


#### 01. Prepare identifiers from hotspot files that match empirical data ####
## Load hotspot and regions file
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


## Merge source files to create hotspot identifier corresponding to VCF file identifier
# note: the VCF file identifier is in the form of JH816222_1:110-279_7

## Identifier build requires hotspot and regions files
head(hotspot.df)
head(regions.df)
# These two files are connectable through the "mname" / "GENE_ID"

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

write.table(x = hotspot_data.df, file = "03_results/hotspot_data.txt", sep = "\t", quote = F, col.names = T, row.names = F)

#### Create position file based on hotspots ####
hotspot_data.df$pos.Locus <- paste0(hotspot_data.df$contig.x, ":", hotspot_data.df$amplicon.start, "-", hotspot_data.df$amplicon.end)
#hotspot_data.df[hotspot_data.df$contig.x=="JH816222.1",]
hotspot_data.df$pos.RefPos <- hotspot_data.df$target.pos.end - hotspot_data.df$amplicon.start
hotspot_data.df$pos.Type <- rep("S", times = nrow(hotspot_data.df))

hotspot_data.df <- separate(data = hotspot_data.df, col = "alleles", into = c("REF", "ALT", "ANCHOR"), sep = ";", remove = F)
head(hotspot_data.df)

hotspot_data.df$pos.ValidAlt <- gsub(pattern = "OBS=", replacement = "", x = hotspot_data.df$ALT)

# Now can make a position file with this info
position_file.df <- hotspot_data.df[, c("pos.Locus", "pos.RefPos", "pos.Type", "pos.ValidAlt")]
head(position_file.df)

colnames(position_file.df) <- gsub(pattern = "pos.", replacement = "", x = colnames(position_file.df))
write.table(x = position_file.df, file = "14_extract_mhap/position_file_hotspots.txt"
            , sep = "\t", row.names = F, col.names = T, quote = F
            )


# Next, go to 01_scripts/hotspot_vs_top_hobs_02_.R
