####
#install.packages("tidyverse")
devtools::install_github("delomast/EFGLmh")

library(tidyverse)
library(EFGLmh)
options(tibble.max_extra_cols = 10)

# Get data into EFGLmh
mhap.FN <- "14_extract_mhap/genos.txt"
mhap.df <- read.delim(file = mhap.FN, header = T)
dim(mhap.df)
head(mhap.df)
mhap.df <- separate(data = mhap.df, col = "Indiv", into = c("dir", "ind"), sep = "/", remove = T)
head(mhap.df)
mhap.df$ind <- gsub(pattern = ".sorted.bam", replacement = "", x = mhap.df$ind)
head(mhap.df)

mhap.df <- mhap.df[, c("ind", "Locus", "Allele1", "Allele2")]
head(mhap.df)

inds <- unique(mhap.df$ind)
head(inds)
loci <- unique(mhap.df$Locus)


results.mat <- matrix(data = NA, nrow = length(inds), ncol = (2+2*length(loci)))
results.mat[1:5,1:5]
colnames(results.mat)  <- c("population", "indiv", sort(c(paste0(loci, ".A1"), paste0(loci, ".A2"))))
results.mat[1:5,1:5]
results.df <- as.data.frame(results.mat)
results.df[1:5,1:6]

results.df$indiv <- inds
results.df[1:10,1:5]


# Fill it in
ioi <- NULL; loi <- NULL
for(i in 1:length(inds)){
  
  ioi <- inds[i]
  
  # per locus
  for(j in 1:length(loci)){
    
    loi <- loci[j]
    
    allele1 <- mhap.df[mhap.df$ind==ioi & mhap.df$Locus==loi, "Allele1"]
    allele2 <- mhap.df[mhap.df$ind==ioi & mhap.df$Locus==loi, "Allele2"]
    
    # Populate table
    results.df[which(results.df$indiv==ioi), which(colnames(results.df)==paste0(loi,".A1"))] <- allele1
    results.df[which(results.df$indiv==ioi), which(colnames(results.df)==paste0(loi,".A2"))] <- allele2
    
    # results.df[which(results.df$indiv==ioi), which(colnames(results.mat)==paste0(loi,".A1"))] <- mhap.df[which(mhap.df$ind==ioi & mhap.df$Locus ==loi), "Allele1"]
    # results.df[which(results.df$indiv==ioi), which(colnames(results.mat)==paste0(loi,".A2"))] <- mhap.df[which(mhap.df$ind==ioi & mhap.df$Locus ==loi), "Allele2"]   
    # 
  }
}

results.df[1:5,1:5]



### NOTE: can make this faster with tidyverse e.g., reshape2 (dcast)

# Then just need to merge it with a population file, and possibly rename samples
annot.FN <- "~/Documents/00_sbio/GBMF_UBC_Pacific_oyster/amplicon_panel_mhap/simple_pop_stats/00_archive/my_data_ind-to-pop_annot.txt"
annot.df <- read.delim2(file = annot.FN, header = T, sep = "\t")
head(annot.df)
annot.df$indiv <- gsub(pattern = ".fastq", replacement = "", x = annot.df$indiv)

annot_results.df <- merge(x = annot.df, y = results.df, by = "indiv", all.y = T)
ncol(annot_results.df)
annot_results.df[1:5,1:10]

annot_results.df <- separate(data = annot_results.df, col = "alt.ID", into = c("run", "barcode", "ID"), sep = "__", remove = T)
annot_results.df[1:5,1:10]

test <- annot_results.df[,c(-1,-3, -4, -6, -7, -8), ]
test[1:5,1:10]
colnames(test)[which(colnames(test)=="pop")] <- "population"
colnames(test)[which(colnames(test)=="ID")] <- "indiv"
test[1:5,1:10]

# Quick fix
test <- test[!duplicated(x = test$indiv), ]


# And finally...
mh.dat <- readInData(input = test) # into EFGLmh format
# , then convert to rubias
parentage_baseline <- exportRubias_baseline(x = mh.dat, pops = c("VIU_F1", "VIU_F2"), repunit = "Pop", collection = "Pop")

# , then amplitools ckmr-sim wrapper
# May come across issues with the naming of .A1 and .A2, but this should be easy to fix, e.g., to .01 and .02


