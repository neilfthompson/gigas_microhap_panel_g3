# Select best replicate from numbers of reads, uses input metadata
# B. Sutherland, 2024-03-21

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Source amplitools initiator

# Read in read counts file
read_counts.FN <- "12_input_mhap/reads_per_sample_table.txt"
read_counts.df <- read.delim2(file = read_counts.FN, sep = " ", header = F)
colnames(read_counts.df) <- c("fastq.name", "num.records")
head(read_counts.df)
# read_counts.df$barcode <- gsub(pattern = ".*_IonCode", replacement = "IonCode", x =read_counts.df$fastq.name)
# read_counts.df$barcode <- gsub(pattern = ".fastq", replacement = "", x = read_counts.df$barcode)
# read_counts.df$run <- gsub(pattern = "_IonCode.*", replacement = "", x =read_counts.df$fastq.name)
# head(read_counts.df)


# Read in metadata
barcode_interp.FN <- "~/Documents/00_sbio/GBMF_UBC_Pacific_oyster/amplicon_panel_mhap/00_archive/GBMF_pilot_sample_ID_barcode_run_from_variantCaller_output_2023-07-27.txt"
sample_interp.FN <- "~/Documents/00_sbio/GBMF_UBC_Pacific_oyster/amplicon_panel_mhap/00_archive/GBMF_pilot_my_data_ind-to-pop_annot_2023-05-01.txt"

barcode_interp.df <- read.delim2(file = barcode_interp.FN)
dim(barcode_interp.df) # 575 lines
barcode_interp.df$matcher <- paste0(barcode_interp.df$Run.Name, "_", barcode_interp.df$Barcode, ".fastq")
head(barcode_interp.df)

sample_interp.df <- read.delim2(file = sample_interp.FN)
dim(sample_interp.df)  # 370 lines
head(sample_interp.df)


### Merge together
head(read_counts.df)
head(barcode_interp.df)

full.df <- merge(x = read_counts.df, y = barcode_interp.df, by.x = "fastq.name", by.y = "matcher"
                    , all.x = T)
head(full.df)
head(sample_interp.df)
dim(full.df)

full.df <- merge(x = full.df, y = sample_interp.df, by.x = "Sample.Name", by.y = "indiv", all.x = T)
dim(full.df)

head(full.df)

# Sort based on number records to keep only one sample per indiv
full.df <- full.df[with(full.df, order(full.df$Sample.Name, full.df$num.records, decreasing = T)), ]
head(full.df)

full.df <- full.df[!duplicated(x = full.df$Sample.Name), ]
head(full.df)
dim(full.df)

table(full.df$pop)


# Optional: limit to specific populations
pop_limited.df <- full.df[grep(pattern = "VIU|PEN", x = full.df$pop), ]
head(pop_limited.df)
dim(pop_limited.df)

write_delim(x = pop_limited.df, file = "00_archive/selected_pops_best_ind.txt")

















