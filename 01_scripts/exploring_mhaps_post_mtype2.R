#### QC steps, Neil F. Thompson original code, adapted as needed by BJGS ####
# Adapted 2024-05-23 (BJGS)
# requires that a datafile produced by mtype2 with the argument --justCount 
# input format: Indiv\t Locus\t Allele\t Count\n

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)

# Load libraries
library(tidyverse)


# User-set variables
input.FN <- "02_input_data/sutherland_viu_read_counts_neil.txt"

# Load data
raw_genos<- read_tsv(file = input.FN)
head(raw_genos)

genos2process<- raw_genos %>%
  mutate(indiv_id = paste0(str_extract(Indiv, "OYR|Oyster"), "_", str_extract(Indiv, "IonCode_[0-9]{1,5}.bam"))) %>%
  dplyr::select(-Indiv) %>%
  arrange(indiv_id, Locus, desc(as.numeric(Count)))

head(genos2process)

# Keep loci that only map to a single location on a chromosome
sam <- read.table("../data/sutherland_snp_panel_bowtie2_map.sam", 
                  header = FALSE, sep = "\t", fill = TRUE)

loci2keep <- 
  sam %>%
  group_by(V1) %>%
  mutate(n_aligns = n()) %>% ungroup() %>%
  filter(n_aligns==1,
         str_detect(V3, "^NC_"))  %>% 
  dplyr::select(V1, V3, V4, V6)
```

tiny bit more data formatting before QA/QC
```{r}
genos2microhaplotopia<- genos2process %>% 
  filter(Locus %in% loci2keep$V1) %>% #keep only single mapping loci
  mutate(count = ifelse(Count=="noReads",0,Count)%>%as.numeric(.),
         haplo = ifelse(Allele=="noReads", NA, Allele)) %>%
  left_join(., loci2keep, by=c("Locus"="V1")) %>%
  mutate(chrom_pos = paste0("NC", str_sub(V3,4,-3),"_",V4)) %>%
  dplyr::select(indiv_id, chrom_pos, haplo, count)
```

#process data with microhaplotpia functions.

filter_raw to get processed genos for CKMRsim
find_contam to identify indivs that have more than 2 alleles passing filter.

input df needs 4 columns:
  indiv.ID, locus, haplo, depth 

```{r, raw-haplo-filter-function, message=FALSE}
filter_raw_microhap_data <- function(hap_raw,
                                     haplotype_depth,
                                     total_depth,
                                     allele_balance,
                                     retain_x_haps = FALSE) {
  
  if (retain_x_haps == FALSE) {
    
    hap_fil1 <- hap_raw %>%
      filter(!str_detect(haplo, "N|X"),
             depth >= haplotype_depth) %>%
      arrange(indiv.ID, locus, desc(depth)) %>%
      group_by(indiv.ID, locus) %>%
      mutate(rank = row_number(),
             allele.balance = depth / depth[1]) %>%
      filter(allele.balance >= allele_balance) %>%
      mutate(gt_type = ifelse(n() > 1, "het", "hom"),
             depth_total = ifelse(gt_type == "het", sum(depth[1], depth[2]), depth)) %>%
      ungroup() %>%
      filter(depth_total >= total_depth)
  } else {
    hap_fil1 <- hap_raw %>%
      filter(!str_detect(haplo, "N"),
             depth >= haplotype_depth) %>%
      arrange(indiv.ID, locus, desc(depth)) %>%
      group_by(indiv.ID, locus) %>%
      mutate(rank = row_number(),
             allele.balance = depth / depth[1]) %>%
      filter(allele.balance >= allele_balance) %>%
      mutate(gt_type = ifelse(n() > 1, "het", "hom"),
             depth_total = ifelse(gt_type == "het", sum(depth[1], depth[2]), depth)) %>%
      ungroup() %>%
      filter(depth_total >= total_depth)
  }
  
  homozygotes2add <- hap_fil1 %>% filter(gt_type == "hom") %>% mutate(rank = 2)
  
  long_genos <- bind_rows(hap_fil1, homozygotes2add) %>%
    group_by(indiv.ID, locus) %>%
    mutate(n_rows = n()) %>%
    ungroup() %>%
    filter(n_rows == 2) %>%
    dplyr::select(-n_rows) %>%
    arrange(indiv.ID, locus)
  
  long_genos
}

find_contaminated_samples <- function(raw_data,
                                      haplotype_depth,
                                      total_depth,
                                      allele_balance) {
  
  hap_fil1 <- raw_data %>%
    filter(!str_detect(haplo, "N|X"),
           depth >= haplotype_depth) %>%
    arrange(indiv.ID, locus, desc(depth)) %>%
    group_by(indiv.ID, locus) %>%
    mutate(rank = row_number(),
           allele.balance = depth / depth[1]) %>%
    filter(allele.balance >= allele_balance) %>%
    mutate(gt_type = ifelse(n() > 1, "het", "hom"),
           depth_total = ifelse(gt_type == "het", sum(depth[1], depth[2]), depth)) %>%
    ungroup() %>%
    filter(depth_total >= total_depth)
  
  homozygotes2add <- hap_fil1 %>% filter(gt_type == "hom") %>% mutate(rank = 2)
  
  long_genos <- bind_rows(hap_fil1, homozygotes2add) %>%
    group_by(indiv.ID, locus) %>%
    mutate(n_rows = n()) %>%
    ungroup() %>%
    filter(n_rows > 2)
  
  long_genos
  
}
```

filter the count data for 
haplotype depth >=5
total depth >=20,
allele.balance >=0.30

```{r, process-raw-count-geno-data}
genos_processed <- filter_raw_microhap_data(
  genos2microhaplotopia%>%rename(indiv.ID = indiv_id, locus = chrom_pos, depth=count),
  5,
  20,
  0.3)
```

remake plots 
total reads per locus
```{r}
genos_processed %>%
  rename(indiv_id=indiv.ID, chrom_pos=locus, count= depth) %>%
  distinct(indiv_id, chrom_pos, .keep_all = TRUE) %>% #get 1 row per locus per indiv
  group_by(chrom_pos) %>%
  summarise(tot_reads = sum(depth_total)) %>% 
  ggplot(., aes(x=fct_reorder(chrom_pos, tot_reads), y=tot_reads)) +
  geom_point() +
  scale_x_discrete(name = "locus", breaks = NULL,expand=c(0, 5)) +
  scale_y_continuous("total reads")+
  theme(axis.text.x = element_blank()) +
  theme_light()+
  ggtitle("reads per locus filtered data-VIU")

#ggsave("../data-microtyper/plots/total-reads-per-locus-filtered.pdf", width = 7, height=4)
```


processed reads per indiv
```{r}
genos_processed %>%
  rename( indiv_id=indiv.ID, chrom_pos=locus ,count= depth) %>%
  distinct(indiv_id, chrom_pos, .keep_all = TRUE) %>% #get 1 row per locus per indiv
  group_by(indiv_id) %>%
  summarise(tot_reads = sum(depth_total),
            mean_reads = mean(depth_total),
            sd_reads = sd(depth_total)) %>%
  ungroup() %>% 
  ggplot(., aes(x = fct_reorder(indiv_id, tot_reads), y=tot_reads)) +
  geom_point() + scale_x_discrete(name = "individuals", breaks = NULL,expand=c(0, 5)) +
  scale_y_continuous(name= "total reads") +
  theme(axis.text.x=element_blank()) +
  theme_light() +
  ggtitle("total reads per indiv filtered VIU")

#ggsave("../data-microtyper/plots/total-reads-per-indiv-filtered.pdf", width=7, height=4)
```

whats the percentage of reads eaten up by whales and hi rollers?
  ```{r}
genos_processed %>%
  rename( indiv_id=indiv.ID, chrom_pos=locus ,count= depth) %>%
  distinct(indiv_id, chrom_pos, .keep_all = TRUE) %>% #get 1 row per locus per indiv
  group_by(indiv_id) %>%
  summarise(tot_reads = sum(depth_total),
            mean_reads = mean(depth_total),
            sd_reads = sd(depth_total)) %>%
  ungroup() %>%
  mutate(indiv_type = case_when(
    tot_reads >=1e6 ~"whale",
    tot_reads>=500000&tot_reads<1e6 ~"hi_roller",
    tot_reads<500000 ~"normal"
  )) %>%
  group_by(indiv_type) %>%
  mutate(tot_read_grp = sum(tot_reads),
         n_indiv = n_distinct(indiv_id)) %>%
  distinct(indiv_type, .keep_all=TRUE) %>%
  dplyr::select(indiv_type, tot_read_grp, n_indiv) %>%
  ungroup() %>%
  mutate(pct_reads = tot_read_grp/sum(tot_read_grp))
```
whales have >1 million reads
hi rollers are 500,000 to 1,000,000 reads
normal are < 500,000 reads


how many loci scored per indiv?
  ```{r}
genos_processed %>%
  distinct(indiv.ID, locus) %>%
  count(indiv.ID) %>%
  ggplot(., aes(x=fct_reorder(indiv.ID, n), y = n)) + 
  geom_point() +
  scale_x_discrete(name = "individuals", breaks = NULL, expand=c(0, 5)) +
  scale_y_continuous(name= "# scored loci") +
  theme(axis.text.x=element_blank()) +
  theme_light()

#ggsave("../data-microtyper/plots/loci-per-indiv-filtered.pdf", width=7, height=4)
```


number of haplos per locus
```{r}
genos_processed  %>%
  group_by(locus) %>%
  summarise(n_haps = n_distinct(haplo)) %>%
  ggplot(., aes(x = fct_reorder(locus, n_haps), y = n_haps)) +
  geom_point() +
  scale_x_discrete(name = "locus") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("number of haplotypes per locus filtered data")

#ggsave("../data-microtyper/plots/number-processed-allele-per-locus.pdf")
```

what is the percent heterozygosity per sample?
  
  ```{r, heterozygosity-per-indiv}
het_tbl <- genos_processed %>% group_by(indiv.ID,locus) %>% slice(1) %>%
  ungroup() %>%
  count(indiv.ID, gt_type) %>%
  group_by(indiv.ID) %>%
  mutate(pct_heterozygous = case_when(
    gt_type=="het" ~ round(n/sum(n),3),
    gt_type=="hom" ~ NA),
    tot_loci = sum(n)) %>%
  filter(!is.na(pct_heterozygous))

library(knitr)
library(kableExtra)
#het_tbl %>%
# kable("html") %>%
#kable_styling("striped")

ggplot(het_tbl, aes(x=tot_loci, y = pct_heterozygous)) +
  geom_point()
```

what is percentage of filtered reads vs total reads per sample?
  ```{r}
indiv_raw_counts <- genos2microhaplotopia %>% group_by(indiv_id) %>%
  summarise(total_raw_reads = sum(count)) %>% ungroup() 

genos_processed %>%
  rename(indiv_id=indiv.ID, chrom_pos=locus ,count= depth) %>%
  distinct(indiv_id, chrom_pos, .keep_all = TRUE) %>% #get 1 row per locus per indiv
  group_by(indiv_id) %>%
  summarise(tot_reads = sum(depth_total),
            mean_reads = mean(depth_total),
            sd_reads = sd(depth_total),
            n_loci = n_distinct(chrom_pos)) %>%
  ungroup()  %>%
  left_join(., indiv_raw_counts, "indiv_id") %>%
  mutate(pct_ontarget = tot_reads/total_raw_reads) %>%
  ggplot(., aes(x=fct_reorder(indiv_id, pct_ontarget), y = pct_ontarget)) +
  geom_point() +
  scale_x_discrete(name = "individuals", breaks = NULL, expand=c(0, 5)) +
  scale_y_continuous(name= "percentage of reads in final genotype", breaks =seq(0,1,0.05)) +
  theme(axis.text.x=element_blank()) +
  theme_light()

#ggsave("../data-microtyper/plots/percentage-reads-final-genotype-indiv.pdf", width=7, height=4)
```


how does percent on target compare to number of loci scored?
  ```{r}
genos_processed %>%
  rename( indiv_id=indiv.ID, chrom_pos=locus ,count= depth) %>%
  distinct(indiv_id, chrom_pos, .keep_all = TRUE) %>% #get 1 row per locus per indiv
  group_by(indiv_id) %>%
  summarise(tot_reads = sum(depth_total),
            mean_reads = mean(depth_total),
            sd_reads = sd(depth_total),
            n_loci = n_distinct(chrom_pos)) %>%
  ungroup()  %>%
  left_join(., indiv_raw_counts, "indiv_id") %>%
  mutate(pct_ontarget = tot_reads/total_raw_reads) %>%
  ggplot(., aes(x=pct_ontarget, y = n_loci)) +
  geom_point() +
  #scale_x_discrete(name = "individuals", breaks = NULL, expand=c(0, 5)) +
  scale_x_continuous(name= "percentage on target reads", breaks =seq(0,1,0.05)) +
  #   theme(axis.text.x=element_blank()) +
  theme_light()
```

lastly what about contaminated samples? These are indivs who have more than 2
alleles pass the microhaplotopia QC filters.

```{r,contaminated-samples-df}
contam_samps <- find_contaminated_samples(
  genos2microhaplotopia%>%rename(indiv.ID = indiv_id, locus = chrom_pos, depth=count), 
  5,20,0.3) %>%
  distinct(indiv.ID, locus) %>%
  count(indiv.ID)  

#write_rds(contam_samps, "../data-microtyper/microtyper-contaminated-samples-viu.rds")
```


#write RDS of processed genotype data for further analysis
```{r}
write_rds(genos_processed, "../data-microtyper/microtyper-processed-haplotypes-viu.rds", compress="xz")
#genos_processed <- readRDS("../data-microtyper/microtyper-processed-haplotypes-USDA.rds")
```


Calc haplotype allele freqs and number of alleles per locus
```{r}
af_mhap <- genos_processed %>%
  count(locus, haplo) %>%
  group_by(locus) %>%
  mutate(freq = n / sum(n),
         n_allele = n(),
         chrom = str_extract(locus, "NC[0-9]{6}"),
         pos = str_extract(locus, "_[0-9]{1,8}")%>%str_replace(., "_","")) %>%
  ungroup() %>%
  select(chrom, pos, locus, haplo, freq, n_allele)

#grab the highest MAF mhap from each locus and plot it.

af_mhap %>%
  arrange(locus, desc(freq)) %>%
  group_by(locus) %>%
  slice(1) %>% ungroup(.) %>%
  #filter(chrom=="NC047559") %>%
  mutate(pos=as.numeric(pos)) %>%
  ggplot(., aes(x=pos, y=freq)) +
  geom_point(aes(size=n_allele))+
  #scale_y_continuous(breaks=seq(0,1,0.1)) +
  theme(axis.text.x = element_text(angle = 70, hjust=1)) +
  facet_grid(rows = vars(chrom), scales = "free_x")

```

do the allele freqs by population
```{r}

genos_processed %>%
  mutate(pop = "VIU") %>%
  filter(!str_detect(indiv.ID, "Sample")) %>%
  count(pop,locus, haplo) %>%
  group_by(pop,locus) %>%
  mutate(freq = n / sum(n),
         n_allele = n(),
         chrom = str_extract(locus, "NC[0-9]{6}"),
         pos = str_extract(locus, "_[0-9]{1,8}")%>%str_replace(., "_","")) %>%
  ungroup() %>%
  select(chrom, pos, pop,locus, haplo, freq, n_allele) %>%
  arrange(locus, desc(freq)) %>%
  group_by(locus) %>%
  slice(1) %>% ungroup(.) %>% 
  #filter(pop=="midori") %>%
  mutate(pos=as.numeric(pos)) %>%
  ggplot(., aes(x=pos, y=freq)) +
  geom_point(aes(size=n_allele, colour=n_allele))+
  #scale_y_continuous(breaks=seq(0,1,0.1)) +
  theme(axis.text.x = element_text(angle = 70, hjust=1)) +
  scale_colour_gradientn(colors=terrain.colors(16))+
  facet_grid(rows = vars(chrom), cols = vars(pop), scales = "free_y")
```

flip it, do n haps as y axis, allele freq as size --> dont like this. i removed
the code to keep this cleaner.

make histo of allele richness by pop
####TO DO HERE#####
```{r}
genos_processed %>%
  mutate(pop = "VIU") %>%
  filter(!str_detect(indiv.ID, "Sample"))%>%
  count(pop,locus, haplo) %>%
  group_by(pop,locus) %>%
  mutate(freq = n / sum(n),
         n_allele = n(),
         chrom = str_extract(locus, "NC[0-9]{6}"),
         pos = str_extract(locus, "_[0-9]{1,8}")%>%str_replace(., "_","")) %>%
  ungroup() %>%
  select(chrom, pos, pop,locus, haplo, freq, n_allele) %>% 
  arrange(locus, desc(freq)) %>%
  group_by(pop,locus) %>%
  slice(1) %>% ungroup(.) %>% 
  group_by(pop) %>%
  count(n_allele) %>%
  mutate(tot_loci = sum(n),
         freq_AR = n/tot_loci %>%round(., 1)) %>%
  ggplot(., aes(x=n_allele, y = freq_AR)) +
  geom_col(width = 0.75) +
  scale_x_continuous(name="Number alleles per locus", breaks=seq(1,17,1)) +
  scale_y_continuous(name = "Proportion of all loci") +
  facet_grid(rows=vars(pop))

#ggsave(filename="../figure3_manuscript.pdf", height = 4.5, width = 6.5)
```

I bet this is skewed by the crappy sample runs....

####TO DO HERE#####
```{r}
af_ar_plots <- function(x) {
  
  chrom = pull(x)
  genos_processed %>%
    mutate(pop="VIU") %>%
    filter(!str_detect(indiv.ID, "Sample"),
           str_detect(locus, chrom)) %>%
    count(pop,locus, haplo) %>%
    group_by(pop,locus) %>%
    mutate(freq = n / sum(n),
           n_allele = n(),
           chrom = str_extract(locus, "NC[0-9]{6}"),
           pos = str_extract(locus, "_[0-9]{1,8}")%>%str_replace(., "_","")%>%as.numeric(.)) %>%
    ungroup() %>%
    arrange(chrom,pos, pop) %>%
    select(chrom, pos, pop,locus, haplo, freq, n_allele) %>%
    ggplot(., aes(x=factor(pos), y=freq)) +
    geom_point(aes( colour=factor(n_allele)))+
    scale_colour_viridis_d(option="D", name="Number of microhaplotype alleles")+
    scale_y_continuous(breaks=seq(0,1,0.1), name = "Microhaplotype allele frequency") +
    theme(axis.text.x = element_text(angle = 70, hjust=1),
          legend.position = "top")+
    xlab("Position (base pair coordinate)")+
    facet_grid(rows = vars(pop), cols = vars(chrom))
  
  ####ggsave(filename = paste0("../figure_supplemental_VIU",x,".pdf"), width = 9, height = 6.5)
  
}

chroms2feed <- genos_processed %>% distinct(locus) %>%
  filter(!is.na(locus))%>%
  mutate(chrom=str_extract(locus, "^NC[0-9]{6}")) %>%
  dplyr::select(chrom) %>% distinct(chrom) %>%
  group_split(chrom)

#map(chroms2feed, af_ar_plots)
```

```{r}
genos_processed %>%
  mutate(pop="VIU") %>%
  count(pop,locus, haplo) %>%
  group_by(pop,locus) %>%
  mutate(freq = n / sum(n),
         n_allele = n(),
         chrom = str_extract(locus, "NC[0-9]{6}"),
         pos = str_extract(locus, "_[0-9]{1,8}")%>%str_replace(., "_","")) %>%
  ungroup() %>%
  select(chrom, pos, pop,locus, haplo, freq, n_allele) %>%
  ggplot(., aes(x=factor(pos), y=freq)) +
  geom_point(aes( colour=factor(n_allele)))+
  scale_colour_viridis_d(option="D", name="Number of microhaplotype alleles")+
  scale_y_continuous(breaks=seq(0,1,0.1), name = "Microhaplotype allele frequency") +
  theme(axis.text.x = element_text(angle = 70, hjust=1),
        legend.position = "top")+
  xlab("Position (base pair coordinate)") + 
  facet_grid(rows = vars(pop), cols = vars(chrom), scales = "free_x")
```
#ok, thats wayyyy to insane to plot 3 pops and all 10 chroms simultaneously,but
I think this idea has some merit. loci within a few hundred kb are too close to
get a good sense of (points overlap eachother), but I think this gets a 
reasonable picture across to the reader.

#le fin.