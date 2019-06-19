###################################### Data 2019
# Mutation: bailey et al, 2019
somatic_mutations = CODE %>% group_by(CODE) %>% summarise(n_mutations=n_distinct(Sample))

# CNV ASCAT results, all samples and QC
somatic_cnv0 = ascat_tcga_ploidy %>% 
  group_by(cancer_type) %>% summarise(n_cnv_notQC=n_distinct(sample)) %>% rename(CODE=cancer_type)
somatic_cnv = ascat_tcga_ploidy[which(ascat_tcga_ploidy$pass==TRUE & ascat_tcga_ploidy$rep==TRUE),] %>% 
  group_by(cancer_type) %>% summarise(n_cnv=n_distinct(sample)) %>% rename(CODE=cancer_type)

# RNASeq in FPKM, downloaded from portal.gdc...
sample_RNA_rosalind$Tumor_code = substring(sample_RNA_rosalind$submitter_id,14,15)
sample_RNA_rosalind$type = ifelse(sample_RNA_rosalind$Tumor_code<10,'Tumor','Normal')
somatic_rnaseq = sample_RNA_rosalind[which(sample_RNA_rosalind$type=='Tumor'),] %>% 
  group_by(CODE) %>% summarise(n_rnaseq=n_distinct(sample))
somatic_rnaseq = somatic_rnaseq[which(is.na(somatic_rnaseq$CODE)==FALSE),]

# All samples that have mutation, cnv, rnaseq data
somatic_mutscnvrnaseq = ascat_tcga_ploidy %>% filter(sample %in% intersect(CODE$Sample,sample_RNA_rosalind$sample))
somatic_mutscnvrnaseq = unique(somatic_mutscnvrnaseq %>% select(cancer_type,sample)) %>% rename(CODE=cancer_type)
somatic_mutscnvrnaseq = somatic_mutscnvrnaseq %>% group_by(CODE) %>% summarise(n_muts_cnv_rna=n_distinct(sample))

somatic_mutscnvrnaseqQC = ascat_tcga_ploidy %>% 
  filter(sample %in% intersect(CODE$Sample,sample_RNA_rosalind$sample)) %>%
  filter(pass==TRUE & rep==TRUE)
somatic_mutscnvrnaseqQC = unique(somatic_mutscnvrnaseqQC %>% select(cancer_type,sample)) %>% rename(CODE=cancer_type)
somatic_mutscnvrnaseqQC = somatic_mutscnvrnaseqQC %>% group_by(CODE) %>% summarise(n_muts_cnv_rna_QC=n_distinct(sample))

###################################### All samples from 2015
# purity
if (FALSE){
library(TCGAbiolinks)
library(dplyr)
purity = Tumor.purity
purity[purity[] == "NaN"] = NA
purity = data.frame(sapply(purity, function(x) gsub(",", ".", x)) )
for (i in 3:7){
  purity[,i] = as.numeric(as.character(purity[,i]))
}
purity = mutate(purity, consensus = NA)

for (i in 1:nrow(purity)){
  if (!is.na(purity$CPE[i])){
    purity$consensus[i] = purity$CPE[i]
    purity$consensus_mean[i] = purity$CPE[i]
  } else {
    purity$consensus[i] = median(c(purity$ESTIMATE[i], purity$ABSOLUTE[i], purity$LUMP[i], purity$IHC[i]), na.rm = TRUE)
    purity$consensus_mean[i] = mean(c(purity$ESTIMATE[i], purity$ABSOLUTE[i], purity$LUMP[i], purity$IHC[i]), na.rm = TRUE)
  }
}
purity$sample = substring(purity$Sample.ID,1,15)
purity$Cancer.type=as.character(purity$Cancer.type)
purity$Cancer.type[which(substring(purity$sample,1,12) %in% oac_ids)]="OAC"
purity$Cancer.type[which(purity$Cancer.type =='ESCA')] ='OSCC'
}
# summarise the 7828 samples annotated in 2015
sample_MutsCNVRNA_2015$CODE = as.character(sample_MutsCNVRNA_2015$CODE)
sample_MutsCNVRNA_2015$CODE[which(substring(sample_MutsCNVRNA_2015$sample,1,12) %in% oac_ids)]="OAC"
sample_MutsCNVRNA_2015$CODE[which(sample_MutsCNVRNA_2015$CODE =='ESCA')] ='OSCC'
somatic_total_2015 = sample_MutsCNVRNA_2015 %>% group_by(CODE) %>% summarise(n_mutcnvrna_2015=n_distinct(sample))
somatic_total_2015_purity = sample_MutsCNVRNA_2015 %>% filter(sample %in% purity$sample) %>% group_by(CODE) %>% summarise(n_mutcnvrna_2015_purity=n_distinct(sample))

#damaging_alterations_2015  =readRDS('~/Rosalind/TCGA_2015_totaltable/total_table_pancan_damaging.rds')
somatic_total_2015_damaging = damaging_alterations_2015 %>% group_by(CODE) %>% summarise(n_mutcnvrna_2015_dam=n_distinct(sample))
somatic_total_2015_damaging_purity = damaging_alterations_2015  %>% filter(sample %in% substring(purity$sample,1,12)) %>% group_by(CODE) %>% summarise(n_mutcnvrna_2015_dam_purity=n_distinct(sample))

#Allcancertypes_pancan_annotation_damaging <- readRDS("~/Trang-mnt/data/Trang_rosalind/Annotation_individual_cancers/Allcancertypes_pancan_annotation_damaging.rds")
somatic_damaging_muts_cnv = Allcancertypes_pancan_annotation_damaging %>% group_by(CODE) %>% summarise(n_damaged_muts_cnv=n_distinct(sample))
#Allcancertypes_pancan_annotation_damaging_withexp <- readRDS("~/Trang-mnt/data/Trang_rosalind/Annotation_individual_cancers/Allcancertypes_pancan_annotation_damaging_withexp.rds")
somatic_damaging_muts_cnv_exp = Allcancertypes_pancan_annotation_damaging_withexp %>% group_by(CODE) %>% summarise(n_damaged_muts_cnv_rna=n_distinct(sample))




somatic = merge(merge(somatic_mutations,merge(somatic_cnv0,somatic_cnv)),somatic_rnaseq)
somatic = merge(somatic, merge(somatic_mutscnvrnaseq,somatic_mutscnvrnaseqQC))
somatic = merge(merge(somatic,somatic_damaging_muts_cnv),somatic_damaging_muts_cnv_exp)
somatic = left_join(somatic,merge(somatic_total_2015, somatic_total_2015_damaging), all.x=TRUE)
somatic = left_join(somatic,merge(somatic_total_2015_purity, somatic_total_2015_damaging_purity), all.x=TRUE)

somatic[is.na(somatic)] <- 0
somatic$diff = somatic$n_damaged_muts_cnv_rna-somatic$n_mutcnvrna_2015_purity
somatic[nrow(somatic)+1,] = c('SUM', colSums(somatic[,-1]))
rm(list = ls(pattern = "somatic_")) 
write.table(somatic, '~/Trang-mnt/data/summary_samples_2019_2015.txt', row.names = F, quote = F, sep='\t')
