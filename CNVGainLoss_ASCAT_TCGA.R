#!/usr/bin/Rscript

## **************************************
#  This script maps CN segments with gene coordinates and classifies them amplified and deleted CNVs 
#  
## **************************************

# use prepare_run_file_parsing_totalTable.R to create submit_CNmap.sh
# bash submit_CNmap.sh
library(plyr)
library(dplyr)
options(stringsAsFactors = F)

code= commandArgs(trailingOnly = TRUE)[1]
ascat_dir = commandArgs(trailingOnly = TRUE)[2]
save_dir = commandArgs(trailingOnly = TRUE)[3]
cnv_annotated = commandArgs(trailingOnly = TRUE)[4]

dir.create(file.path(save_dir, "CNV"), showWarnings = TRUE)
save_dir=paste0(save_dir,"/CNV/")

annotate_genes = function(x, samplename, cnv_path,gene_coord_fn="/mnt/lustre/users/k1623452/Software/NCG/NCG6/ncg6_blat_exon_coordinates_hg19.RData",save_dir=NULL){
  
  if(is.null(gene_coord_fn) | is.null(save_dir)){
    stop("annotate_genes: parameters missing")
  }
  
  load(gene_coord_fn)
  coord_file <- paste0(save_dir, samplename,".tmp.blat.genes.bed")
  
  write.table(unique(coord[,c("chrom","start","end","symbol","Entrez")]),file=coord_file,row.names=F,col.names = F,quote = F,sep="\t")
  write.table(x[,c("Chromosome","Start","End", "sample")], file=paste0(save_dir, samplename,".bed"), row.names=F,col.names = F,quote = F,sep="\t")
  
  #system("intersectBed -a tmp.blat.genes.bed -b tmp.bed -f 0.25 -wa -wb > tmp.genes.bed", wait = TRUE)
  system(paste0("intersectBed -a ", coord_file, " -b ", paste0(save_dir, samplename,".bed"), " -f 0.25 -wa -wb > ", paste0(save_dir, samplename,".genes.bed")), wait = TRUE)
  tmp = read.table(paste0(save_dir, samplename,".genes.bed"), h=F)
  #tmp = read.table(cnv_path,h=F)
  tmp$key=paste0(tmp$V6,".",tmp$V7,'.',tmp$V8,'.',tmp$V9)
  key=paste0(x$Chromosome,".",x$Start,'.',x$End,'.',x$sample)
  y = cbind(tmp[,1:5], x[match(tmp$key, key),])
  colnames(y)[1:5] = c("chrom","start","end","symbol_19549","entrez_19549")
  y
}

path=paste0(ascat_dir,'/segments/')

# Load ASCAT summary and QC file
ascat_tcga_ploidy <- read.delim(paste0(ascat_dir,"/summary.ascatTCGA.penalty70.txt"))
PassSample = unique(ascat_tcga_ploidy$name[which((ascat_tcga_ploidy$rep==TRUE)&(ascat_tcga_ploidy$pass==TRUE))])
ascat_tcga_ploidy = ascat_tcga_ploidy[which(ascat_tcga_ploidy$name%in%PassSample),]

if (TRUE){
  tmp=ascat_tcga_ploidy[which(ascat_tcga_ploidy$cancer_type==code),]$name
  xx=NULL
  for (i in tmp){
    ascat_tcga_segments = read.delim(paste0(path,i,".segments.txt"))
    x = ascat_tcga_segments 
    colnames(x)[1]="name"
    x$Total_CN = x$nMajor + x$nMinor
    x = left_join(x, ascat_tcga_ploidy %>% select(name,patient, purity, ploidy,cancer_type, pass,rep), by="name")
    x = x %>% mutate(CNV_type_corrected=ifelse(Total_CN >= 2 * ploidy, "Gain", ifelse(Total_CN < 2, "Loss", NA)))
    x = x %>% rename(Chromosome = chr, Start = startpos, End = endpos, sample = name)
    x$Chromosome = paste0("chr",x$Chromosome)
    xx = rbind(xx,x)
  }
  xx = annotate_genes(xx, code,cnv_path,save_dir =save_dir)
  message("Done intersecting CNV with genes coord")
  
  saveRDS(xx, file = paste0(save_dir,"CNVGainLoss.rds"))
}
