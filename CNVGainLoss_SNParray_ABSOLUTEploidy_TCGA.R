#!/usr/bin/Rscript

## **************************************
#  Start with CN data mapped with gene coordinates, filtered for CN matching >25% gene length
#  This was done in local terminal with command `intersectBed -a <hg19.genes.coord.bed> -b <CN.coord.bed> -f 0.25 -wa -wb > <tmp.genes.bed>`
#  This script calculates CNVGain and CNVLoss in 3 ways (hard threshold 0.3, CNVGain>=2*ploidy, CNVGain>ploidy)
#  Change `if (FALSE)` to `if (TRUE)` to run
## **************************************


setwd("/mnt/lustre/users/k1801782/dataset/CNV/scripts")

## Command line arguments
filepath = commandArgs(trailingOnly = TRUE)[1]

print(filepath)

filepath = sub("CN_intersectBed","CN_intersectBed_filtered", filepath)
filepath = sub(".bed",".txt", filepath)
cnv_ann = read.table(filepath)
cnv_ann$Copy_number = round(cnv_ann$CN)#cnv_ann$CN#round(cnv_ann$CN)
cnv_ann$Sample = cnv_ann$sampleID
cnv_ann$sampleID =  substring(cnv_ann$sampleID,1,19)

if (FALSE){
# WAY 1
cnv_ann1 = subset(cnv_ann, abs(segment_mean)>0.3 & segment_mean<1.5)
cnv_ann1$CNV_type <- ifelse(cnv_ann1$segment_mean>0.3, "Gain","Loss")

cat(paste("segment abs(0.3):", "Gains:", nrow(subset(cnv_ann1, CNV_type=="Gain")), sep = " "), "\n")
cat(paste("segment abs(0.3):", "Loss:", nrow(subset(cnv_ann1, CNV_type=="Loss")), sep = " "), "\n")
fileout = sub(".txt","_gainloss_segmentmean.rds", filepath)
saveRDS(cnv_ann1,fileout)

print(fileout)
cat('Done annotating gain/loss based on segment mean',"\n")
}


# WAY 2
if (FALSE){
ploidy = read.delim('/mnt/lustre/users/k1801782/dataset/CNV/TCGA_mastercalls.broad.ploidy.txt',header=T)
ploidy$sampleID = substring(ploidy$sample,1,19)
cat(paste("number of samples with both CN and ploidy:", length(intersect(unique(cnv_ann$sampleID),ploidy$sampleID)), sep=" "),"\n")
cnv_ann = merge(cnv_ann, ploidy, by='sampleID')
cnv_ann$CNV_type <- ifelse(cnv_ann$Copy_number>=2*cnv_ann$ploidy,"Gain", ifelse(cnv_ann$Copy_number<2,"Loss", NA))

fileout = sub(".txt","_gainloss_segmentmean.rds", filepath)
fileout = sub("_gainloss_segmentmean.rds", "_gainloss_ploidy_CNrounded_19charID.rds", fileout)
saveRDS(cnv_ann,fileout)
cat('Done annotating gain/loss based on 2*ploidy',"\n")
}


# WAY 3
if (FALSE){
ploidy = read.delim('/mnt/lustre/users/k1801782/dataset/CNV/TCGA_mastercalls.broad.ploidy.txt',header=T)
ploidy$sampleID = substring(ploidy$sample,1,19)
cat(paste("number of samples with both CN and ploidy:", length(intersect(unique(cnv_ann$sampleID),ploidy$sampleID)), sep=" "),"\n")
cnv_ann = merge(cnv_ann, ploidy, by='sampleID')
#cnv_ann$CNV_type <- ifelse(cnv_ann$Copy_number>=2*cnv_ann$ploidy,"Gain", ifelse(cnv_ann$Copy_number<2,"Loss", NA))
cnv_ann$CNV_type <- ifelse(cnv_ann$Copy_number>cnv_ann$ploidy,"Gain", ifelse(cnv_ann$Copy_number<2,"Loss",NA))

fileout = sub(".txt","_gainloss_segmentmean.rds", filepath)
fileout = sub("_gainloss_segmentmean.rds", "_gainloss_ploidy_CNrounded_19charID.rds", fileout)
saveRDS(cnv_ann,fileout)
cat('Done annotating gain/loss based on ploidy',"\n")
}
cat(paste("number of samples with both CN and ploidy:", length(unique(cnv_ann$sampleID)), sep=" "),"\n")
cat(paste("ploidy:", "Gains:", nrow(subset(cnv_ann, CNV_type=="Gain")), sep = " "), "\n")
cat(paste("ploidy:", "Loss:", nrow(subset(cnv_ann, CNV_type=="Loss")), sep = " "), "\n")


