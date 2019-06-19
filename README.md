# sysSVM-Annotation

This repo is dedicated to annotation of molecular properties in TCGA and multiple OAC datasets as input to sysSVM

## Notes on the latest results
TCGA mutation files: `/mnt/lustre/users/k1801782/dataset/mc3_bailey.bed` or `/home/camp/let/working/Trang/data/mc3.v0.2.8.PUBLIC.code.filtered.maf` <br />
TCGA_ASCAT segment files : `/mnt/lustre/users/k1801782/dataset/CNV/ASCATforTCGA/segments` <br />
TCGA_ASCAT ploidy file : `/home/camp/let/working/Trang/data/ascat_tcga_ploidy.txt`<br />
TCGA expression files: `/home/camp/let/working/Trang/data/TCGA-RNASeq-download16042019/`. In which `TCGA-RNASeq` is the raw RNASeq file downloaded from TCGA portal, `TCGA-RNASeq-CODE` is the RNASeq for each cancer type (per .rds file)

Script to plot expression of damaged genes for all cancer types: `/home/camp/let/working/Trang/scripts/TCGA_AmpHomoLoss.R`

Annovar, OncodriveClust and CNV script (for the above result, also in Github): `/mnt/lustre/users/k1801782/dataset/Annotation_individual_cancers`<br />
Annovar, OncodriveClust and CNV results for each cancer type: `/mnt/lustre/users/k1801782/dataset/Annotation_individual_cancers`<br />
Inside each folder (for each cancer type) you will find the following main files:<br />
`ANNOVAR/muts_ann_damaging.rds` : annotated mutations, and damaging classification according to function, splicing and conservation <br />
`OncodriveClust/muts_ann_onco_damaging.rds`: annotated mutations (including damaging classification) and OncodriveClust results for individual cancer type <br />
`OncodriveClust/muts_ann_onco_damaging_pancan.rds`:  annotated mutations (including damaging classification) and OncodriveClust results for all cancer types <br />
`CNV/CNVGainLoss.rds`: results of the intersection between ASCAT segments and hg19 gene coordinates, with `CNV_type_corrected` <br />
`total_table.rds`: table counting number of truncating, non truncating damaging mutation, gain of function, CNVGain and Loss for each gene in each patient, based on Annovar and OncodriveClust (individual cancer) annotation <br />
`total_table_pancan.rds`: table counting number of truncating, non truncating damaging mutation, gain of function, CNVGain and Loss for each gene in each patient, based on Annovar and OncodriveClust (pan cancer) annotation <br />

Combination of files and filtering (for 34 cancer types):
`Allcancertypes_individual_annotation_damaging.rds`: total_table of 34 cancer types with only damaging mutations (truncating, non truncating damaging mutation, gain of function, CNVGain and Loss), based on ANNOVAR, OncodriveClust (individual cancer) and CNV annotation <br />
`Allcancertypes_pancan_annotation_damaging.rds`: total_table of 34 cancer types with only damaging mutations (truncating, non truncating damaging mutation, gain of function, CNVGain and Loss), based on ANNOVAR, OncodriveClust (pan cancer) and CNV annotation  <br />
`Allcancertypes_pancan_annotation_damaging_falseHomotoHeteroDel.rds`: total_table (annovar, oncodriveclust and CNV) with damaging mutations after removing false positive Homozygous deletion (genes with CN = 0 but expression (RNA_seq) > 1/cancertype_purity), and move those with a damaging mutation to Heterozygous deletion (change CN ==1) (This is the final file to use, but if you want to filter for different expression threshold, then use the above 2 files)
 
Other files that might be important/interesting:
`/mnt/lustre/users/k1801782/dataset/TCGA_2015_totaltable`: Total table for 31 TCGA cancer types annotated by Thanos in 2015 (based on individual files in `/mnt/lustre/users/k1469280/mourikisa/data/TCGA/01_03_2015`) <br />
`/mnt/lustre/users/k1801782/dataset/TCGA_2015_totaltable/total_table_pancan_damaging.rds`: Total table with only damaging mutations for 31 TCGA cancer types annotated by Thanos in 2015 <br />
`/mnt/lustre/users/k1801782/dataset/TCGA_2015_totaltable/total_table_pancan_damaging_withpurity.rds`: Total table with only damaging mutations for 31 TCGA cancer types annotated by Thanos in 2015, removing samples without purity values

Softwares: `/mnt/lustre/users/k1801782/softwares` or `/home/camp/let/working/Trang/softwares`
