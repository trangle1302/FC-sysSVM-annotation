###################################################################
#### This is a file with random functions that accummulates 
###################################################################


MC3_bailey = MC3_bailey_filtered
intersect(intersect(grep("LOW", MC3_bailey$IMPACT),
                    grep("benign", MC3_bailey$PolyPhen)),
          grep("tolerated|tolerated_low_confidence", MC3_bailey$SIFT)) -> id
MC3_bailey_filtered = MC3_bailey

MC3_bailey = MC3_bailey_filtered
library(dplyr)
summary_cohort = MC3_bailey %>%
  group_by(CODE) %>%
  summarise(num_samples = n_distinct(Tumor_Sample_Barcode),
            num_alterations = n())
summary_cohort$CODE = factor(summary_cohort$CODE,levels = summary_cohort$CODE[order(summary_cohort$num_samples)])
library(dplyr)
library(ggplot2)
library(ggrepel)
library(forcats)
library(scales)

df=melt(summary_cohort)
p<-ggplot(data=df, aes(x=CODE, y=value)) +
  geom_col()+
  facet_grid( variable ~ ., scales = "free")+
  #  ggtitle("Number of samples in each cancer type") +
  ylab("Number of samples and alterations") + xlab("Cancer types") +
  theme(axis.text.x=element_text(angle=30,hjust=0.75,vjust=0.75, size = 12, face = "bold"))
p 
ggsave("/Users/let/Desktop/a.png",width = 12, height = 9)

df = data.frame(sort(summary(MC3_bailey_filtered$Variant_Classification))) 
df$MutClass = rownames(df)
colnames(df)[1] = "value"
df %>%
  mutate(prop = round(df$value*100 / sum(df$value), 2)) -> df
df$position = cumsum(df$prop) - df$prop/2
pie <- ggplot(df, aes(x ="", y = prop, fill =  fct_reorder(MutClass,prop, desc))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") +
  geom_label_repel(aes(y=position, label = value), size=5, show.legend = F, nudge_x = 1) +
  guides(fill = guide_legend(title = "Mutation Class"))
pie

library(plotrix)
pie3D(df$prop, labels = df$MutClass, main = "Mutation Class", 
      explode=0.1, radius=.9, labelcex = 1.2,  start=0.7)
CGC_types = Census_allMon.Oct.29.14_00_10.2018$Tumour.Types.Somatic.
uniqueCGCtypes = unique(unlist(strsplit(CGC_types, ", ", fixed = TRUE, perl = FALSE, useBytes = FALSE)))

CGC_genelist = list()
CGC_genelist$ACC = CGC$Gene.Symbol[grep('adrenocortical', CGC$Tumour.Types.Somatic.)]
CGC_genelist$BLCA = CGC$Gene.Symbol[setdiff(grep('bladder|urothelial', CGC$Tumour.Types.Somatic.),
                                            grep('gall', CGC$Tumour.Types.Somatic.))]
CGC_genelist$BRCA = CGC$Gene.Symbol[grep('breast', CGC$Tumour.Types.Somatic.)] #breast_cancer and triple_negative_breast_cancer
CGC_genelist$CESC = CGC$Gene.Symbol[grep('cervical', CGC$Tumour.Types.Somatic.)]
CGC_genelist$CHOL = CGC$Gene.Symbol[grep('cholangiocarcinoma|biliary tract', CGC$Tumour.Types.Somatic.)]
CGC_genelist$COAD = CGC$Gene.Symbol[grep('colon|colorectal', CGC$Tumour.Types.Somatic.)] #colorectal = colon + rectum => combine or not
CGC_genelist$DLBC = CGC$Gene.Symbol[grep('DLBCL', CGC$Tumour.Types.Somatic.)] #or anything with B cell carcinoma ("B-")?
CGC_genelist$ESCA = CGC$Gene.Symbol[grep('oesophageal', CGC$Tumour.Types.Somatic.)]
CGC_genelist$GBM = CGC$Gene.Symbol[grep('glioblastoma|GBM', CGC$Tumour.Types.Somatic.)]
CGC_genelist$HNSC = CGC$Gene.Symbol[grep('HNSCC|head|neck|oral squamous', CGC$Tumour.Types.Somatic.)]
CGC_genelist$KICH = CGC$Gene.Symbol[grep('kidney', CGC$Tumour.Types.Somatic.)] # no chromophobe mensioned 
CGC_genelist$KIRC = CGC$Gene.Symbol[grep('kidney|clear cell renal|RCC', CGC$Tumour.Types.Somatic.)] #CGC only have kidney cancer
CGC_genelist$KIRP = CGC$Gene.Symbol[grep('kidney|papillary renal', CGC$Tumour.Types.Somatic.)] # 90% of kidney cancers are clear cell
CGC_genelist$LAML = CGC$Gene.Symbol[grep('AML', CGC$Tumour.Types.Somatic.)] 
CGC_genelist$LGG = CGC$Gene.Symbol[setdiff(grep('glioma', CGC$Tumour.Types.Somatic.),
                                           grep('paraganglioma', CGC$Tumour.Types.Somatic.))]
CGC_genelist$LIHC = CGC$Gene.Symbol[grep('hepatocellular', CGC$Tumour.Types.Somatic.)]
CGC_genelist$LUAD = CGC$Gene.Symbol[grep('lung adenocarcinoma', CGC$Tumour.Types.Somatic.)] #lung, lung cancer , lung carcinoma, NSCLC can be both
CGC_genelist$LUSC = CGC$Gene.Symbol[grep('lung SCC', CGC$Tumour.Types.Somatic.)]
CGC_genelist$MESO = CGC$Gene.Symbol[grep('mesothelioma', CGC$Tumour.Types.Somatic.)]
CGC_genelist$OV = CGC$Gene.Symbol[setdiff(grep('ovarian', CGC$Tumour.Types.Somatic.),
                                          grep('mixed germ cell tumour', CGC$Tumour.Types.Somatic.))]
CGC_genelist$PAAD = CGC$Gene.Symbol[setdiff(grep('pancrea', CGC$Tumour.Types.Somatic.),
                                            grep("neuroendocrine", CGC$Tumour.Types.Somatic.))]
CGC_genelist$PCPG = CGC$Gene.Symbol[grep('pheochromocytoma|paraganglioma', CGC$Tumour.Types.Somatic.)]
CGC_genelist$PRAD = CGC$Gene.Symbol[grep('prostate', CGC$Tumour.Types.Somatic.)]
CGC_genelist$READ = CGC$Gene.Symbol[grep('colorectal', CGC$Tumour.Types.Somatic.)]
CGC_genelist$SARC = CGC$Gene.Symbol[setdiff(grep('sarcoma', CGC$Tumour.Types.Somatic.),
                                            grep('uterine|endometrial', CGC$Tumour.Types.Somatic.))] 
CGC_genelist$SKCM = CGC$Gene.Symbol[setdiff(grep('melanoma', CGC$Tumour.Types.Somatic.),
                                            grep('mucosal|soft|uveal', CGC$Tumour.Types.Somatic.))]
CGC_genelist$STAD = CGC$Gene.Symbol[grep('stomach|gastric', CGC$Tumour.Types.Somatic.)]
CGC_genelist$TGCT = CGC$Gene.Symbol[grep('testicular', CGC$Tumour.Types.Somatic.)]
CGC_genelist$THCA = CGC$Gene.Symbol[setdiff(grep('thyroid', CGC$Tumour.Types.Somatic.),
                                            grep('parathyroid', CGC$Tumour.Types.Somatic.))]
CGC_genelist$THYM = CGC$Gene.Symbol[grep('thymoma', CGC$Tumour.Types.Somatic.)] #there is no thymic or thymoma
CGC_genelist$UCEC = CGC$Gene.Symbol[setdiff(grep('uterine serous|endometri', CGC$Tumour.Types.Somatic.),
                                            grep('stromal', CGC$Tumour.Types.Somatic.))]
CGC_genelist$UCS = CGC$Gene.Symbol[grep('uterine carcinosarcoma', CGC$Tumour.Types.Somatic.)]
CGC_genelist$UVM = CGC$Gene.Symbol[grep('uveal', CGC$Tumour.Types.Somatic.)]

library(reshape2)
CGC_genes = melt(CGC_genelist)
CGC_genes_CODEsum = as.data.frame(table(CGC_genes$L1))
p<-ggplot(data=CGC_genes_CODEsum, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1))
p

# Over lapping the number of samples and alterations
trainsetCGC = data.frame(matrix(ncol=115))
colnames(trainsetCGC) = colnames(MC3_bailey_filtered)
for (i in 1:length(CGC_genelist)){
  trainsetCGC = rbind(trainsetCGC,
                      MC3_bailey_filtered[intersect(which(MC3_bailey_filtered$CODE == names(CGC_genelist)[i]),
                                                    which(MC3_bailey_filtered$Hugo_Symbol %in% unlist(CGC_genelist[i]))),])
}
trainsetCGC = trainsetCGC[-1,]
df = trainsetCGC %>%
  group_by(CODE) %>%
  summarise(num_samples = n_distinct(Tumor_Sample_Barcode),
            num_alterations = n())

##############################################################
#### Search in NCG6
##############################################################
NCGgene = NCG6_systemslevelproperties$symbol[grep("cgc|vog", NCG6_systemslevelproperties$cancer_type)]
#NCG6_cancergenes = NCG6_cancergenes.tsv
#NCG6_cancergenes.tsv = NCG6_cancergenes[which(NCG6_cancergenes$symbol %in% NCGgene),]
NCG_genelist = list()
#common = NCG6_cancergenes.tsv$symbol[grep('multiple', NCG6_cancergenes.tsv$primary_site)]
NCG_genelist$ACC = NCG6_cancergenes.tsv$symbol[grep('adrenocortical', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$BLCA = NCG6_cancergenes.tsv$symbol[grep('bladder_cancer', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$BRCA = NCG6_cancergenes.tsv$symbol[grep('breast', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$CESC = NCG6_cancergenes.tsv$symbol[grep('cervical', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$CHOL = NCG6_cancergenes.tsv$symbol[grep('cholangio|biliary_tract_cancer', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$COAD = NCG6_cancergenes.tsv$symbol[grep('colorectal', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$DLBC = NCG6_cancergenes.tsv$symbol[grep('B-cell|cutaneous_DLBCL', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$ESCA = NCG6_cancergenes.tsv$symbol[grep('esophageal', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$GBM = NCG6_cancergenes.tsv$symbol[grep('glioblastoma', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$HNSC = NCG6_cancergenes.tsv$symbol[grep('squamous_head_and_neck', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$KICH = NCG6_cancergenes.tsv$symbol[grep('chromophobe', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$KIRC = NCG6_cancergenes.tsv$symbol[grep('clear_cell_renal_cancer', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$KIRP = NCG6_cancergenes.tsv$symbol[grep('papillary_renal_cell_carcinoma', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$LAML = NCG6_cancergenes.tsv$symbol[grep('acute_myeloid_leukemia', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$LGG = NCG6_cancergenes.tsv$symbol[grep('low_grade_glioma', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$LIHC = NCG6_cancergenes.tsv$symbol[grep('hepatocellular_carcinoma|pan-liver', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$LUAD = NCG6_cancergenes.tsv$symbol[grep('lung_adenocarcinoma|lung_cancer', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$LUSC = NCG6_cancergenes.tsv$symbol[grep('lung_squamous|non-small_cell_lung_cancer', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$MESO = NCG6_cancergenes.tsv$symbol[grep('meso', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$OV = NCG6_cancergenes.tsv$symbol[grep('ovarian_serous|ovarian_cancer', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$PAAD = NCG6_cancergenes.tsv$symbol[grep('pancreatic', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$PCPG = NCG6_cancergenes.tsv$symbol[grep('pheochromocytoma,_paraganglioma', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$PRAD = NCG6_cancergenes.tsv$symbol[grep('prostate_cancer', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$READ = NCG6_cancergenes.tsv$symbol[grep('colorectal', NCG6_cancergenes.tsv$cancer_type)] # colon and rectal are bound
NCG_genelist$SARC = NCG6_cancergenes.tsv$symbol[setdiff(grep('sarcoma', NCG6_cancergenes.tsv$cancer_type),
                                                        grep('uterine', NCG6_cancergenes.tsv$cancer_type))] 
NCG_genelist$SKCM = NCG6_cancergenes.tsv$symbol[setdiff(grep('melanoma', NCG6_cancergenes.tsv$cancer_type),
                                                        union(grep('mucosal', NCG6_cancergenes.tsv$cancer_type),
                                                              grep('uvea', NCG6_cancergenes.tsv$primary_site)))]
NCG_genelist$STAD = NCG6_cancergenes.tsv$symbol[grep('gastric', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$TGCT = NCG6_cancergenes.tsv$symbol[grep('testicular', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$THCA = NCG6_cancergenes.tsv$symbol[grep('thyroid', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$THYM = NCG6_cancergenes.tsv$symbol[grep('thymic', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$UCEC = NCG6_cancergenes.tsv$symbol[grep('endometrial', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$UCS = NCG6_cancergenes.tsv$symbol[grep('uterine', NCG6_cancergenes.tsv$cancer_type)]
NCG_genelist$UVM = NCG6_cancergenes.tsv$symbol[grep('uvea', NCG6_cancergenes.tsv$primary_site)]

'''
for (i in 1:length(NCG_genelist)){
NCG_genelist[i] = list(union(unlist(NCG_genelist[i]), common))
}
'''
NCG_genes = melt(NCG_genelist)
NCG_genes_CODEsum = as.data.frame(table(NCG_genes$L1))
p<-ggplot(data=NCG_genes_CODEsum, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1, size = 12))
p + geom_text(aes(label = Freq), size = 4, hjust = 0.5, vjust = -1)

##############################################################
#### Overlaps between TCGA and NCG
##############################################################
trainsetNCG = data.frame(matrix(ncol=115))
colnames(trainsetNCG) = colnames(MC3_bailey_filtered)
for (i in 1:length(NCG_genelist)){
  trainsetNCG = rbind(trainsetNCG,
                      MC3_bailey_filtered[intersect(which(MC3_bailey_filtered$CODE == names(NCG_genelist)[i]),
                                                    which(MC3_bailey_filtered$Hugo_Symbol %in% unlist(NCG_genelist[i]))),])
}
trainsetNCG = trainsetNCG[-1,]

df = trainsetNCG %>%
  group_by(CODE) %>%
  summarise(num_samples = n_distinct(Tumor_Sample_Barcode),
            num_alterations = n())
"""k
alterations_all = as.data.frame(table(MC3_bailey_filtered$CODE))
alterations_train = as.data.frame(table(trainsetNCG$CODE))
alterations_all = merge(alterations_all,alterations_train, by = "Var1")
colnames(alterations_all) = c("CancerTypes","all_alterations","train_alteration")
alterations_all$test_alteration = alterations_all$all_alterations - alterations_all$train_alteration
alterations_all$train_percentage = alterations_all$train_alteration/alterations_all$all_alterations*100

df = melt(alterations_all)
"""
#df$variable = as.character(df$variable)
p = ggplot(df[which(!df$variable %in% c("all_alterations","train_percentage")),], 
           aes(x=CancerTypes, y = value, fill=variable))+ 
  geom_bar(stat = "identity") + 
  ylab("Number of alterations") + xlab("Cancer types") +
  theme(axis.text.x=element_text(angle=30,hjust=0.75,vjust=0.75, size = 12, face = "bold"))

p + geom_text(data=df[which(df$variable=="train_percentage"),],aes(y=df[which(df$variable=="all_alterations"),]$value,label = round(value)), 
              size = 4, hjust = 0.5, vjust = -1)

# Number of samples left after minus training set
for (i in 1:nrow(alterations_all)){
  MC3_bailey_filtered$CODE == 
    length(unique(MC3_bailey_filtered$Tumor_Sample_Barcode))
}


##############################################################
#### Overlaps between TCGA and Vogelstein
##############################################################
trainsetVogelstein = data.frame(matrix(ncol=115))
colnames(trainsetVogelstein) = colnames(MC3_bailey_filtered)
for (i in 1:nrow(Vogelstein)){
  trainsetVogelstein = rbind(trainsetVogelstein,
                             MC3_bailey_filtered[which(MC3_bailey_filtered$Hugo_Symbol %in% unlist(Vogelstein$`Gene Symbol`[i])),])
}
trainsetVogelstein = trainsetVogelstein[-1,]

df = trainsetVogelstein %>%
  group_by(CODE) %>%
  summarise(num_samples = n_distinct(Tumor_Sample_Barcode),
            num_alterations = n())

df = merge(df, summary_cohort, by="CODE")
colnames(df) = c("CancerTypes","n.sample.train","n.alterations.train","n.sample.all","n.alterations.all")
df$n.sample.test = df$n.sample.all - df$n.sample.train
df$percent.sample.train = df$n.sample.train/df$n.sample.all
df$n.alterations.test = df$n.alterations.all - df$n.alterations.train
df$percent.alterations.train = df$n.alterations.train/df$n.alterations.all

plotthis <- melt(df)
p = ggplot(plotthis[which(plotthis$variable %in% c("n.sample.train","n.sample.test")),], 
           aes(x=CancerTypes, y = value, fill=variable))+ 
  geom_bar(stat = "identity") + 
  ylab("Number of samples") + xlab("Cancer types") +
  theme(axis.text.x=element_text(angle=30,hjust=0.75,vjust=0.75, size = 12, face = "bold"))

p + geom_text(data=plotthis[which(plotthis$variable=="percent.sample.train"),],
              aes(y=plotthis[which(plotthis$variable %in% c("n.sample.all")),]$value,label = round(value*100)), 
              size = 4, hjust = 0.5, vjust = -1)


###################################################################
#### Intersection of mc3_bailey with ncg6 genes
###################################################################

# Intersecting Bailey and NCG6 file
a = c()
for (i in 1:nrow(NCG6_systemslevelproperties)){
  k = grep(NCG6_systemslevelproperties$symbol[i],MC3_bailey$Hugo_Symbol)
  a = c(a,k)
  #  overlapNCG = rbind(overlapNCG, MC3_bailey[k,])
}
MC3_bailey_overlapNCG = MC3_bailey[unique(a),]


# Intersecting Bailey and NCG6 711 cancer genes
a = c()
for (i in 1:711){
  k = grep(NCGgene[i],MC3_bailey$Hugo_Symbol)
  a = c(a,k)
  #  overlapNCG = rbind(overlapNCG, MC3_bailey[k,])
}
MC3_bailey_overlap711NCG = MC3_bailey[unique(a),]

#################

a = data.frame(table(MC3_bailey[,c("Tumor_Sample_Barcode","CODE")])) #freq = number of mutated genes in the original driver discovery file
a[which(a$Freq)]
for (i in 1:3){
  k = MC3_bailey$Gene[which(MC3_bailey$Tumor_Sample_Barcode == a[i,1])]
}

# counting number of discovered driver genes per sample
a = data.frame(table(MC3_bailey[,c("Tumor_Sample_Barcode","CODE")])) #freq = number of mutated genes in the original driver discovery file
a = a[which(a$Freq!=0)] #1 sample per tumor type
a$drivergene = 0

MC3_bailey_overlapNCG_damaged=muts1_exon_ncg_onco_cnv_damaging#muts1_exon #muts1_exon_ncg_onco_damaging #muts1_exon #
cancertype = unique(Bailey_finalgenelist$Cancer)
MC3_bailey_overlapNCG_damaged$Driver_status = 0
MC3_bailey_overlapNCG_damaged$Driver_status_Pancan = 0
cancertype=unique(MC3_bailey_overlapNCG_damaged$Cancer_type)
colnames(MC3_bailey_overlapNCG_damaged)[17]="CODE"
for (i in cancertype){
  # For each cancer type in the consensus list, grab all the driver gene names
  drivers = Bailey_finalgenelist$Gene[grep(i, Bailey_finalgenelist$Cancer)]
  
  # Find in a specific cancer, which sample contain driver gene
  # Accounting for special cases of COADREAD and PANCAN
  if (i == "COADREAD"){
    k = intersect(which(MC3_bailey_overlapNCG_damaged$CODE %in% c("COAD","READ")), which(MC3_bailey_overlapNCG_damaged$symbol_19549 %in% drivers))
    # add driver status to the specific gene & tumor sample & cancer type combination
    MC3_bailey_overlapNCG_damaged$Driver_status[k] = MC3_bailey_overlapNCG_damaged$Driver_status[k] + 1 
  } else if (i == "ESCA"){
    k = intersect(which(MC3_bailey_overlapNCG_damaged$CODE %in% c("OAC","OSCC")), which(MC3_bailey_overlapNCG_damaged$symbol_19549 %in% drivers))
    # add driver status to the specific gene & tumor sample & cancer type combination
    MC3_bailey_overlapNCG_damaged$Driver_status[k] = MC3_bailey_overlapNCG_damaged$Driver_status[k] + 1 
  } else if (i == "PANCAN"){
    k = which(MC3_bailey_overlapNCG_damaged$symbol_19549 %in% drivers)
    MC3_bailey_overlapNCG_damaged$Driver_status_Pancan[k] = MC3_bailey_overlapNCG_damaged$Driver_status_Pancan[k] + 1 
  } else {
    k = intersect(which(MC3_bailey_overlapNCG_damaged$CODE == i), which(MC3_bailey_overlapNCG_damaged$symbol_19549 %in% drivers))
    # add driver status to the specific gene & tumor sample & cancer type combination
    MC3_bailey_overlapNCG_damaged$Driver_status[k] = MC3_bailey_overlapNCG_damaged$Driver_status[k] + 1 
  }
}
table(MC3_bailey_overlapNCG_damaged$Driver_status)
length(unique(MC3_bailey_overlapNCG_damaged$Sample[which(MC3_bailey_overlapNCG_damaged$Driver_status==1)]))


length(intersect(substring(MC3_bailey$Tumor_Sample_Barcode,1,19), substring(CN_broad$Sample, 1, 19))) #8775 samples matched

length(intersect(substring(MC3_bailey$Tumor_Sample_Barcode,1,19), 
                 substring(CN_broad$Sample, 1, 19))) #8775 samples matched

###################################################################
#### Filtering for damaging mutation
###################################################################

truncating_alt = which(MC3_bailey_overlapNCG$Consequence %in% c("frameshift_variant","stop_gained","stop_lost",
                                                                "incomplete_terminal_codon_variant","start_lost")) #different to stop lost?
nonframeshift = which(MC3_bailey_overlapNCG$Consequence %in% c("inframe_deletion","inframe_insertion"))
nonsynonymous = which(MC3_bailey_overlapNCG$Consequence %in% c("missense_variant"))
splicing = which(MC3_bailey_overlapNCG$Consequence %in% c("splice_acceptor_variant","splice_donor_variant","splice_region_variant"))

notsure = c("protein_altering_alteration","transcript_ablation")

SIFT_Polyphen2 = union(grep("deleterious", MC3_bailey_overlapNCG$SIFT),# this includes also deleterious_low_confidence
                       grep("damaging", MC3_bailey_overlapNCG$PolyPhen)) # includes also probably_damaging
VEP = which(MC3_bailey_overlapNCG$IMPACT %in% c("HIGH","MODERATE"))
dam2 = intersect(union(nonframeshift, nonsynonymous), SIFT_Polyphen2)
dam_splice = intersect(splicing, VEP)

dam_all = union(union(truncating_alt, dam2), dam_splice)
MC3_bailey_overlapNCG_damaged = MC3_bailey_overlapNCG[dam_all,]

'''
ascat_acf_ploidy$Sample<-gsub("[.]","-",ascat_acf_ploidy$Sample)
length(intersect(substring(ascat_acf_ploidy$Sample,1,19),
substring(NC_broad$Sample,1,19)))
'''
###################################################################
#### Copy Number data mapping to genes, Gain Loss Distribution
###################################################################
df = merge(Thanos[,c("Sample","Copy_number","Entrez")],
           sum_seg,
           by.x=c("Sample","Entrez"), by.y =c("Sample","entrez"))
colnames(df)[c(3,9)] = c("Copy_number_Thanos","Copy_number_Trang")
tmp = melt(table(df[,c(3,9)]))
ggplot(tmp, aes(x=Copy_number_Thanos, y =Copy_number_Trang, fill=value))+
  geom_tile()+
  scale_fill_continuous(high = "red", low = "white")
k=which(df$Copy_number_Thanos != df$Copy_number_Trang) #the difference between thanos CN and mine CN

# mapping done in HPC
path= "~/Rosalind/CNV/CN_intersectBed_filtered"
code = c("GBM","OV","LUAD","LUSC", "PRAD", "UCEC", "BLCA", "TGCT", "ESCA", "PAAD", "KIRP", "LIHC", "SARC", "BRCA", "THYM",
         "MESO", "COAD", "STAD", "SKCM", "CHOL", "KIRC", "THCA", "CESC", "HNSC", "LAML", "READ", "LGG",  "DLBC", "KICH", "UCS", 
         "ACC", "PCPG", "UVM")
GainLoss_summary <- function(CODE){
  tmp1 = readRDS(sprintf("%s/%s_gainloss_2ploidy.rds",path, CODE))
  tmp2 = readRDS(sprintf("%s/%s_gainloss_segmentmean.rds",path ,CODE))
  
  df = base::merge(tmp1 %>% group_by(sampleID) %>% summarise(CNVGain_ploidy = sum(CNV_type=="Gain",na.rm = TRUE), CNVLoss_ploidy = sum(CNV_type=="Loss",na.rm = TRUE)),
                   tmp2 %>% group_by(sampleID) %>% summarise(CNVGain_seg = sum(CNV_type=="Gain"), CNVLoss_seg = sum(CNV_type=="Loss")),
                   all=TRUE)
  df$CODE = CODE
  df
}
df = data.frame()
for (i in unique(MC3_bailey$CODE)) {
  if (i =="ESCA"){
    tmp1 = readRDS(sprintf("%s/%s_gainloss_2ploidy.rds", path,i))
    tmp2 = readRDS(sprintf("%/%s_gainloss_segmentmean.rds", i))
    
    k=which(tmp1$sampleID%in% substring(oac,1,15))
    saveRDS(tmp1[k,],sprintf("%s/OAC_gainloss_2ploidy.rds", path))
    saveRDS(tmp1[-k,],sprintf("%s/OSCC_gainloss_2ploidy.rds", path))
    
    k=which(tmp2$sampleID%in% substring(oac,1,15))
    saveRDS(tmp2[k,],sprintf("%s/OAC_gainloss_OAC_gainloss_segmentmean.rds", path))
    saveRDS(tmp2[-k,],sprintf("%s/OAC_gainloss_OSCC_gainloss_segmentmean.rds", path))
    
    tmp = GainLoss_summary("OAC")
    df = rbind(df,tmp)
    tmp = GainLoss_summary("OSCC")
    df = rbind(df,tmp)
  } else {
    tmp = GainLoss_summary(i)
    df = rbind(df,tmp)
  }
}

saveRDS(df, sprintf("%s/summary.rds", path))




for (i in code){
  cnv_ann = read.table(sprintf('%s/%s.txt', path, i))
  sprintf('%s has %s samples 15 char and %s samples 19 char', i, length(unique(cnv_ann$sampleID)), length(unique(cnv_ann$Sample)))
  id = unique(cnv_ann$sampleID)
  saveRDS(id, sprintf("%s/%s_sample.rds",path,i))
}

df <- readRDS("~/Rosalind/CNV/CN_intersectBed_filtered/summary.rds")
CN_mc3 = df[which(df$sampleID %in% substring(code$Sample,1,15)),]
plotthis = melt(CN_mc3)
library(ggplot2)
ggplot(plotthis, aes(CODE,value, color=variable))+
  geom_boxplot()+
  theme_boss()

ggplot(plotthis, aes(variable,value))+
  geom_boxplot()+
  #  stat_summary(aes(label=round(..y..)), fun.y=mean, geom="text", size=4,position = position_nudge(y = 50)) +
  stat_summary(aes(label=..y..), fun.y=min, geom="text", size=4,position = position_nudge(y = -100)) +
  stat_summary(aes(label=..y..), fun.y=max, geom="text", size=4,position = position_nudge(y = 100)) +
  theme_boss()

tmp = plotthis %>% group_by(variable) %>% summarise(n_sample = n_distinct(sampleID),
                                                    CNV = sum(value, na.rm = T))

theme_boss <- function(base_size = 12, base_family = "sans"){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border=element_blank(),
      axis.line.x = element_line(color="black", size = 0.5),
      axis.line.y = element_line(color="black", size = 0.5),
      axis.text=element_text(size=16),
      axis.title=element_text(size=14),
      axis.text.x=element_text(angle=30,hjust=0.75,vjust=0.75, size = 16, face = "bold")
    )
}

###################################################################
#### Distribution of cancer genes per sample
###################################################################
# this status only flag all driver genes in the samples. This means that it can't differentiate which mutation is the driver in that gene
# MC3_bailey$Driver_status
# MC3_bailey$Driver_status_Pancan

#All exonic mutations
muts1=MC3_bailey
MC3_bailey_exonic = muts1 %>% subset(!grepl("downstream_gene_variant", muts1$Consequence) &
                                       !grepl("upstream_gene_variant", muts1$Consequence) &
                                       #                                       !grepl("intergenic", muts1$Func.refGene) &
                                       !grepl("non_coding_transcript_variant",muts1$Consequence) &
                                       #                                       !grepl("non_coding_transcript_exon_variant", muts1$Consequence) &
                                       !grepl("intron_variant", muts1$Consequence) &
                                       !grepl("5_prime_UTR_variant", muts1$Consequence) &
                                       !grepl("3_prime_UTR_variant", muts1$Consequence))


library(dplyr)
MC3_bailey_overlapNCG_damaged$Sample = substring(MC3_bailey_overlapNCG_damaged$sample,1,12) 
df = MC3_bailey_overlapNCG_damaged[which(MC3_bailey_overlapNCG_damaged$Driver_status_Pancan==1),] %>%
  group_by(Sample) %>%
  summarise(num_driver = n_distinct(symbol_19549))
nodriver = setdiff(unique(CODE$Sample), df$Sample)
df[(nrow(df)+1):(nrow(df)+length(nodriver)),1] = nodriver
df[(9079-length(nodriver)+1):9079,2] = 0

#code = unique(MC3_bailey_overlapNCG_damaged[,c("Sample","CODE")])
df = merge(df, CODE, by='Sample')

#driver_per_type = df
#driver_per_type$Analysis_Type = "individual cancer"
df$Analysis_Type = 'PANCAN'
df2 = rbind(df, driver_per_type)

#Driver_exon2= df2
Driver_damaging2=df2
library(ggplot2)
library(reshape2)
plotthis = df2 #data.frame(table(df$CODE,df$num_driver))
plotthis %>% mutate(label=replace(label, "PANCAN (9072)")) %>%
  +     as.data.frame()
### Driver distribution (bailey mutation VEP)
tmp_sum = driver_per_type %>% group_by(CODE) %>% summarise(n_all = length(unique(Sample)))
tmp_sum$label = paste0(tmp_sum$CODE, " (",tmp_sum$n_all,")")
#tmp_sum[35,] = c("PANCAN",9072,"PANCAN (9072)")
plotthis = merge(plotthis,tmp_sum)
plotthis$category = "individual"
plotthis = rbind(plotthis, plotthis %>% mutate(label=replace(label, TRUE,"PANCAN (9072)"), category="pancan") %>% as.data.frame())

order = 
p <- ggplot(plotthis, aes(label, num_driver),dodge=Type) + 
  geom_boxplot(aes(fill = Analysis_Type), alpha=0.2) + 
  geom_hline(yintercept=3, linetype="dashed", color = "red")+
  #stat_summary(aes(label=..y..), fun.y=median, geom="label", size=4) +
  #stat_summary(aes(label=..y..), fun.y=min, geom="label", size=4) +
  #stat_summary(aes(label=..y..), fun.y=max, geom="label", size=4) +
  xlab('Cancer types') + ylab('number of driver genes per sample')+
  facet_grid(.~category,space = "free_x", scales = "free_x")+
  theme_boss_xtilted()
p
ggsave(filename=paste0("~/Rosalind/Plots/Bailey_Driver_per_sample_damaging.png"), 
       plot=p,width = 25, height = 15, dpi = 300)


### Plotthing the number of samples with 1,2,3 drivers etc
df = rbind(Driver_exon2%>% mutate(Type="exonic"),
           Driver_damaging2 %>% mutate(Type="damaging"))
df$Type = factor(df$Type,levels = c("exonic", "damaging"))
df = data.frame(table(df[,c(2,5,6)]))
df = df %>% group_by(Analysis_Type, Type) %>% mutate(percent = Freq/sum(Freq),
                                               label = ifelse(as.numeric(as.character(num_driver))>10,">10", as.character(num_driver)))
df$label=factor(df$label,levels=c(0:10,">10"))
df = df %>% group_by(Analysis_Type, Type,label) %>% summarise(percent_corrected = sum(percent),
                                                      Freq_corrected = sum(Freq))

df$percent_corrected_char = paste0(as.character(round(df$percent_corrected*100)),"%")
##  Bar plot
p <- ggplot(df, aes(label, Freq_corrected,fill=Analysis_Type)) + 
  geom_bar(stat='identity',position=position_dodge(),alpha=0.8) + 
  ylab('Number of samples') + xlab('Number of cancer genes per sample')+
  facet_wrap(.~Type,strip.position ="top")+
  geom_text(aes(y=Freq_corrected, label=percent_corrected_char), 
            vjust=1.5, color="black", position = position_dodge(0.9), size=7)+
  #  theme(axis.text.x=element_text(angle=30,hjust=0.75,vjust=0.75, size = 12, face = "bold"))
  theme_boss(base_size = 40)
p

pdf(file = "~/Trang/Plots/Bailey_NumberDriver_per_sample_exonic_damaging_19022019.pdf", width = 25, height = 10)
p
dev.off()

ggsave(filename=paste0("~/Trang/Plots/Bailey_NumberDriver_per_sample_exonic_damaging.pdf"), 
       plot=p,width = 25, height = 10, dpi = 300)

 ##  Piechart
p <- ggplot(df[which(df$Type="exonic"),]) + 
  geom_bar(aes(x="", percent_corrected,fill=label), stat='identity',width = 1,alpha=0.8) + 
  coord_polar("y", start=0)+
  theme_void()+
  geom_text(aes(x=1, y = cumsum(percent_corrected) - percent_corrected/2, label=percent_corrected))+ 
  ylab('Number of samples') + xlab('Number of driver genes per sample')+
  facet_wrap(.~Analysis_Type)+
  theme_boss(base_size = 40)
p

ggsave(filename=paste0("~/Trang/Plots/Bailey_NumberDriver_per_sample_exonic_damaging.png"), 
       plot=p,width = 25, height = 10, dpi = 300)

#### sample ID: [1:19] characters are portion for TCGA
k=intersect(unique(substring(CN_broad$Sample,1,19)), unique(substring(MC3_bailey$Tumor_Sample_Barcode,1,19)))
CN_broad$SampleID = substring(CN_broad$Sample,1,19)
CN_9079 = CN_broad[which(CN_broad$SampleID %in% k),]
CN_9079 = merge(CN_9079, code, by='SampleID', all.x =TRUE)
#CN = CN[,c(2,3,4,5,6,1,8,10)]
for (i in unique(code$CODE)){
  tmp = CN_9079[which(CN$CODE == i),]
  write.table(tmp, file=paste0('/Users/let/Rosalind/CN_segmentmean/',i,".bed"),
              row.names = F, col.names = F, sep = "\t", quote = F)
}


###################################################################
#### Distribution of damaging mutation per cohort
###################################################################
theme_boss <- function(base_size = 12, base_family = "sans"){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border=element_blank(),
      axis.line.x = element_line(color="black", size = 0.5),
      axis.line.y = element_line(color="black", size = 0.5),
      axis.text=element_text(size=base_size),
      axis.title=element_text(size=base_size, face = "bold")
    )
}

theme_boss_xtilted <- function(base_size = 12, base_family = "sans"){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border=element_blank(),
      axis.line.x = element_line(color="black", size = 0.5),
      axis.line.y = element_line(color="black", size = 0.5),
      axis.text=element_text(size=16),
      axis.title=element_text(size=14),
      axis.text.x=element_text(angle=30,hjust=0.75,vjust=0.75, size = 12, face = "bold")
    )
}
library(dplyr)


df = muts1_exon_ncg %>%
  group_by(Sample) %>%
  summarise(num_alterations = n(),
            num_gene = n_distinct(Gene.refGene))



dist_bailey = merge(df, code, by='Tumor_Sample_Barcode')
dist_bailey$Type = "all"
dist_overlappNCG = merge(df, code, by='Tumor_Sample_Barcode')
dist_overlappNCG$Type = "all_ncg6"
dist_damaging = merge(df, code, by.x = "Sample", by.y='Tumor_Sample_Barcode')
dist_damaging$Type = "damaging"

all_dist = rbind(dist_bailey[,-1],dist_overlappNCG[,-1],dist_damaging[,-1])
plotthis = all_dist[which(all_dist$CODE=="ESCA"),]
plotthis = plotthis[which(plotthis$SampleID %in% substring(oacs,1,19)),]

library(ggplot2)
total_table_CODE <- readRDS("~/Trang/tabl") #readRDS("~/Rosalind/table_muts_snvindelGOF_CODE.rds")
total_table_CODE$CODE[which(total_table_CODE$CODE=="ESCA")] = "OSCC"
code = unique(total_table_CODE[,c("CODE", "Sample")])
plotthis = melt(total_table_CODE)
for (i in unique(code$CODE)) {
  tmp = plotthis[which(plotthis$CODE==i),]
  tmp_sum = tmp %>% group_by(Type) %>% summarise(n_all = length(unique(Sample)))
  stat = tmp %>% group_by(variable,Type) %>% summarise(med = median(value),maxim = max(value),minim = min(value))
  p=ggplot(tmp %>% mutate(value = replace(value, value==0, 0.9999)), aes(Type, value)) + 
    geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(Type)), alpha=0.1) +
    scale_y_log10() + 
    geom_boxplot(aes(fill = Type), alpha=0.2) + 
    xlab('Cancer types') + ylab('number per sample')+
    stat_summary(aes(label=..y..), fun.y=median, geom="label", size=4) +
    #stat_summary(aes(label=..y..), fun.y=min, geom="label", size=4) +
    #stat_summary(aes(label=..y..), fun.y=max, geom="label", size=4) +
    facet_grid(.~variable)+
    annotate("label", x=stat$Type, y=replace(stat$med, stat$med==0, 0.9999), label=stat$med) +
    annotate("label", x=1:nrow(stat), y=replace(stat$maxim, stat$maxim==0, 0.9999), label=stat$maxim) + 
    annotate("label", x=1:nrow(stat), y=replace(stat$minim, stat$minim==0, 0.9999), label=stat$minim) +
    xlab("")+
    ggtitle(paste0(i, sprintf(" (n_all=%s,n_ncg6=%s,n_damaging=%s)", tmp_sum$n_all[1], tmp_sum$n_all[2], tmp_sum$n_all[3]))) +
    theme(axis.text.x=element_text(angle=30,hjust=0.75,vjust=0.75, size = 12, face = "bold"))
  ggsave(filename=paste0("~/Rosalind/Plots/damaging_snvindel_Onco/",i,"_withOutliers.png"), 
         plot=p,width = 16, height = 9, dpi = 100)
}

##########################################
###### CN & segment mean distribution
##########################################
plotthis=rbind(rbind(bailey_0,bailey_1),rbind(bailey_2,bailey_3))
plotthis = merge(plotthis, CODE)
tmp=melt(plotthis)
tmp = tmp[which(tmp$variable=="mutated_genes"),]
#tmp_sum = tmp[!which(tmp$Type=="damaging"),] %>% group_by(CODE) %>% summarise(n_all = length(unique(Sample)))
#tmp_sum = tmp %>% group_by(Type) %>% summarise(n_all = length(unique(Sample)))
#tmp_sum$label = paste0("PANCAN", " (",tmp_sum$n_all,")")

tmp$label = "PANCAN (9079)"
stat = tmp %>% group_by(label,Type) %>% summarise(med = median(value),maxim = max(value),minim = min(value))# %>% arrange(med) #arrange(desc(med))

p=ggplot(tmp %>% mutate(value = replace(value, value==0, 0.9999)), aes(x=Type, y=value)) +
  geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(Type)), alpha=0.1) +
  scale_color_manual(values = c("dodgerblue1","deepskyblue2","blue3")) +
  geom_boxplot(alpha=0.2) + #aes(fill=factor(category))
  scale_y_log10() + 
  theme_boss(base_size = 40)+
  #facet_wrap(~category, scales = 'free')+
  annotate("label", x=1:nrow(stat), y=replace(stat$med, stat$med==0, 0.9999), label=stat$med, size = 13) +
  annotate("label", x=1:nrow(stat), y=replace(stat$maxim, stat$maxim==0, 0.9999), label=stat$maxim, size = 13) + 
  annotate("label", x=1:nrow(stat), y=replace(stat$minim, stat$minim==0, 0.9999), label=stat$minim, size = 13) +
  ylab("Number of mutated genes") + xlab("Category")
p
ggsave(filename=paste0("~/Trang//Plots/damaging_snvindel_Onco/MutatedGenes_all_distribution.png"), 
       plot=p,width = 20, height = 10, dpi = 300)


### plot the mutations and mutated genes distribution across cancers
total_table_CODE <- readRDS("~/Trang/table_muts_snvindelGOF_CODE.rds")
total_table_CODE$CODE[which(total_table_CODE$CODE=="ESCA")] = "OSCC"
code = unique(total_table_CODE[,c("CODE", "Sample")])
plotthis = melt(total_table_CODE)

tmp=melt(plotthis)
tmp_sum = tmp[which(tmp$Type=="exonic_ncg6"),] %>% group_by(CODE) %>% summarise(n_all = length(unique(Sample)))
tmp_sum$label = paste0(tmp_sum$CODE, " (",tmp_sum$n_all,")")

plotthis2 = merge(tmp[intersect(which(tmp$variable=="mutations"), which(tmp$Type == "exonic_ncg6")),], tmp_sum) #merge(tmp,tmp_sum)
plotthis2 = plotthis2 %>% group_by(CODE) %>% mutate(median = median(value))
plotthis2$category = "cancertype"

plotthis_pancan = plotthis2
plotthis_pancan$label="PANCAN (9072)"
plotthis_pancan$category = "PANCAN"
plotthis_pancan$median = median(plotthis_pancan$value)

stat = plotthis2 %>% group_by(label,category) %>% summarise(med = median(value),maxim = max(value),minim = min(value)) %>% arrange(med) #arrange(desc(med))
stat = rbind(stat, plotthis_pancan %>% group_by(label,category) %>% summarise(med = median(value),maxim = max(value),minim = min(value)))

plotthis2 = rbind(plotthis2, plotthis_pancan)

order = c(levels(reorder(plotthis2$label, plotthis2$median))[-19],"PANCAN (9072)")
p=ggplot(plotthis2 %>% mutate(value = replace(value, value==0, 0.9999)), aes(x=factor(label, order), y=value)) +
  geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(category)), alpha=0.1) +
  geom_boxplot(alpha=0.2) + #aes(fill=factor(category))
  scale_y_log10() + 
  theme_boss_xtilted()+
  #facet_wrap(~category, scales = 'free') +#, space = "free_x", scales = "free_x")+
  annotate("label", x=1:nrow(stat), y=replace(stat$med, stat$med==0, 0.9999), label=stat$med, na.rm = TRUE) +
  annotate("label", x=1:nrow(stat), y=replace(stat$maxim, stat$maxim==0, 0.9999), label=stat$maxim) + 
  annotate("label", x=1:nrow(stat), y=replace(stat$minim, stat$minim==0, 0.9999), label=stat$minim) +
  ylab("Number of exonic mutations") + xlab("Cohort(#samples)")
p
ggsave(filename=paste0("~/Rosalind/Plots/damaging_snvindel_Onco/ExonicMutations_ncg6_distribution_withOutliers.png"), 
       plot=p,width = 25, height = 15, dpi = 300)
dev.off()

ggsave(filename=paste0("~/Rosalind/Plots/damaging_snvindel_Onco/Mutation_damaging_PANCAN_distribution_withOutliers.png"), 
       plot=p,width = 4, height = 15, dpi = 300)

###### CN & segment mean distribution
code = unique(mc3.v0.2.8.PUBLIC.code.filtered[,c("CODE","Tumor_Sample_Barcode")])
code$sampleID=substring(code$Tumor_Sample_Barcode,1,15)
code$CODE=as.character(code$CODE)
code[which(code$sampleID %in% oac),]$CODE ="OAC"
code[which(code$CODE =="ESCA"),]$CODE ="OSCC"

broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted$sampleID=substring(broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted$Sample,1,15)
CN = merge(broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted,code, by="sampleID")

summarySample = cbind(table(code$CODE),
                      CN %>% group_by(CODE) %>% summarise(count=n_distinct(sampleID)))
summarySample = summarySample[,c(3,2,4)]
colnames(summarySample)[2:3] = c("n_muts","n_cnv")

#Density plots of segment means in 34 cancer types
p <- ggplot(CN, aes(Segment_Mean, colour=CODE, fill=CODE))+ 
  geom_density(alpha=0.55)
ggsave(filename=paste0("~/Rosalind/Plots/CN_segmentmean/CNsegmentmean_distribution.png"), 
       plot=p,width = 25, height = 15, dpi = 300)


df = CN %>% group_by(sampleID) %>% summarise(CN_mean = mean(Segment_Mean))
df = merge(df,code, by="sampleID")
df_pancan = df 
df_pancan$CODE ="PANCAN"
df = rbind(df,df_pancan)
df = merge(df,
           df %>% group_by(CODE) %>% summarise(n_sample = n_distinct(sampleID)),
           by="CODE")
df$label = paste0(df$CODE," (",df$n_sample,")")
ggplot(df , aes(reorder(label,CN_mean),CN_mean))+geom_boxplot()+
  stat_summary(aes(label=round(..y..,2)), fun.y=median, geom="text", size=4) +
  stat_summary(aes(label=round(..y..,2)), fun.y=min, geom="text", size=4) +
  stat_summary(aes(label=round(..y..,2)), fun.y=max, geom="text", size=4) +
  ylab('Mean of CN Segment_Mean per sample')+
  xlab('Cancer type')+
  theme(axis.text.x=element_text(angle=30,hjust=0.75,vjust=0.75, size = 15, face = "bold"))
ggsave(filename=paste0("~/Rosalind/Plots/CN_segmentmean/CNsegmentmean_meanperSample_distribution.png"), width = 25, height = 15, dpi = 300)


TCGA_mastercalls.abs_tables_JSedit.fixed <- read.delim("~/Downloads/TCGA_mastercalls.abs_tables_JSedit.fixed.txt")
ploidy = TCGA_mastercalls.abs_tables_JSedit.fixed[,c(1,2,4,5)]
ploidy$sampleID = substring(ploidy$sample,1,15)
ploidy2 = merge(ploidy,code,by.x = "array", by.y = "sampleID")
tmp=ploidy2 %>% group_by(CODE) %>% summarise(n_ploidy = n_distinct(array))
summarySample = merge(summarySample,tmp)
summarySample$label = paste0(summarySample$CODE," (",summarySample$n_ploidy,")")
ploidy2 = merge(ploidy2, summarySample[,c("CODE","label")], all = F)
ploidy2 %>% group_by(CODE) %>% summarise(Mean = mean(ploidy))
geom_boxplot()+
  theme(axis.text.x=element_text(angle=30,hjust=0.75,vjust=0.75, size = 15, face = "bold"))

###################################################################
#### intersection of snv cnv and indels
###################################################################

library(VennDiagram)
draw.triple.venn(area1 = 9079, area2 = 10965, area3 = 10786, n12 = 8774, n23 = 8025, n13 = 6349, n123 = 6184, 
                 category = c("snv_indels (9079)", "cnv (10965)", "ploidy (8524)"), 
                 lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid") , 
                 cex=2, cat.cex=2, 
                 cat.fontfamily = rep("serif", 3))
dev.off()

length(unique(substring(CN_broad$Sample,1,16)))
draw.pairwise.venn(area1 = 9079, area2 = 10965, cross.area = 8774, 
                   category = c("snv_indels", "cnv"), 
                   lty = "blank", 
                   fill = c("skyblue", "pink1") , 
                   cex=2, cat.cex=2)

###########################################################
## This function combines all types of data to single table 
## and prepares them for the extraction of drivers and prediction
createTotalTable = function(muts=NULL, cnvs=NULL, svs=NULL, exclude_samples=NULL){
  
  ns = c("nonsynonymous","stopgain","frameshift deletion","splicing","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","stoploss","frameshift substitution")
  dam = c("nonsynonymous","frameshift deletion","frameshift insertion","frameshift substitution","splicing","stopgain","stoploss")
  trunc = c("frameshift deletion","frameshift insertion","frameshift substitution","stopgain","stoploss") ## Always damaging==TRUE
  non_trunc = c("nonsynonymous","splicing")
  
  ## Make the lists
  message("Integrating SNVs...")
  df_mut = muts
  
  if(!is.null(exclude_samples)){
    df_mut = df_mut %>% subset(!sample%in%exclude_samples)
    message(paste0("Samples excluded in mutation data: ", paste0(exclude_samples, collapse=",")))
  }
  
  ## Fix nonsilent here
  df_mut = df_mut %>% select(-nonsilent)
  df_mut=is_nonsilent(df_mut)
  
  ## In order to get the number of all mutations per gene and because
  ## I have WGS data, I exclude mutations that fall in the following categories
  
  ## Exclude genes that are not in 19014
  df_mut = df_mut %>% subset(!is.na(entrez_19549))
  
  ## Add oncodriveClust
  df_mut = df_mut %>% left_join(onco_out %>% select(IGV, oncodriveClust), by = 'IGV')
  
  ## Create the total table
  total_muts = ddply(df_mut, .(sample, symbol_19549, entrez_19549), dplyr::summarise,
                     no_ALL_muts=n(),
                     no_NSI_muts=sum(nonsilent),
                     no_TRUNC_muts = sum(ExonicFunc.refGene %in% trunc),
                     no_NTDam_muts = sum(ExonicFunc.refGene %in% non_trunc & damaging),
                     no_GOF_muts = sum(oncodriveClust), .progress = 'text'
  )
  
  ## Add protein position if needed
  #aa = df_mut %>% select(sample, entrez_19549, symbol_19014, AAChange.refGene) %>% group_by(sample, entrez_19549, symbol_19014) %>% summarise(AAChange=paste(unique(AAChange.refGene), collapse=","))
  #total_muts = total_muts %>% left_join(aa)
  
  # ## Check that you see a difference in the number of total mutations
  # geneInfo_fn="/Volumes/mourikisa/data/geneInfoNCG5.Rdata"
  # cancerGenes_fn="/Volumes/mourikisa/data/cancerGenesNCG5.Rdata"
  # load(geneInfo_fn)
  # load(cancerGenes_fn)
  # ## Fix gene info table from NCG
  # geneInfo = geneInfo %>% select(entrez, cancer_type) %>% unique
  # ## Get a cancer gene with all the associated primary sites and cancer sites
  # cancerGenes = cancerGenes %>% select(entrez, primary_site, cancer_site) %>%
  #     group_by(entrez) %>% summarise(primary_site=paste(unique(primary_site), collapse=","),
  #                                    cancer_site=paste(unique(cancer_site), collapse=",")) %>%
  #     ungroup
  # geneInfo = geneInfo %>% left_join(cancerGenes, by=c("entrez"))
  # test = total_muts %>% left_join(geneInfo, by=c("entrez_19549"="entrez"))
  # 
  # test = test %>% mutate(sumDrivers = rowSums(.[6:8]))
  # test = test %>% mutate(dVa=sumDrivers/no_ALL_muts)
  # test %>% subset(dVa==1) %>% head
  # wilcox.test(test%>%subset(cancer_type=="cgc")%>%.$dVa, test%>%subset(cancer_type=="can")%>%.$dVa)
  # summary(test%>%subset(cancer_type=="cgc")%>%.$dVa)
  # summary(test%>%subset(cancer_type=="can")%>%.$dVa)
  
  
  ## Bring in the total_table the SVs
  message("Integrating SVs...")
  
  if(!is.null(exclude_samples)){
    svs = svs %>% subset(!sample%in%exclude_samples)
    message(paste0("Samples excluded in SV data: ", paste0(exclude_samples, collapse=",")))
  }
  
  if(!is.null(svs)){
    ## There are 2 genes duplicated in 4 samples due to aliases
    svs = svs %>% subset(!is.na(entrez_19549)) %>% mutate(key=paste(sample, entrez_19549, sep=".")) %>% subset(!duplicated(key)) %>% select(-key)
    svs = svs %>% subset(BND>0 | INS>0 | INV>0) %>% select(-symbol, -cancer_type, -primary_site, -cancer_site)
    ## And put them in the total table as well
    total_table = total_muts %>% full_join(svs%>%select(sample, BND, INS, INV, entrez_19549)%>%subset(!is.na(entrez_19549)))
    total_table = total_table %>% mutate(key=paste(sample, entrez_19549, sep="."))
  }else{
    total_table = total_muts %>% mutate(BND=0, INS=0, INV=0)
    total_table = total_table %>% mutate(key=paste(sample, entrez_19549, sep="."))
  }
  
  
  
  
  message("Integrating CNVs...")
  
  ## Add also the genes that are in muts and SVs to get their ploidy and Copy number
  cnvs = cnvs %>% mutate(key=paste(sample, entrez_19549, sep="."))
  df_cnv = cnvs %>% subset(key%in%total_table$key | !is.na(CNV_type_corrected)) ## Also get the real CNVs
  
  
  if(!is.null(exclude_samples)){
    df_cnv = df_cnv %>% subset(!sample%in%exclude_samples)
    message(paste0("Samples excluded in CNV data: ", paste0(exclude_samples, collapse=",")))
  }
  
  ## define Gains and Losses - this was done in previous step
  ## Deduplicate the CNV data, because some genes may fall into two regions (sometimes it can be gain and loss)
  message("Resolving duplicated entries in CNVs...")
  dups = df_cnv %>% group_by(key) %>% mutate(n=n(), types=paste(unique(CNV_type_corrected), collapse=","), ntypes=length(unique(CNV_type_corrected))) %>% subset(n>1)
  dups = dups %>% ungroup()
  ## In here you will find two kinds of duplications those that are duplicates but associated with one type of CNV (i.e Gain/Loss)
  ## And those that are associated with two types of CNVs
  ## I didn't use overlap function in the end
  overlap <- function(start1, end1, start2, end2){
    res = pmin(end1, end2) - pmax(start2, start1)
    if(res>=0){
      return(res) 
    }else if(res<0){
      res=0
      return(res) 
    } 
  }
  dups_refined = NULL
  for (t in unique(dups$types)){
    if (grepl(",NA|NA,", t)){ ## When arrange always on top will be Gain/Loss and those will be selected when deduplicate
      d = dups %>% subset(types==t) %>% arrange(key, desc(CNV_type_corrected)) %>% subset(!duplicated(key))
      dups_refined = rbind(dups_refined, d)
    }else if (t=="Gain" | t=="NA" | t=="Loss"){ ## Choose the one with the highest overlap
      d = dups %>% subset(types==t) %>% mutate(overlap=(overlap(start, end, Start, End)/(end-start))*100) %>% arrange(key, desc(overlap)) %>% subset(!duplicated(key)) %>% select(-overlap)
      dups_refined = rbind(dups_refined, d)
    }else if (grepl("Gain", t) & grepl("Loss", t)){ ## For those genes we have both gain and loss, I set CNV_type to NA because we cannot distinguish between the two
      d = dups %>% subset(types==t) %>% mutate(CNV_type_corrected=NA) %>% arrange(key) %>% subset(!duplicated(key))
      dups_refined = rbind(dups_refined, d)
    }
  }
  
  ## For now deduplicte them and keep as CNV_type_corrected the concatenation of both types to see how many they are in the drivers
  ## Take them out first from the df_cnv
  df_cnv = df_cnv %>% subset(!key%in%dups$key)
  df_cnv$n = 1
  ## Fix dups
  dups_refined = dups_refined %>% select(-types, -ntypes)
  ## Put the back in the df_cnv
  df_cnv = rbind(df_cnv, dups_refined)
  
  ## Create total table from mutations and CNVs
  total_table = total_table %>% subset(!is.na(entrez_19549)) %>% 
    full_join(df_cnv%>%select(sample, entrez_19549, Total_CN, CNV_type_corrected, ploidy, n)%>%rename(CNV_entries=n)%>%subset(!is.na(entrez_19549)))
  
  
  total_table$na_19549 = apply(total_table[,c("symbol_19549", "entrez_19549")], 1, function(x) length(x[is.na(x)]))
  total_table = total_table %>% mutate(in_19549=ifelse(na_19549<2, TRUE, FALSE)) %>% select(-na_19549) %>% subset(in_19549==TRUE)
  
  return(total_table)
  
}
for (i in 1:nrow(CN_9079)){
  CN_9079[i,"gain_loss"] = 1
}

dist_damaging %>%
  group_by(CODE) %>%
  summarise(n_sample = n_distinct(Sample),
            n_mutation = sum(num_alterations)) -> df

snvindel[which(snvindel$damaging==TRUE),] %>%
  group_by(CODE) %>%
  summarise(n_mutation = n(),
            n_gene = n_distinct(symbol_19549)) ->df2


######################################################################################
### separate snvindels to 34 files for each cancer (easier processed to total table)
######################################################################################
snvindel = set_IGV_code(mc3_bailey_ann_damaging_oncodriveclust)
code = unique(snvindel$CODE)
for (i in code){
  tmp = snvindel[which(snvindel$CODE==i),]
  saveRDS(tmp, sprintf("/Users/let/Trang/OncodriveClust/result_12022019/34cancer_types/%s.rds", i))
}


######################################################################################
### create total table
######################################################################################
library(plyr)
library(dplyr)
library(tidyr)

for (i in code){
  df_mut = readRDS(sprintf("/Users/let/Rosalind/OncodriveClust/result_30112018/34cancer_types/%s.rds", i))
  df_mut$sample = substring(df_mut$Sample,1,15)
  
  ## Fix nonsilent here
  df_mut = df_mut %>% select(-nonsilent)
  df_mut=is_nonsilent(df_mut)
  
  # Exclude all non-exonic mutations
  df_mut = df_mut %>% subset(Func.refGene!="" &
                               !grepl("downstream", df_mut$Func.refGene) &
                               !grepl("upstream", df_mut$Func.refGene) &
                               !grepl("intergenic", df_mut$Func.refGene) &
                               !grepl("ncRNA", df_mut$Func.refGene) &
                               !grepl("intronic", df_mut$Func.refGene) &    
                               !grepl("UTR", df_mut$Func.refGene))
  
  df_mut = df_mut %>% subset(!is.na(entrez_19549))
  
  
  ## Create the total table
  total_muts = ddply(df_mut, .(sample, symbol_19549, entrez_19549), dplyr::summarise,
                     no_ALL_muts=n(),
                     no_NSI_muts=sum(nonsilent),
                     no_TRUNC_muts = sum(ExonicFunc.refGene %in% trunc),
                     no_NTDam_muts = sum(ExonicFunc.refGene %in% non_trunc & damaging),
                     no_GOF_muts = sum(oncodriveClust), .progress = 'text'
  )
  
  total_table = total_muts %>% mutate(BND=0, INS=0, INV=0)
  total_table = total_table %>% mutate(key=paste(sample, entrez_19549, sep="."))
  
  message("Integrating CNVs...")
  
  cnvs = readRDS(sprinf("/Users/let/Rosalind/CNV/CN_intersectBed_filtered/%s_gainloss_2ploidy.rds",i))
  colnames(cnvs)[c(1,13)] = c("sample","entrez_19549")
  ## Add also the genes that are in muts and SVs to get their ploidy and Copy number
  cnvs$key =paste(cnvs$sample, cnvs$entrez_19549, sep=".")
  
  
  #This will unfortunately remove all non-exonic regions from cnv file (so no Copy_number = 0 anymore)
  df_cnv = cnvs %>% subset(key%in%total_table$key| !is.na(CNV_type)) ## Also get the real CNVs
  
  ## define Gains and Losses - this was done in previous step
  ## Deduplicate the CNV data, because some genes may fall into two regions (sometimes it can be gain and loss)
  message("Resolving duplicated entries in CNVs...")
  dups = df_cnv %>% group_by(key) %>% mutate(n=n(), types=paste(unique(CNV_type), collapse=","), ntypes=length(unique(CNV_type))) %>% subset(n>1)
  dups = dups %>% ungroup()
  ## In here you will find two kinds of duplications those that are duplicates but associated with one type of CNV (i.e Gain/Loss)
  ## And those that are associated with two types of CNVs
  ## I didn't use overlap function in the end
  overlap <- function(start1, end1, start2, end2){
    res = pmin(end1, end2) - pmax(start2, start1)
    if(res>=0){
      return(res) 
    }else if(res<0){
      res=0
      return(res) 
    } 
  }
  dups_refined = NULL
  for (t in unique(dups$types)){
    if (t=="NA"){ ## When arrange always on top will be Gain/Loss and those will be selected when deduplicate
      d = dups %>% subset(types==t) %>% arrange(key, desc(CNV_type)) %>% subset(!duplicated(key))
      dups_refined = rbind(dups_refined, d)
    }else if (t=="Loss,NA" | t=="Gain,NA" | t=="NA,Loss" | t=="Loss"){ ## Choose the one with the highest overlap
      d = dups %>% subset(types==t) %>% mutate(overlap=(overlap(Start, End, gstart, gend)/(gend-gstart))*100) %>% arrange(key, desc(overlap)) %>% subset(!duplicated(key)) %>% select(-overlap)
      dups_refined = rbind(dups_refined, d)
    }else if (grepl("Gain", t) & grepl("Loss", t)){ ## For those genes we have both gain and loss, I set CNV_type to NA because we cannot distinguish between the two
      d = dups %>% subset(types==t) %>% mutate(CNV_type=NA) %>% arrange(key) %>% subset(!duplicated(key))
      dups_refined = rbind(dups_refined, d)
    }
  }
  
  ## For now deduplicte them and keep as CNV_type_corrected the concatenation of both types to see how many they are in the drivers
  ## Take them out first from the df_cnv
  df_cnv = df_cnv %>% subset(!key%in%dups$key)
  df_cnv$n = 1
  ## Fix dups
  dups_refined = dups_refined %>% select(-types, -ntypes)
  ## Put the back in the df_cnv
  df_cnv = rbind(df_cnv, dups_refined)
  
  ## Create total table from mutations and CNVs
  total_table = total_table %>% subset(!is.na(entrez_19549)) %>% 
    full_join(df_cnv%>%select(sample, entrez_19549, Copy_number, CNV_type, ploidy, n)%>%rename(CNV_entries=n)%>%subset(!is.na(entrez_19549)))
  
  saveRDS(total_table,sprintf("/Users/let/Rosalind/OncodriveClust/result_30112018/34cancer_types/%s.rds", i))
}



muts1_exon_ncg_onco_damaging = muts2_totalTable %>% subset(no_TRUNC_muts!=0 |
                                                             no_NTDam_muts!=0 |
                                                             no_GOF_muts!=0
                                                           (CNVGain==1) | 
                                                             (Copy_number==0) |
                                                             ((Copy_number==1) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0))
                                                           # BND!=0 | 
                                                           # INS!=0 | 
                                                           # INV!=0)
)



###############################################################################


library(plyr)
library(dplyr)
library(tidyr)


ns = c("nonsynonymous","stopgain","frameshift deletion","splicing","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","stoploss","frameshift substitution")
dam = c("nonsynonymous","frameshift deletion","frameshift insertion","frameshift substitution","splicing","stopgain","stoploss")
trunc = c("frameshift deletion","frameshift insertion","frameshift substitution","stopgain","stoploss") ## Always damaging==TRUE
non_trunc = c("nonsynonymous","splicing")

is_nonsilent=function(x){
  x$nonsilent=x$ExonicFunc.refGene%in%ns
  x
}

muts_path = "/mnt/lustre/users/k1801782/dataset/OncodriveClust/result_30112018/34cancer_types"
cnv_path = "/mnt/lustre/users/k1801782/dataset/CNV/CN_intersectBed_filtered"

code = c("LUAD","GBM","OV","LUSC","TGCT","OAC","LIHC","SARC","SKCM","UCEC","COAD","BRCA","LAML","READ",
         "KIRC","KIRP","THCA","BLCA","STAD","CESC","HNSC","LGG","PRAD","PAAD","OSCC","UCS","MESO","THYM",
         "UVM","DLBC","ACC","PCPG","CHOL","KICH")

for (i in code){
  if(FALSE){
    df_mut = readRDS(sprintf("%s/%s.rds", muts_path, i))
    df_mut$sample = substring(df_mut$Sample,1,15)
    
    ## Fix nonsilent here
    df_mut = df_mut %>% select(-nonsilent)
    df_mut=is_nonsilent(df_mut)
    
    # Exclude all genes not in ncg6
    df_mut = df_mut %>% subset(!is.na(entrez_19549))
    
    # Exclude all non-exonic mutations
    df_mut2 = df_mut %>% subset(Func.refGene!="." &
                                  !grepl("downstream", df_mut$Func.refGene) &
                                  !grepl("upstream", df_mut$Func.refGene) &
                                  !grepl("intergenic", df_mut$Func.refGene) &
                                  !grepl("ncRNA", df_mut$Func.refGene) &
                                  !grepl("intronic", df_mut$Func.refGene) &    
                                  !grepl("UTR", df_mut$Func.refGene))
    
    df_mut = df_mut %>% mutate(key=paste(sample, entrez_19549, sep="."))
    
    total_table = total_muts %>% mutate(key=paste(sample, entrez_19549, sep="."))
    total_table2 = total_muts2 %>% mutate(key=paste(sample, entrez_19549, sep="."))
    
    message("Integrating CNVs...")
  }
  #cnvs = readRDS(sprintf("%s/%s_gainloss_2ploidy.rds",cnv_path,i))
  cnvs = read.table(sprintf("%s/%s.txt",cnv_path,i))
  
  cnvs$Copy_number = round(cnvs$CN)
  colnames(cnvs)[13] = "entrez_19549"#c("sample","entrez_19549")
  ## Add also the genes that are in muts and SVs to get their ploidy and Copy number
  cnvs$sample = substring(cnvs$sampleID,1,15)
  cnvs$key =paste(cnvs$sample, cnvs$entrez_19549, sep=".")
  
  
  #This will unfortunately remove all non-exonic regions from cnv file (so no Copy_number = 0 anymore)
  df_cnv0 = cnvs %>% subset(Copy_number==0)
  if (FALSE){
    df_cnv1 = cnvs %>% subset(Copy_number==1)
    df_cnv1 = df_cnv1 %>% subset(key%in%total_table$key)
    
    k = which(cnvs$key%in% total_table$key)
    k_exonic = which(cnvs$key%in% total_table2$key)
    cnvs$group[-k] = "cnvs_map_no_muts"
    cnvs$group[setdiff(k,k_exonic)] = "cnvs_map_nonexonic_muts"
    cnvs$group[k_exonic] = "cnvs_map_exonic_muts"
    
    df = cnvs %>% group_by(Copy_number, group) %>% summarise(n_gene = n_distinct(gname),
                                                             n_alts = n())
    df$log_n_alts = log(df$n_alts)
    df$CODE = i
    
    
    p=ggplot(df, aes(x=Copy_number,y=n_gene, fill=group))+
      geom_bar(stat="identity")+
      scale_fill_grey(start = 0, end = .9)+
      xlab('Copy Number') + ylab('number of genes')+
      ggtitle(sprintf("%s (n=%s,n_map_exonic=%s)", i, length(unique(cnvs$sampleID)), length(unique(cnvs$sampleID)))) +
      theme_boss()
    ggsave(filename=paste0("~/Rosalind/Plots/CN_map_muts/",i,"_CNmappingmuts.png"), 
           plot=p,width = 20, height = 12, dpi = 100)
  }
  write.table(df_cn0, "/mnt/lustre/users/k1801782/dataset/CNV/CN0_all.txt", col.names=F, row.names = F,append = T)
  #  write.table(df_cn1, "/mnt/lustre/users/k1801782/dataset/CNV/CN1_map_exonic_nonexonic.txt", col.names=F, append = T)
}

library(ggplot2)
library(reshape2)
df = cnvs %>% group_by(Copy_number) %>% summarise(n_gene = n_distinct(gname),
                                                  n_alts = n())
df$log_n_alts = log(df$n_alts)
p=ggplot(df, aes(x=Copy_number,y=n_alts))+
  geom_bar(stat="identity")+
  xlab('Copy Number') + ylab('Number of alterations')+
  #  geom_text(aes(label=n_alts),position=position_dodge(width=0.9), vjust=-0.25)+
  theme_boss()
p
ggsave(filename=paste0("~/Desktop/OAC_alts2_CNmappingmuts.png"), 
       plot=p,width = 20, height = 12, dpi = 100)

p=ggplot(df, aes(x=Copy_number,y=n_gene))+
  geom_bar(stat="identity")+
  xlab('Copy Number') + ylab('Number of genes')+
  geom_text(aes(label=n_gene),position=position_dodge(width=0.9), vjust=-0.25)+
  theme_boss()
p
ggsave(filename=paste0("~/Desktop/OAC_uniquegene_CNmappingmuts.png"), 
       plot=p,width = 20, height = 12, dpi = 100)

df = cnvs %>% 
  group_by(Copy_number, group) %>% 
  summarise(n_gene = n_distinct(gname),n_alts = n())

#df = data.frame(table(cnvs$group,cnvs$Copy_number))
#colnames(df) = c("group","Copy_number","n_alterations")
df$log_n_alts = log(df$n_alts)

ggplot(df, aes(x=Copy_number,y=log_n_alts))+
  geom_bar(stat="identity")+
  theme_boss()


p=ggplot(df, aes(x=Copy_number,y=n_alts, fill=group))+
  geom_bar(stat="identity")+
  scale_fill_grey(start = 0, end = .9)+
  xlab('Copy Number') + ylab('number of alterations')+
  #  geom_text(aes(label=n_gene),position="stack", vjust=-0.25)+
  #  ggtitle(sprintf("%s (n=%s,n_map_exonic=%s)", i, length(unique(cnvs$sampleID)), length(unique(cnvs$sampleID)))) +
  theme_boss()
p
ggsave(filename=paste0("~/Desktop/OAC_CNmappingmuts.png"), 
       plot=p,width = 20, height = 12, dpi = 100)


p=ggplot(df, aes(x=Copy_number,y=log_n_alts, fill=group))+
  geom_bar(stat="identity")+
  scale_fill_grey(start = 0, end = .9)+
  xlab('Copy Number') + ylab('number of genes')+
  geom_text(aes(label=n_gene),position="stack", vjust=-0.25)+
  #  ggtitle(sprintf("%s (n=%s,n_map_exonic=%s)", i, length(unique(cnvs$sampleID)), length(unique(cnvs$sampleID)))) +
  theme_boss()
p
ggsave(filename=paste0("~/Desktop/OAC_CNmappingmuts.png"), 
       plot=p,width = 20, height = 12, dpi = 100)


summary_mlinput = data.frame()
summary_damaging = data.frame()
for (i in code){
  mlinput = readRDS(sprintf("~/Rosalind/OncodriveClust/result_30112018/34cancer_types/%s_total_table_ploidy.rds",i))
  damaging = readRDS(sprintf("~/Rosalind/OncodriveClust/result_30112018/34cancer_types/%s_damaging_snvindelgofCNV_ploidy.rds",i))
  tmp = mlinput %>% group_by(sample) %>% summarise(CNVgain = sum(CNV_type=="Gain",na.rm = T), CNVLoss= sum(CNV_type=="Loss",na.rm = T))
  tmp$CODE =i
  tmp$label = paste0(i, " (",length(tmp$sample),")")
  summary_mlinput = rbind(summary_mlinput,tmp)
  
  tmp = damaging %>% group_by(sample) %>% summarise(CNVgain = sum(CNVGain), CNVLoss= sum(CNVLoss))
  tmp$CODE=i
  tmp$label = paste0(i, " (",length(tmp$sample),")")
  summary_damaging = rbind(summary_damaging,tmp)
}

tmp=melt(summary_damaging[,-3], id=c("sample","CODE","label"))
tmp = tmp%>%group_by(label)%>%mutate(mean=median(value))

p=ggplot(tmp, aes(reorder(label,mean),value))+
  geom_boxplot(aes(fill = variable)) + 
  xlab('Cancer types') + ylab('number of CNVGain per sample')+
  stat_summary(aes(label=round(..y..)), fun.y=min, geom="label", size=4, position = position_nudge(y=-10)) +
  stat_summary(aes(label=round(..y..)), fun.y=median, geom="label", size=4) +
  stat_summary(aes(label=round(..y..)), fun.y=max, geom="label", size=4) +
  #  ylim(c(0,500))+ #for CNVGain
  #  ylim(c(0,200))+ #for CNVLoss
  theme_boss()

p

ggsave(filename=paste0("~/Desktop/CNVGain_all_damaging.png"), 
       plot=p,width = 20, height = 12, dpi = 100)


##################################################
cnv_path="/mnt/lustre/users/k1801782/dataset/CNV/CNV_intersectBed_filtered"
muts_path ="~/Rosalind/OncodriveClust/result_30112018/34cancer_types"
cnv_path="~/Rosalind/CNV/CN_intersectBed_filtered"
muts_path ="~/Rosalind/OncodriveClust/result_30112018/34cancer_types"

library(dplyr)
library(reshape2)

summary_table_damaging = data.frame()
for (i in code){
  round = readRDS(sprintf("%s/%s_damaging_snvindelgofCNV_ploidy_CNrounded_19charID.rds",muts_path,i))
  round_not = readRDS(sprintf("%s/%s_damaging_snvindelgofCNV_ploidy_notround_19charID.rds",muts_path,i))
  tmp = merge(round %>% group_by(sample) %>% summarise(CNVGain_CNround = sum(CNVGain),CNVLoss_CNround = sum(CNVLoss)),
              round_not %>% group_by(sample) %>% summarise(CNVGain_CNdecimal = sum(CNVGain),CNVLoss_CNdecimal = sum(CNVLoss)),
              by="sample") %>% mutate(CODE =i)
  
  #  round_not = readRDS(sprintf("%s/%s_gainloss_2ploidyi_notround_19charID.rds",cnv_path,i))
  #  round = readRDS(sprintf("%s/%s_gainloss_2ploidy_CNrounded_19charID.rds",cnv_path,i))
  #  tmp = merge(round %>% group_by(sampleID) %>% summarise(CNVGain_CNround = sum(CNV_type=="Gain", na.rm = T),CNVLoss_CNround = sum(CNV_type=="Loss", na.rm = T)),
  #              round_not %>% group_by(sampleID) %>% summarise(CNVGain_CNdecimal = sum(CNV_type=="Gain", na.rm = T),CNVLoss_CNdecimal = sum(CNV_type=="Loss", na.rm = T)),
  #              by="sampleID") %>% mutate(CODE =i)
  summary_table_damaging = rbind(summary_table_damaging,tmp)
  cat(i," done","\n")
}
colnames(summary_table) = c("sampleID","CNVGain_CNround","CNVLoss_CNround","CNVGain_CNdecimal","CNVLoss_CNdecimal","CODE","n","label")

saveRDS(summary_table,sprintf("%s/summary_table_CNVGainLoss.rds",cnv_path))

summary_table_damaging=summary_table_damaging %>% group_by(CODE) %>% mutate(n=n())
summary_table_damaging$label = paste0(summary_table_damaging$CODE, " (",summary_table_damaging$n, ")")

plotthis = melt(summary_table_damaging[,-c(4,5)], id=c("sample","CODE","n","label"))
plotthis$logCNVGainLoss = log(plotthis$value)
plotthis = plotthis %>% group_by(CODE,variable) %>% mutate(mean = mean(value))


plotthis$label = factor(plotthis$label,
                        levels=(plotthis[which(plotthis$variable=="CNVGain_CNround"),] %>%
                                  subset(is.finite(logCNVGainLoss)) %>%
                                  group_by(label,variable) %>% summarise(order = median(logCNVGainLoss)) %>% 
                                  arrange(order))$label,
                        ordered = T)
p=ggplot(plotthis, aes(label,logCNVGainLoss))+
  geom_boxplot(aes(fill = variable)) + 
  xlab('Cancer types') + ylab('log of number of CNVGain/Loss per sample')+
  #  geom_text(aes(x=label, y=order, label=mean))+
  #  stat_summary(aes(label=round(..y..)), fun.y=median, geom="label", size=4) +#, position = position_nudge(y=-10)) +
  #  facet_grid(variable~.)+
  #  ylim(c(0,500))+ #for CNVGain
  #  ylim(c(0,200))+ #for CNVLoss
  theme_boss_xtilted()
p

ggplot(summary_table, aes(CNVLoss_CNround, CNVLoss_CNdecimal))+
  geom_point()+
  theme_boss()

#PANCAN
plotthis = melt(summary_table, id=c("sampleID","CODE","n","label"))
plotthis$logCNVGainLoss = log(plotthis$value)
ggplot(plotthis,aes(variable, logCNVGainLoss))+
  geom_boxplot()+
  xlab('') + ylab('log of number of CNVGain/CNVLoss per sample')+
  stat_summary(aes(label=round(..y..,1)), fun.y=min, geom="label", size=4) +#, position = position_nudge(y=-10)) +
  stat_summary(aes(label=round(..y..,1)), fun.y=median, geom="label", size=4) +
  stat_summary(aes(label=round(..y..,1)), fun.y=max, geom="label", size=4) +
  theme_boss()

ggplot(ploidy_8611, aes(reorder(CODE,ploidy),ploidy))+
  geom_boxplot()+
  theme_boss_xtilted()

###################################################################
#### Plotting of damaging mutations across cohorts based on snvindel gof and cnvs (2*ploidy)
###################################################################

# All exonic mutations
muts1 = snvindel
muts1_exon = muts1 %>% subset(Func.refGene!="." &
                                !grepl("downstream", muts1$Func.refGene) &
                                !grepl("upstream", muts1$Func.refGene) &
                                !grepl("intergenic", muts1$Func.refGene) &
                                !grepl("ncRNA", muts1$Func.refGene) &
                                !grepl("intronic", muts1$Func.refGene) &    
                                !grepl("UTR", muts1$Func.refGene))

bailey_1 = left_join(muts1_exon %>% group_by(Sample) %>% dplyr::summarise(mutations=n()), 
                     muts1_exon %>% group_by(Sample) %>% dplyr::summarise(mutated_genes=n_distinct(Gene.refGene)))
bailey_1$Type = "all_snvindel"
bailey_1$Sample = substring(bailey_1$Sample,1,19)

# All exonic mutations in ncg
muts1_exon_ncg = muts1_exon %>% subset(!is.na(entrez_19549))
bailey_2 = left_join(muts1_exon_ncg %>% group_by(Sample) %>% dplyr::summarise(mutations=n()), 
                     muts1_exon_ncg %>% group_by(Sample) %>% dplyr::summarise(mutated_genes=n_distinct(Gene.refGene)))
bailey_2$Type = "all_snvindel_ncg6"
bailey_2$Sample = substring(bailey_2$Sample,1,19)

# All damaging exonic mutation based on snvindel gof and cnvs
muts_path ="~/Trang/OncodriveClust/result_30112018/34cancer_types"
total_table_snvindelgofcnv = data.frame()
for (i in code){
  tmp = readRDS(sprintf("%s/%s_total_table_ploidy_CNrounded_19charID.rds",muts_path,i))
  #  tmp=tmp[!is.na(tmp$no_ALL_muts),]
  total_table_snvindelgofcnv=rbind(total_table_snvindelgofcnv,tmp)
  cat(i,"has",sum(is.na(tmp$no_TRUNC_muts)),"CNV alterations without mutations","\n")
}

#Counting damaging mutations
total_table_snvindelgofcnv = getMLinput(total_table_snvindelgofcnv)

total_table_snvindelgofcnv = total_table_snvindelgofcnv %>% 
  mutate(damaging_snvindel_cnv = no_TRUNC_muts + no_NTDam_muts + no_GOF_muts + (CNVGain==1) + (Copy_number==0) + ((Copy_number==1) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)))

muts1_exon_ncg_onco_cnv_damaging = total_table_snvindelgofcnv  %>% subset(no_TRUNC_muts!=0 |
                                                                            no_NTDam_muts!=0 |
                                                                            no_GOF_muts!=0 |
                                                                            (CNVGain==1) | 
                                                                            (Copy_number==0) |
                                                                            ((Copy_number==1) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0))
                                                                          # BND!=0 |9 
                                                                          # INS!=0 | 
                                                                          # INV!=0)
)

bailey_3 = left_join(muts1_exon_ncg_onco_cnv_damaging %>% group_by(sample) %>% dplyr::summarise(mutations=sum(no_TRUNC_muts + no_NTDam_muts + no_GOF_muts + 
                                                                                                                (CNVGain==1) + (Copy_number==0) + ((Copy_number==1) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0)))),
                     muts1_exon_ncg_onco_cnv_damaging %>% group_by(sample) %>% dplyr::summarise(mutated_genes=n_distinct(entrez_19549)))

colnames(bailey_3)[1] = "Sample"
bailey_3$Type = "damaging_snvindelgofcnv"
plotthis=rbind(bailey_1,rbind(bailey_2,bailey_3))
plotthis=merge(plotthis,CODE, by.x ="Sample", by.y = "sampleID" )
plotthis2=melt(plotthis[,-c(5)],id=c("Sample","CODE","Type"))

plotthis$logMutations = log(plotthis$mutations)
plotthis$logMutatedGenes = log(plotthis$mutated_genes)


library(ggplot2)
df = plotthis[which(plotthis$Type=="damaging_snvindelgofcnv"),]
df = df %>% group_by(CODE) %>% mutate(n=n())
df$label = paste0(df$CODE, " (",df$n, ")")
df = df %>% group_by(CODE) %>% mutate(order=median(mutations))

p <- ggplot(df, aes(reorder(label, order), logMutatedGenes)) + 
  geom_boxplot(aes(fill = Type)) + 
  xlab('Cancer types') + ylab('log of number of damaged genes per sample')+
  stat_summary(aes(label=round(..y..)), fun.y=median, geom="text", size=4, vjust = -0.5) +
  stat_summary(aes(label=round(..y..)), fun.y=min, geom="text", size=4, vjust = 1) +
  stat_summary(aes(label=round(..y..)), fun.y=max, geom="text", size=4, vjust = -0.5)+
  theme_boss_xtilted()
p

ggsave(filename=paste0("~/Rosalind/Plots/damaging_snvindel_Onco_cnv/LogMutatedGenes_damaging_snvindelgofcnv_34cancers.png"), 
       plot=p,width = 16, height = 9, dpi = 100)

df$label = "PANCAN (8577)"
p <- ggplot(df, aes(label, logMutatedGenes)) + 
  geom_boxplot(aes(fill = Type)) + 
  xlab('Cancer types') + ylab('log of number of damaged genes per sample')+
  stat_summary(aes(label=round(..y..)), fun.y=median, geom="text", size=4, vjust = -0.5) +
  stat_summary(aes(label=round(..y..)), fun.y=min, geom="text", size=4, vjust = 1) +
  stat_summary(aes(label=round(..y..)), fun.y=max, geom="text", size=4, vjust = -0.5)+
  theme_boss_xtilted()
p

ggsave(filename=paste0("~/Rosalind/Plots/damaging_snvindel_Onco_cnv/LogMutatedGenes_damaging_snvindelgofcnv_pancan.png"), 
       plot=p,width = 4, height = 12, dpi = 100)


for (i in code) {
  tmp = plotthis2[which(plotthis2$CODE==i),]
  tmp_sum = tmp %>% group_by(Type) %>% summarise(n_all = length(unique(Sample)))
  p=ggplot(tmp, aes(Type, value)) + 
    geom_boxplot(aes(fill = Type)) + 
    xlab('Cancer types') + ylab('number per sample')+
    stat_summary(aes(label=..y..), fun.y=median, geom="label", size=4) +
    stat_summary(aes(label=..y..), fun.y=min, geom="label", size=4) +
    stat_summary(aes(label=..y..), fun.y=max, geom="label", size=4) +
    #    ylim(c(0,600))+
    facet_grid(.~variable)+
    xlab("")+
    ggtitle(paste0(i, sprintf(" (n_all=%s,n_ncg6=%s,n_damaging=%s)", tmp_sum$n_all[1], tmp_sum$n_all[2], tmp_sum$n_all[3]))) +
    theme(axis.text.x=element_text(angle=30,hjust=0.75,vjust=0.75, size = 12, face = "bold"))
  ggsave(filename=paste0("~/Rosalind/Plots/damaging_snvindel_Onco_cnv/",i,"_withOutliers.png"), 
         plot=p,width = 16, height = 9, dpi = 100)
}



bailey_3 = left_join(muts1_exon_ncg_onco_damaging %>% group_by(sample) %>% dplyr::summarise(mutations=sum(no_TRUNC_muts+no_NTDam_muts+no_GOF_muts+CNVGain+CNVLoss)),
                     muts1_exon_ncg_onco_damaging %>% group_by(sample) %>% dplyr::summarise(mutated_genes=n_distinct(entrez_19549)))




###################################################################
#### getMLinput funciton
###################################################################

getMLinput <- function(df){
  df = df %>% mutate(Cancer_type="OAC") %>% rename(Entrez=entrez_19549, Copy_number=Copy_number, CNV_type=CNV_type) #%>% select(-symbol_19549, -in_19549, -ploidy, -CNV_entries, -key)
  ## Replace numbers with names here
  ## We assume every gene with no mutation data that it's not mutated
  message("Fixing mutations...")
  df[,c("no_ALL_muts", "no_NSI_muts", "no_TRUNC_muts", 
        "no_NTDam_muts",
        "no_GOF_muts")][is.na(df[,c("no_ALL_muts", "no_NSI_muts", "no_TRUNC_muts",
                                    "no_NTDam_muts",
                                    "no_GOF_muts")])] <- 0
  
  message("Fixing CNVs...")
  ## Copy number (where copy number is NA, put copy number equal to 2)
  ## I integrated copy number data by selecting segment mean >|0.3| therefore I took only gains and losses
  ## But at the same time the unique number of genes in the CNV data is quite high, therefore whatever is left with NA is probably 2
  df$Copy_number[is.na(df$Copy_number)] <- 2
  ## Convert categorical features to multiple factors
  message("Performing cleaning of categorical variables...")
  ## CNV type
  df <- df %>% 
    mutate(CNVGain=ifelse(is.na(CNV_type), 0, ifelse(CNV_type=="Gain",1, 0)), 
           CNVLoss=ifelse(Copy_number==0 | Copy_number==1, 1, 0)) %>%
    select(-CNV_type)
  ## Before you add (change NAs in these columns to 0)
  #df[,c("High", "Low", "Medium", "NotExpressed")][is.na(df[,c("High", "Low", "Medium", "NotExpressed")])] <- 0
  #df <- df %>% ungroup %>% mutate(tot.tissues=High+Low+Medium)
  # df <- data.frame(df)
  return(df)
}

#for OAC
OAC = left_join(OAC_total_table_1ploidy_CNrounded_19charID %>% group_by(sample) %>% dplyr::summarise(CNVGain_1ploidy=sum(CNV_type=="Gain",na.rm = T),
                                                                                                     CNVLoss_1ploidy=sum(CNV_type=="Loss",na.rm = T)),
                OAC_total_table_ploidy_CNrounded_19charID %>% group_by(sample) %>% dplyr::summarise(CNVGain_2ploidy=sum(CNV_type=="Gain",na.rm = T), 
                                                                                                    CNVLoss_2ploidy=sum(CNV_type=="Loss",na.rm = T)))

OAC = left_join(OAC_gainloss_ploidy_CNrounded_19charID %>% group_by(sample) %>% dplyr::summarise(CNVGain_1ploidy=sum(CNV_type=="Gain",na.rm = T),
                                                                                                 CNVLoss_1ploidy=sum(CNV_type=="Loss",na.rm = T)),
                OAC_gainloss_2ploidy_CNrounded_19charID %>% group_by(sample) %>% dplyr::summarise(CNVGain_2ploidy=sum(CNV_type=="Gain",na.rm = T), 
                                                                                                  CNVLoss_2ploidy=sum(CNV_type=="Loss",na.rm = T)))

plotthis = melt(OAC[,-c(3,5)])
plotthis$logValue=log(plotthis$value)
p=ggplot(plotthis,aes(variable, logValue))+
  geom_boxplot()+
  stat_summary(aes(label=round(..y..)), fun.y=min, geom="label", size=4) +
  stat_summary(aes(label=round(..y..)), fun.y=max, geom="label", size=4) +
  stat_summary(aes(label=round(..y..)), fun.y=median, geom="label", size=4) +
  theme_boss()
p

p=ggplot(plotthis,aes(variable, value))+
  geom_boxplot()+
  stat_summary(aes(label=..y..), fun.y=min, geom="label", size=4) +
  stat_summary(aes(label=..y..), fun.y=max, geom="label", size=4) +
  stat_summary(aes(label=..y..), fun.y=median, geom="label", size=4) +
  theme_boss()
p

ggsave(filename=paste0("~/Rosalind/Plots/CNVGain_1ploidyvs2ploidy.png"), 
       plot=p,width = 12, height = 9, dpi = 100)


###################################################################
#### CNV Gain/Loss landscape
###################################################################
# gain = CNV>=2*ploidy
# loss <2
library(dplyr)
# cancertype = unique(muts1_exon_ncg_onco_damaging$CODE)
# path = '~/Trang/CNV/ASCATforTCGA/CNVGainLoss'
cancertype = c("GBM","OV","LUSC","TGCT","ESCA", "SARC","LIHC","LUAD","BRCA","COAD","UCEC","KIRC","STAD","CESC","HNSC","SKCM","LGG","BLCA","PRAD","THCA",
               "KIRP", "MESO", "READ","DLBC", "ACC","LAML", "UCS","PAAD", "THYM", "CHOL", "KICH", "PCPG", "UVM")
path = '/mnt/lustre/users/k1801782/dataset/CNV/ASCATforTCGA/CNVGainLoss'
path = '~/Trang/CNV/ASCATforTCGA/CNVGainLoss'
CNVGainLoss = NULL
for (code in cancertype){
  tmp = readRDS(sprintf('%s/%s.rds', path, code))
  tmp = tmp %>% group_by(sample) %>% summarise(CNVGain = sum(CNV_type_corrected=="Gain",na.rm=T),
                                               CNVLoss = sum(CNV_type_corrected=="Loss",na.rm=T),
                                               CN0 = sum(Total_CN==0), CN1=sum(Total_CN==1))
  CNVGainLoss = rbind(CNVGainLoss,tmp)
  print(paste0("done ",code))
}
saveRDS(CNVGainLoss,'~/Trang/CNV/ASCATforTCGA/CNVGainLoss/countCNV.rds')
k = which(ascat_tcga_ploidy$patient %in% oac)
ascat_tcga_ploidy[,"CODE"] = as.character(ascat_tcga_ploidy$cancer_type)
ascat_tcga_ploidy[k,"CODE"] = "OAC"
ascat_tcga_ploidy[which(ascat_tcga_ploidy$CODE=="ESCA"),"CODE"] = "OSCC"

# local
countCNV <- readRDS("~/Trang/CNV/ASCATforTCGA/CNVGainLoss/countCNV.rds")
colnames(countCNV)[1] = "patient"
countCNV$patient = substring(countCNV$patient,1,12)
countCNV = merge(countCNV, ascat_tcga_ploidy %>% select(patient, CODE, pass, rep))


library(ggplot2)
library(reshape2)
plotthis = melt(countCNV %>% select(patient, CODE,CNVGain, CN0, CN1))
plotthis = plotthis %>% group_by(CODE) %>% mutate(n=n_distinct(patient))
plotthis$label = paste0(plotthis$CODE, " (",plotthis$n, ")")
plotthis$Analysis_Type="individual cancers"
tmp = rbind(plotthis, plotthis%>%mutate(label="PANCAN (9873)", Analysis_Type="PANCAN"))
tmp = tmp[which(tmp$variable=="CN1"),]
stat = plotthis[which(plotthis$variable=="CN1"),] %>% group_by(label) %>% summarise(med = median(value),maxim = max(value),minim = min(value)) %>% arrange(med) #arrange(desc(med))
stat = rbind(stat, plotthis[which(plotthis$variable=="CN1"),] %>%mutate(label="PANCAN (9873)") %>% group_by(label) %>% summarise(med = median(value),maxim = max(value),minim = min(value)))

order = c(levels(reorder(stat$label, stat$med))[-18],"PANCAN (9873)")
tmp$label = factor(tmp$label,order)

p=ggplot(tmp %>% mutate(value = replace(value, value==0, 0.9999)), aes(x=label, y=value)) +
  geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(Analysis_Type)), alpha=0.1) +
  #scale_color_manual(values = c("dodgerblue1","deepskyblue2","blue3")) +
  geom_boxplot(alpha=0.2)+#, aes(fill=factor(Analysis_Type)))+
  scale_y_log10() + 
  theme_boss_xtilted()+
  #theme_boss(base_size = 10)+
  #facet_wrap(~category, scales = 'free')+
  annotate("label", x=1:nrow(stat), y=replace(stat$med, stat$med==0, 0.9999), label=stat$med,size=6) +
  annotate("label", x=1:nrow(stat), y=replace(stat$maxim, stat$maxim==0, 0.9999), label=stat$maxim,size=6) +
  annotate("label", x=1:nrow(stat), y=replace(stat$minim, stat$minim==0, 0.9999), label=stat$minim, size=6) +
  ylab("Number of heterozygous CN loss per sample") + xlab("Category")
p

ggsave(filename=paste0("~/Trang/Plots/CN1_all_ascat.png"), 
       plot=p,width = 25, height = 15, dpi = 300)


countCNV = countCNV %>% rename(CN_gain = CNVGain, heterozygous_loss=CN1, homozygous_loss=CN0)
tmp = melt(countCNV[,-c(3,6,7,8)])
stat = tmp %>% group_by(variable) %>% summarise(med = median(value),maxim = max(value),minim = min(value))
#tmp$variable = factor(tmp$variable, levels = c("CN_gain","heterozygous_loss", "homozygous_loss"))

p=ggplot(tmp %>% mutate(value = replace(value, value==0, 0.9999)), aes(x=variable, y=value)) +
  geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(variable)), alpha=0.1) +
  scale_color_manual(values = c("olivedrab2","tomato3","red")) +
  geom_boxplot(alpha=0.2)+#, aes(fill=factor(Analysis_Type)))+
  scale_y_log10() + 
  theme_boss(base_size = 30)+
  annotate("label", x=1:nrow(stat), y=replace(stat$med, stat$med==0, 0.9999), label=stat$med,size=10) +
  ylab("Number of copy number alterations per sample") + xlab("Category")
p

ggsave(filename=paste0("~/Trang/Plots/CNVGainLoss_all_ascat.png"), 
       plot=p,width = 25, height = 15, dpi = 300)

###################################################################
#### Damaging mutations by snvs, indels, gof, and cnvs
###################################################################
cancertype = unique(muts1_exon_ncg_onco_damaging$CODE) # unique(ascat_tcga_ploidy$cancer_type)
#path = '~/Trang/CNV/Total_table_26022019'
total_table_pancanrun = NULL
for (code in cancertype){ #
  tmp = readRDS(sprintf('~/Rosalind/Annotation_individual_cancers/%s/total_table_pancan.rds', code))#readRDS(sprintf('%s/%s.rds', path, code))
  tmp$CODE=code
  total_table_pancanrun = rbind(total_table_pancanrun,tmp)
  print(paste0("done ",code))
}
# Removing duplicate OAC samples that were wrongly labeled as OSCC
# k = intersect(which(total_table$Cancer_type=="OSCC"),which(total_table$sample %in% oac))
# total_table_filtered = total_table[-k,]
saveRDS(total_table,'~/Trang/CNV/Total_table_26022019/PANCAN_unfiltered.rds')

saveRDS(total_table_pancanrun,'~/Rosalind/Annotation_individual_cancers/Allcancertypes_pancan_annotation_unfiltered.rds')
saveRDS(total_table_individualrun,'~/Rosalind/Annotation_individual_cancers/Allcancertypes_individual_annotation_unfiltered.rds')




PassSample = ascat_tcga_ploidy$name[which((ascat_tcga_ploidy$pass==TRUE)&(ascat_tcga_ploidy$rep==TRUE))]
PassSample_mc3 = intersect(substring(PassSample,1,12), substring(unique(muts1_exon_ncg$Sample),1,12))

damaging_muts_cnvs = function(total_table, PassSample_mc3,ascat_tcga_ploidy){
  library(dplyr)
  total_table = total_table %>% filter((sample %in% PassSample_mc3))
  # Fill in ploidy for the (21) samples that has no cnv alterations but has mutations 
  #ascat_tcga_ploidy = ascat_tcga_ploidy %>% filter(pass==TRUE)
  #total_table = total_table %>% left_join(ascat_tcga_ploidy %>% rename(sample=patient) %>% select(sample, ploidy))
  
  #a = total_table_filtered2 %>% group_by(Cancer_type) %>% summarise(n = n_distinct(sample))
  #saveRDS(total_table_filtered2,'~/Trang/CNV/Total_table_13022019/PANCAN_bailey.rds')
  damaging_alterations = total_table %>% subset(no_TRUNC_muts!=0 |
                                                  no_NTDam_muts!=0 |
                                                  no_GOF_muts!=0 |
                                                  (CNVGain==1) | 
                                                  (Copy_number==0) |
                                                  ((Copy_number==1) & (no_TRUNC_muts!=0 | no_NTDam_muts!=0))
                                                   # BND!=0 | 
                                                   # INS!=0 | 
                                                   # INV!=0)
  )
  #damaging_alterations = damaging_alterations[order(damaging_alterations$key),]
  return(damaging_alterations)
}


muts1_exon_ncg_onco_cnv_damaging = damaging_muts_cnvs(total_table_pancanrun,PassSample_mc3,ascat_tcga_ploidy)
muts1_exon_ncg_onco_cnv_damaging = muts1_exon_ncg_onco_cnv_damaging %>% filter(!sample %in% c('TCGA-KN-8419','TCGA-DW-5561','TCGA-SR-A6MT')) #TCGA-SR-A6MT,TCGA-B1-7332
saveRDS(muts1_exon_ncg_onco_cnv_damaging,'~/Rosalind/Annotation_individual_cancers/Allcancertypes_pancan_annotation_damaging.rds')

muts1_exon_ncg_onco_cnv_damaging_ind = damaging_muts_cnvs(total_table_individualrun,PassSample_mc3,ascat_tcga_ploidy)
muts1_exon_ncg_onco_cnv_damaging_ind = muts1_exon_ncg_onco_cnv_damaging_ind %>% filter(!sample %in% c('TCGA-KN-8419','TCGA-DW-5561','TCGA-SR-A6MT'))
saveRDS(muts1_exon_ncg_onco_cnv_damaging_ind,'~/Rosalind/Annotation_individual_cancers/Allcancertypes_individual_annotation_damaging.rds')



stat_dam = left_join(muts1_exon_ncg_onco_cnv_damaging %>% group_by(sample, Cancer_type) %>% dplyr::summarise(damaging_alterations=sum(no_TRUNC_muts+no_NTDam_muts+no_GOF_muts+CNVGain+CNVLoss)),
                     muts1_exon_ncg_onco_cnv_damaging %>% group_by(sample, Cancer_type) %>% dplyr::summarise(damaged_genes=n_distinct(Entrez)))

stat_dam= stat_dam %>% group_by(Cancer_type)%>% mutate(n=n_distinct(sample))
stat_dam$label = paste0(stat_dam$Cancer_type, " (",stat_dam$n, ")")

# Plotting damaging alterations and damaged genes distribution in the order of median
plotthis = melt(stat_dam[,-c(3,5)])  
plotthis$Analysis_Type="individual cancers"
tmp = rbind(plotthis, plotthis%>%mutate(label="PANCAN (8113)", Analysis_Type="PANCAN"))
stat = plotthis %>% group_by(label) %>% summarise(med = median(value),maxim = max(value),minim = min(value)) %>% arrange(med) #arrange(desc(med))
stat = rbind(stat, plotthis %>% mutate(label="PANCAN (8113)") %>% group_by(label) %>% summarise(med = median(value),maxim = max(value),minim = min(value)))

order = c(levels(reorder(stat$label, stat$med))[-17],"PANCAN (8113)") 
tmp$label = factor(tmp$label,order)

p=ggplot(tmp %>% mutate(value = replace(value, value==0, 0.9999)), aes(x=label, y=value)) +
  geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(Analysis_Type)), alpha=0.1) +
  geom_boxplot(alpha=0.2)+#, aes(fill=factor(Analysis_Type)))+
  scale_y_log10() + 
  theme_boss_xtilted()+
  #theme_boss(base_size = 10)+
  #facet_wrap(~category, scales = 'free')+
  annotate("label", x=1:nrow(stat), y=replace(stat$med, stat$med==0, 0.9999), label=stat$med,size=6) +
  annotate("label", x=1:nrow(stat), y=replace(stat$maxim, stat$maxim==0, 0.9999), label=stat$maxim,size=6) +
  annotate("label", x=1:nrow(stat), y=replace(stat$minim, stat$minim==0, 0.9999), label=stat$minim,size=6) +
  ylab("Number of damaged genes per sample") + xlab("Cohort")
p


pdf(file = "~/Rosalind/Plots/damaging_snvindel_Onco_cnv/DamagedGenes_all_cohorts.pdf", width = 25, height = 15)
p
dev.off()
ggsave(filename=paste0("~/Trang/Plots/damaging_snvindel_Onco_cnv/DamagingMutations_all_cohorts.png"), 
       plot=p,width = 25, height = 15, dpi = 300)





tmp_damjan = THCA_damjan %>% mutate(key2 = paste0(sample,Entrez)) 
tmp_damjan = damaging_muts_cnvs(tmp_damjan,PassSample_mc3,ascat_tcga_ploidy)
tmp_damjan = tmp_damjan[order(tmp_damjan$key2),]
THCA_damjan = left_join(tmp_damjan %>% group_by(sample, Cancer_type) %>% dplyr::summarise(damaging_alterations=sum(no_TRUNC_muts+no_NTDam_muts+no_GOF_muts+CNVGain+CNVLoss)),
                       tmp_damjan %>% group_by(sample, Cancer_type) %>% dplyr::summarise(damaged_genes=n_distinct(Entrez)))

tmp_trang = THCA_trang %>% mutate(key2 = paste0(sample,Entrez))
tmp_trang = damaging_muts_cnvs(tmp_trang, PassSample_mc3, ascat_tcga_ploidy)
tmp_trang = tmp_trang[order(tmp_trang$key2),]
THCA_trang = left_join(tmp_trang %>% group_by(sample, Cancer_type) %>% dplyr::summarise(damaging_alterations=sum(no_TRUNC_muts+no_NTDam_muts+no_GOF_muts+CNVGain+CNVLoss)),
                       tmp_trang %>% group_by(sample, Cancer_type) %>% dplyr::summarise(damaged_genes=n_distinct(Entrez)))

summary(THCA_damjan)
summary(THCA_trang)

tmp = rbind(THCA_damjan,PCPG_damjan,LUSC_damjan)
tmp$Analyst="Damjan"
tmp2 = rbind(THCA_trang,PCPG_trang, LUSC_trang)
tmp2$Analyst ="Trang"
tmp=rbind(tmp,tmp2)
plotthis = melt(tmp) 
stat = plotthis %>% group_by(label) %>% summarise(med = median(value),maxim = max(value),minim = min(value)) %>% arrange(med) #arrange(desc(med))

ggplot(plotthis %>% mutate(value = replace(value, value==0, 0.9999)), aes(x=Analyst, y=value)) +
  geom_boxplot(alpha=0.2)+
  scale_y_log10() + 
  theme_boss(base_size = 10)+
  facet_grid(variable~Cancer_type)+
  ylab("Number of damaged genes/damaging alterations per sample") + xlab("Cohort")




Genes_m025 = list()
for (code in names(Genes_f025)){
  tmp = readRDS(sprintf('~/Trang/CNV/ASCATforTCGA/CNV_intersectBed_max0.25/%s.rds',code))
  Genes_m025[[code]] =unique(tmp$symbol_19549)
}

higher0.25=unique(unlist(Genes_f025))
lower0.25= unique(unlist(Genes_m025))

# summarise
df = muts1_exon_ncg_onco_cnv_damaging %>% group_by(CODE) %>%
  summarise(n_sample=n_distinct(sample),
            n_genes_avg = n_distinct(symbol_19549)/n_distinct(sample),
            n_alts_avg = n()/n_distinct(sample),
            NSI_sampleavg = sum(no_NSI_muts)/n_distinct(sample),
            TRUNC_sampleavg = sum(no_TRUNC_muts)/n_distinct(sample),
            NTDam_sampleavg = sum(no_NTDam_muts)/n_distinct(sample),
            GOF_sampleavg = sum(no_GOF_muts)/n_distinct(sample),
            CNVGain_sampleavg = sum(CNVGain)/n_distinct(sample),
            HomozygousCNVLoss = sum(Copy_number==0)/n_distinct(sample),
            HeterozygousCNVLoss = sum(Copy_number==1)/n_distinct(sample))
df[,3:11]=round(df[,3:11])

df_0 = muts1_exon_ncg_onco_cnv_damaging %>% group_by(Cancer_type) %>%
  summarise(n_sample=n_distinct(sample),
            n_dam_genes = n_distinct(symbol_19549),
            n_dam_alts = sum(no_TRUNC_muts+no_NTDam_muts+no_GOF_muts+CNVGain+CNVLoss),
#            all_alt = sum(no_ALL_muts),
            NSI_all = sum(no_NSI_muts),
            TRUNC_all = sum(no_TRUNC_muts),
            NTDam_all = sum(no_NTDam_muts),
            GOF_all = sum(no_GOF_muts),
            CNVGain_all = sum(CNVGain),
            HomozygousCNVLoss = sum(Copy_number==0),
            HeterozygousCNVLoss = sum(Copy_number==1))

#################################################################
###### Percentage of alterations, damaging alterations breakdown
#################################################################
damaging_2019 = Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKMover1 %>% group_by(Cancer_type) %>%
  summarise(n_sample=n_distinct(sample),
            n_dam_genes = n_distinct(symbol_19549),
            n_dam_alts = sum(no_TRUNC_muts+no_NTDam_muts+no_GOF_muts+CNVGain+CNVLoss),
            #            all_alt = sum(no_ALL_muts),
            NSI_all = sum(no_NSI_muts),
            TRUNC_all = sum(no_TRUNC_muts),
            NTDam_all = sum(no_NTDam_muts),
            GOF_all = sum(no_GOF_muts),
            CNVGain_all = sum(CNVGain),
            HomozygousCNVLoss = sum(Copy_number==0),
            HeterozygousCNVLoss = sum(Copy_number==1))

#damaging_alterations_2015 all samples with MutsCNVRNA in 2015
tmp = damaging_alterations_2015 %>% filter(sample %in% purity$sample)
damaging_2015 = tmp %>% group_by(Cancer_type) %>%
  summarise(n_sample=n_distinct(sample),
            n_dam_genes = n_distinct(symbol_19014),
            n_dam_alts = sum(no_TRUNC_muts+no_NTDam_muts+no_GOF_muts+CNVGain+CNVLoss),
            #            all_alt = sum(no_ALL_muts),
            NSI_all = sum(no_NSI_muts),
            TRUNC_all = sum(no_TRUNC_muts),
            NTDam_all = sum(no_NTDam_muts),
            GOF_all = sum(no_GOF_muts),
            CNVGain_all = sum(CNVGain),
            HomozygousCNVLoss = sum(Copy_number==0),
            HeterozygousCNVLoss = sum(Copy_number==1))

breakdown_portion_persample=function(df){
  df = df %>% mutate(portion_TRUNC = TRUNC_all/n_dam_alts,
                                portion_NTDam = NTDam_all/n_dam_alts,
                                portion_GOF = GOF_all/n_dam_alts,
                                portion_CNVGain = CNVGain_all/n_dam_alts,
                                portion_HomoLoss=HomozygousCNVLoss/n_dam_alts,
                                portion_HeteLoss=HeterozygousCNVLoss/n_dam_alts,
                                
                                persample_n_dam_alts = n_dam_alts/n_sample,
                                persample_TRUNC = TRUNC_all/n_sample,
                                persample_NTDam = NTDam_all/n_sample,
                                persample_GOF = GOF_all/n_sample,
                                persample_CNVGain = CNVGain_all/n_sample,
                                persample_HomoLoss=HomozygousCNVLoss/n_sample,
                                persample_HeteLoss=HeterozygousCNVLoss/n_sample,
                                label = paste0(Cancer_type," (",n_sample,")"))
  df
}

damaging_2015 = breakdown_portion_persample(damaging_2015)
damaging_2019 = breakdown_portion_persample(damaging_2019)
#damaging_2019 = damaging_2019[which(damaging_2019$Cancer_type %in% damaging_2015$Cancer_type),]

#Number of mutations per sample in TCGA15 vs 19
tmp = cbind(data.frame(colMeans(damaging_2015[,c(18:24)])),data.frame(colMeans(damaging_2019[,c(18:24)])))
colnames(tmp) = c('Ann_2015','Ann_2019')

stat = merge(damaging_2015 %>% select(Cancer_type, n_sample, n_dam_alts),
             damaging_2019 %>% select(Cancer_type, n_sample, n_dam_alts), by='Cancer_type', all.y = T)

stat[is.na(stat)] =0
stat = stat %>% mutate(label = paste0(Cancer_type," (", n_sample.x,"vs",n_sample.y,")"))
tmp = rbind(melt(damaging_2015[,c(1,19:24)], id="Cancer_type") %>% mutate(ann_year=2015),
            melt(damaging_2019[,c(1,19:24)], id="Cancer_type") %>% mutate(ann_year=2019))
tmp = merge(tmp, stat %>% select(Cancer_type, label), by='Cancer_type', all.x = T)
#tmp$variable =  as.factor(gsub("persample_","",tmp$variable))
tmp = tmp %>% filter(Cancer_type=='ACC')
p <- ggplot(aes(y = value, x = factor(ann_year), fill = variable), 
            data=tmp)+
            #data =tmp%>% mutate(value = replace(value, value==0, 0.0001))) + 
            #data =tmp[which(!tmp$variable %in% c('persample_CNVGain','persample_HomoLoss','persample_HeteLoss')),])+ # %>% mutate(value = replace(value, value==0, 0.0001))) + 
            #data =tmp[which(!tmp$variable %in% c('CNVGain_all','HomozygousCNVLoss','HeterozygousCNVLoss')),] %>% mutate(value = replace(value, value==0, 0.0001))) + 
  geom_bar(stat="identity",position = 'stack') +
  #geom_text(data=tmp, aes(x = label, y = value, label = paste0(round(value),"%")), size=4) +
  xlab('Year')+ ylab('Number of damaging mutation by type per sample')+#ylab('Number of damaging mutation by type (log10)')+
  #scale_y_log10()+
  #facet_grid(~label,labeller = label_wrap_gen(width=4))+
  scale_fill_brewer(palette="Set3")+
  theme_boss_xtilted(base_size = 12)
p

tmp = cbind(data.frame(colMeans(damaging_2015[,c(18:24)])),data.frame(colMeans(damaging_2019[,c(18:24)])))
colnames(tmp) = c('TCGA_2015','TCGA_2019')

length(intersect(sample_Bailey,sample_Ascat_notQC))

length(intersect(sample_Bailey,sample_AscatQC))
length(intersect(sample_Bailey,sample_Ascat))
length(intersect(intersect(sample_Bailey,sample_Ascat),sample_Thorsson))
length(setdiff(intersect(sample_Bailey,sample_Ascat),sample_Thorsson))
length(setdiff(intersect(sample_Bailey,sample_Thorsson),sample_Ascat))
length(setdiff(intersect(sample_Ascat,sample_Thorsson),sample_Bailey))
length(intersect(sample_Bailey,sample_Thorsson))
length(intersect(sample_Bailey,sample_AscatQC))



###################################################################
#### Homozygous loss and Amplification
###################################################################

tmp = RNASeq %>% group_by(Cancer_type,CNVGain,HomoLoss) %>% 
  summarise(n=n(),min_exp = min(value, na.rm = TRUE), mean_exp = mean(value,na.rm=TRUE),median_exp = median(value,na.rm=TRUE), max_exp= max(value,na.rm=TRUE))
CNVLoss=RNASeq %>% group_by(Copy_number) %>% summarise(min_exp = min(value, na.rm = TRUE), mean_exp = mean(value,na.rm=TRUE), max_exp= max(value, na.rm=TRUE))
CNVGain=RNASeq_damaged %>% group_by(Cancer_type.x,CNVGain.x) %>% summarise(n=n(),
                                                             n_genes=n_distinct(symbol_19549.x),
                                                             n_sample=n_distinct(sample),
                                                             min_exp = min(value, na.rm = TRUE), median_exp = median(value,na.rm=TRUE),mean_exp = mean(value,na.rm=TRUE), max_exp= max(value, na.rm=TRUE))
gene_amp = unique(RNASeq[which(RNASeq$CNVGain==1),"symbol_19549"])
gene_del = unique(RNASeq[which(RNASeq$Copy_number==0),"symbol_19549"])

CNVLoss = RNASeq_damaged[which(RNASeq_damaged$symbol_19549.x %in% gene_del),] %>% group_by(Copy_number.x) %>%
  summarise(n=n(),min_exp = min(value), median_exp = median(value),mean_exp = mean(value), max_exp= max(value))

ggplot(CNVLoss, aes(x=Copy_number.x, y=mean_exp)) + 
  geom_point(aes(size=log(n))) +
  scale_y_log10()+
  xlab('Copy Number') + ylab('Log10(Average expression)')+
  theme_boss(base_size = 15)

tmp = RNASeq_damaged[which(RNASeq_damaged$symbol_19549.x %in% gene_del),]

ggplot(data=RNASeq_damaged_del %>% mutate(value = replace(FPKM, value==0, 0.9999)), aes(x=as.factor(HomoLoss),y=FPKM))+
      #geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(HomoLoss)), alpha=0.1)+
       geom_boxplot(aes(fill = as.factor(HomoLoss)), alpha=0.2) +
       scale_y_log10()+
       #facet_grid(~Cancer_type.x)+ 
       xlab('Homozygously deleted genes') + ylab('FPKM (log10)')+
       theme_boss(base_size = 15)+theme(legend.position="none")


KICH=RNASeq[which(RNASeq$Cancer_type=="KICH"),]
print(nrow(KICH))
KICH=merge(KICH, PANCAN_total_table_altered %>% mutate(key=paste0(sample,"-",symbol_19549)), by="key")
print(nrow(KICH))
KICH %>% group_by(HomoLoss) %>% 
  summarise(n=n(),min_exp = min(value, na.rm = TRUE), mean_exp = mean(value,na.rm=TRUE), max_exp= max(value, na.rm=TRUE))

RNASeq_damaged= RNASeq
RNASeq_damaged[which(is.na(RNASeq_damaged$value)),'value']=0
RNASeq_damaged[which(RNASeq_damaged$value<0),'value']=0

#RNASeq_damaged= RNASeq_damaged[which(!is.na(RNASeq_damaged$value)),]
RNASeq_damaged = merge(RNASeq_damaged, muts1_exon_ncg_onco_cnv_damaging%>% mutate(key=paste0(sample,"-",symbol_19549)), by="key")


 
RNASeq_damaged_amp=NULL
RNASeq_damaged_del=NULL
RNASeq_rosalind =NULL
for (i in cancertype){
  #tmp = RNASeq_damaged[which(RNASeq_damaged$Cancer_type.x==i),]
  tmp = muts1_exon_ncg_onco_cnv_damaging %>% filter(Cancer_type == i)%>% mutate(key=paste0(sample,"-",symbol_19549))
  if (i %in% c('OAC','OSCC')){
    protmp = readRDS(sprintf('~/Trang/CNV/RNASeq_rosalind/RNASeq_rosalind/%s.rds','ESCA'))
  } else {
    protmp = readRDS(sprintf('~/Trang/CNV/RNASeq_rosalind/RNASeq_rosalind/%s.rds',i))
  }
  gene_amp = unique(tmp$symbol_19549[which(tmp$CNVGain==1)])
  gene_del = unique(tmp$symbol_19549[which(tmp$Copy_number==0)])
  protmp = protmp %>% mutate(key = paste(substring(patient,1,12), symbol, sep = '-'))
  protmp = protmp %>% filter(is.na(symbol)==FALSE)
  if (FALSE){
  tmp_amp =tmp[which(tmp$CNVGain.x==1),]
  
  tmp_notamp =tmp[which(tmp$CNVGain.x==1),]
  tmp_geneamp = unique(tmp_amp$symbol_19549.x)
  tmp_geneamp2=NULL
  for (ii in tmp_geneamp) {
    if (tmp_geneamp>0){tmp_geneamp2 = tmp_gen}
    }
  }
  tmp_geneamp = unique(tmp[which(tmp$CNVGain==1),"symbol_19549"])
  tmp_genedel = unique(tmp[which(tmp$Copy_number==0),"symbol_19549"])
  #tmp = merge(tmp,protmp, by='key')
  tmp_geneamp = protmp[which(protmp$symbol %in% tmp_geneamp),]
  tmp_geneamp$CNVGain = ifelse(tmp_geneamp$symbol %in% gene_amp, 1,0)
  tmp_genedel = protmp[which(protmp$symbol %in% tmp_genedel),]
  tmp_genedel$Homoloss = ifelse(tmp_genedel$symbol %in% gene_del, 1,0)
    
  RNASeq_damaged_amp=rbind(RNASeq_damaged_amp,tmp_geneamp)
  RNASeq_damaged_del=rbind(RNASeq_damaged_del,tmp_genedel)
  RNASeq_rosalind = rbind(RNASeq_rosalind, tmp)
  message(sprintf('done %s', i))
}
ggplot(data=RNASeq_damaged_amp %>% mutate(FPKM = replace(FPKM, FPKM==0, 1)), aes(x=as.factor(CNVGain),y=FPKM))+
  #geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(CNVGain.x)), alpha=0.1)+
  geom_boxplot(aes(fill = as.factor(CNVGain)), alpha=0.2) +
  scale_y_log10()+
  #facet_grid(~Cancer_type)+ 
  #stat_summary(aes(label=round(..y..,2), y=0), fun.y=mean, geom="label", size=2) +
  #annotate("label", x=stat$Type, y=replace(stat$med, stat$med==0, 0.9999), label=stat$med) +
  #annotate("label", x=1:nrow(CNVGain), y=replace(CNVGain$median_exp, CNVGain$median_exp==0, 0.9999), label=CNVGain$median_exp) + 
  xlab('Gene amplification') + ylab('FPKM (log10)')+ #ylab('GeneExpression (upper quartile normalized RSEM)')+
  theme_boss(base_size = 15)+theme(legend.position="none")
tmp=RNASeq_damaged_amp %>% group_by(Cancer_type,CNVGain) %>%summarise(n=n(),n_gene= n_distinct(symbol_19549),nsample=n_distinct(sample),
                                                                          #min_RSEMestimate = min(value), median_RSEMestimate = median(value),mean_RSEMestimate = mean(value), max_RSEMestimate= max(value),
                                                                          min_FPKM = min(FPKM), median_FPKM = median(FPKM),mean_FPKM = mean(FPKM), max_FPKM= max(FPKM))
tmp_alldam=muts1_exon_ncg_onco_cnv_damaging %>% group_by(Cancer_type,CNVGain) %>%summarise(n=n(),n_gene= n_distinct(symbol_19549),nsample=n_distinct(sample))


write.table(tmp,'~/Desktop/TCGA_geneAmp_exp.txt')


ggplot(data=RNASeq_rosalind %>% mutate(FPKM = replace(FPKM, FPKM==0, 0.9999)), aes(x=as.factor(HomoLoss),y=FPKM))+
  #geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(CNVGain.x)), alpha=0.1)+
  geom_boxplot(aes(fill = as.factor(HomoLoss)), alpha=0.2) +
  scale_y_log10()+
  #facet_grid(~Cancer_type)+ 
  #stat_summary(aes(label=round(..y..,2), y=0), fun.y=mean, geom="label", size=2) +
  #annotate("label", x=stat$Type, y=replace(stat$med, stat$med==0, 0.9999), label=stat$med) +
  #annotate("label", x=1:nrow(CNVGain), y=replace(CNVGain$median_exp, CNVGain$median_exp==0, 0.9999), label=CNVGain$median_exp) + 
  xlab('Gene homozygous deletion') + ylab('FPKM (log10)')+ #ylab('GeneExpression (upper quartile normalized RSEM)')+
  theme_boss(base_size = 15)+theme(legend.position="none")+
  geom_hline(aes(yintercept=log10(1)))
tmp=RNASeq_damaged_del %>% group_by(Cancer_type,HomoLoss) %>%summarise(n=n(),n_gene= n_distinct(symbol_19549),nsample=n_distinct(sample),
                                                                         #min_RSEMestimate = min(value), median_RSEMestimate = median(value),mean_RSEMestimate = mean(value), max_RSEMestimate= max(value),
                                                                          min_FPKM = min(FPKM), median_FPKM = median(FPKM),mean_FPKM = mean(FPKM), max_FPKM= max(FPKM))
 
  
write.table(tmp,'~/Desktop/TCGA_geneDel_exp.txt')

get_19549_entrez = function(x, geneSymbols){
  require(dplyr)
  geneSymbols$Entrez=as.character(geneSymbols$Entrez)
  geneSymbols= unique(geneSymbols[,c('NCG_symbol','Entrez')])
  x$Entrez=as.character(x$Entrez)
  x = merge(x, geneSymbols, by='Entrez', all.x = T) %>% rename(symbol_19549=NCG_symbol, entrez_19549=Entrez)
  x =x %>%select(-Symbol)
}


tmp = RNAseq_all[which(RNAseq_all$patient_barcode %in% allRNAseq_ID_patient),]
tmp = tmp %>% group_by(cancer.type) %>% summarise(n_samples = n_distinct(patient_barcode))

tmp = RNASeq_damaged[which(RNASeq_damaged$symbol_19549.x %in% gene_del),]

Total_table_2015 = NULL
CNV_2015 = NULL
for (i in cancertype) {
  load(sprintf("~/Thanos/%s/Tumor/Somatic_Mutations/TCGA_NSI/%s_concatenated_deduplicated_MutExpCNV_1SamplePerPatient_MMF_NonSi_InsDel_NoHM_MMFfilters_Annovar_dbNSFP_gAnn.Rdata",
               i,i))
  load(sprintf('~/Thanos/%s/Tumor/CNV/%s_concatenated_deduplicated_MutExpCNV_1SamplePerPatient_GainLoss_BLAT_overlaps_segmentsUnique.Rdata', i,i))
  
  for (ii in names(somatic_mutations)){
    tmp = data.frame(somatic_mutations[[ii]])
    Total_table_2015 = rbind(Total_table_2015, tmp)
  }
  for (ii in names(somatic_cnv)){
    tmp = data.frame(somatic_cnv[[ii]])
    CNV_2015 = rbind(CNV_2015, tmp)
  }
}


KICH_9= muts1_exon_ncg_onco_cnv_damaging %>% filter(Cancer_type=='KICH') %>% filter(Copy_number==0)


#### Check samples that have 25% segments with CN0, or 50% segments with CN0
for (i in ascat_tcga_ploidy$name){
  tmp = read.table(sprintf('~/Trang/CNV/ASCATforTCGA/segments/%s.txt',i))
  
}



muts_damaging_pancan= muts1_exon_ncg_onco_cnv_damaging[intersect(which(muts1_exon_ncg_onco_cnv_damaging$sample %in% sample_muts_cnv_exp),
                                           which(muts1_exon_ncg_onco_cnv_damaging$symbol_19549 %in% genes_muts_exp)),]

muts_damaging_ind = muts1_exon_ncg_onco_cnv_damaging_ind[intersect(which(muts1_exon_ncg_onco_cnv_damaging_ind$sample %in% sample_muts_cnv_exp),
                                                                 which(muts1_exon_ncg_onco_cnv_damaging_ind$symbol_19549 %in% genes_muts_exp)),]

tmp2 = merge(tmp_del[which(tmp_del$Homoloss==0),], tmp_del[which(tmp_del$Homoloss==1),], by=c("symbol","CODE"))
tmp3=tmp2[which(tmp2$CODE=='COAD'),]
tmp3=tmp3[which(t)]
plot(log10(tmp3$median_exp.x), log10(tmp3$median_exp.y),ylim = c(-4,4), 
     xlab="Not deleted", ylab = "Deleted",
     xlim = c(-4,4),cex=0.2)
abline(coef = c(0,1))

tmp = TCGA_RNASeq_damaged_homodel[setdiff(which(TCGA_RNASeq_damaged_homodel$CODE=='KICH'),
                                          which(TCGA_RNASeq_damaged_homodel$patient=='TCGA-KN-8419')),] %>% 
  group_by(Homoloss) %>% summarise(mean = mean(FPKM), median=median(FPKM))


tmp2 = merge(tmp_amp[which(tmp_amp$CNVGain==0),], tmp_amp[which(tmp_amp$CNVGain==1),], by=c("symbol","CODE"))
tmp3=tmp2[which(tmp2$CODE=='LAML'),]
plot(log10(tmp3$median_exp.x), log10(tmp3$median_exp.y),ylim = c(-4,4), 
     xlab="Not amplified", ylab = "Amplified",
     xlim = c(-4,4),cex=0.2)
abline(coef = c(0,1))

###########################################################
### Amplified and deleted genes TCGA
###########################################################
library(ggplot2)
library(reshape2)
GeneAmp = read.table('~/Trang-mnt/data/TCGA_geneAmp_exp_statsperGene.txt')
GeneAmp = merge(GeneAmp %>% filter(CNVGain==0), GeneAmp %>% filter(CNVGain==1), by=c('CODE','label','symbol'))
GeneAmp$foldchange = GeneAmp$mean_FPKM.y/GeneAmp$mean_FPKM.x # amplififed/notamplified
message(sprintf('Percentage of amplified genes with amplified expression signal is %s %s',
                round(length(which(GeneAmp$foldchange>1))*100 / nrow(GeneAmp)),'%'))
GeneAmp[is.na(GeneAmp)]<-0 


GeneDel = read.table('~/Trang-mnt/data/TCGA_geneDel_exp_statsperGene.txt')
tmp = GeneDel%>% group_by(CODE,Homoloss) %>% summarise(ngene)

GeneDel = merge(GeneDel %>% filter(Homoloss==0), GeneDel %>% filter(Homoloss==1), by=c('CODE','label','symbol'))
message(sprintf('Percentage of homozygously deleted genes with expression < 1FPKM is %s %s',
                round(length(which(GeneDel$mean_FPKM.y<1))*100 / nrow(GeneDel)),'%'))
tmp = merge(GeneAmp%>% group_by(label) %>% summarise(ngene_amp= n_distinct(symbol),percentAmpExp =sum(foldchange>1)*100/n()),
            GeneDel %>% group_by(label) %>% summarise(ngene_del= n_distinct(symbol),percentNoExp =(sum(mean_FPKM.y<1)*100/n())))
tmp$percentAmpExp_lab = paste0(round(tmp$percentAmpExp),"% (/",tmp$ngene_amp,")")
tmp$percentNoExp_lab = paste0(round(tmp$percentNoExp),"% (/",tmp$ngene_del,")")

ggplot(tmp[,c(1,2,3)], aes(label,percentAmpExp))+
  geom_bar(stat="identity")+
  xlab('Percentage of genes ')+
  #geom_text(aes(label=len), vjust=1.6, color="white", size=3.5)+
  theme_boss_xtilted()

t=cbind(melt(tmp[,c(1,3,5)],id='label'), melt(tmp[,c(1,6,7)],id='label'))
#colnames(t) = 
ggplot(melt(tmp[,c(1,3,5)],id='label'), aes(label,value,fill=variable))+
  geom_bar(position="dodge", stat="identity")+
  ylab('Percentage of proteins with increased exp \n and homozygously deleted proteins with low exp')+
  #facet_grid(~label,labeller = label_wrap_gen(width=4))+ 
  #geom_text(aes(label=len), vjust=1.6, color="white", size=3.5)+
  theme_boss_xtilted()

total_table = NULL
for (i in cancer){
  tmp = readRDS(sprintf('~/Rosalind/TCGA_2015_totaltable/%s_total_table.rds',i))
  if (ncol(tmp)==19){
    tmp = tmp %>% select(-gname)
  }
  
  total_table=rbind(total_table,tmp)
  message(sprintf('%s done', i))
}

View(GeneAmp_summary[whic])


#########################################################
#### Homozygously deleted genes expression
#########################################################
tmp = merge(Allcancertypes_pancan_annotation_damaging_withexp %>% group_by(CODE) %>% summarise(n_gene_dam = n_distinct(symbol_19549)),
            Allcancertypes_pancan_annotation_damaging_withexp %>% filter(Copy_number==0) %>% group_by(CODE) %>% 
  summarise(n_gene_del = n_distinct(symbol_19549),
            n_del_frag = sum(Copy_number==0)),
  by='CODE')


tmp2 = GeneDel %>% filter(Homoloss==1) %>% group_by(CODE,label) %>% 
  summarise(n_gene_del = n(), # all deleted genes,
            FPKMunder1 = sum(median_FPKM<1),
            FPKMhigher1 = sum(median_FPKM>1),
            FPKM1_keep = sum(median_FPKM<1)/n(),
            FPKMunder2 = sum(median_FPKM<2),
            FPKMhigher2 = sum(median_FPKM>2),
            FPKM2_keep = sum(median_FPKM<2)/n())
tmp = merge(tmp,tmp2, by=c('CODE','n_gene_del'))
write.table(tmp,'~/Desktop/Homodel_gene_FPKM12.txt', sep = '\t')


####Overall non-redundant damaged gene change:
Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKM12purity = Allcancertypes_pancan_annotation_damaging_withexp %>% mutate(HeteLossfromHomo1 = ifelse(key %in% TCGA_geneDel_toRMFPKM1purity$key,1,0),
                                                                                                                                       HeteLossfromHomo2 = ifelse(key %in% TCGA_geneDel_toRMFPKM2purity$key,1,0))

damaging_2019 = Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKM12purity %>% group_by(Cancer_type) %>%
  summarise(n_sample=n_distinct(sample),
            n_dam_genes = n_distinct(symbol_19549),
            n_dam_genes_redundant = n(),
            n_dam_alts = sum(no_TRUNC_muts+no_NTDam_muts+no_GOF_muts+CNVGain+CNVLoss),
            #            all_alt = sum(no_ALL_muts),
            NSI_all = sum(no_NSI_muts),
            TRUNC_all = sum(no_TRUNC_muts),
            NTDam_all = sum(no_NTDam_muts),
            GOF_all = sum(no_GOF_muts),
            CNVGain_all = sum(CNVGain),
            HomozygousCNVLoss1 = sum(Copy_number==0)-sum(HeteLossfromHomo1),
            HeterozygousCNVLoss1 = sum(Copy_number==1)+sum(HeteLossfromHomo1),
            HomozygousCNVLoss2 = sum(Copy_number==0)-sum(HeteLossfromHomo2),
            HeterozygousCNVLoss2 = sum(Copy_number==1)+sum(HeteLossfromHomo2))

Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKMover2purity = Allcancertypes_pancan_annotation_damaging_withexp %>% filter(!key %in% TCGA_geneDel_toRMFPKM2purity$key)


tmp = Allcancertypes_pancan_annotation_damaging_withexp %>% group_by(CODE) %>% 
  summarise(n_sample = n_distinct(sample),redundant_dam_genes = n(), homoloss = sum(Copy_number==0),unique_dam_gene = n_distinct(symbol_19549))
tmp = merge(tmp,
            merge(Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKMover1purity %>% group_by(CODE) %>% 
               summarise(redundant_dam_genes_1FPKMpurity = n(), homoloss_1FPKMpurity = sum(Copy_number==0),unique_dam_gene_1FPKMpurity = n_distinct(symbol_19549)),
             Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKMover2purity %>% group_by(CODE) %>% 
               summarise(redundant_dam_genes_2FPKMpurity = n(), homoloss_2FPKMpurity = sum(Copy_number==0),unique_dam_gene_2FPKMpurity = n_distinct(symbol_19549))))

# Purity consensus 
lm_eqn <- function(df){
  m <- lm(purity ~ consensus, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2),
                        b = format(coef(m)[2], digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


ascat_tcga_ploidy=ascat_tcga_ploidy %>% mutate(Sample.ID = substring(barcodeTumour,1,16))
purity_tmp = merge(ascat_tcga_ploidy %>% filter(rep==T&pass==T) %>% select(sample, Sample.ID, cancer_type, patient, purity),
                   purity %>% filter(Tumor_code<10) %>% select(sample, Sample.ID,Cancer.type, patient, ABSOLUTE,LUMP,IHC,CPE, consensus_mean), 
                   by.x=c("sample","Sample.ID","cancer_type","patient"), by.y=c('sample',"Sample.ID",'Cancer.type','patient'), all = T)
colnames(purity_tmp)[c(5,10)] = c('ASCAT','Consensus')
library(corrplot)
corrplot.mixed(cor(purity_tmp[,c(5:10)], use = "pairwise.complete.obs"))
sapply(purity_tmp, function(x) sum(is.na(x)==FALSE))

summary(lm(ASCAT~Consensus,purity_tmp)) # R2 = 0.42
cor.test(purity_tmp$ASCAT, purity_tmp$Consensus) # cor_coeff = 0.65
ggplot(purity_tmp,aes(y=purity, x =consensus))+
  geom_point(size=1) +
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_smooth(method="lm",se=FALSE) +
  #geom_text(x = 0.8,y = 0.2, label = lm_eqn(purity_tmp %>% select(purity, consensus)), parse = TRUE)+
  ylab('ASCAT purity') + xlab('TCGA consensus purity (ABSOLUTE,LUMP,IHC,CPE)')+
  theme_boss()


hist(purity_tmp$ASCAT)
hist(purity_tmp$Consensus)  
tmp = melt(purity_tmp[,c(2,5,10)])
ggplot(tmp, aes(x = value,  stat(count),fill =variable)) + geom_density(alpha = 0.5)+ xlab('Purity value') + theme_boss()
mid = intersect(PassSample_mc3, purity$Sample.ID)
mid_damaging = intersect(mid, purity$Sample.ID)
mid_damaged = intersect(mid, substring(sample_damaged,1,16))

tmp = merge(merge(purity_tmp[,c(2,3,5,10)], 
                  purity_tmp[which(purity_tmp$Sample.ID %in% mid),c(2,3,5,10)],by=c('Sample.ID','cancer_type'), all = T),
            purity_tmp[which(purity_tmp$Sample.ID %in% mid_damaged),c(2,3,5,10)],by=c('Sample.ID','cancer_type'), all = T)

colnames(tmp)[3:8] = c('ASCAT (9873)', 'TCGAbiolink (9364)','ASCAT (8158)','TCGAbiolink (8158)', 'ASCAT (6260)', 'TCGAbiolink (6260)')
ggplot(melt(tmp, id=c('Sample.ID','cancer_type')), aes(x=variable, y=value))+
  geom_boxplot()+
  stat_summary(aes(label=round(..y..,2)), fun.y=min, geom="text", size=4,position = position_nudge(y = -0.02)) +
  stat_summary(aes(label=round(..y..,2)), fun.y=median, geom="text", size=4,position = position_nudge(y = 0.02)) +
  stat_summary(aes(label=round(..y..,2)), fun.y=max, geom="text", size=4,position = position_nudge(y = 0.02)) +
  ylab('Purity value')+
  theme_boss()



##### plot damaging alterations distribution across cohorts
damaging_2019 = Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKM12purity_homotohetero %>% group_by(Cancer_type, sample) %>%
  summarise(n_sample=n_distinct(sample),
            n_dam_genes = n_distinct(symbol_19549),
            n_dam_genes_redundant = n(),
            n_dam_alts = sum(no_TRUNC_muts+no_NTDam_muts+no_GOF_muts+CNVGain+CNVLoss),
            #            all_alt = sum(no_ALL_muts),
            NSI_all = sum(no_NSI_muts),
            TRUNC_all = sum(no_TRUNC_muts),
            NTDam_all = sum(no_NTDam_muts),
            GOF_all = sum(no_GOF_muts),
            CNVGain_all = sum(CNVGain),
            HomozygousCNVLoss = sum(Copy_number==0),
            HeterozygousCNVLoss = sum(Copy_number==1),
            HomozygousCNVLoss1 = sum(Copy_number==0)-sum(HeteLossfromHomo1),
            HeterozygousCNVLoss1 = sum(Copy_number==1)+sum(HeteLossfromHomo1))
damaging_2019 = damaging_2019 %>% group_by(Cancer_type) %>% mutate(n_sample = n())
damaging_2019$label = paste0(damaging_2019$Cancer_type,' (',damaging_2019$n_sample,')')

damaging_2019_1 = Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKM12purity_homotohetero %>% group_by(Cancer_type) %>%
  summarise(n_sample=n_distinct(sample),
            n_dam_genes = n_distinct(symbol_19549),
            n_dam_genes_redundant = n(),
            n_dam_alts = sum(no_TRUNC_muts+no_NTDam_muts+no_GOF_muts+CNVGain+CNVLoss),
            #            all_alt = sum(no_ALL_muts),
            NSI_all = sum(no_NSI_muts),
            TRUNC_all = sum(no_TRUNC_muts),
            NTDam_all = sum(no_NTDam_muts),
            GOF_all = sum(no_GOF_muts),
            CNVGain_all = sum(CNVGain),
            HomozygousCNVLoss = sum(Copy_number==0),
            HeterozygousCNVLoss = sum(Copy_number==1),
            HomozygousCNVLoss1 = sum(Copy_number==0)-sum(FalsePosHomo1),
            HeterozygousCNVLoss1 = sum(Copy_number==1)+sum(HeteLossfromHomo1))

damaging_2019_1$label = paste0(damaging_2019_1$Cancer_type,' (',damaging_2019_1$n_sample,')')
damaging_2019_1$HomozygousCNVLoss-damaging_2019_1$HomozygousCNVLoss1

damaging_2019_1$FalseHomoCNVLoss=damaging_2019_1$HomozygousCNVLoss-damaging_2019_1$HomozygousCNVLoss1
damaging_2019_1$SavedHeteroCNVLoss=damaging_2019_1$HeterozygousCNVLoss1-damaging_2019_1$HeterozygousCNVLoss
plotthis = melt(damaging_2019_1[,c(13,12,15,16,17)], id= 'label')
plotthis$cat = ifelse(plotthis$variable %in% c('HomozygousCNVLoss1','FalseHomoCNVLoss'),'Homozygous loss','Heterozygous loss')
plotthis$cat = factor(plotthis$cat, levels = c('Homozygous loss','Heterozygous loss'))
#melt(damaging_2019_1[,c(15,16,17)], id= 'label')
#levels(plotthis$variable) = c('False positive homozygous loss', 'Saved Heterozygous loss (False Homozygous + 1 damaging mutation)')
ggplot(plotthis, aes(label, value)) + 
  geom_bar(stat = 'identity') + geom_text(aes(label=value), vjust=-0.5)+
  facet_wrap(variable~., nrow = 2,scales = 'free') +#, space = "free_x", scales = "free_x")+
  theme_boss_xtilted()


ggplot(plotthis, aes(label, value, fill=variable)) + 
  geom_bar(stat = 'identity') + geom_text(aes(label=value),position = position_stack(vjust = 0.75))+
  facet_wrap(cat~., nrow = 2,scales = 'free') +#, space = "free_x", scales = "free_x")+
  #scale_y_log10()+
  theme_boss_xtilted()

#Homozygous loss
tmp = data.frame(type = c('All homozygous deletion','HomoDel with FPKM<1', 'HomoDel with FPKM<1/purity', 'Homodel with FPKM<2','Homodel with FPKM<2/purity'),
                 tsg = c(2477,1304,1588,1626,1875),
                 tsg_per = c(1,1-0.47,1-0.36,1-0.34,1-0.24),
                 og = c(832,345,429,446,539),
                 og_per = c(1, 1-0.59,1-0.48,1-0.46,1-0.35))
plotthis = cbind(melt(tmp[,c(1,3,5)]), melt(tmp[,c(1,2,4)]))
colnames(plotthis) = c('type','var_por','portion','type1','var_var','value')
ggplot(plotthis, aes(type, value, fill=var_var)) + geom_bar(position ='dodge', stat = 'identity') + 
  geom_text(aes(label=paste0(value,'\n (',portion,')')), position=position_dodge(width=0.9), vjust=0.5) +theme_boss()

########### Plot the number of damaging alterations per sample and number of damaged genes per sample across cohorts
tmp = melt(damaging_2019[,c(16,2,4,6)], id=c('label', 'sample'))

plotthis2 = tmp[which(tmp$variable=="CNVGain_all"),]
plotthis2 = plotthis2 %>% group_by(label) %>% mutate(median = median(value))
plotthis2$category = "cancertype"

plotthis_pancan = plotthis2
plotthis_pancan$label="PANCAN (7630)"
plotthis_pancan$category = "PANCAN"
plotthis_pancan$median = median(plotthis_pancan$value)
stat = plotthis2 %>% group_by(label,category) %>% summarise(med = median(value),maxim = max(value),minim = min(value)) %>% arrange(med) #arrange(desc(med))
stat = rbind(stat, plotthis_pancan %>% group_by(label,category) %>% summarise(med = median(value),maxim = max(value),minim = min(value)))
plotthis2 = rbind(plotthis2, plotthis_pancan)
order = c(levels(reorder(plotthis2$label, plotthis2$median))[-18],"PANCAN (7630)")

ggplot(plotthis2 , aes(x=factor(label, order), y=value))+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(category)), alpha=0.1) +
  geom_boxplot(alpha=0.2)+
  #scale_y_log10()+
  #facet_wrap(~category, scales = 'free') +#, space = "free_x", scales = "free_x")+
  annotate("label", x=1:nrow(stat), y=replace(stat$med, stat$med==0, 0.9999), label=stat$med, na.rm = TRUE) +
  annotate("label", x=1:nrow(stat), y=replace(stat$maxim, stat$maxim==0, 0.9999), label=stat$maxim) + 
  annotate("label", x=1:nrow(stat), y=replace(stat$minim, stat$minim==0, 0.9999), label=stat$minim) +
  ylab("Number of damaging alterations") + xlab("Cancer type (#samples)")+
  theme_boss_xtilted()

########### Plot breakdown of the  damaged genes per sample across cohorts
tmp = melt(damaging_2019[,c(16,2,8:13)], id=c('label', 'sample'))

plotthis2 = tmp %>% group_by(label) %>% mutate(median = median(value))
plotthis2$category = "cancertype"

plotthis_pancan = plotthis2
plotthis_pancan$label="PANCAN (7630)"
plotthis_pancan$category = "PANCAN"
plotthis_pancan$median = median(plotthis_pancan$value)

stat = plotthis2 %>% group_by(label,category) %>% summarise(med = median(value),maxim = max(value),minim = min(value)) %>% arrange(med) #arrange(desc(med))
stat = rbind(stat, plotthis_pancan %>% group_by(label,category) %>% summarise(med = median(value),maxim = max(value),minim = min(value)))
plotthis2 = rbind(plotthis2, plotthis_pancan)
order = c(levels(reorder(plotthis2$label, plotthis2$median))[-16],"PANCAN (7630)")


ggplot(plotthis2 , aes(x=label, y=value, fill=variable))+
  geom_boxplot()+
  ylim(c(0,500))+
  #geom_bar(stat="identity",position = 'stack') +
  facet_grid(~category, scales = 'free', space = 'free') +#, space = "free_x", scales = "free_x")+
  ylab("Number of damaging alterations") + xlab("Cancer type (#samples)")+
  theme_boss_xtilted()

k=NULL
for (i in which(Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKM12purity$HeteLossfromHomo1==1)){
  if (Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKM12purity$no_TRUNC_muts[i]==1 |
      Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKM12purity$no_NTDam_muts[i]==1 ){
    k = c(k,i)
  }
}
#all k to remove (false positive homozygous deletion, not to be moved to heterozygous)

Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKM12purity_homotohetero = Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKM12purity
Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKM12purity_homotohetero$Copy_number[k]=1

k = setdiff(which(Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKM12purity$HeteLossfromHomo1==1), k) 
Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKM12purity_homotohetero$FalsePosHomo1=Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKM12purity_homotohetero$HeteLossfromHomo1
Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKM12purity_homotohetero$HeteLossfromHomo1[k]=0

Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKM12purity_homotohetero=Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKM12purity_homotohetero[-k,]
damaging_2019_1 = Allcancertypes_pancan_annotation_damaging_withexp_rmHomoDelFPKM12purity_homotohetero %>% group_by(Cancer_type) %>%
  summarise(n_sample=n_distinct(sample),
            n_dam_genes = n_distinct(symbol_19549),
            n_dam_genes_redundant = n(),
            n_dam_alts = sum(no_TRUNC_muts+no_NTDam_muts+no_GOF_muts+CNVGain+CNVLoss),
            NSI_all = sum(no_NSI_muts),
            TRUNC_all = sum(no_TRUNC_muts),
            NTDam_all = sum(no_NTDam_muts),
            GOF_all = sum(no_GOF_muts),
            CNVGain_all = sum(CNVGain),
            HomozygousCNVLoss = sum(Copy_number==0),
            HeterozygousCNVLoss = sum(Copy_number==1))
damaging_2019_1$label = paste0(damaging_2019_1$Cancer_type,' (',damaging_2019_1$n_sample,')')
write.table(damaging_2019_1, '~/Desktop/tmp.txt')



