library(ggplot2)
###################################################################
#### Imputation of system level properties
###################################################################
CGC = read.delim("~/Desktop/Annotation/SNV_annotation/CGC_phen.tsv")
NCG6_systemslevelproperties$cancer = ifelse(NCG6_systemslevelproperties$entrez %in% CGC$Entrez.GeneId,TRUE,FALSE)
# properties from column 25-36

# Check which properties have missing data
sapply(NCG6_systemslevelproperties, function(x) sum(is.na(x))/length(x))
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
# Impute only those properties

NCG6_systemslevelproperties_imputed = NCG6_systemslevelproperties %>%
  #mutate(cancer_type_status = ifelse(cancer_type %in% c("can","cgc","vog"),"can","rst")) %>%
  mutate(cancer = ifelse(entrez %in% CGC$Entrez.GeneId,1,0)) %>%
  mutate(origin=as.character(origin))

NCG6_systemslevelproperties_imputed=NCG6_systemslevelproperties_imputed %>% 
  group_by(cancer) %>%
  mutate(expressed_tissues_protein=ifelse(is.na(expressed_tissues_protein), median(expressed_tissues_protein, na.rm=TRUE), expressed_tissues_protein),
         essentiality_percentage=ifelse(is.na(essentiality_percentage), median(essentiality_percentage, na.rm=TRUE), essentiality_percentage),
         expressed_tissues_rnaseq=ifelse(is.na(expressed_tissues_rnaseq), median(expressed_tissues_rnaseq, na.rm=TRUE), expressed_tissues_rnaseq),
         ppin_clustering=ifelse(is.na(ppin_clustering), median(ppin_clustering, na.rm=TRUE), ppin_clustering),
         origin=ifelse(origin=="-", getmode(origin), origin),
         age=ifelse(origin %in% c("Eukaryotes","Opisthokonts","Last Universal Common Ancestor"),"old","young"), #Last Universal Common Ancestor=LUCA, missing data?
         hub=ifelse(ppin_degree>quantile(ppin_degree)["75%"],1,0),
         duplicability=ifelse(duplicability>0,1,0),
         central = ifelse(ppin_betweenness>median(ppin_betweenness),1,0)
  )

# lacking: young/old, exp.breadth = expressed_tissues_protein+expressed_tissues_rnaseq ?
# adding essentiality_precentage
# what to do with functional annotation




# Not used but potentially useful functions
continuous = NCG6_systemslevelproperties[ ,!sapply(NCG6_systemslevelproperties, is.factor)]
categorical = NCG6_systemslevelproperties[ ,sapply(NCG6_systemslevelproperties, is.factor)]

imputeMedian<-function(x) apply(x, 2, function(x){x[is.na(x)]<-median(x, na.rm=T); x})
imputeMode<-function(x) apply(x, 2, function(x){x[is.na(x)]<-mode(x, na.rm=T); x})
tmp = data.frame(lapply(NCG6_systemslevelproperties[,c(17,25:36)],function(x) {ifelse(is.numeric(x), imputeMedian(x), imputeMode(x))}))

for (i in colnames(NCG6_systemslevelproperties)){
  print(paste0(colnames(NCG6_systemslevelproperties)[i]," ",is.factor(NCG6_systemslevelproperties%>% select(i))))
  print()
}
NCG6_systemslevelproperties%>% group_by(cancer_type) %>% summarise(med=median(i))


###################################################################
#### Plotting
###################################################################
tmp = NCG6_systemslevelproperties_imputed[which(NCG6_systemslevelproperties_imputed$cancer==1),]
# Betweenness
ggplot()+
  geom_histogram(data=tmp, aes(x=log10(ppin_betweenness+1)), alpha=0.3)+
  geom_vline(data=tmp,
             aes(xintercept=log10(median(ppin_betweenness)+1)), 
             colour="#BB0000", linetype="dashed") + 
  annotate("text", x = log10(median(tmp$ppin_betweenness)+1)+1, y = 4000, colour="red", size=5,
           label = paste0("Median: ", round(median(tmp$ppin_betweenness))))+
  theme_boss()
  
geom_vline(data=tmp,
             aes(xintercept=log10(mean(ppin_betweenness)+1)), 
             colour="#BB0000", linetype="dashed")


# Hub
ggplot()+
  geom_histogram(data=tmp, aes(x=ppin_degree), alpha=0.3)+
  xlim(c(0,500))+
  ylim(c(0,200))+
  geom_vline(data=tmp,
             aes(xintercept=quantile(ppin_degree)["75%"]), 
             colour="#BB0000", linetype="dashed") + 
  annotate("text", x = quantile(tmp$ppin_degree)["75%"]+50, y = 150, colour="red", size=5,
           label = paste0(".75 quantile: ", round(quantile(tmp$ppin_degree)["75%"])))+
  theme_boss()
quantile(ppin_degree)["75%"]