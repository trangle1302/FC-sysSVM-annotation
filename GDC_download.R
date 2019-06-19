devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
cancer_barcode = CODE$sampleID[which(CODE$CODE =='GBM')]
query.exp.hg19 <- GDCquery(project = "TCGA-ESCA",
                           data.category = "Gene expression",
                           data.type = "Gene expression quantification",
                           platform = "Illumina HiSeq", 
                           file.type  = "normalized_results",
                           #barcode = cancer_barcode,
                           experimental.strategy = "RNA-Seq",
                           legacy = TRUE)
GDCdownload(query.exp.hg19, method = "api", files.per.chunk = 10)
data <- GDCprepare(query.exp.hg19)
datatable(as.data.frame(colData(data)), 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
rowRanges(data)
