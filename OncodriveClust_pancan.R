library(plyr)
library(dplyr)
library(tidyr)
options(stringsAsFactors = F)

###################################################################
#### Running OncodriveClust for pancancer
###################################################################

annovar_path='/mnt/lustre/users/k1801782/softwares/annovar'
onco_path='/mnt/lustre/users/k1623452/Software/Anaconda5.3-Python3.7/bin/oncodriveclust'
gene_coord_fn="/mnt/lustre/users/k1623452/Software/NCG/NCG6/ncg6_blat_exon_coordinates_hg19.RData"
gene_symbols_fn="/mnt/lustre/users/k1623452/Software/NCG/NCG6/alias.long.RData"
cgc_phen_path="/mnt/lustre/users/k1623452/ICGC_OAC/data/annotation.pipeline/ncg6/somatic/CGC_phen_v84.tsv"
OAC_ID='/mnt/lustre/users/k1801782/dataset/CNV/OAC_ids.rds'

save_dir= "/mnt/lustre/users/k1801782/dataset/Annotation_individual_cancers/"

runOncodriveClust = function(muts=NULL,
                             annovar_path=annovar_path,
                             ucsc_refgene_path=NULL,
                             onco_path = onco_path,
                             cgc_phen_path = cgc_phen_path,
                             save_dir = save_dir){
  
  #dir.create() does not crash in the directory already exist, it will just give a warning
  dir.create(file.path(save_dir, "OncodriveClust"), showWarnings = TRUE)
  save_dir=paste0(save_dir,"/OncodriveClust")
  
  if(is.null(muts)){
    stop("No mutations provided")
  }
  
  # === ONCOCLUSTER =====
  ## For OncodriveClust I need 1: Synonymous mutations, 2: Non-Synonymous mutations, 3: Transcript lengths
  message("Running OncodriveClust...")
  # Prepare INPUT files (if no path provided)
  if(is.null(ucsc_refgene_path)){
    ucsc_refgene_path <- paste0(save_dir, "/transcript_length.tsv")
    
    ucsc_refgene <- read.table(paste0(annovar_path, "/humandb/hg19_refGene.txt"), header = F, stringsAsFactors = F)
    colnames(ucsc_refgene) <- c("bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames")
    ## Replace the exon start with the cdsStart to exclude UTRs 5&3 and same for the cdsEnd
    pb <- txtProgressBar(min = 0, max = nrow(ucsc_refgene), style = 3)
    for (row in 1:nrow(ucsc_refgene)){
      ucsc_refgene$exonStarts[row] <- paste(paste(c(ucsc_refgene$cdsStart[row],unlist(strsplit(ucsc_refgene$exonStarts[row], "\\,"))[which(unlist(strsplit(ucsc_refgene$exonStarts[row], "\\,")) > ucsc_refgene$cdsStart[row] & unlist(strsplit(ucsc_refgene$exonStarts[row], "\\,")) < ucsc_refgene$cdsEnd[row])]), collapse=","), ",", sep = "")
      ucsc_refgene$exonEnds[row] <- paste(paste(c(unlist(strsplit(ucsc_refgene$exonEnds[row], "\\,"))[which(unlist(strsplit(ucsc_refgene$exonEnds[row], "\\,")) < ucsc_refgene$cdsEnd[row] & unlist(strsplit(ucsc_refgene$exonEnds[row], "\\,")) > ucsc_refgene$cdsStart[row])], ucsc_refgene$cdsEnd[row]), collapse=","), ",", sep = "")
      setTxtProgressBar(pb, row)
    }
    close(pb)
    ucsc_refgene <- ucsc_refgene %>% select(name2, name, chrom, exonStarts, exonEnds)
    ucsc_refgene <- do.call(rbind, lapply(split(ucsc_refgene, rownames(ucsc_refgene)), function(x) cbind(x[,1],x[,2],x[,3],  unlist(strsplit(x[,4],"\\,")), unlist(strsplit(x[,5],"\\,"))))) %>%
      data.frame(stringsAsFactors=F)
    colnames(ucsc_refgene) <- c("Symbol", "Transcript.id", "chrom", "exon_start", "exon_end")
    ucsc_refgene$exon_start <- as.numeric(ucsc_refgene$exon_start)
    ucsc_refgene$exon_end <- as.numeric(ucsc_refgene$exon_end)
    #ucsc_refgene = ddply(ucsc_refgene, .(Symbol, Transcript.id), mutate, n=1:length(exon_start), .progress = 'text')
    ucsc_refgene <- ucsc_refgene %>% mutate(length=exon_end-exon_start) %>% group_by(Symbol, Transcript.id) %>% summarize(CDS.length=sum(length)) %>% ungroup %>% subset(CDS.length!=0)
    write.table(ucsc_refgene, file=ucsc_refgene_path, quote = F, row.names = F, sep = "\t")
  }
  
  muts = muts #readRDS(muts) 
  message('loaded Annovar results')
  df_mut = as.data.frame(muts, stringsAsFactors=F)
  df_mut = unrowname(df_mut)
  nsyn = subset(df_mut, ExonicFunc.refGene%in%c("nonsynonymous SNV","nonsynonymous"))
  syn = subset(df_mut, ExonicFunc.refGene%in%c("synonymous SNV","synonymous"))
  
  
  tmp = do.call('rbind', strsplit(sapply(strsplit(as.character(nsyn$AAChange.refGene),"\\,"), function(x) x[1] ),"\\:"))
  tmp = data.frame(
    symbol =  tmp[,1],
    Transcript.id =  tmp[,2],
    exon =  tmp[,3],
    nChange =  tmp[,4],
    aaChange =  tmp[,5],
    key = with(nsyn, paste0("chr",Chr,'.',Start,'.',End,'.',Ref,'.',Alt, '.', Sample)),
    aa.position = gsub("[^0-9]", "", tmp[,5])
  )
  nsyn_fn = paste0(save_dir, "/nsyn.onco")
  write.table(tmp, file=nsyn_fn, row.names = F, quote = F, sep = '\t')
  
  tmp = do.call('rbind', strsplit(sapply(strsplit(syn$AAChange.refGene,"\\,"), function(x) x[1] ),"\\:"))
  tmp = data.frame(
    symbol =  tmp[,1],
    Transcript.id =  tmp[,2],
    exon =  tmp[,3],
    nChange =  tmp[,4],
    aaChange =  tmp[,5],
    key = with(syn, paste0("chr",Chr,'.',Start,'.',End,'.',Ref,'.',Alt, '.', Sample)),
    aa.position = gsub("[^0-9]", "", tmp[,5])
  )
  syn_fn = paste0(save_dir, "/syn.onco")
  write.table(tmp, file=syn_fn, col.names=T, row.names = F, quote = F, sep = '\t')
  
  onco_out_fn = paste0(save_dir, "/oncodriveclust-results.tsv")
  message(nsyn_fn)
  message(syn_fn)
  message(onco_out_fn)
  
  cmd = paste0(onco_path,' -c -m 5 --cgc ', cgc_phen_path," -o ", onco_out_fn, " ", nsyn_fn, " ", syn_fn, " ", ucsc_refgene_path)
  message(noquote(cmd))
  system(noquote(cmd))
  
  message("Parsing the output of OncodriverClust...")
  df_mut = tryCatch(
    {  
      onco_out <- read.table(onco_out_fn, header = T, sep = "\t", stringsAsFactors = F)
      onco_out <- onco_out %>% subset(QVALUE<=0.1)
      onco_in <- read.table(nsyn_fn, header = T, sep = "\t", stringsAsFactors = F)
      onco_in$Sample <- apply(onco_in, 1, function(x) unlist(strsplit(x[6], "\\."))[6])
      
      ## Each cluster a separate row - it will be easier to gather patients and check if clusters are costant across cancer types
      onco_out <- onco_out %>% mutate(CLUST_COORDS=strsplit(CLUST_COORDS, "\\,\\[")) %>%
        unnest(CLUST_COORDS) %>% mutate(CLUST_COORDS=gsub("\\[|\\]", "", CLUST_COORDS)) %>%
        separate(CLUST_COORDS, into=c("CLUST_COORDS_START", "CLUST_COORDS_END"), sep="\\,") %>%
        separate(CLUST_COORDS_END, into=c("CLUST_COORDS_END", "NUMBER_OF_MUTS_IN_CLUST"), sep="\\:")
      
      ## Find which samples have mutations in the clusters
      find_samples <- function(gene, start, end, onco_in){
        keys <- onco_in %>% subset(symbol==gene & as.numeric(aa.position) >= start & as.numeric(aa.position) <= end) %>% select(key)
        keys <- paste(keys$key, collapse=",")
        return(keys)
      }
      onco_out$KEY <- apply(onco_out, 1, function(x) find_samples(x[1], as.numeric(x["CLUST_COORDS_START"]), as.numeric(x["CLUST_COORDS_END"]), onco_in))
      
      ## Expand the key components in the table
      onco_out <- onco_out %>% mutate(KEY=strsplit(KEY, "\\,")) %>% unnest(KEY)
      onco_out$SAMPLE <- apply(onco_out, 1, function(x) unlist(strsplit(x[14], "\\."))[6])
      
      df_mut$key = with(df_mut, paste0("chr",Chr,'.',Start,'.',End,'.',Ref,'.',Alt, '.', Sample))
      df_mut$oncodriveClust = as.character(df_mut$key)%in%onco_out$KEY
      return(df_mut)
    },
    error=function(e){
      message('OncodriveClust identifies no cluster, might be caused by lack of coding silent entries for background model')
      print(e)
      df_mut$key = with(df_mut, paste0("chr",Chr,'.',Start,'.',End,'.',Ref,'.',Alt, '.', Sample))
      df_mut[,'oncodriveClust'] = FALSE
      return(df_mut)
    },
    warning=function(w){
      message('OncodriveClust identifies no cluster, might be caused by lack of coding silent entries for background model')
      print(w)
      df_mut$key = with(df_mut, paste0("chr",Chr,'.',Start,'.',End,'.',Ref,'.',Alt, '.', Sample))
      df_mut[,'oncodriveClust'] = FALSE
      return(df_mut)
    }
  )
  
  return(df_mut)
}



cancer = list.dirs(path = save_dir, full.names = FALSE, recursive = FALSE)


muts=NULL
for (i in cancer[which(cancer != 'OncodriveClust')]){
  annovar_output = readRDS(sprintf('/mnt/lustre/users/k1801782/dataset/Annotation_individual_cancers/%s/ANNOVAR/muts_ann_damaging.rds',i))
  annovar_output$CODE = i
  muts=rbind(muts, annovar_output)
}


df_mut = runOncodriveClust(muts=muts,
                           annovar_path=annovar_path,
                           ucsc_refgene_path=NULL,
                           onco_path = onco_path,
                           cgc_phen_path = cgc_phen_path,
                           save_dir = save_dir)

for (i in cancer[which(cancer != 'OncodriveClust')]){
  tmp = df_mut[which(df_mut$CODE==i),] # %>% filter(CODE == i)
  saveRDS(df_mut,file = paste0(save_dir,i,"/OncodriveClust/muts_ann_onco_damaging_pancan.rds"))
  message(sprintf('%s done', i))
}

###################################################################
#### TOTAL TABLE
###################################################################
saving_path = "/mnt/lustre/users/k1801782/dataset/Annotation_individual_cancers/"
cnv_annotated = paste0(saving_path,"/CNV/CNVGainLoss.rds")
mc3_annotated = paste0(saving_path,"/OncodriveClust/muts_ann_onco_damaging_pancan.rds")

get_19549 = function(x, geneSymbols){
  require(dplyr)
  
  ix = which(colnames(x)=="Gene.refGene")
  x$symbol_19549 = sapply(as.character(x[,ix]), function(x) {
    spl = unlist(strsplit(x, ";"))
    unique(geneSymbols[as.character(geneSymbols$Symbol)%in%spl,'NCG_symbol'])[1]
  } )
  x = x %>% left_join(geneSymbols%>%select(symbol_19549=NCG_symbol, entrez_19549=Entrez)%>%unique(), by='symbol_19549')
  x
}

is_nonsilent = function(x){
  x$nonsilent=x$ExonicFunc.refGene%in%ns
  x
}

remove_dups = function(df_cnv){
  
  dups = df_cnv %>% group_by(key) %>% mutate(n=n(), types=paste(unique(CNV_type_corrected), collapse=","), ntypes=length(unique(CNV_type_corrected))) %>% subset(n>1)
  if (nrow(dups)>0){
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
  } else {
    df_cnv$n = 1
  }
  return(df_cnv)
}

ns = c("nonsynonymous","stopgain","frameshift deletion","splicing","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","stoploss","frameshift substitution")
dam = c("nonsynonymous","frameshift deletion","frameshift insertion","frameshift substitution","splicing","stopgain","stoploss")
trunc = c("frameshift deletion","frameshift insertion","frameshift substitution","stopgain","stoploss") ## Always damaging==TRUE
non_trunc = c("nonsynonymous","splicing")

createTotalTable = function(muts=NULL, cnvs=NULL, svs=NULL, exclude_samples=NULL){
  
  ns = c("nonsynonymous","stopgain","frameshift deletion","splicing","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","stoploss","frameshift substitution")
  dam = c("nonsynonymous","frameshift deletion","frameshift insertion","frameshift substitution","splicing","stopgain","stoploss")
  trunc = c("frameshift deletion","frameshift insertion","frameshift substitution","stopgain","stoploss") ## Always damaging==TRUE
  non_trunc = c("nonsynonymous","splicing")
  
  ## Make the lists
  message("Integrating SNVs...")
  df_mut = muts
  df_mut$sample = substring(df_mut$Sample,1,12)
  # rm(muts)
  
  if(!is.null(exclude_samples)){
    df_mut = df_mut %>% subset(!sample%in%exclude_samples)
    message(paste0("Samples excluded in mutation data: ", paste0(exclude_samples, collapse=",")))
  }
  
  ## Fix nonsilent here
  df_mut = df_mut %>% select(-nonsilent)
  df_mut=is_nonsilent(df_mut)
  
  ## In order to get the number of all mutations per gene and because
  ## I have WGS data, I exclude mutations that fall in the following categories
  df_mut = df_mut %>% subset(Func.refGene!="" &
                               !grepl("downstream", df_mut$Func.refGene) &
                               !grepl("upstream", df_mut$Func.refGene) &
                               !grepl("intergenic", df_mut$Func.refGene) &
                               !grepl("ncRNA", df_mut$Func.refGene) &
                               !grepl("intronic", df_mut$Func.refGene) &    
                               !grepl("UTR", df_mut$Func.refGene))
  
  ## Exclude genes that are not in 19014
  df_mut = df_mut %>% subset(!is.na(entrez_19549))
  
  ## Create the total table
  total_muts = ddply(df_mut, .(sample, symbol_19549, entrez_19549), dplyr::summarise,
                     no_ALL_muts=n(),
                     no_NSI_muts=sum(nonsilent),
                     no_TRUNC_muts = sum(ExonicFunc.refGene %in% trunc),
                     no_NTDam_muts = sum(ExonicFunc.refGene %in% non_trunc & damaging),
                     no_GOF_muts = sum(oncodriveClust), .progress = 'text'
  )
  
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
  cnvs$sample = substring(cnvs$sample,1,12)
  message(sprintf("There are %s samples with snvs and cnvs for %s", as.character(length(intersect(unique(df_mut$sample),unique(cnvs$sample)))), code))
  cnvs = cnvs %>% mutate(key=paste(sample, entrez_19549, sep="."))
  df_cnv = cnvs %>% subset(key%in%total_table$key | !is.na(CNV_type_corrected)) ## Also get the real CNVs
  
  
  if(!is.null(exclude_samples)){
    df_cnv = df_cnv %>% subset(!sample%in%exclude_samples)
    message(paste0("Samples excluded in CNV data: ", paste0(exclude_samples, collapse=",")))
  }
  
  ## define Gains and Losses - this was done in previous step
  ## Deduplicate the CNV data, because some genes may fall into two regions (sometimes it can be gain and loss)
  message("Resolving duplicated entries in CNVs...")
  
  df_cnv = remove_dups(df_cnv)
  
  
  ## Create total table from mutations and CNVs
  total_table = total_table %>% subset(!is.na(entrez_19549)) %>% 
    full_join(df_cnv%>%select(sample, symbol_19549, entrez_19549, Total_CN, CNV_type_corrected, ploidy, n)%>%rename(CNV_entries=n)%>%subset(!is.na(entrez_19549)))
  
  
  total_table$na_19549 = apply(total_table[,c("symbol_19549", "entrez_19549")], 1, function(x) length(x[is.na(x)]))
  total_table = total_table %>% mutate(in_19549=ifelse(na_19549<2, TRUE, FALSE)) %>% select(-na_19549) %>% subset(in_19549==TRUE)
  message('Done create total table')  
  return(total_table)
  
}

getMLinput <- function(df, code, geneProperties_dir="~/Thanos/mourikisa/data/geneProperties_final_mmImputed.Rdata"){
  df = df %>% mutate(Cancer_type=code) %>% rename(Entrez=entrez_19549, Copy_number=Total_CN, CNV_type=CNV_type_corrected) #%>% select(-symbol_19549, -in_19549, -ploidy, -CNV_entries, -key)
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
  
  message("Fixing SVs...")
  df$BND[is.na(df$BND)] = 0
  df$INS[is.na(df$INS)] = 0
  df$INV[is.na(df$INV)] = 0
  
  ## Convert categorical features to multiple factors
  message("Performing cleaning of categorical variables...")
  ## CNV type
  df <- df %>% 
    mutate(CNVGain=ifelse(is.na(CNV_type), 0, ifelse(CNV_type=="Gain",1, 0)), 
           CNVLoss=ifelse(Copy_number==0 | Copy_number==1, 1, 0)) %>%
    select(-CNV_type)
  
  
  if (FALSE){ #adding system-level properties in NCG6 - ongoing...
    message("Joining table with systems-level properties...")
    load(geneProperties_dir)
    geneProperties = tmp#geneProperties_mmImputed
    geneProperties = geneProperties %>% select(-symbol, -cancer_type, -cancer_dom, -cancer_rec,-cancer_type_status) %>% rename(Entrez=entrez)
    
    df <- df %>% left_join(geneProperties, by=c("Entrez"))
    
    ## age
    df <- df %>% 
      mutate(old=ifelse(is.na(age), NA, ifelse(age=="old",1, 0)),
             young=ifelse(is.na(age), NA, ifelse(age=="young",1, 0))) %>%
      select(-age)
    
    ## origin
    df <- df %>% 
      mutate(luca=ifelse(is.na(origin), NA, ifelse(origin=="LUCA",1, 0)), 
             eukaryotes=ifelse(is.na(origin), NA, ifelse(origin=="Eukaryotes",1, 0)),
             metazoans=ifelse(is.na(origin), NA, ifelse(origin=="Metazoans",1, 0)),
             vertebrates=ifelse(is.na(origin), NA, ifelse(origin=="Vertebrates",1, 0)),
             opisthokonts=ifelse(is.na(origin), NA, ifelse(origin=="Opisthokonts",1, 0)),
             mammals=ifelse(is.na(origin), NA, ifelse(origin=="Mammals",1, 0)),
             primates=ifelse(is.na(origin), NA, ifelse(origin=="Primates", 1, 0))) %>% 
      select(-origin)
    
    ## Essentiality
    df <- df %>%
      mutate(essentiality_percentage = )
    
    ## exp.breadth.class
    df <- df %>% 
      mutate(selective=ifelse(is.na(exp.breadth), NA, ifelse(exp.breadth=="Selective",1, 0)), 
             always.expressed=ifelse(is.na(exp.breadth), NA, ifelse(exp.breadth=="AlwaysExpressed",1, 0)),
             middle=ifelse(is.na(exp.breadth), NA, ifelse(exp.breadth=="Middle",1, 0)),
             one.tissue=ifelse(is.na(exp.breadth), NA, ifelse(exp.breadth=="OneTissue",1, 0)),
             never.expressed=ifelse(is.na(exp.breadth), NA, ifelse(exp.breadth=="Neverexpressed",1, 0))) %>%
      select(-exp.breadth)
    
    df <- data.frame(df)
    
    message("Converting features to factors...")
    fcols <- c("duplicated",
               "WGD", "hub", "central", "CNVGain", "CNVLoss",
               "ExpT_ME", "ExpT_HE", "ExpT_LE", "ExpT_NE",
               "ExpT_NET", "old", "young", "luca", "eukaryotes",
               "metazoans", "vertebrates", "opisthokonts",
               "mammals", "primates", "selective", "always.expressed",
               "middle", "one.tissue", "never.expressed")
    cols <- which(colnames(df) %in% fcols)
    for(i in cols){
      df[,i] = factor(df[,i], levels = c(0,1))
    }
    
    ## Reorder columns
    df = df[,c(12, 1:2, 3:11, 13:length(df))]
  }
  return(df)
}



