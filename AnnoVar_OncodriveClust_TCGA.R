library(plyr)
library(dplyr)
library(tidyr) #tidyr is not installed in rosalind

options(stringsAsFactors = F)

print(commandArgs)

code = commandArgs(trailingOnly = TRUE)[1]
annovar_path = commandArgs(trailingOnly = TRUE)[2]
onco_path = commandArgs(trailingOnly = TRUE)[3]
mc3 = commandArgs(trailingOnly = TRUE)[4]
cgc_phen_path = commandArgs(trailingOnly = TRUE)[5]
gene_symbols_fn = commandArgs(trailingOnly = TRUE)[6]
save_dir=commandArgs(trailingOnly = TRUE)[7]


ns = c("nonsynonymous","stopgain","frameshift deletion","splicing","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonframeshift substitution","stoploss","frameshift substitution")
dam = c("nonsynonymous","frameshift deletion","frameshift insertion","frameshift substitution","splicing","stopgain","stoploss")
trunc = c("frameshift deletion","frameshift insertion","frameshift substitution","stopgain","stoploss") ## Always damaging==TRUE
non_trunc = c("nonsynonymous","splicing")

fix_splicing=function(x){
  x[which(x$Func.refGene=="splicing"),"ExonicFunc.refGene"]="splicing"
  x
}

fix_exonic=function(x){
  x[which(x$Func.refGene=="exonic;splicing"),"Func.refGene"]="exonic"
  x
}

fix_exonicFunc=function(x){
  if(length(which(x$ExonicFunc.refGene=="nonsynonymous SNV"))>0) x[which(x$ExonicFunc.refGene=="nonsynonymous SNV"),"ExonicFunc.refGene"]="nonsynonymous"
  if(length(which(x$ExonicFunc.refGene=="synonymous SNV"))>0) x[which(x$ExonicFunc.refGene=="synonymous SNV"),"ExonicFunc.refGene"]="synonymous"
  if(length(which(x$ExonicFunc.refGene=="stopgain SNV"))>0) x[which(x$ExonicFunc.refGene=="stopgain SNV"),"ExonicFunc.refGene"]="stopgain"
  if(length(which(x$ExonicFunc.refGene=="stoploss SNV"))>0) x[which(x$ExonicFunc.refGene=="stoploss SNV"),"ExonicFunc.refGene"]="stoploss"
  x
}

is_nonsilent=function(x){
  x$nonsilent=as.character(x$ExonicFunc.refGene)%in%ns
  x
}


set_IGV_code = function(x){
  x$IGV = paste(x$chr,":",x$end,sep="",coll="")
  x
}


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


## Careful here with the version of human genome (cmds below use hg19 add a parameter fi you want that changed or change it directly to the commands)
annotateMutations = function(mutation = mc3,
                             gene_symbols_fn=gene_symbols_fn,
                             save_dir=save_dir){
    
  #dir.create() does not crash in the directory already exist, it will just give a warning
  dir.create(file.path(save_dir, "ANNOVAR"), showWarnings = TRUE)
  save_dir=paste0(save_dir,"/ANNOVAR")
  annovar_input = paste0(save_dir, '/muts_mc3_original.txt')
  annovar_output = paste0(save_dir, '/muts_mc3_ann.txt')
  
  options(stringsAsFactors = F)
  
  ## Running ANNOVAR
  oldwd = getwd()
  setwd(annovar_path)
  mc3_muts = read.delim(mc3)
  mc3_muts = mc3_muts[which(mc3_muts$CODE==code),]
  mc3_muts = mc3_muts[,c("Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode","CODE")]
  
  # Fix end position for Annovar insertion annotation
  for (i in 1:nrow(mc3_muts)){
    if (mc3_muts$Reference_Allele[i]=='-'){
      mc3_muts$End_Position[i]=mc3_muts$Start_Position[i]
    }
  }
  
  write.table(mc3_muts, file=annovar_input, col.names = F, row.names = F, quote = F)
  
  cmd = paste0('perl table_annovar.pl ', annovar_input, ' humandb/ -buildver hg19 -out ', annovar_output, ' -protocol refGene,esp6500siv2_all,esp6500siv2_ea,1000g2015aug_all,1000g2015aug_eur,avsnp150,cytoBand,exac03,clinvar_20180603,cosmic87_coding,cosmic87_noncoding,dbnsfp35a,dbscsnv11 -operation g,f,f,f,f,f,r,f,f,f,f,f,f -otherinfo -nastring "." ')
  system(noquote(cmd))
  setwd(oldwd)
  
  ## Parse ANNOVAR's output
  message("Parsing ANNOVAR's results...")
  a = read.delim(file=paste0(annovar_output, ".hg19_multianno.txt"))
  a$Sample = sapply(strsplit(as.character(a$Otherinfo)," "),function(x) x[1])
  a$CODE = sapply(strsplit(as.character(a$Otherinfo)," "),function(x) x[2])
  a$key = with(a,paste0("chr",Chr,'.',Start,'.',End,'.',Ref,'.',Alt, '.', Sample))
  muts=a
  
  ns = subset(muts, as.character(ExonicFunc.refGene)%in%c("nonsynonymous SNV","nonsynonymous"))
  muts=fix_splicing(muts)
  muts=fix_exonic(muts)
  muts=fix_exonicFunc(muts)
  muts=is_nonsilent(muts)
  muts=set_IGV_code(muts)
  
  saveRDS(muts, file=paste0(save_dir, "/muts_annovar.rds"))
  
  ## Vectors of Booleans defining which variants are damaging according to function, splicing and conservation
  damFunc <- apply(muts[,c("SIFT_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred","LRT_pred","MutationTaster_pred","MutationAssessor_pred","FATHMM_pred")],
                   1, function(x){
                     sift = !is.na(x[1]) & x[1]=='D'
                     polyphen_hdiv = !is.na(x[2]) & x[2]=='D'
                     polyphen_hvar = !is.na(x[3]) & x[3]=='D'
                     lrt = !is.na(x[4]) & x[4]=='D'
                     mutTast = !is.na(x[5]) & (x[5]=='D' | x[5]=='A')
                     mutAss = !is.na(x[6]) & x[6]=='H'
                     fathmm = !is.na(x[7]) & x[7]=='D'
                     sift + polyphen_hdiv + polyphen_hvar + lrt + mutTast + mutAss + fathmm >= 5
                   })
  damSC <- apply(muts[,c("dbscSNV_ADA_SCORE","dbscSNV_RF_SCORE")], 1, function(x){ 
    x = suppressWarnings(as.numeric(x))
    any(!is.na(x) & x>=0.6) 
  })
  damCons <- apply(muts[,c("phyloP100way_vertebrate","GERP.._RS","SiPhy_29way_logOdds")], 1, function(x){
    x = suppressWarnings(as.numeric(x))
    phylop = !is.na(x[1]) & x[1]>1.6
    gerp = !is.na(x[2]) & x[2]>4.4
    siphy = !is.na(x[3]) & x[3]>12.17
    phylop + gerp + siphy >= 2 
  })
  
  
  ## Finally define damaging as TRUE or FALSE regardless if it is nonsynonymous or splicing
  muts$damagingFunc = muts$nonsilent & damFunc
  muts$damagingCons = muts$nonsilent & damCons
  muts$damagingNS = muts$damagingFunc | muts$damagingCons
  muts$damagingSC = damSC
  muts$truncating = as.character(muts$ExonicFunc.refGene) %in% c("stopgain","stoploss","stopgain","stoploss", "frameshift deletion", "frameshift insertion", "frameshift substitution")
  muts$damaging = muts$damagingNS | muts$damagingSC | muts$truncating
  
  saveRDS(muts, file=paste0(save_dir, "/muts_annovar_dbnsfp.rds"))
  
  # adding 19549
  load(gene_symbols_fn) ## object name symbol2entrez
  ## Check which of the gene symbol are in the symbols of our 19549
  ## First get index of Gene.refGene
  message("Annotating 19,549...")
  muts = get_19549(muts, geneSymbols=symbol2entrez)
  return(muts)
  
}

## Feed the mutations from above to the oncodrive clust function
## Here I need to run all the samples together because OncodriveClust does the prediction using recurrent mutations
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
  
  muts = readRDS(muts) 
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

if (TRUE){
  muts = annotateMutations(mutation = mc3,
                           gene_symbols_fn=gene_symbols_fn,
                           save_dir=save_dir)
  
  saveRDS(muts, paste0(save_dir,"/ANNOVAR/muts_ann_damaging.rds"))
}

df_mut = runOncodriveClust(muts=paste0(save_dir,"/ANNOVAR/muts_ann_damaging.rds"),
                           annovar_path=annovar_path,
                             ucsc_refgene_path=NULL,
                             onco_path = onco_path,
                             cgc_phen_path = cgc_phen_path,
                             save_dir = save_dir)

saveRDS(df_mut,file = paste0(save_dir,"/OncodriveClust/muts_ann_onco_damaging.rds"))


