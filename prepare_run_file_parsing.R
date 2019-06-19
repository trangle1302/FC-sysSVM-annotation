# This was run in my local computer to create the submit_CNmap.sh (through sshfs - mounded drive)
lustre_mnt = "/Users/let/Trang/SKCM_control/annotation_scripts_test2/" # mounted path to Rosalind
lustre_wd = "/mnt/lustre/users/k1801782/dataset/" # working directory in Rosalind

# setwd("~/Documents/rosalind-mnt-group-ciccarelli/OAC/data/scripts/ncg6_OncodriveClust_1.0.0_annovar_2018Apr16")
CN_segment_DIR <- sub("SKCM_control/annotation_scripts_test2/","CNV/CN_segmentmean/",lustre_mnt)
submit_fn <- paste0(lustre_mnt,"submit_CNmap.sh")
file.create(submit_fn)

CN_cancer_paths <- list.files(CN_segment_DIR, recursive = F, full.names = F)[-6]
script_path = paste0(lustre_wd,'SKCM_control/annotation_scripts_test2/run.sh')
x = lapply(CN_cancer_paths, FUN = function(CN_cancer_path){
  path_to_logfiles = paste0(lustre_wd,'CNV/log/log_',CN_cancer_path)
  qsub_command <- paste0(
    "qsub -N CNmap", CN_cancer_path,
    " -q LowMemShortterm.q -l h_vmem=16G ",
    " -o ", path_to_logfiles, ".out",
    " -e ", path_to_logfiles, ".err ",
    script_path,
    " -c ", sub(".bed","",CN_cancer_path), 
    " -s ", lustre_wd, "Annotation_individual_cancers/", sub(".bed","",CN_cancer_path)
  )
  print(paste0("Writing ", CN_cancer_path, " into ", submit_fn))
  write(qsub_command, file = submit_fn, append =T)
})
 

