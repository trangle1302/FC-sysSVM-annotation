#!/bin/sh

# On Rosalind STDERR and STDOUT go to the same file by default.
# n Rosalind STDERR and STDOUT go to the same file by default.
# If you want them to go to separate files add the following line
#$ -j n

# To explicitly specify where you want STDOUT and STDERR to go
#$ -o /mnt/lustre/users/k1801782/dataset/CNV/log
#$ -e /mnt/lustre/users/k1801782/dataset/CNV/log

# Alternatively you can have these go the current directory with the following.
#  # $ -cwd

# How much memory do you need?
#$ -l h_vmem=18G

# Get email alerts when your job begins and ends
#$ -M trang.le@kcl.ac.uk
#$ -m be

# Which queues would you like to submit to?
#$ -q HighMemShortterm.q

# Todo: create $CODE as command line arguments, generalizing the rest of config.conf with <cancer_type>, and replace all following <cancer_type> in Rscripts with CODE
source /mnt/lustre/users/k1801782/dataset/SKCM_control/annotation_scripts_test2/config.conf
echo "Done sourcing configuration file"

while getopts ':c:s:' flag; do
  case "${flag}" in
    c) CODE="${OPTARG}" ;;
    s) SAVE_DIR="${OPTARG}" ;;
  esac
done

echo ${SAVE_DIR}
mkdir -p ${SAVE_DIR}

## CONFIGURATION
## --------------
module load general/R/3.4.1 #3.5.0 #tidyr is not installed in 3.5!
module load bioinformatics/bedtools2/2.25.0

echo "Done loading R and bedtools"

## Rscripts to be executed
## --------------
Rscript1='/mnt/lustre/users/k1801782/dataset/SKCM_control/annotation_scripts_test2/CNVGainLoss_ASCAT_TCGA.R'
Rscript2='/mnt/lustre/users/k1801782/dataset/SKCM_control/annotation_scripts_test2/AnnoVar_OncodriveClust_TCGA.R'
Rscript3='/mnt/lustre/users/k1801782/dataset/SKCM_control/annotation_scripts_test2/create_totaltable.R'

## Rscripts to be executed
## --------------

# CNV
# Rscript $Rscript1 ${CODE} ${ASCAT_DIR} ${SAVE_DIR}
# echo "Done annotating CNVs"

# Mutation
# Rscript $Rscript2 ${CODE} ${ANNOVAR_PATH} ${ONCO} ${MC3} ${cgc_phen_path} ${gene_symbols_fn} ${SAVE_DIR}
# echo "Done annotating snvs and indels (AnnoVAR and OncodriveClust)"

# Total_table
Rscript $Rscript3 ${CODE} ${SAVE_DIR}
echo "Total_table is created"

echo $HOSTNAME

