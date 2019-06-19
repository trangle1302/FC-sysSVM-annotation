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

module load general/R/3.4.1 #3.5.0 #tidyr is not installed in 3.5!
module load bioinformatics/bedtools2/2.25.0

Rscript '/mnt/lustre/users/k1801782/dataset/SKCM_control/annotation_scripts_test2/OncodriveClust_pancan.R'
