#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=sv264@kent.ac.uk
#SBATCH --output=/home/sv264/%j.out

# Testing parallelisation of GATk HaplotypeCaller - may crash. (It did not! Resulted in 2x speedup)
# NOTE: this is a haploid organism. For diploid organism, change "ploidy" argument to 2.
# Changes required in the script:
# VARIABLES
# Reference - the genome reference used in read mapping.
# INSIDE THE GATK command:
# To specify which BAM mapping files (Out1 from pre_SNP_calling_cleanup.sh, RefName ending with "_rg" -> that is, with
# read group added) are to be used in SNP calling, use the -I argument with full path to each file following after that.
# Each new BAM file has to be specified after a separate -I

Reference=$1
Isolate=$2

# Project=/home/sv264
Project=/home/sv264
# OutDir=analysis/popgen/SNP_calling
# OutDir=analysis/popgen/SNP_calling
OutDir=analysis/popgen/SNP_calling/Y-11545_v2/vs_${Isolate}
mkdir $OutDir
# Reference=$(ls /home/sv264/assembly/misc_publications/p.stipitis/Y-11545_v2/pichia.fasta)


RefName=$(basename "$Reference")
Out1=$OutDir/"${RefName%.*}_temp.vcf"
Out2=$OutDir/"${RefName%.*}.vcf"

# ProgDir=/home/sv264/git_repos/tools/seq_tools/gatk-4.2.6.1/

# java -jar $ProgDir/GenomeAnalysisTK.jar \
#      -R $Project/$Reference \
#      -T HaplotypeCaller \
#      -ploidy 1 \
#      -nct 24 \
#      --allow_potentially_misencoded_quality_scores \
#      -I $Project/analysis/popgen/vs_Y-11545_v2/s.stipitis/ab918/ab918_vs_Y-11545_v2_sorted_nomulti_proper_sorted_nodup_rg.bam \
#      -I $Project/analysis/popgen/vs_Y-11545_v2/s.stipitis/ab920/ab920_vs_Y-11545_v2_sorted_nomulti_proper_sorted_nodup_rg.bam \
#      -o $Out1

# conda activate gatk4

$ProgDir/gatk --java-options HaplotypeCaller \
     -R $Project/$Reference \
     -ploidy 1 \
     -I $Project/analysis/popgen/vs_Y-11545_v2/s.stipitis/ab918/ab918_vs_Y-11545_v2_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/vs_Y-11545_v2/s.stipitis/ab920/ab920_vs_Y-11545_v2_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -O $Out1

#Break down complex SNPs into primitive ones with VariantsToAllelicPrimitives
#This tool will take an MNP (e.g. ACCCA -> TCCCG) and break it up into separate records for each component part (A-T and A->G).

# gatk VariantsToAllelicPrimitives \
#    -R $Reference \
#    -V $Out1 \
#    -o $Out2

#####################################
# Notes on GATK parallelisation
#####################################

# http://gatkforums.broadinstitute.org/gatk/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster
