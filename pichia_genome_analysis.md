Pichia Stipitis Assemblies analysis. 

QC of MiSeq Data for Pichia stipitis.

programs needed: Install with Conda: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:

```bash
for RawData in $(ls projects/niab_assemblies/p.stipitis/raw_dna/*/*/*/*.fq.gz); do
    ProgDir=~/git_repos/tools/seq_tools/dna_qc
    echo $RawData;
    sbatch $ProgDir/run_fastqc_slurm.sh $RawData
  done
```
```bash
for StrainPath in $(ls -d projects/niab_assemblies/p.stipitis/raw_dna/paired/ab920/F); do
    ProgDir=~/git_repos/tools/seq_tools/rna_qc
    IlluminaAdapters=/home/sv264/git_repos/tools/seq_tools/ncbi_adapters.fa
    ReadsF=$(ls $StrainPath/*.fq*)
    ReadsR=$(ls $StrainPath/*.fq*)
    echo $ReadsF
    echo $ReadsR
    sbatch $ProgDir/rna_qc_fastq-mcf_slurm.sh $ReadsF $ReadsR $IlluminaAdapters DNA
  done
```
Busco was used to assess completeness of the genomes. 

```bash
for Assembly in $(ls niab_assemblies/*/920/*.fasta); do
echo $Assembly
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=~/git_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d ~/git_repos/tools/gene_prediction/busco/saccharomycetes_odb10)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
sbatch $ProgDir/sub_busco3_slurm.sh $Assembly $BuscoDB $OutDir
done
```
```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/p*/*/*/*/short_summary*.txt); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

```
BUSCO Results                                                   Comp.   Dup.   Frag.  Miss.   Total 

short_summary.specific.saccharomycetes_odb10.918_Y-11543.txt    2125    0        6     6       2137
short_summary.specific.saccharomycetes_odb10.920_Y-11543.txt    2126    2        5     6       2137

List of missing BUSCOS in 918:
# Busco id
16391at4891
20594at4891
24627at4891
33944at4891
38279at4891
41181at4891

List of missing BUSCOS in 920:
24627at4891
25236at4891
33944at4891
38279at4891
41181at4891
6800at4891
```

GENOME ALIGNMENTS

First we did the genome alignments with Mummer for all the possible combinations required:

All versus Illumina of Y-11545.

```bash
Reference=$(ls assembly/misc_publications/p.stipitis/*/pichia.fasta)
for Query in $(ls niab_assemblies/pichia/*/*.fasta); do
Strain=$(echo $Query | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_Y-11545
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=git_repos/tools/seq_tools/genome_alignment/MUMmer
sbatch $ProgDir/sub_nucmer_slurm.sh $Reference $Query $Prefix $OutDir
done
```
```bash
Reference=$(ls assembly/misc_publications/p.stipitis/*/pichia.fasta)
for Query in $(ls niab_assemblies/pichia/*/*/*.fasta); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_Y-11545_split_chromosomes
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=git_repos/tools/seq_tools/genome_alignment/MUMmer
sbatch $ProgDir/sub_nucmer_slurm.sh $Reference $Query $Prefix $OutDir
done
```
918 vs 920

```bash
Reference=$(ls niab_assemblies/pichia/918/*.fasta)
for Query in $(ls niab_assemblies/pichia/920/*.fasta); do
Strain=$(echo $Query | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_918
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=git_repos/tools/seq_tools/genome_alignment/MUMmer
sbatch $ProgDir/sub_nucmer_slurm.sh $Reference $Query $Prefix $OutDir
done
```
918 vs 1082

```bash
Reference=$(ls niab_assemblies/pichia/918/*.fasta)
for Query in $(ls niab_assemblies/pichia/1082/*.fasta); do
Strain=$(echo $Query | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_918
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=git_repos/tools/seq_tools/genome_alignment/MUMmer
sbatch $ProgDir/sub_nucmer_slurm.sh $Reference $Query $Prefix $OutDir
done
```
918 and 920 vs Y-7124

```bash
Reference=$(ls assembly/misc_publications/p.stipitis/*/*.fasta)
for Query in $(ls niab_assemblies/pichia/9*/*.fasta); do
Strain=$(echo $Query | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_Y-7124
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=git_repos/tools/seq_tools/genome_alignment/MUMmer
sbatch $ProgDir/sub_nucmer_slurm.sh $Reference $Query $Prefix $OutDir
done
```
With Mummer alignments circos plots were generated (see specific folder)

A second round of genome alignments was done with bwa to plot read coverage over the assembly of Y-11545.

```bash
for Assembly in $(ls niab_assemblies/pichia/1082/*.fasta); do
Reference=$(ls assembly/misc_publications/p.stipitis/Y-11545_v2/pichia.fasta)
Strain=$(echo $Assembly | rev | cut -f2 -d '/'| rev)
Organism=$(echo $Reference | rev | cut -f3 -d '/' | rev)
Reads=$(ls niab_assemblies/pichia/1082/*.fastq.gz)
Prefix="${Organism}_${Strain}"
OutDir=analysis/genome_alignment/minimap/$Organism/vs_${Strain}
mkdir -p $OutDir
ProgDir=git_repos/tools/seq_tools/genome_alignment
sbatch $ProgDir/minimap/slurm_minimap2.sh $Reference $Reads $OutDir
done
```

(To see how the circos plots were generated see specific folder)


GENE MODELS 

Gene prediction was performed using fungap.

First FunGap was installed following the install instructions and the trial run was tested succesfully: https://github.com/CompSynBioLab-KoreaUniv/FunGAP/blob/master/INSTALL.md


Data download: RNAseq reads were downloaded for the reference genome

```bash
conda create -n sra-tools
conda activate sra-tools
conda install -c bioconda sra-tools
```
```bash
ProjDir=
cd $ProjDir
OutDir=raw_rna/paired/P.stipitis/Y-7124
mkdir -p $OutDir
fastq-dump --split-files --gzip --outdir $OutDir SRR8420582
```
This file is .gz but for Fungap needs to be .fastq, so gzip it before running FunGAP and then delete after using it to avoid using too much space.

```bash
gunzip -cf SRR8420582_1.fastq.gz > SRR8420582_1.fastq
gunzip -cf SRR8420582_2.fastq.gz > SRR8420582_2.fastq
```

Proteins sequences of related species were dowloaded
```bash
conda activate fungap 
FUNGAP_DIR=~/git_repos/tools/seq_tools/FunGAP
$FUNGAP_DIR/download_sister_orgs.py \
  --taxon "Pichia stipitis" \
  --email_address sv264@kent.ac.uk \
  --num_sisters 1
zcat sister_orgs/GCF_019049425.1_protein.faa.gz > prot_db_pichia.faa
```

Augustus species were obtained

```bash
conda activate fungap  # if you didn't do it already
$FUNGAP_DIR/get_augustus_species.py \
  --genus_name "Pichia stipitis" \
  --email_address sv264@kent.ac.uk

#Pichia stipitis 
```

Run Fungap

Note: This is a long-running script. As such, these commands were run using 'screen' to allow jobs to be submitted and monitored in the background. This allows the session to be disconnected and reconnected over time. Moreover, the program has to be run from the FunGAP folder, when ran from main folder an error looking for trinity outputs shows up. Probably could be solved by changing core script, but did no try.


```bash
screen -a

conda activate fungap

for Assembly in $(ls ~/niab_assemblies/*/920/*.fasta); do
echo $Assembly
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Reads1=$(ls ~/raw_rna/*/*/*/*_1.fastq)
Reads2=$(ls ~/raw_rna/*/*/*/*_2.fastq)
python fungap.py \
  --genome_assembly $Assembly \
  --trans_read_1 $Reads1 \
  --trans_read_2 $Reads2 \
  --augustus_species pichia_stipitis \
  --busco_dataset saccharomycetes_odb10 \
  --sister_proteome prot_db_pichia.faa \
  --num_cores 8 
done

Ctrl a + d
```
To run the same program with Slurm:

```bash
conda activate fungap

for Assembly in $(ls ~/niab_assemblies/*/1082/AB1082.fasta); do
echo $Assembly
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Reads1=$(ls ~/raw_rna/*/*/*/*_1.fastq)
Reads2=$(ls ~/raw_rna/*/*/*/*_2.fastq)
sbatch fungap_slurm_3.sh
done
```

SNP Calling.

First we need to do the alignment of the sequencing reads versus the genomes of interest with bwa. We will align both 918 and 920 to Y-11545 and then 918 to 1082.

1. AB918 and AB920 against Y-11545.

```bash
for Reference in $(ls assembly/misc_publications/p.stipitis/Y-11545_v2/*/GCF_000209165.1_ASM20916v1_genomic.fna); do
for StrainPath in $(ls -d qc_dna/paired/s.stipitis/ab9*/*); do
Strain=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f3 -d '/' | rev)
Ref=$(echo $Reference | rev | cut -f3 -d '/' | rev)
F_Read=$(ls $StrainPath/*_trim.fq.gz)
R_Read=$(ls $StrainPath/*_trim.fq.gz)
echo $F_Read
echo $R_Read
echo $Ref
Prefix="${Organism}_${Strain}"
OutDir=analysis/genome_alignment/bwa/$Organism/$Strain/vs_${Ref}
ProgDir=~/git_repos/tools/seq_tools/genome_alignment/bwa
sbatch $ProgDir/sub_bwa_slurm.sh $Prefix $Reference $F_Read $R_Read $OutDir
done
done
```

1.1 Pre SNP calling clean up.

1.1.1

Rename input mapping files in each folder by prefixing with the strain ID.

Alignment of the evolved strains were made against the parental strain.

```bash
  for File in $(ls analysis/genome_alignment/bwa/*/*/vs_Y-11545_v2/*_sorted.bam); do
    Strain=$(echo $File | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f4 -d '/' | rev)
    Prefix="vs_Y-11545_v2"
    echo "$Organism - $Strain"
    OutDir=analysis/popgen/${Prefix}/$Organism/$Strain
    CurDir=$PWD
    echo $OutDir
    mkdir -p $OutDir
    cp $CurDir/$File $OutDir/${Strain}_${Prefix}_sorted.bam
  done
```

1.1.2

Remove multimapping reads, discordant reads. PCR and optical duplicates, and add read group and sample name to each mapped read (preferably, the shortest ID possible)

Convention used: qsub $ProgDir/sub_pre_snp_calling.sh <SAMPLE_ID>

```bash
Reference=$(ls assembly/misc_publications/p.stipitis/Y-11545_v2/*/GCF_000209165.1_ASM20916v1_genomic.fna)
for Sam in $(ls analysis/popgen/vs_Y-11545_v2/*/*/*_sorted.bam); do
Strain=$(echo $Sam | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Sam | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=~/git_repos/scripts/pichia_analysis/popgen
sbatch $ProgDir/slurm_pre_snp_calling_2.sh $Sam $Strain $Reference
done
```
Prepare genome reference indexes required by GATK. Prepare for 918 and 920.

```bash
for Reference in $(ls assembly/misc_publications/p.stipitis/Y-11545_v2/*/GCF_000209165.1_ASM20916v1_genomic.fna); do
OutName=$(echo $Reference | sed 's/.fna/.dict/g')
OutDir=$(dirname $Reference)
mkdir -p $OutDir
ProgDir=/home/sv264/local/bin/picard-2.27.3
java -jar $ProgDir/picard.jar CreateSequenceDictionary R=$Reference O=$OutName
samtools faidx $Reference
done
```
SNP calling of 918 vs the reference genome

```bash
Isolate="ab918"
Reference=$(ls assembly/misc_publications/p.stipitis/Y-11545_v2/*/GCF_000209165.1_ASM20916v1_genomic.fna)
Ref=$(echo $Reference | rev | cut -f1 -d '/' | rev)
CurDir=$PWD
OutDir=analysis/popgen/SNP_calling/Y-11545_v2/vs_${Isolate}
mkdir -p $OutDir
ProgDir=/home/sv264/git_repos/scripts/pichia_analysis/popgen
sbatch $ProgDir/slurm_SNP_calling_918_vs_ref_2.sh $Reference $Isolate
```
918 vs ref Filter SNPs

```bash
Vcf=$(ls analysis/popgen/SNP_calling/Y-11545_v2/vs_ab918/_temp.vcf)
vcftools=/home/sv264/local/bin/vcftools_0.1.13/bin
vcflib=/home/sv264/miniconda3/pkgs/vcflib-1.0.0_rc2-h56106d0_2/bin
mq=40
qual=30
dp=10
gq=30
na=1.00
removeindel=N
echo "count prefilter"
cat ${Vcf} | grep -v '#' | wc -l

export LD_LIBRARY_PATH=/home/sv264/miniconda3/pkgs/bzip2-1.0.8-h7b6447c_0/lib
export LD_LIBRARY_PATH=/home/sv264/miniconda3/lib
export PYTHONPATH=/usr/lib64/python2.7
$vcflib/vcffilter -f "QUAL > $qual & MQ > $mq" $Vcf \
| $vcflib/vcffilter -g "DP > $dp & GQ > $gq" > ${Vcf%.vcf}_qfiltered.vcf

echo "count qfilter"
cat ${Vcf%.vcf}_qfiltered.vcf | grep -v '#' | wc -l
$vcftools/vcftools --vcf ${Vcf%.vcf}_qfiltered.vcf --max-missing $na --remove-indels --recode --out ${Vcf%.vcf}_qfiltered_presence
$vcftools/vcftools --vcf ${Vcf%.vcf}_qfiltered.vcf --max-missing $na --keep-only-indels --recode --out ${Vcf%.vcf}_indels_qfiltered_presence
```
```
count prefilter
50253

count qfilter
49481

After filtering, kept 2 out of 2 Individuals
Outputting VCF file...
After filtering, kept 45319 out of a possible 49481 Sites
Run Time = 1.00 seconds

After filtering, kept 2 out of 2 Individuals
Outputting VCF file...
After filtering, kept 3630 out of a possible 49481 Sites
Run Time = 0.00 seconds
```

Collect VCF stats

General VCF stats (remember that vcftools needs to have the PERL library exported)

```bash
  Isolate="Y-11545_v2"
  VcfTools=/home/sv264/local/bin/vcftools_0.1.13/perl
  export PERL5LIB="$VcfTools:$PERL5LIB"
  VcfFiltered=$(ls analysis/popgen/SNP_calling/${Isolate}/*/*_qfiltered_presence*.vcf | grep -v 'indels')
  Stats=$(echo $VcfFiltered | sed 's/.vcf/.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
  VcfFiltered=$(ls analysis/popgen/SNP_calling/${Isolate}/*/*_qfiltered_presence*.vcf | grep 'indels')
  Stats=$(echo $VcfFiltered | sed 's/.vcf/_indels.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
```
Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
Isolate="Y-11545_v2"
  for Vcf in $(ls analysis/popgen/SNP_calling/${Isolate}/*/*_qfiltered_presence*.vcf | grep -v 'indels'); do
      ProgDir=~/git_repos/scripts/popgen/snp
      $ProgDir/similarity_percentage.py $Vcf
  done
```

Create custom SnpEff genome database
```bash
SnpEff=/home/sv264/local/bin/snpEff
nano $SnpEff/snpEff.config
```
```
Add the following lines to the section with databases:

#---
# EMR Databases
#----
# Fus2 genome
Fus2v1.0.genome : Fus2
# Bc16 genome
Bc16v1.0.genome: BC-16
# P414 genome
P414v1.0.genome: 414
# 62471 genome
62471v1.0.genome: 62471
# R36_14 genome
R36_14v1.0.genome: R36_14
# SCRP371 genome
SCRP371v1.0.genome: SCRP371
# P. stipis
Ps589v1.0.genome: 589
PsY-11545v1.0.genome: Y-11545
PsCBS6054v1.0.genome: CBS6054
```

Collect input files

```bash
Organism="P.stipitis"
Strain="CBS6054"
DbName="PsCBS6054v1.0"
ProjDir=$PWD
Reference=$(ls $ProjDir/assembly/misc_publications/p.stipitis/Y-11545_v2/*/GCF_000209165.1_ASM20916v1_genomic.fasta)
Gff=$(ls $ProjDir/assembly/misc_publications/p.stipitis/Y-11545_v2/GCF_000209165.1_ASM20916v1_genomic.gff)
SnpEff=/home/sv264/local/bin/snpEff
mkdir $SnpEff/data/${DbName}
cp $Reference $SnpEff/data/${DbName}/sequences.fa
cp $Gff $SnpEff/data/${DbName}/genes.gff

#Build database using GFF3 annotation
java -jar $SnpEff/snpEff.jar build -gff3 -v ${DbName}
```

Annotate VCF files

```bash
Organism="P.stipitis"
Isolate="CBS6054"
Strain="Y-11545_v2"
DbName="PsCBS6054v1.0"
CurDir=/home/sv264
cd $CurDir
  for Vcf in $(ls analysis/popgen/SNP_calling/${Strain}/*/*_qfiltered_presence.recode.vcf | grep -v 'indels'); do
    echo $Vcf
    filename=$(basename "$Vcf")
    Prefix=${filename%.vcf}
    OutDir=$(dirname $Vcf)
    SnpEff=/home/sv264/local/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 ${DbName} $Vcf > $OutDir/"$Prefix"_annotated.vcf
    mv snpEff_genes.txt $OutDir/snpEff_genes_"$Prefix".txt
    mv snpEff_summary.html $OutDir/snpEff_summary_"$Prefix".html
    # mv 414_v2_contigs_unmasked_filtered* $OutDir/.
    #-
    #Create subsamples of SNPs containing those in a given category
    #-
    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'synonymous_variant') || (ANN[0].EFFECT has 'intron_variant') || (ANN[*].EFFECT has 'splice_region_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_gene.vcf
    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_nonsyn.vcf
    #synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_syn.vcf
    #Four-fold degenrate sites (output file suffix: 4fd)
    ProgDir=/home/sv264/local/bin/popgen/summary_stats
    python $ProgDir/parse_snpeff_synonymous.py $OutDir/"$Prefix"_syn.vcf
    AllSnps=$(cat $OutDir/"$Prefix"_annotated.vcf | grep -v '#' | wc -l)
    GeneSnps=$(cat $OutDir/"$Prefix"_gene.vcf | grep -v '#' | wc -l)
    CdsSnps=$(cat $OutDir/"$Prefix"_coding.vcf | grep -v '#' | wc -l)
    NonsynSnps=$(cat $OutDir/"$Prefix"_nonsyn.vcf | grep -v '#' | wc -l)
    SynSnps=$(cat $OutDir/"$Prefix"_syn.vcf | grep -v '#' | wc -l)
    printf "$filename\t$AllSnps\t$GeneSnps\t$CdsSnps\t$NonsynSnps\t$SynSnps\n"
done
```
```bash
_temp_qfiltered_presence.recode.vcf     45319   22162   19594   4791    14803
```
Perform interproscan annotation of the reference protein sequences. Interproscan needs java 11, which is installed in an individual environment called interproscan

```bash
export PATH=$PATH:/home/sv264/local/bin/my_interproscan/interproscan-5.56-89.0
Fasta=$(ls assembly/misc_publications/p.stipitis/Y-11545_v2/GCF_000209165.1_ASM20916v1_protein.faa)
OutDir=$(dirname $Fasta | sed 's&assembly/misc_publications&gene_pred/interproscan&g')
ProgDir=/home/sv264/git_repos/tools/seq_tools/feature_annotation/interproscan
sbatch $ProgDir/slurm_interproscan_2.sh $Fasta $OutDir
```

```bash
ProgDir=/home/sv264/git_repos/scripts/pichia_analysis/popgen
$ProgDir/extract_ref_annotations_2.py --genes_gff assembly/misc_publications/p.stipitis/Y-11545_v2/GCF_000209165.1_ASM20916v1_genomic.gff --refseq assembly/misc_publications/p.stipitis/Y-11545_v2/GCF_000209165.1_ASM20916v1_feature_table.txt --snp_vcf analysis/popgen/SNP_calling/Y-11545_v2/vs_ab918/_temp_qfiltered_presence.recode_nonsyn.vcf --InterPro gene_pred/interproscan/p.stipitis/Y-11545_v2/GCF_000209165.1_ASM20916v1_protein.faa.tsv > analysis/popgen/SNP_calling/Y-11545_v2/vs_ab918/ab918_annotation_table.tsv
```
