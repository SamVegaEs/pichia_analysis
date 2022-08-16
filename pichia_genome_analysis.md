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

To generate the coverage plots with circos we need to generate the .tsv files. This can be also used for coverage position comparison when viewed in IGV

```bash
for Sam in $(ls analysis/genome_alignment/minimap/p.stipitis/vs_*/*_aligned_sorted.bam); do
  Target=$(echo $Sam | rev | cut -f3 -d '/' | rev)
  Strain=$(echo $Sam | rev | cut -f2 -d '/' | rev)
  echo "$Strain-$Target"
  OutDir=$(dirname $Sam)
  samtools depth -aa $Sam > $OutDir/${Strain}_${Target}_depth.tsv
done

for Strain in 580 918 1082; do
  for Cov in $(ls analysis/genome_alignment/minimap/p.stipitis/vs_*/*_depth.tsv); do
    echo ${Cov} | cut -f4,5,6 -d '/' --output-delimiter " - "
    cat $Cov | cut -f3 | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'
  done
done > analysis/genome_alignment/minimap/read_coverage.txt
```
```bash
  for Cov in $(ls analysis/genome_alignment/minimap/p.stipitis/vs_*/*_depth.tsv); do
    Strain=$(echo $Cov | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Cov | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=$(dirname $Cov)
    ProgDir=git_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Strain}_${Organism}_depth.tsv > $OutDir/${Organism}_${Strain}_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_depth_10kb.tsv
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
1.2 Create databases needed for the runs

We created a database for Y-11545 used for the SNP calling between the strains and the reference genome (Y-11545) and another one with the genome and gene models prediction of AB1082 for the comparison between ab918 and ab1082.

    1.2.1 Database for Y-11545.
    
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
    1.2.2 Database for ab1082
    
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
PsAB108.genome: ab1082
```

Collect input files

```bash
Organism="P.stipitis"
Strain="ab1082"
DbName="PsAB1082"
ProjDir=$PWD
Reference=$(ls $ProjDir/niab_assemblies/pichia/1082/AB1082.fasta)
Gff=$(ls git_repos/tools/seq_tools/FunGAP/fungap_out/fungap_out/*.gff3)
SnpEff=/home/sv264/local/bin/snpEff
mkdir $SnpEff/data/${DbName}
cp $Reference $SnpEff/data/${DbName}/sequences.fa
cp $Gff $SnpEff/data/${DbName}/genes.gff

#Build database using GFF3 annotation
java -jar $SnpEff/snpEff.jar build -gff3 -v ${DbName}
```

1.3 SNP calling of 918 vs the reference genome

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
    1.3.1 918 vs ref Filter SNPs

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

    1.3.2 Collect VCF stats

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
    1.3.3 Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
Isolate="Y-11545_v2"
  for Vcf in $(ls analysis/popgen/SNP_calling/${Isolate}/*/*_qfiltered_presence*.vcf | grep -v 'indels'); do
      ProgDir=~/git_repos/scripts/popgen/snp
      $ProgDir/similarity_percentage.py $Vcf
  done
```

    1.3.4 Annotate VCF files using the specific database.

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
    1.3.5 Perform interproscan annotation of the reference protein sequences. 

Interproscan needs java 11, which is installed in an individual environment called interproscan

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

1.4 SNP calling of 920 vs the reference genome

```bash
Isolate="ab920"
Reference=$(ls assembly/misc_publications/p.stipitis/CBS6054/GCF_000209165.1_ASM20916v1_genomic.fna)
Ref=$(echo $Reference | rev | cut -f1 -d '/' | rev)
CurDir=$PWD
OutDir=analysis/popgen/SNP_calling/CBS6054/vs_${Isolate}
mkdir -p $OutDir
ProgDir=/home/sv264/git_repos/scripts/pichia_analysis/popgen
sbatch $ProgDir/slurm_SNP_calling_920_vs_ref_2.sh $Reference $Isolate
```

    1.4.1 920 vs ref Filter SNPs

```bash
Vcf=$(ls analysis/popgen/SNP_calling/CBS6054/vs_ab920/_temp.vcf)
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
    1.4.2 Collect VCF stats

General VCF stats (remember that vcftools needs to have the PERL library exported)

```bash
  Isolate="CBS6054"
  VcfTools=/home/sv264/local/bin/vcftools_0.1.13/perl
  export PERL5LIB="$VcfTools:$PERL5LIB"
  VcfFiltered=$(ls analysis/popgen/SNP_calling/${Isolate}/*/*_qfiltered_presence*.vcf | grep -v 'indels')
  Stats=$(echo $VcfFiltered | sed 's/.vcf/.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
  VcfFiltered=$(ls analysis/popgen/SNP_calling/${Isolate}/*/*_qfiltered_presence*.vcf | grep 'indels')
  Stats=$(echo $VcfFiltered | sed 's/.vcf/_indels.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
```

    1.4.3 Annotate VCF files

```bash
Organism="P.stipitis"
Isolate="CBS6054"
Strain="CBS6054"
DbName="PsCBS6054v1.0"
CurDir=/home/sv264
cd $CurDir
  for Vcf in $(ls analysis/popgen/SNP_calling/${Strain}/vs_ab920/*_qfiltered_presence.recode.vcf | grep -v 'indels'); do
    echo $Vcf
    filename=$(basename "$Vcf")
    Prefix=${filename%.vcf}
    OutDir=$(dirname $Vcf)
    SnpEff=/home/sv264/local/bin/snpEff_v4_2_core/snpEff
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
2. AB1082 vs AB918

2.0 Alignment of AB918 and AB1082 with bwa. We do not have the Illumina reads of AB1082, so I will use the assembly for AB1082 as reference and then call the SNPs with the illumina reads of AB918

```bash
for Reference in $(ls niab_assemblies/pichia/1082/AB1082.fasta); do
for StrainPath in $(ls -d qc_dna/paired/s.stipitis/ab918/*); do
Strain=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f3 -d '/' | rev)
Ref=$(echo $Reference | rev | cut -f2 -d '/' | rev)
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

2.1 Pre SNP calling clean up.

2.1.1

Rename input mapping files in each folder by prefixing with the strain ID.

Alignment of the evolved strains were made against the parental strain.

```bash
for File in $(ls analysis/genome_alignment/bwa/*/*/vs_1082/*_sorted.bam); do
Strain=$(echo $File | rev | cut -f3 -d '/' | rev)
Organism=$(echo $File | rev | cut -f4 -d '/' | rev)
Prefix="vs_1082"
echo "$Organism - $Strain"
OutDir=analysis/popgen/${Prefix}/$Organism/$Strain
CurDir=$PWD
echo $OutDir
cp $CurDir/$File $OutDir/${Strain}_${Prefix}_sorted.bam
done
```

2.1.2

Remove multimapping reads, discordant reads. PCR and optical duplicates, and add read group and sample name to each mapped read (preferably, the shortest ID possible)

Convention used: qsub $ProgDir/sub_pre_snp_calling.sh <SAMPLE_ID>

```bash
Reference=$(ls niab_assemblies/pichia/1082/AB1082.fasta)
for Sam in $(ls analysis/popgen/vs_1082/*/*/*_sorted.bam); do
Strain=$(echo $Sam | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Sam | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=~/git_repos/scripts/pichia_analysis/popgen
sbatch $ProgDir/slurm_pre_snp_calling_2.sh $Sam $Strain $Reference
done
```
Prepare genome reference indexes required by GATK.

```bash
for Reference in $(ls niab_assemblies/pichia/1082/AB1082.fasta); do
OutName=$(echo $Reference | sed 's/.fasta/.dict/g')
OutDir=$(dirname $Reference)
ProgDir=~/local/bin/picard-2.27.3
java -jar $ProgDir/picard.jar CreateSequenceDictionary R=$Reference O=$OutName
samtools faidx $Reference
done
```
2.2 SNP calling of 1082 vs the 918 genome
 
```bash
Isolate="ab918"
Reference=$(ls niab_assemblies/pichia/1082/AB1082.fasta)
Ref=$(echo $Reference | rev | cut -f1 -d '/' | rev)
CurDir=$PWD
OutDir=analysis/popgen/SNP_calling/ab1082/vs_${Isolate}
mkdir -p $OutDir
ProgDir=~/git_repos/scripts/pichia_analysis/popgen
sbatch $ProgDir/slurm_SNP_calling_1082_vs_918.sh $Reference $Isolate
```

    2.2.1 918 vs 1082 Filter SNPs

```bash
Vcf=$(ls analysis/popgen/SNP_calling/ab1082/vs_ab918/_temp.vcf)
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
266

count qfilter
259

After filtering, kept 1 out of 1 Individuals
Outputting VCF file...
After filtering, kept 104 out of a possible 259 Sites
Run Time = 0.00 seconds

After filtering, kept 1 out of 1 Individuals
Outputting VCF file...
After filtering, kept 154 out of a possible 259 Sites
Run Time = 0.00 seconds
```

    2.2.2 Collect VCF stats

General VCF stats (remember that vcftools needs to have the PERL library exported)

```bash
  Isolate="ab918"
  VcfTools=/home/sv264/local/bin/vcftools_0.1.13/perl
  export PERL5LIB="$VcfTools:$PERL5LIB"
  VcfFiltered=$(ls analysis/popgen/SNP_calling/ab1082/vs_${Isolate}/*_qfiltered_presence*.vcf | grep -v 'indels')
  Stats=$(echo $VcfFiltered | sed 's/.vcf/.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
  VcfFiltered=$(ls analysis/popgen/SNP_calling/ab1082/vs_${Isolate}/*_qfiltered_presence*.vcf | grep 'indels')
  Stats=$(echo $VcfFiltered | sed 's/.vcf/_indels.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
```
    2.2.3 Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
Isolate="ab918"
  for Vcf in $(ls analysis/popgen/SNP_calling/ab1082/vs_${Isolate}/*_qfiltered_presence*.vcf | grep -v 'indels'); do
      ProgDir=~/git_repos/scripts/popgen/snp
      $ProgDir/similarity_percentage.py $Vcf
  done
```
    2.2.4 Annotate VCF files

```bash
Organism="P.stipitis"
Isolate="ab918"
Strain="ab918"
DbName="PsAB1082" 
CurDir=/home/sv264
cd $CurDir
  for Vcf in $(ls analysis/popgen/SNP_calling/ab1082/vs_${Strain}/*_qfiltered_presence.recode.vcf | grep -v 'indels'); do
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
_temp_qfiltered_presence.recode.vcf     104     95      86      32      54
```
Perform interproscan annotation of the reference protein sequences. Interproscan needs java 11, which is installed in an individual environment called interproscan

```bash
export PATH=$PATH:/home/sv264/local/bin/my_interproscan/interproscan-5.56-89.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PAT:~/miniconda3/lib
Fasta=$(ls gene_pred/FunGAP/pichia/1082/*.faa)
OutDir=$(dirname $Fasta | sed 's&assembly/misc_publications&gene_pred/interproscan&g')
ProgDir=/home/sv264/git_repos/tools/seq_tools/feature_annotation/interproscan
sbatch $ProgDir/slurm_interproscan_2.sh $Fasta $OutDir
```

```bash
ProgDir=/projects/oldhome/armita/git_repos/emr_repos/scripts/pichia/popgen

ProgDir=/home/sv264/git_repos/scripts/pichia_analysis/popgen
$ProgDir/extract_ref_annotations.py --genes_gff git_repos/tools/seq_tools/FunGAP/fungap_out/fungap_out/*.gff3 --refseq assembly/misc_publications/p.stipitis/Y-11545_v2/GCF_000209165.1_ASM20916v1_feature_table.txt --snp_vcf analysis/popgen/SNP_calling/ab1082/vs_ab918/_temp_qfiltered_presence.recode_nonsyn.vcf --InterPro gene_pred/interproscan/p.stipitis/Y-11545_v2/GCF_000209165.1_ASM20916v1_protein.faa.tsv > analysis/popgen/SNP_calling/ab1082/vs_ab918/ab918_annotation_table.tsv
```
