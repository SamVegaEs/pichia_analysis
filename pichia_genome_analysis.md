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


#Busco was used to assess completeness of the genomes. 

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
