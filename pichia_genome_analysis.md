#Gene models 

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
