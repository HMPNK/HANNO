

# HANNO: efficient High-throughput ANNOtation of protein coding genes in eukaryote genomes
* high speed, can annotate a large vertebrate genome in below 1h  
* transcript models include UTRs  
* can work completely comparatively (using NCBI Refseq annotations), no need for RNAseq from your organism (but still beneficial of course)  
* completely evidence based (no weird gene models from *ab initio* prediction)
* includes functional annotation using eggNog and protein homology


### You will need MAMBA for installation!
```sh
#if you have conda already installed do
conda create -n MAMBA -c conda-forge mamba
conda activate MAMBA
```

### Installation (tested on: Ubuntu 16.04.7 LTS, Ubuntu 22.04.3 LTS)


```sh
git clone https://github.com/HMPNK/HANNO.git
cd HANNO
chmod -R 750 scripts/*
mamba create -n HANNO -c conda-forge -c bioconda bedtools=2.27.1 seqtk samtools minimap2 miniprot last stringtie eggnog-mapper transdecoder ucsc-gtftogenepred ucsc-genepredtobed busco perl-bioperl
#it is important to install bedtools=2.27.1 , I had massive issues using bedtools v2.31.1 ("bedtools intersect" showing strange behaviour!)
mamba activate HANNO

## TACO is not compatible with the environment create its own:
mamba create -n TACO -c conda-forge -c bioconda taco

## get EGGNOG DATA
mkdir EGGNOGG-DBs
export EGGNOG_DATA_DIR=$PWD/EGGNOGG-DBs/
download_eggnog_data.py

## get BUSCO databases from https://busco-data.ezlab.org/v5/data/lineages/
## For example fishes:
wget https://busco-data.ezlab.org/v5/data/lineages/actinopterygii_odb10.2021-02-19.tar.gz
tar xvf actinopterygii_odb10.2021-02-19.tar.gz
rm actinopterygii_odb10.2021-02-19.tar.gz

##change paths in main script "HANNO.v0.4.pl" to fit to your system:

my $scr = "/path_to/HANNO/scripts";
my $mamba="/path_to/miniconda2/envs/MAMBA/bin";
my $eggnog="/path_to/HANNO/EGGNOGG-DBs";

#IMPORTANT sometimes the script run_BUSCO.py is not available (check by "which run_BUSCO.py" !). If it is not available, add it like this to the HANNO environment bin dir:
find $HOME/ | grep HANNO | grep run_BUSCO.py$
#this should output /path/to/run_BUSCO.py
#change permissions for this file:
chmod 770 /path/to/run_BUSCO.py
#copy file to bin directory of HANNO environment, similar like this example:
cp /home/user/miniconda2/envs/MAMBA/envs/HANNO/lib/python3.7/site-packages/busco/run_BUSCO.py /home/user/miniconda2/envs/MAMBA/envs/HANNO/bin/run_BUSCO.py
#check again
which run_BUSCO.py

## now it should be ready to run, try to execute:
scripts/HANNO.v0.4.pl

HANNO version 0.4 (High-throughput ANNOtation for eukaryote genomes)
Author: Heiner Kuhl, Phd (heiner.kuhl@igb-berlin.de)

THIS SCRIPT CREATES THE PIPELINE AS A BASH script

INPUT DATA: genome assembly, proteins and mRNAs from related organism or denovo transcriptome assemblies, gtf from stringtie reference guided transcriptome assembly, BUSCO database

ALWAYS USE RELATIVE PATH (e.g. "../../assembly/asm.fasta"), IF INPUT data is not in current directory!

      Options:
		-a your genome assembly (fasta, fasta.gz)
		-d output directory to be created
		-p proteins used for gene-modeling (fasta, fasta.gz)
		-r mRNAs to be used for gene modeling (fasta, fasta.gz)
		-g stringtie assembled transcripts (gtf, be sure the gtf was created using the genome assembly provided with \"-a\" )
		-b path to busco lineage database (e.g. /home/user/eukaryota_odb9)
		-B BUSCO only on all Models (0=on; 1=off, default=0)
		-P PROTEIN DB for functional annotation (fasta)
		-t number of threads to use (default 8)
		-E skip EGGNOG functional annotation (0 or 1=skip, default=0)

This script generates a bash script for running the pipeline! Write script to file and run by: nohup bash <script> & !
ERROR: NEED an assembly to annotate!!!

## test HANNO INSTALLATION:
unzip TTN-TEST-RUNS.zip
mamba activate HANNO
bash TTN-TEST-RUNS.sh

##This will run all use cases of the pipeline on the largest known vertebrate gene TTN
##CHECK the logs for early stops
##If the pipeline finishes successfully all logs will start and end with a date
```

### Output files
ALLMODELS-FINAL.bedDB -> Database of all transcript models and all assigned functional information   
BESTMODELS-FINAL.bedDB -> putative best scoring transcript model of a gene from above   
BESTMODELS-FINAL.gtf -> same as above as gtf format   
BESTMODELS-FINAL.mRNA.fa -> corresponding mRNA sequences as fasta   
BESTMODELS-FINAL.CDS.fa -> corresponding CDS sequences as fasta   
BESTMODELS-FINAL.AA.faa -> corresponding aminoacid sequences as fasta   
