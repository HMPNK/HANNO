

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
mamba create -n HANNO -c conda-forge -c bioconda bedtools samtools minimap2 miniprot last stringtie eggnog-mapper transdecoder ucsc-gtftogenepred ucsc-genepredtobed busco perl-bioperl
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

##change paths in main scripts starting with "HANNO.v0.3.pl" to fit to your system:

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

## now it should be ready to run ##
## test HANNO INSTALLATION:
unzip TTN-TEST-RUNS.zip
mamba activate HANNO
bash TTN-TEST-RUNS.sh

##This will run all use cases of the pipeline on the largest known vertebrate gene TTN
##CHECK the logs for early stops
##If the pipeline finishes successfully all logs will start and end with a date
```

# Output files
Output will bed12 format, you may convert to gtf using the scripts:
```sh
cut -f 1-12 ALLMODELS.bed12 | awk -f bed12ToGTF.awk > ALLMODELS.gtf
cut -f 1-12 ALLMODELS.bed12 | awk -f bed12ToGTF_addscore.awk > ALLMODELS.gtf # here the score field of CDS will be the total length of the ORF
#only the putative "best" model of a cluster of Models
cut -f 1-12 ALLMODELS.bed12 | grep -wFf BESTofCDScluster.list | awk -f bed12ToGTF.awk > ALLMODELS.gtf
cut -f 1-12 ALLMODELS.bed12 | grep -wFf BESTofCDScluster.list | > ALLMODELS.gtf
```

### Functional annotations for all transcripts with assigned CDS are in tables:
ALLMODELS.lastp.description.txt 
ALLMODELS.eggnog.description.txt
