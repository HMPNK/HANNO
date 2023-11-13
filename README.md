

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
mamba create -n HANNO -c conda-forge -c bioconda bedtools samtools minimap2 miniprot last eggnog-mapper transdecoder ucsc-gtftogenepred ucsc-genepredtobed busco perl-bioperl
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

##change paths in main scripts starting with "HANNO.v0.1-NCBI"  to fit to your system:

SCRIPTS="/path_to/HANNO/scripts";
MAMBA="/path_to/miniconda2/envs/MAMBA/bin";
EGGNOG="/path_to/HANNO/EGGNOGG-DBs";
THREADS=24;

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

# put your data in the current directory! You will need:

# genome.fasta = the assembly to be annotated
# workingdir = just a name, this directory will be created
# XXX_protein.faa = a protein fasta file from a NCBI Refseq annotation of a resonably close species (the closer the better, but diverged species work!)
# Alternatively you may use denovo transcriptome translated ORFs here, if you have RNAseq. You can concatenate multiple Proteomes into the file.
# XXX_rna_from_genomic.fna = a RNA fasta file from a NCBI Refseq annotation of a resonably close species (the closer the better, but diverged species work!)
# Alternatively, you may use denovo assembled transcripts here if you have RNAseq. You can concatenate multiple transcriptomes in the file.
# busco_lineage_dir = full path to directory ( /your/path/to/actinopterygii_odb10 )  where the BUSCO lineage data ist stored

##run HANNO:
mamba activate HANNO
bash ./scripts/HANNO.v0.1-NCBI.sh genome.fasta workingdir XXX_protein.faa XXX_rna_from_genomic.fna busco_lineage_dir > workingdir.log 2>&1

##You may also include stringtie transcriptome assemblies as gtf file in the annotation (make sure you have used the same genome reference for hisat"/stringtie as you are using here("genome.fasta")):

mamba activate HANNO
bash ./scripts/HANNO.v0.1-NCBI+GTF.sh genome.fasta workingdir XXX_protein.faa XXX_rna_from_genomic.fna busco_lineage_dir StringTie.gtf > workingdir.log 2>&1

##You may include a second protein database, which is ONLY used for final functional assignment.

mamba activate HANNO
scripts/HANNO.v0.1-NCBI+REFPROTs.sh genome.fasta workingdir XXX_protein.faa XXX_rna_from_genomic.fna busco_lineage_dir REFPROTDB.faa > workingdir.log 2>&1
##or with stringtie GTF
scripts/HANNO.v0.1-NCBI+GTF+REFPROTs.sh genome.fasta workingdir XXX_protein.faa XXX_rna_from_genomic.fna busco_lineage_dir StringTie.gtf REFPROTDB.faa > workingdir.log 2>&1
```

# Output files
Output will bed12 format, you may convert to gtf using the scripts:
```sh
bed12ToGTF.awk file.bed12 > file.gtf
bed12ToGTF_addscore.awk file.bed12 > file.gtf # here the score field of CDS will be the total length of the ORF
```

### Functional annotations for all transcripts with assigned CDS are in tables:
**out.emapper.annotations** -> all EGGNOG annotations  
**out.emapper.annotations.xlsx** -> same as above in Excel format  
**final_description.txt** -> gathered information from EGGNOG and LAST against protein database to populate the bed12 name fields  

### Annotated bed12 output:
**final.all.eggnog.bed12** -> all transcripts including eggnog annotations (eggnog results are more reliable than LAST, but also fewer annotatted transcripts)  
**final.all.last.bed12** ->  all transcripts including annotations from protein database (no score/evalue cut-off values below score 200 to be treated with care!), gene symbols are taken from eggnog  
**final.best.eggnog.bed12** -> one transcripts per gene (longest Orf) including eggnog annotations  
**final.best.last.bed12** -> one transcripts per gene (longest Orf) including annotations from protein database (no score/evalue cut-off values below score 200 to be treated with care!), gene symbols are taken from eggnog  
