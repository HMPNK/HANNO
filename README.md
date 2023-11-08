

# **HANNO**: efficient **H**igh-throughput **ANNO**tation of eukaryote genomes

### You will need MAMBA for installation!

## Installation

```sh
git clone https://github.com/HMPNK/HANNO.git
cd HANNO
mamba create -n HANNO bedtools samtools minimap2 miniprot last eggnog-mapper transdecoder ucsc-gtftogenepred ucsc-genepredtobed
mamba activate HANNO
## for some reason I had to install BUSCO after creating the environment:
mamba install busco
## Also need bioperl
mamba install -c bioconda perl-bioperl
## TACO is not compatible with the environment create its own:
mamba create -n TACO taco

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

SCRIPTS="path_to/HANNO/scripts";
MAMBA="path_to/miniconda2/envs/MAMBA/bin";

## ready to run ##

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
scripts/HANNO.v0.1-NCBI+GTF+REFPROTs.sh genome.fasta workingdir XXX_protein.faa XXX_rna_from_genomic.fna busco_lineage_dir StringTie.gtf REFPROTDB.faa > workingdir.log 2>&1```

