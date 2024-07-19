[![DOI](https://zenodo.org/badge/716269109.svg)](https://zenodo.org/doi/10.5281/zenodo.11532369)

### HANNO: efficient High-throughput ANNOtation of protein coding genes in eukaryote genomes
* Outperforms most methods in speed by a factor of at least 10-20X
* Can annotate a large vertebrate genome in below 1h.
* Includes functional annotation using eggNog, protein homology and BUSCO.
* Transcript models include UTRs
* Works well with medium diverged input evidence (i.e. using NCBI Refseq annotations of species that diverged less than 50 MYA).
* No need for RNAseq from your organism (but still beneficial of course)
* Completely evidence based.
* No repeat masking necessary
* Tested on: Ubuntu 16.04.7 LTS, Ubuntu 22.04.3 LTS, Fedora 8.8

  
<p align="center">
<img src="https://github.com/HMPNK/HANNO/assets/51913753/1ffd4d02-f148-4214-90ad-628ce828d050" width="75%" height="75%">
</p>


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
mamba create -n HANNO -c conda-forge -c bioconda bedtools=2.27.1 seqtk samtools minimap2 miniprot last stringtie eggnog-mapper transdecoder ucsc-gtftogenepred ucsc-genepredtobed busco=3.0.2 perl-bioperl mawk
#it is important to install bedtools=2.27.1 , I had massive issues using bedtools v2.31.1 ("bedtools intersect" showing strange behaviour!)
#it is important to install busco=3.0.2 as some important outputs changed in newer version!
mamba activate HANNO

## TACO is not compatible with the environment create its own:
mamba create -n TACO -c conda-forge -c bioconda taco

## get EGGNOG DATA
mkdir EGGNOGG-DBs
export EGGNOG_DATA_DIR=$PWD/EGGNOGG-DBs/
download_eggnog_data.py

## get BUSCO databases from https://busco-data.ezlab.org/v5/data/lineages/
## For example fishes:
wget https://busco-data.ezlab.org/v5/data/lineages/actinopterygii_odb10.2024-01-08.tar.gz
tar xvf actinopterygii_odb10.2024-01-08.tar.gz
rm actinopterygii_odb10.2024-01-08.tar.gz

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
mamba activate HANNO
bash TTN-TEST-RUNS.sh

##This will run all use cases of the pipeline on the largest known vertebrate gene TTN
##CHECK the logs for early stops
##If the pipeline finishes successfully all logs will start and end with a date
##Do not wonder that BUSCO results are 0%, the TTN gene tested is not present in the BUSCO db used. 
```

### Input files
* Proteins ("-p" / "-P") should be formatted like they come from NCBI to make use of it in functional annotation:
  >\>XP_032605298.2 titin [Taeniopygia guttata]
  MTTKAPTFTQPLQSVVALEGSAATFEAHISGFPVPEVNWFRDGQVLSAATLPGVQISFSDGRARLVIP...
* If no input via "-P", protein input via "-p" will be used for gene-modeling and for functional annotation.
* If input via "-p" and "-P", "-p" will be used for modelling and -P will be used for functional annotation.
* long RNA sequence (Refseq mRNA, ISOSEQ from your organism) should be input via "-r"
* If using NCBI Refseq \*rna_from_genomic\* files, <ins>**consider removing "miscrna", "ncrna", "precursorrna", "rrna" and "trna" to increase gene-level specificity!**</ins> This is especially true, if the human genome reference geneset is used as it contains much more of these non-coding RNAs than other annotations, which will induce HANNO to build more spurious gene-models! Similarly, if too many genes are predicted when using RNAseq or Helixer inputs (often occurs, if _de novo_ transcript assemblies are used!), one may consider removing transcripts which have no functional annotations by eggNog, BUSCO or protein homology from the output files.
```sh
#remove non-coding RNAs from NCBI RefSeq "*rna_from_genomic*" input:
seqtk comp GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz | grep  -vE 'miscrna|ncrna|precursorrna|rrna|trna' | cut -f 1 | seqtk subseq GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz /dev/stdin | gzip -c > GCF_009914755.1_T2T-CHM13v2.0_mRNA_from_genomic.fna.gz

#Alternatively, screen functionally annotated gene models after HANNO run (if de novo or reference guided assembled RNAseq or Helixer annotation was used as input):
awk 'BEGIN{OFS="\t";FS="\t"} {if($26!="-" || $39!="-" || ($21>60 || $17>300)){print}}' BESTMODELS-FINAL.bedDB > BESTMODELS-FINAL-with-DB-hits.bedDB
```
* short read RNAseq should be assembled reference guided by stringtie and the resulting gtf should be input by "-g"
* Alternatively, but time consuming, denovo assemblies of short read RNAseq can be input via "-r".
* Currently, HANNO will run only, if a protein input is provided. To override this to check performance on pure transcript data, one may input a single mappable protein via "-p".
* see use cases below
  
### Output files
* ALLMODELS-FINAL.bedDB -> bed12 like Database of all transcript models and all assigned functional information   
* BESTMODELS-FINAL.bedDB -> putative best scoring transcript models for each transcript cluster from above
* BESTMODELS-FINAL.gtf -> same as above as gtf format   
* BESTMODELS-FINAL.mRNA.fa -> corresponding mRNA sequences as fasta
* BESTMODELS-FINAL.CDS.fa -> corresponding CDS sequences as fasta   
* BESTMODELS-FINAL.AA.faa -> corresponding aminoacid sequences as fasta
  
  Note that "\*.bedDB" files can be viewed in IGV, after renaming them to "\*.bed". Due to the one-line format, they can also be browsed and edited in a table calculation software. For instance one can add fields of functional annotations (fields 13 to xy) to the bed name field (field 4) to add more information to the IGV visualization.
```sh
#Example for adding functional annotation from ".bedDB" to ".bed" for viewing in IGV (using field 26=gene symbol by eggNog and field 15=description from best protein hit):
awk 'BEGIN{FS="\t";OFS="\t"} {$4=$26" | "$15" | "$4;print}' BESTMODELS-FINAL.bedDB | cut -f1-12 > BESTMODELS-FINAL.for-IGV.bed
```     

### BENCHMARKING HANNO BY TEST-RUNS ON VERTEBRATE GENOMES

HANNO benchmark runtimes on different vertebrate clades (HPC-server: Intel(R) Xeon(R) CPU E7-8890 v4 @ 2.20GHz; 96-threads; 1TB RAM). Note that 6 (Birds, Amphibians), 8 (Fish) and 11 (Mammals) genome annotations were run in parallel.

![image](https://github.com/HMPNK/HANNO/assets/51913753/eeb80ea5-dc7b-4a29-9b87-8bada1712076)

Note that speed of HANNO can be increased significantly by switching off BUSCO and/or eggNog steps (no "-b" provided and/or "-E 1" ) on the cost of reduced efficency of best-model selection.

HANNO was tested on fish (n=8), amphibian (n=6), bird (n=6) and mammal genomes (n=11), by transferring Refseq Annotations (from _P. flavescens_, _B. bufo_, _T. guttata_ and _H. sapiens_) to closely related and diverged species. Results below show the trends for "Complete", "Fragmented" and "Missing" BUSCOs for species ordered by divergence time from reference species (confidence intervalls from www.timetree.org) as well as total recovery of CDS and UTR compared to the 0 MYA diverged reference genome. HANNOtations of genomes by reference proteins and mRNA that diverged less than 50 MYA are typically yielding good results in terms of total recovered CDS and UTR sequences and BUSCO statistics for all vertebrate clades tested.

![image](https://github.com/HMPNK/HANNO/assets/51913753/55905ae2-0c2d-4468-9881-5609a66b525d)

![image](https://github.com/HMPNK/HANNO/assets/51913753/61d148f2-8a46-42c3-aa27-62796ec124ae)

![image](https://github.com/HMPNK/HANNO/assets/51913753/d3d0a892-7a5d-4659-886f-b2ee32db788f)

![image](https://github.com/HMPNK/HANNO/assets/51913753/b628356e-a8ce-42a1-bd2e-a96bdfd7844c)

![image](https://github.com/HMPNK/HANNO/assets/51913753/c130f4ee-1461-476c-9075-79bc7c3cefd0)

```sh

#test-runs on diverged vertebrate genomes
mamba activate HANNO

#remove non-coding RNAs from input:
seqtk comp GCF_905171765.1_aBufBuf1.1_rna_from_genomic.fna.gz | grep  -vE 'miscrna|ncrna|precursorrna|rrna|trna' | cut -f 1 | seqtk subseq GCF_905171765.1_aBufBuf1.1_rna_from_genomic.fna.gz /dev/stdin | gzip -c > GCF_905171765.1_aBufBuf1.1_mRNA_from_genomic.fna.gz
#AMPHIBIANS                 
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-BUFVIR-v0.4-mRNA -a GCA_037900795.1_ASM3790079v1_genomic.fna.gz -p GCF_905171765.1_aBufBuf1.1_protein.faa.gz -r GCF_905171765.1_aBufBuf1.1_mRNA_from_genomic.fna.gz -b vertebrata_odb9 | bash > HANNO-BUFVIR-v0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-LEPFUS-v0.4-mRNA -a GCA_031893025.1_aLepFus1.hap1_genomic.fna.gz -p GCF_905171765.1_aBufBuf1.1_protein.faa.gz -r GCF_905171765.1_aBufBuf1.1_mRNA_from_genomic.fna.gz -b vertebrata_odb9 | bash > HANNO-LEPFUS-v0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-XENTRO-v0.4-mRNA -a GCF_000004195.4_UCB_Xtro_10.0_genomic.fna.gz -p GCF_905171765.1_aBufBuf1.1_protein.faa.gz -r GCF_905171765.1_aBufBuf1.1_mRNA_from_genomic.fna.gz -b vertebrata_odb9 | bash > HANNO-XENTROP-v0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-SPEBOM-v0.4-mRNA -a GCF_027358695.1_aSpeBom1.2.pri_genomic.fna.gz -p GCF_905171765.1_aBufBuf1.1_protein.faa.gz -r GCF_905171765.1_aBufBuf1.1_mRNA_from_genomic.fna.gz -b vertebrata_odb9 | bash > HANNO-SPEBOM-v0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-BOMBOM-v0.4-mRNA -a GCF_027579735.1_aBomBom1.pri_genomic.fna.gz -p GCF_905171765.1_aBufBuf1.1_protein.faa.gz -r GCF_905171765.1_aBufBuf1.1_mRNA_from_genomic.fna.gz -b vertebrata_odb9 | bash > HANNO-BOMBOM-v0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-BUFBUF-v0.4-mRNA -a GCF_905171765.1_aBufBuf1.1_genomic.fna.gz -p GCF_905171765.1_aBufBuf1.1_protein.faa.gz -r GCF_905171765.1_aBufBuf1.1_mRNA_from_genomic.fna.gz -b vertebrata_odb9 | bash > HANNO-BUFBUF-v0.4-mRNA.log 2>&1
                 
#BIRDS
#remove non-coding RNAs from input:
seqtk comp GCF_003957565.2_bTaeGut1.4.pri_rna_from_genomic.fna.gz | grep  -vE 'miscrna|ncrna|precursorrna|rrna|trna' | cut -f 1 | seqtk subseq GCF_003957565.2_bTaeGut1.4.pri_rna_from_genomic.fna.gz /dev/stdin | gzip -c > GCF_003957565.2_bTaeGut1.4.pri_mRNA_from_genomic.fna.gz
#                 
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-TAEGUT-V0.4-mRNA -a GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz -p GCF_003957565.2_bTaeGut1.4.pri_protein.faa.gz -r GCF_003957565.2_bTaeGut1.4.pri_mRNA_from_genomic.fna.gz -b aves_odb9 | bash > HANNO-TAEGUT-V0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-SERCAN-V0.4-mRNA -a GCF_022539315.1_serCan2020_genomic.fna.gz -p GCF_003957565.2_bTaeGut1.4.pri_protein.faa.gz -r GCF_003957565.2_bTaeGut1.4.pri_mRNA_from_genomic.fna.gz -b aves_odb9 | bash > HANNO-SERCAN-V0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-MYICAY-V0.4-mRNA -a GCF_022539395.1_myiCay2020_genomic.fna.gz -p GCF_003957565.2_bTaeGut1.4.pri_protein.faa.gz -r GCF_003957565.2_bTaeGut1.4.pri_mRNA_from_genomic.fna.gz -b aves_odb9 | bash > HANNO-MYICAY-V0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-PHASUP-V0.4-mRNA -a GCA_023637945.1_DD_ASM_B1_genomic.fna.gz -p GCF_003957565.2_bTaeGut1.4.pri_protein.faa.gz -r GCF_003957565.2_bTaeGut1.4.pri_mRNA_from_genomic.fna.gz -b aves_odb9 | bash > HANNO-PHASUP-V0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-GALGAL-V0.4-mRNA -a GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz -p GCF_003957565.2_bTaeGut1.4.pri_protein.faa.gz -r GCF_003957565.2_bTaeGut1.4.pri_mRNA_from_genomic.fna.gz -b aves_odb9 | bash > HANNO-GALGAL-V0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-DROHOL-V0.4-mRNA -a GCA_016128335.2_ZJU2.0_genomic.fna.gz -p GCF_003957565.2_bTaeGut1.4.pri_protein.faa.gz -r GCF_003957565.2_bTaeGut1.4.pri_mRNA_from_genomic.fna.gz -b aves_odb9 | bash > HANNO-DROHOL-V0.4-mRNA.log 2>&1
                 
#FISHES                 
#remove non-coding RNAs from input:
seqtk comp GCF_004354835.1_PFLA_1.0_rna_from_genomic.fna.gz | grep  -vE 'miscrna|ncrna|precursorrna|rrna|trna' | cut -f 1 | seqtk subseq GCF_004354835.1_PFLA_1.0_rna_from_genomic.fna.gz /dev/stdin | gzip -c > GCF_004354835.1_PFLA_1.0_mRNA_from_genomic.fna.gz
#
../scripts/HANNO.v0.4.pl -t 8 -d HANNO-PERFLA-V0.4-mRNA -a GCF_004354835.1_PFLA_1.0_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_mRNA_from_genomic.fna.gz -b actinopterygii_odb9 | bash > HANNO-PERFLA-V0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 8 -d HANNO-SANLUC-V0.4-mRNA -a GCF_008315115.2_SLUC_FBN_1.2_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_mRNA_from_genomic.fna.gz -b actinopterygii_odb9 | bash > HANNO-SANLUC-V0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 8 -d HANNO-EPIFUS-V0.4-mRNA -a GCF_011397635.1_E.fuscoguttatus.final_Chr_v1_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_mRNA_from_genomic.fna.gz -b actinopterygii_odb9 | bash > HANNO-EPIFUS-V0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 8 -d HANNO-DICLAB-V0.4-mRNA -a GCF_905237075.1_dlabrax2021_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_mRNA_from_genomic.fna.gz -b actinopterygii_odb9 | bash > HANNO-DICLAB-V0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 8 -d HANNO-ESOLUC-V0.4-mRNA -a GCF_011004845.1_fEsoLuc1.pri_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_mRNA_from_genomic.fna.gz -b actinopterygii_odb9 | bash > HANNO-ESOLUC-V0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 8 -d HANNO-CLAGAR-V0.4-mRNA -a GCF_024256425.1_CGAR_prim_01v2_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_mRNA_from_genomic.fna.gz -b actinopterygii_odb9 | bash > HANNO-CLAGAR-V0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 8 -d HANNO-DANRER-V0.4-mRNA -a GCA_033170195.1_ASM3317019v1_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_mRNA_from_genomic.fna.gz -b actinopterygii_odb9 | bash > HANNO-DANRER-V0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 8 -d HANNO-ORENIL2-V0.4-mRNA -a GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_mRNA_from_genomic.fna.gz -b actinopterygii_odb9 | bash > HANNO-ORENIL2-V0.4-mRNA.log 2>&1
                 
#MAMMALS
#remove non-coding RNAs from input:
seqtk comp GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz | grep  -vE 'miscrna|ncrna|precursorrna|rrna|trna' | cut -f 1 | seqtk subseq GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz /dev/stdin | gzip -c > GCF_009914755.1_T2T-CHM13v2.0_mRNA_from_genomic.fna.gz                 
#run HANNO without ncRNA in -r input
../scripts/HANNO.v0.4.pl -t 6 -d T2T-CHM13v2.0.0.4-mRNA -a GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_mRNA_from_genomic.fna.gz -b mammalia_odb10 | bash > T2T-CHM13v2.0.0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d PanPan.0.4-mRNA -a GCF_029289425.1_NHGRI_mPanPan1-v1.1-0.1.freeze_pri_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_mRNA_from_genomic.fna.gz -b mammalia_odb10 | bash > PanPan.0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d PonAbe1.0.4-mRNA -a GCF_028885655.1_NHGRI_mPonAbe1-v1.1-hic.freeze_pri_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_mRNA_from_genomic.fna.gz -b mammalia_odb10 | bash > PonAbe1.0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d CalJa1.0.4-mRNA -a GCF_011100555.1_mCalJa1.2.pat.X_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_mRNA_from_genomic.fna.gz -b mammalia_odb10 | bash > CalJa1.0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d GRCm39.0.4-mRNA -a GCF_000001635.27_GRCm39_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_mRNA_from_genomic.fna.gz -b mammalia_odb10 | bash > GRCm39.0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d ARS-UCD2.0.0.4-mRNA -a GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_mRNA_from_genomic.fna.gz -b mammalia_odb10 | bash > ARS-UCD2.0.0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d ARS1.2.0.4-mRNA -a GCF_001704415.2_ARS1.2_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_mRNA_from_genomic.fna.gz -b mammalia_odb10 | bash > ARS1.2.0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d OrcOrc.0.4-mRNA -a GCF_937001465.1_mOrcOrc1.1_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_mRNA_from_genomic.fna.gz -b mammalia_odb10 | bash > OrcOrc.0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d RhiFer1.0.4-mRNA -a GCF_004115265.2_mRhiFer1_v1.p_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_mRNA_from_genomic.fna.gz -b mammalia_odb10 | bash > RhiFer1.0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d EleMax1.0.4-mRNA -a GCF_024166365.1_mEleMax1_primary_haplotype_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_mRNA_from_genomic.fna.gz -b mammalia_odb10 | bash > EleMax1.0.4-mRNA.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d phaCin.0.4-mRNA -a GCF_002099425.1_phaCin_unsw_v4.1_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_mRNA_from_genomic.fna.gz -b mammalia_odb10 | bash > phaCin.0.4-mRNA.log 2>&1

```

### IMPROVING ANNOTATION BY ADDING SPECIES-LEVEL TRANSCRIPTOME DATA (genus-level should work, too)
To benchmark improvement of annotation by adding species-level RNAseq short reads, the _E. lucius_ run from above was supported by _E. lucius_ brain, ovary and testis data. Results show improvement in BUSCO statistics is due to reduced number of fragmented BUSCOs. Annotated CDS-length is improved by 18% and UTR-length shows massive improvement. A run of Helixer AI gene prediction was added for comparison.

### INCLUDING HELIXER ANNOTATION IN HANNO

Including the Helixer annotation into HANNO via the "-g" parameter can improve BUSCO scoring, but may cost a bit of gene-level precision. It seems to be a good choice, if only diverged reference proteins/mRNAs are available for your species.

![image](https://github.com/HMPNK/HANNO/assets/51913753/b7d4c916-eda1-4eee-bdbb-2521bb897e94)


![image](https://github.com/HMPNK/HANNO/assets/51913753/7a019f27-83ff-4bff-a03d-8968f9e3fca6)




```sh
#Build hisat2 index for genome to be annotated
hisat2-build -p 24 GCF_011004845.1_fEsoLuc1.pri_genomic.fna ESOLUC
#Map RNAseq data
hisat2 -p 16 --dta --summary-file ESOLUC.BRAIN.log -x ESOLUC -1 SRR1533652_1.fastq.gz -2 SRR1533652_2.fastq.gz -S ESOLUC.BRAIN.sam
hisat2 -p 16 --dta --summary-file ESOLUC.OVARY.log -x ESOLUC -1 SRR1533651_1.fastq.gz -2 SRR1533651_2.fastq.gz -S ESOLUC.OVARY.sam
hisat2 -p 16 --dta --summary-file ESOLUC.TESTIS.log -x ESOLUC -1 SRR1533661_1.fastq.gz -2 SRR1533661_2.fastq.gz -S ESOLUC.TESTIS.sam
#Convert SAM to BAM and sort
samtools sort -@16 -m 10G -o ESOLUC.BRAIN.srt.bam ESOLUC.BRAIN.sam
samtools sort -@16 -m 10G -o ESOLUC.OVARY.srt.bam ESOLUC.OVARY.sam
samtools sort -@16 -m 10G -o ESOLUC.TESTIS.srt.bam ESOLUC.TESTIS.sam

#remove SAM files to free disk:
rm *sam

#assemble transcripts from BAM-files:
stringtie -p 16 -l BRAIN -o ESOLUC.BRAIN.gtf ESOLUC.BRAIN.srt.bam
stringtie -p 16 -l OVARY -o ESOLUC.OVARY.gtf ESOLUC.OVARY.srt.bam
stringtie -p 16 -l TESTIS -o ESOLUC.TESTIS.gtf ESOLUC.TESTIS.srt.bam
#This took about 20 minutes (if running the 3 hisat2 jobs, the 3 samtools and the 3 stringtie jobs in parallel)

#just concatenate assembled trancript gtf-files prior to feeding it into HANNO:
cat ESOLUC.OVARY.gtf ESOLUC.TESTIS.gtf ESOLUC.BRAIN.gtf > ESOLUC.O+T+B.gtf

#TEST HANNO with transcriptome only (use dummy file to override HANNOS need for input "-p")
#create "dummy.faa" first that contains only a single protein sequence (must be a mappable one otherwise HANNO will stop, because miniprot gtf is empty)
seqtk seq -l 0 GCF_004354835.1_PFLA_1.0_protein.faa.gz | head -2 > dummy.faa
#Run HANNO (use all proteins via "-P" in functional annotation only)
../scripts/HANNO.v0.4.pl -t 80 -d HANNO-ESOLUC-V0.4-TRANS-only -a GCF_011004845.1_fEsoLuc1.pri_genomic.fna.gz -p dummy.faa -b actinopterygii_odb9 -g ESOLUC.O+T+B.gtf -P GCF_004354835.1_PFLA_1.0_protein.faa.gz | bash > HANNO-ESOLUC-V0.4-TRANS-only.log 2>&1
#Using 80 threads this took 34 minutes (mainly due to BUSCO and eggNog)!

#Run HANNO with proteins and mRNA from a diverged species and reference guided transcript assemblies from same species:
../scripts/HANNO.v0.4.pl -t 80 -d HANNO-ESOLUC-V0.4-TRANS -a GCF_011004845.1_fEsoLuc1.pri_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_mRNA_from_genomic.fna.gz -b actinopterygii_odb9 -g ESOLUC.O+T+B.gtf | bash > HANNO-ESOLUC-V0.4-TRANS.log 2>&1

#Using 80 threads HANNO finished in 48 minutes!
#It was also tested to input denovo assembled instead of reference guided transcripts from the same species, but it delivers similar results and is much less efficient in terms of computing time (denovo transcript assembly takes too much time!)

```


```sh
#convert Helixer gff3 to stringtie-like gtf first:
#create helixer gtf as input for HANNO "-g"
awk -f ../scripts/helixergff3Togtf.awk HELIXER-GCF_011004845.1_fEsoLuc1.pri_genomic.fna.gff3 |gtfToGenePred stdin stdout | genePredToBed stdin stdout | awk -f ../scripts/bed12ToGTF_addscore-tacoFake.awk > helixer-input.gtf

#run HANNO with Helixer-annotation input via "-g"
../scripts/HANNO.v0.4.pl -t 80 -d HANNO-ESOLUC-V0.4+helixer -a GCF_011004845.1_fEsoLuc1.pri_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_mRNA_from_genomic.fna.gz -b actinopterygii_odb9 -g helixer-input.gtf | bash > HANNO-ESOLUC-V0.4+helixer.log 2>&1

#to also add RNAseq, just concatenate helixer and stringtie gtfs
cat ESOLUC.O+T+B.gtf helixer-input.gtf > ESOLUC.O+T+B+Helixer.gtf
../scripts/HANNO.v0.4.pl -t 80 -d HANNO-ESOLUC-V0.4-trans+helixer -a GCF_011004845.1_fEsoLuc1.pri_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_mRNA_from_genomic.fna.gz -b actinopterygii_odb9 -g ESOLUC.O+T+B+Helixer.gtf | bash > HANNO-ESOLUC-V0.4-trans+helixer.log 2>&1
```

### Merging gene-models from protein (miniprot) and mRNA(minimap2) evidence improves recall and precision
A short glimpse how HANNO performs with **protein only** input compared to **protein + corresponding mRNA** input (using 25 vertebrate species from above)

![image](https://github.com/HMPNK/HANNO/assets/51913753/1f4495d0-202d-40be-9686-17cc35bf4303)

Alignments of mRNAs strongly support recall of small terminal cds exons (start and stop codons next to or interrupted by splice sites), these are hard to find with protein alignment alone. This helping effect is limited by the divergence of the input dataset, as it is dependent on UTR alignability.

<p align="center">
<img src="https://github.com/HMPNK/HANNO/assets/51913753/62c4c214-7c7a-4fb3-8026-0d4752618053" width="75%" height="75%">
</p>

### HANNO compared to other pipelines
The _Gallus gallus_ dataset (proteins and transcriptomes, RefSeq annotation GCF_000002315.6_GRCg6a_genomic.gff for comparison) used in the [Braker3 manuscript](https://genome.cshlp.org/content/early/2024/05/28/gr.278090.123.abstract) was used with HANNO. HANNO was tested with different evidences (P = protein, R = mRNA corresponding to P, T = transcriptome and combinations thereof). Braker1/2/3 (including GeneMarkET/EP/ETP and using A = _ab initio_ predictions, P and T) were re-run for evaluation. A Helixer annotation was done on a Nvidia RTX6000 Ada GPU. The Braker3 ("braker.gtf") and GeneMarkETP ("genemark_supported.gtf") results were also piped through HANNO to accomodate for its gene selection by functional annotation. Values for Maker2 and Funnotate were taken directly from the publication. HANNO outperforms all methods in speed by a factor of at least 10-20X (same number of CPU threads, Intel(R) Xeon(R) CPU E7-8890 v4 @ 2.20GHz from 2017), while being in the top field of gene predictors regarding precision on the cds gene-level. Regarding annotation completeness, HANNO achieves the highest BUSCO scores and precision is inbetween GeneMarkETP and Braker3. It should be mentioned that HANNO produces respectable results on "same order mRNA-only" input (R), which no other method seems to be capable of. This makes HANNO promising for gene annotation with mRNA-like ISOseq data (to be tested).

![image](https://github.com/user-attachments/assets/49fcbffd-3c4b-446d-b017-15c6cc38ad20)

### INDEPENDENT TESTING BY COLLABORATORS
Martin Racoupeau and Christophe Klopp (INRAE, Toulouse, France) tested HANNO on a newly assembled Plant genome (245 Mbp, contig N50 22.25 Mbp) and compared results to Braker3 and Helixer. They used RNAseq data from the same species and protein data from related species.
To make HANNO run on their cluster, they had to change "source path/to/activate" with "conda activate" in the bash script.

![image](https://github.com/user-attachments/assets/221ba526-aff5-4089-9005-fb3c01003f28)



### UNDER DEVELOPMENT

* converter for bedDB to NCBI Refseq-like gtf and gff3 files
* filter some strange gene models that appear (introns spanning multiple gene models)
* improve annotation of tandem gene copies

### HISTORY AND ACKNOWLEDGEMENTS
This tool has been reaping in my mind over more than a decade and developed with every genome project I have been involved with. I thank all those colleagues who worked with me in those genome projects. I thank Martin Racoupeau and Christophe Klopp from INRAE, Toulouse for independent testing.
Special thanks go out to Heng Li for his ground breaking work in bioinformatic tools. His MINIPROT tool has replaced SPALN2 in prior versions of HANNO and made the whole pipeline much more efficient and easy-to-use.
Parts of this work were supported by my by DFG grant KU 3596/1-1; project number:
324050651).
