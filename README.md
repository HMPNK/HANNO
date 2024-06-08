

# HANNO: efficient High-throughput ANNOtation of protein coding genes in eukaryote genomes
* high speed, can annotate a large vertebrate genome in below 1h  
* works well with medium diverged input evidence (i.e. using NCBI Refseq annotations of species that diverged less than 50 MYA ), no need for RNAseq from your organism (but still beneficial of course)  
* transcript models include UTRs  
* completely evidence based (no weird gene models from *ab initio* prediction)
* includes functional annotation using eggNog and protein homology

![image](https://github.com/HMPNK/HANNO/assets/51913753/1ffd4d02-f148-4214-90ad-628ce828d050)

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
mamba activate HANNO
bash TTN-TEST-RUNS.sh

##This will run all use cases of the pipeline on the largest known vertebrate gene TTN
##CHECK the logs for early stops
##If the pipeline finishes successfully all logs will start and end with a date
##Do not wonder that BUSCO results are 0%, the TTN gene tested is not present in the BUSCO db used. 
```

### Output files
ALLMODELS-FINAL.bedDB -> bed12 like Database of all transcript models and all assigned functional information   
BESTMODELS-FINAL.bedDB -> putative best scoring transcript models for each transcript cluster from above
BESTMODELS-FINAL.gtf -> same as above as gtf format   
BESTMODELS-FINAL.mRNA.fa -> corresponding mRNA sequences as fasta
BESTMODELS-FINAL.CDS.fa -> corresponding CDS sequences as fasta   
BESTMODELS-FINAL.AA.faa -> corresponding aminoacid sequences as fasta   

### BENCHMARKING HANNO BY TEST-RUNS ON VERTEBRATE GENOMES

HANNO benchmark runtimes on different vertebrate clades (HPC-server: Intel(R) Xeon(R) CPU E7-8890 v4 @ 2.20GHz; 96-threads; 1TB RAM). Note that 6 (Birds, Amphibians), 8 (Fish) and 11 (Mammals) genome annotations were run in parallel.

![image](https://github.com/HMPNK/HANNO/assets/51913753/eeb80ea5-dc7b-4a29-9b87-8bada1712076)

HANNO was tested on fish (n=8), amphibian (n=6), bird (n=6) and mammal genomes (n=11), by transferring Refseq Annotations (from _P. flavescens_, _B. bufo_, _T. guttata_ and _H. sapiens_) to closely related and diverged species. Results below show the trends for "Complete", "Fragmented" and "Missing" BUSCOs for species ordered by divergence time from reference species (confidence intervalls from www.timetree.org) as well as total recovery of CDS and UTR compared to the 0 MYA diverged reference genome. HANNOtations of genomes by reference proteins and mRNA that diverged less than 50 MYA are typically yielding good results in terms of total recovered CDS and UTR sequences and BUSCO statistics for all vertebrate clades tested.

![image](https://github.com/HMPNK/HANNO/assets/51913753/97d72096-64d2-4624-bda3-5f6f1f338772)

![image](https://github.com/HMPNK/HANNO/assets/51913753/784a1528-f595-4a7e-8b55-30d35e4780e0)

![image](https://github.com/HMPNK/HANNO/assets/51913753/0f3ebb3b-51c7-47d3-b1ab-b5c390d5f16e)

![image](https://github.com/HMPNK/HANNO/assets/51913753/4bbe958d-e508-4239-8442-464fbfdf4446)


```sh

#test-runs on diverged vertebrate genomes

#AMPHIBIANS                 
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-BUFVIR-v0.4 -a GCA_037900795.1_ASM3790079v1_genomic.fna.gz -p GCF_905171765.1_aBufBuf1.1_protein.faa.gz -r GCF_905171765.1_aBufBuf1.1_rna_from_genomic.fna.gz -b vertebrata_odb9 | bash > HANNO-BUFVIR-v0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-LEPFUS-v0.4 -a GCA_031893025.1_aLepFus1.hap1_genomic.fna.gz -p GCF_905171765.1_aBufBuf1.1_protein.faa.gz -r GCF_905171765.1_aBufBuf1.1_rna_from_genomic.fna.gz -b vertebrata_odb9 | bash > HANNO-LEPFUS-v0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-XENTRO-v0.4 -a GCF_000004195.4_UCB_Xtro_10.0_genomic.fna.gz -p GCF_905171765.1_aBufBuf1.1_protein.faa.gz -r GCF_905171765.1_aBufBuf1.1_rna_from_genomic.fna.gz -b vertebrata_odb9 | bash > HANNO-XENTROP-v0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-BUFGAR-v0.4 -a GCF_014858855.1_ASM1485885v1_genomic.fna.gz -p GCF_905171765.1_aBufBuf1.1_protein.faa.gz -r GCF_905171765.1_aBufBuf1.1_rna_from_genomic.fna.gz -b vertebrata_odb9 | bash > HANNO-BUFGAR-v0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-SPEBOM-v0.4 -a GCF_027358695.1_aSpeBom1.2.pri_genomic.fna.gz -p GCF_905171765.1_aBufBuf1.1_protein.faa.gz -r GCF_905171765.1_aBufBuf1.1_rna_from_genomic.fna.gz -b vertebrata_odb9 | bash > HANNO-SPEBOM-v0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-BOMBOM-v0.4 -a GCF_027579735.1_aBomBom1.pri_genomic.fna.gz -p GCF_905171765.1_aBufBuf1.1_protein.faa.gz -r GCF_905171765.1_aBufBuf1.1_rna_from_genomic.fna.gz -b vertebrata_odb9 | bash > HANNO-BOMBOM-v0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-BUFBUF-v0.4 -a GCF_905171765.1_aBufBuf1.1_genomic.fna.gz -p GCF_905171765.1_aBufBuf1.1_protein.faa.gz -r GCF_905171765.1_aBufBuf1.1_rna_from_genomic.fna.gz -b vertebrata_odb9 | bash > HANNO-BUFBUF-v0.4.log 2>&1
                 
#BIRDS                 
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-TAEGUT-V0.4 -a GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz -p GCF_003957565.2_bTaeGut1.4.pri_protein.faa.gz -r GCF_003957565.2_bTaeGut1.4.pri_rna_from_genomic.fna.gz -b aves_odb9 | bash > HANNO-TAEGUT-V0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-SERCAN-V0.4 -a GCF_022539315.1_serCan2020_genomic.fna.gz -p GCF_003957565.2_bTaeGut1.4.pri_protein.faa.gz -r GCF_003957565.2_bTaeGut1.4.pri_rna_from_genomic.fna.gz -b aves_odb9 | bash > HANNO-SERCAN-V0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-MYICAY-V0.4 -a GCF_022539395.1_myiCay2020_genomic.fna.gz -p GCF_003957565.2_bTaeGut1.4.pri_protein.faa.gz -r GCF_003957565.2_bTaeGut1.4.pri_rna_from_genomic.fna.gz -b aves_odb9 | bash > HANNO-MYICAY-V0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-PHASUP-V0.4 -a GCA_023637945.1_DD_ASM_B1_genomic.fna.gz -p GCF_003957565.2_bTaeGut1.4.pri_protein.faa.gz -r GCF_003957565.2_bTaeGut1.4.pri_rna_from_genomic.fna.gz -b aves_odb9 | bash > HANNO-PHASUP-V0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-GALGAL-V0.4 -a GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz -p GCF_003957565.2_bTaeGut1.4.pri_protein.faa.gz -r GCF_003957565.2_bTaeGut1.4.pri_rna_from_genomic.fna.gz -b aves_odb9 | bash > HANNO-GALGAL-V0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 12 -d HANNO-DROHOL-V0.4 -a GCA_016128335.2_ZJU2.0_genomic.fna.gz -p GCF_003957565.2_bTaeGut1.4.pri_protein.faa.gz -r GCF_003957565.2_bTaeGut1.4.pri_rna_from_genomic.fna.gz -b aves_odb9 | bash > HANNO-DROHOL-V0.4.log 2>&1
                 
#FISHES                 
../scripts/HANNO.v0.4.pl -t 8 -d HANNO-PERFLA-V0.4 -a GCF_004354835.1_PFLA_1.0_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_rna_from_genomic.fna.gz -b actinopterygii_odb9 | bash > HANNO-PERFLA-V0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 8 -d HANNO-SANVIT-V0.4 -a GCA_031162955.1_sanVit1_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_rna_from_genomic.fna.gz -b actinopterygii_odb9 | bash > HANNO-SANVIT-V0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 8 -d HANNO-EPIFUS-V0.4 -a GCF_011397635.1_E.fuscoguttatus.final_Chr_v1_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_rna_from_genomic.fna.gz -b actinopterygii_odb9 | bash > HANNO-EPIFUS-V0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 8 -d HANNO-DICLAB-V0.4 -a GCF_905237075.1_dlabrax2021_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_rna_from_genomic.fna.gz -b actinopterygii_odb9 | bash > HANNO-DICLAB-V0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 8 -d HANNO-ORENIL-V0.4 -a GCA_922820385.1_Nile_Tilapia_GIFT_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_rna_from_genomic.fna.gz -b actinopterygii_odb9 | bash > HANNO-ORENIL-V0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 8 -d HANNO-ESOLUC-V0.4 -a GCF_011004845.1_fEsoLuc1.pri_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_rna_from_genomic.fna.gz -b actinopterygii_odb9 | bash > HANNO-ESOLUC-V0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 8 -d HANNO-CLAGAR-V0.4 -a GCF_024256425.1_CGAR_prim_01v2_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_rna_from_genomic.fna.gz -b actinopterygii_odb9 | bash > HANNO-CLAGAR-V0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 8 -d HANNO-DANRER-V0.4 -a GCA_033170195.1_ASM3317019v1_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa.gz -r GCF_004354835.1_PFLA_1.0_rna_from_genomic.fna.gz -b actinopterygii_odb9 | bash > HANNO-DANRER-V0.4.log 2>&1
                 
#MAMMALS                 
../scripts/HANNO.v0.4.pl -t 6 -d GRCm39.0.4 -a GCF_000001635.27_GRCm39_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz -b mammalia_odb10 | bash > GRCm39.0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d ARS1.2.0.4 -a GCF_001704415.2_ARS1.2_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz -b mammalia_odb10 | bash > ARS1.2.0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d phaCin.0.4 -a GCF_002099425.1_phaCin_unsw_v4.1_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz -b mammalia_odb10 | bash > phaCin.0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d ARS-UCD2.0.0.4 -a GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz -b mammalia_odb10 | bash > ARS-UCD2.0.0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d RhiFer1.0.4 -a GCF_004115265.2_mRhiFer1_v1.p_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz -b mammalia_odb10 | bash > RhiFer1.0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d T2T-CHM13v2.0.0.4 -a GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz -b mammalia_odb10 | bash > T2T-CHM13v2.0.0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d CalJa1.0.4 -a GCF_011100555.1_mCalJa1.2.pat.X_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz -b mammalia_odb10 | bash > CalJa1.0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d talOcc4v2.1.0.4 -a GCF_014898055.3_MPIMG_talOcc4v2.1_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz -b mammalia_odb10 | bash > talOcc4v2.1.0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d EleMax1.0.4 -a GCF_024166365.1_mEleMax1_primary_haplotype_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz -b mammalia_odb10 | bash > EleMax1.0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d PonAbe1.0.4 -a GCF_028885655.1_NHGRI_mPonAbe1-v1.1-hic.freeze_pri_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz -b mammalia_odb10 | bash > PonAbe1.0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d PanPan.0.4 -a GCF_029289425.1_NHGRI_mPanPan1-v1.1-0.1.freeze_pri_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz -b mammalia_odb10 | bash > PanPan.0.4.log 2>&1
../scripts/HANNO.v0.4.pl -t 6 -d OrcOrc.0.4 -a GCF_937001465.1_mOrcOrc1.1_genomic.fna.gz -p GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz -r GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz -b mammalia_odb10 | bash > OrcOrc.0.4.log 2>&1

```

### IMPROVING ANNOTATION BY ADDING TRANSCRIPTOME DATA
To benchmark improvement of annotation by adding RNAseq the _E. lucius_ run from above, was supported by _E. lucius_ Brain, Ovary and Testis RNAseq data:

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

#Merge assembled trancript gtf-files priot feeding it into HANNO:
cat ESOLUC.OVARY.gtf ESOLUC.TESTIS.gtf ESOLUC.BRAIN.gtf > ESOLUC.O+T+B.gtf

#Run HANNO with diverged proteins and mRNA and assembled transcript from same species:
../scripts/HANNO.v0.4.pl -t 80 -d HANNO-ESOLUC-V0.4-TRANS -a GCF_011004845.1_fEsoLuc1.pri_genomic.fna.gz -p GCF_004354835.1_PFLA_1.0_protein.faa -r GCF_004354835.1_PFLA_1.0_rna_from_genomic.fna -b ../../../home/osboxes/BUSCOVM/lineages/actinopterygii_odb9 -g ESOLUC.O+T+B.gtf | bash > HANNO-ESOLUC-V0.4-TRANS.log 2>&1

```
