#!/usr/bin/env bash
##High-throughput Annotation Pipeline - HANNO
set -e
set -o pipefail

export SCRIPTS="/data2/HANNO/scripts";
MAMBA="$HOME/miniconda2/envs/MAMBA/bin";
EGGNOG="/data2/HANNO/EGGNOGG-DBs"
THREADS=24;

##VARIABLES: $1=genome.fasta, $2=workingdir, $3=protein.faa (from NCBI refSeq), $4=_rna_from_genomic.fna (from NCBI refSeq) $5=BUSCO-LINEAGE-DIR


##create working dir
mkdir -p $2
cd $2

##load reference from $1
ln -s ../$1 .
samtools faidx $1

##map proteins
date
ln -s ../$3 .
miniprot -t $THREADS -d $1.mpi $1
miniprot -Iu -t $THREADS --gtf $1.mpi $3 | gtfToGenePred stdin stdout | genePredToBed stdin stdout | awk -f $SCRIPTS/bed12ToGTF_addscore-tacoFake.awk > $3.stringtie-like.gtf

##merge protein mappings using TACO
date
source $MAMBA/activate TACO
ls $3.stringtie-like.gtf > $3.gtf.list
taco_run --ref-genome-fasta $1 -p $THREADS -o output $3.gtf.list
source $MAMBA/deactivate

##TEST output
date
bedtools getfasta -nameOnly -s -split -bed output/assembly.bed -fi $1 -fo stdout | $SCRIPTS/TRANSLATE.sh > assembly_CDS.faa
run_BUSCO.py -i assembly_CDS.faa -o BUSCO1 -l $5 -m proteins -c $THREADS
#rm -rf tmp run_BUSCO1

##map denovo TRANSCRIPTOME same genus
date
ln -s ../$4 .
minimap2 -I 100G -t $THREADS -x splice -a --junc-bed output/assembly.bed $1  $4 > $4.sam
samtools sort -o $4.bam $4.sam
rm -f $4.sam

##convert to bed12 and correct strands and check for bad transcripts ("exon length 0")
samtools view $4.bam | awk 'BEGIN{OFS="\t";FS="\t";} {counter[$1]++;$1=$1"_"counter[$1];n=split($0,a,"\t");printf $1"\t";for(x=1;x<=n;x++){if(substr(a[x],1,3)=="ts:"){printf a[x]}};printf "\n"}' | grep "ts:A:" | sed "s/ts:A://g" > $4.strands
bamToBed -split -bed12 -i $4.bam > $4.bed12
awk -v infile=$4.strands 'BEGIN{OFS="\t";FS="\t";while(getline l < infile){split(l,a,"\t");h[a[1]]=a[2]}} {if($10==1){frame="."} else{frame=$6};counter[$4]++;$4=$4"_"counter[$4];if(h[$4]=="-" && h[$4]!=""){if($6=="+"){$6="-"} else{$6="+"}};if(frame=="."){$6="."};print}' $4.bed12 | awk '{n=split($11,a,",");for(x=1;x<n;x++){if(a[x]==0){t=1} else{t=0}};if(t==0){print}}' > $4.stranded.bed12

##merge protein mappings + denovo transcritome using TACO
date

source $MAMBA/activate TACO
cat output/assembly.bed $4.stranded.bed12 | awk -f $SCRIPTS/bed12ToGTF_addscore-tacoFake.awk > prot-denovo.stringtie-like.gtf
ls prot-denovo.stringtie-like.gtf > $3.gtf.list2
taco_run --ref-genome-fasta $1 -p $THREADS -o output2 $3.gtf.list2
source $MAMBA/deactivate

##TRANSDECODER
date
gtf_to_alignment_gff3.pl output2/assembly.gtf > Merged.gff
gtf_genome_to_cdna_fasta.pl output2/assembly.gtf $1 > Merged.fa
TransDecoder.LongOrfs -S -m 80 -t Merged.fa

##
lastdb -P $THREADS -p PROTDB $3
lastal -P $THREADS -p BL80 -m100 -K 1 PROTDB Merged.fa.transdecoder_dir/longest_orfs.pep > ALLPROT.maf
maf-convert blasttab ALLPROT.maf| sort -k1,1V -k12,12rn| awk '{n=split($1,a,"\.");name=a[1];for(x=2;x<n;x++){name=name"."a[x]};print name"\t"$0}'| awk '{if(o!=$1){print};o=$1}' > ALLPROT.bestORF
awk '{if($13>=200){print $2 > "ALLPROT.bestORF.keeplist";print $1".p" > "ALLPROT.bestORF.remove"}}' ALLPROT.bestORF

##
awk '{if($3=="CDS") {split($0,d,"=");print $1"\t"$5-$4+1"\t"d[3]}}' Merged.fa.transdecoder_dir/longest_orfs.gff3 | sort -k1,1 -k2,2rn | awk '{if(oldid!=$1) {split($3,a,";");print a[1]};oldid=$1}' > transdecoder.keeplist
grep -vFf ALLPROT.bestORF.remove transdecoder.keeplist | cat - ALLPROT.bestORF.keeplist > transdecoder.keeplist2
grep -w -F -f transdecoder.keeplist2 Merged.fa.transdecoder_dir/longest_orfs.gff3 > transdecoder.gff3_BEST2

##shift CDS to ATG codons
grep -w CDS transdecoder.gff3_BEST2 | awk '{print $1"\t"$4-1"\t"$5}' > start.bed
bedtools getfasta -split -s -bed start.bed -fi Merged.fa -fo start.fa
awk '{if(substr($1,1,1)==">"){split($0,a,/[>:]/)} else{seq=toupper($1);for(x=1;x<=length(seq);x=x+3){if(substr(seq,x,3)=="ATG"){n=x;break}};print a[2]"\t"n-1; }}' start.fa > start.shift
awk -v infile=start.shift -v maxshiftpercent=20 'BEGIN{OFS="\t";FS="\t";while(getline l < infile){split(l,a,"\t");h[a[1]]=a[2]}} {if($3=="CDS" && 100*h[$1]/($5-$4-1)<=maxshiftpercent &&  h[$1]>0){$4=$4+h[$1];shift++};print} END{print "shifted "shift" start-codons!" > "/dev/stderr"}' transdecoder.gff3_BEST2 > transdecoder.gff3_BEST2.corrATG

##
cdna_alignment_orf_to_genome_orf.pl transdecoder.gff3_BEST2.corrATG Merged.gff Merged.fa > Merged.CDS2.gff3
awk -f  $SCRIPTS/TRANSDECODER-newversion-GFF3toGTF.awk Merged.CDS2.gff3 > Merged.CDS2.gtf
gtfToGenePred Merged.CDS2.gtf Merged.CDS2.gp
genePredToBed Merged.CDS2.gp Merged.CDS2.bed12
bash $SCRIPTS/CDS_gtfToBed12 Merged.CDS2.gtf > Merged.CDS2only.bed12
bedtools getfasta -split -nameOnly -s -bed Merged.CDS2only.bed12 -fi $1 -fo /dev/stdout | $SCRIPTS/TRANSLATE.sh > Merged.CDS2.faa
run_BUSCO.py -i Merged.CDS2.faa -o Merged.CDS2.BUSCO -l $5 -m proteins -c 40 
#rm tmp run_Merged.CDS* -rf

##functional annotation with EggNog use before best of cluster selection
date
#cp -r $EGGNOG /dev/shm/EGGNOGG-DBs
emapper.py --cpu $THREADS --mp_start_method forkserver --data_dir $EGGNOG -o out --output_dir ./ --temp_dir ./ --override -m diamond --dmnd_ignore_warnings -i Merged.CDS2.faa --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 --itype proteins --tax_scope auto --target_orthologs all --go_evidence non-electronic --pfam_realign none --report_orthologs --decorate_gff yes --excel  > emapper.out  2>emapper.err
#rm -rf /dev/shm/EGGNOGG-DBs

##functional annotation using custom PROTDB
date
lastal -P $THREADS -p BL80 -m 100 -K1 PROTDB Merged.CDS2.faa > Merged.CDS2_vs_PROTDB.maf
grep ">" $3 | awk '{l=length($1);print substr($0,2,l-1)"\t"substr($0,l+2)}' > PROT-DB.desc.txt
maf-convert tab Merged.CDS2_vs_PROTDB.maf| grep -v '^#'| sed "s/_frame.\t/\t/g" | sed "s/EG2=//g" | sort --buffer-size=128G --parallel=16 -k7,7V -k13,13g | awk '{if($7!=o && $13<=1){print};o=$7}'| awk -v infile=PROT-DB.desc.txt 'BEGIN{while(getline l<infile){split(l,a,"\t");h[a[1]]=a[2]}} {print $7"\t"h[$2]" | "$2"\t"$1" e-val: "$13}' > Merged.CDS2_vs_PROTDB.description.txt

##merge functional annotations
awk -v eggnog=out.emapper.annotations -v last=Merged.CDS2_vs_PROTDB.description.txt 'BEGIN{OFS="\t";FS="\t";while(getline l < eggnog){n=split(l,a,"\t");if(substr(l,1,1)!="#"){split(a[1],b,".");i++;nlist[i]=b[1];t[b[1]]=1;egg3[b[1]]=a[3];egg4[b[1]]=a[4];egg8[b[1]]=a[8];egg9[b[1]]=a[9];;egg11[b[1]]=a[11];}};while(getline l < last){n=split(l,a,"\t");if(substr(l,1,1)!="#"){split(a[1],b,".");if(t[b[1]]==""){i++;nlist[i]=b[1];};last2[b[1]]=a[2];last3[b[1]]=a[3]}};exit} END{for(x=1;x<=i;x++){y=nlist[x];$0=y"\t"egg9[y]"\t"egg8[y]"\t"egg11[y]"\t"egg3[y]"\t"egg4[y]"\t"last2[y]"\t"last3[y];print}}' | sort -k1,1V > final_description.txt

##BEST MODEL OF GENE
date
bash $SCRIPTS/BESTCLUSTERGTF-SEQOUT-e.sh Merged.CDS2.gtf 0.1 $1
run_BUSCO.py -i Merged.CDS2.gtf.clustered.cds.faa -o Merged.CDS3.BUSCO -l $5 -m proteins -c 40 
#rm tmp run_Merged.CDS* -rf

awk -v description=final_description.txt 'BEGIN{OFS="\t";FS="\t";while(getline l < description){split(l,a,"\t");egs[a[1]]=a[6];ege[a[1]]=a[5];egg[a[1]]=a[2];egd[a[1]]=substr(a[3],1,140)}} {split($4,a,".");$4=a[1]; if(egg[$4]!="" && egg[$4]!="-"){$4=$4" | "egg[$4]" | "egd[$4]" | eggnog score: "egs[$4]" evalue: "ege[$4]};print}' Merged.CDS2.bed12  | sed "s/ /_/g" >  final.all.eggnog.bed12
awk -v description=final_description.txt 'BEGIN{OFS="\t";FS="\t";while(getline l < description){split(l,a,"\t");gene[a[1]]=a[2];lastd[a[1]]=a[7];lasts[a[1]]=a[8]};laste[a[1]]=a[9]} {split($4,a,".");$4=a[1]; if(gene[$4]==""){gene[$4]="-"}; if(lasts[$4]!=""){$4=$4" | "gene[$4]" | "lastd[$4]" | last score: "lasts[$4]" "laste[$4]};print}' Merged.CDS2.bed12  | sed "s/ /_/g" >  final.all.last.bed12
awk -v description=final_description.txt 'BEGIN{OFS="\t";FS="\t";while(getline l < description){split(l,a,"\t");egs[a[1]]=a[6];ege[a[1]]=a[5];egg[a[1]]=a[2];egd[a[1]]=substr(a[3],1,140)}} {split($4,a,".");$4=a[1]; if(egg[$4]!="" && egg[$4]!="-"){$4=$4" | "egg[$4]" | "egd[$4]" | eggnog score: "egs[$4]" evalue: "ege[$4]};print}' Merged.CDS2.gtf.clustered.bed12  | sed "s/ /_/g" >  final.best.eggnog.bed12
awk -v description=final_description.txt 'BEGIN{OFS="\t";FS="\t";while(getline l < description){split(l,a,"\t");gene[a[1]]=a[2];lastd[a[1]]=a[7];lasts[a[1]]=a[8]};laste[a[1]]=a[9]} {split($4,a,".");$4=a[1];if(gene[$4]==""){gene[$4]="-"};if(lasts[$4]!=""){$4=$4" | "gene[$4]" | "lastd[$4]" | last score: "lasts[$4]" "laste[$4]};print}' Merged.CDS2.gtf.clustered.bed12  | sed "s/ /_/g" >  final.best.last.bed12

##clean-up
ls | grep -vE 'final|out.emapper.annotations' | xargs rm -rf

date
##END
