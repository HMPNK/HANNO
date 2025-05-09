#!/usr/bin/env perl

##PERL script to run HANNO-PIPE with different parameters

  use strict;
  use Getopt::Std;

##define path to HANNO scripts and MAMBA installation (depend on your system, change before first run!)
my $scr = "/data2/HANNO/scripts";
my $mamba="/home/$USER/miniconda3/envs/MAMBA/bin";
my $eggnog="/data2/HANNO/EGGNOGG-DBs";

##VARIABLES
my $asm = "";
my $dir = "hanno-out";
my $pro = "";
my $rna = "";
my $gtf = "";
my $odb = "";
my $pdb = "";
my $cpu = 8;
my $files = "";
my $gtf2 = "";
my $egg = 0;
my $ort = 0;
my $utr = "true";

##GET OPTIONS
my %options=();
getopts("a:d:p:P:m:b:g:t:r:E:B:U:", \%options);

if(! $ARGV[0] && ! $options{a})       {
       print STDERR "
HANNO version 0.5 (High-throughput ANNOtation for eukaryote genomes)
Author: Heiner Kuhl, Phd (heiner.kuhl\@igb-berlin.de)

THIS SCRIPT CREATES THE PIPELINE AS A BASH script

INPUT DATA: genome assembly, proteins and mRNAs from related organism or denovo transcriptome assemblies, gtf from stringtie reference guided transcriptome assembly, BUSCO database

ALWAYS USE RELATIVE PATH (e.g. \"../../assembly/asm.fasta\"), IF INPUT data is not in current directory!

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
		-U polish UTRs (default=true, no=false), recommended if using diverged mRNA evidence (will remove surplus UTR exons far fom CDS) 

This script generates a bash script for running the pipeline! Write script to file and run by: nohup bash <script> & !
";       };
die "ERROR: NEED an assembly to annotate!!!\n" if(! $options{a});
die "ERROR: PROVIDE AT LEAST A PROTEIN FASTA FILE FOR ANNOTATION BY HOMOLOGY!!!\n" if(! $options{p});
die "ERROR: unknown option $ARGV[0]" if $ARGV[0];

$asm = $options{a} if($options{a});
$dir = $options{d} if($options{d});
$pro = $options{p} if($options{p});
$rna = $options{r} if($options{r});
$gtf = $options{g} if($options{g});
$odb = $options{b} if($options{b});
$pdb = $options{P} if($options{P});
$cpu = $options{t} if($options{t});
$egg = $options{E} if($options{E});
$ort = $options{B} if($options{B});
$utr = $options{U} if($options{U});

##GET OPTIONS END
##CHECK FILE PATHs
die "\nERROR: Invalid (-a) assembly file $asm ... Check PATH ...Exiting.\n\n" if(! -e $asm);
die "\nERROR: Invalid (-p) protein file $pro ... Check PATH ...Exiting.\n\n" if(! -e $pro);
die "\nERROR: Invalid (-b) ODB files $odb ... Check PATH ...Exiting.\n\n" if(! -e $odb && $odb ne "");
die "\nERROR: Invalid (-r) mRNA file $rna ... Check PATH ...Exiting.\n\n" if(! -e $rna && $rna ne "");
die "\nERROR: Invalid (-g) GTF file $gtf ... Check PATH ...Exiting.\n\n" if(! -e $gtf && $gtf ne "");
die "\nERROR: Invalid (-P) protein file $pdb ... Check PATH ...Exiting.\n\n" if(! -e $pdb && $pdb ne "");
##CREATE HANNO assembly pipeline:
my $COMMAND = "#!/usr/bash\nset -e\n#set -o pipefail\n
export SCRIPTS=$scr

##RUN HANNO assembly pipeline:\n
mkdir $dir
cd $dir
##MAP PROTEINS in $pro
date
seqtk seq ../$asm > asm.fa
seqtk seq ../$pro > protein.input.fa
miniprot --outn=10 -p 0.5 --outc=0.5 --outs=0.5 -Iu -t $cpu --gtf asm.fa protein.input.fa | gtfToGenePred stdin stdout | genePredToBed stdin stdout | awk -f $scr/bed12ToGTF_addscore-tacoFake.awk > asm.prot.gtf
##MERGE USING STRINGTIE
stringtie --merge -o prot.merge.gtf asm.prot.gtf
##MERGE USING TACO
source $mamba/activate TACO
ls asm.prot.gtf > TACO.gtf.list
taco_run --max-isoforms 7 -p $cpu -o output1 TACO.gtf.list
source $mamba/deactivate
##GET BESTMODELS FROM CLUSTER OF PREDICTIONS
awk \'BEGIN{OFS=\"\\t\";FS=\"\\t\"} {if(\$3==\"exon\"){\$3=\"CDS\";print}}\' prot.merge.gtf > prot.merge1.cds.gtf
bash $scr/BESTCLUSTERGTF-SEQOUT-e.sh prot.merge1.cds.gtf 0.1 asm.fa
awk \'BEGIN{OFS=\"\\t\";FS=\"\\t\"} {if(\$3==\"exon\"){\$3=\"CDS\";print}}\' output1/assembly.gtf > prot.merge2.cds.gtf
bash $scr/BESTCLUSTERGTF-SEQOUT-e.sh prot.merge2.cds.gtf 0.1 asm.fa
cat prot.merge1.cds.gtf.clustered.gtf prot.merge2.cds.gtf.clustered.gtf > MODELS1.gtf 
cat prot.merge1.cds.gtf.clustered.bed12 prot.merge2.cds.gtf.clustered.bed12 > MODELS1.bed12
cat prot.merge1.cds.gtf.clustered.cds.faa prot.merge2.cds.gtf.clustered.cds.faa > MODELS1.faa
rm -f prot.merge1* prot.merge2*
";

##optional BUSCO
if($odb ne "" && $ort == 1)	{
$COMMAND = "$COMMAND
##TEST BY BUSCO using $odb
run_BUSCO.py -i MODELS1.faa -o BUSCO1 -l ../$odb -m proteins -c $cpu -f
";
		}

##MAP mRNAs
if($rna ne "")  {
$COMMAND = "$COMMAND
##MAP mRNAs in $rna		
seqtk seq ../$rna > rna.input.fa
awk \'BEGIN{OFS=\"\\t\";FS=OFS;}{\$3=\"exon\";print;}\' MODELS1.gtf | paftools.js gff2bed -j - > splicejunctions.bed
minimap2 -I 100G -t $cpu -x splice -a --junc-bed splicejunctions.bed asm.fa rna.input.fa > rna.input.sam
samtools sort -o rna.input.bam rna.input.sam
rm -f rna.input.sam
##convert to bed12 and correct strands and check for bad transcripts (\"exon length 0\")
samtools view rna.input.bam | awk \'BEGIN{OFS=\"\\t\";FS=\"\\t\";} {counter[\$1]++;\$1=\$1\"_\"counter[\$1];n=split(\$0,a,\"\\t\");printf \$1\"\\t\";for(x=1;x<=n;x++){if(substr(a[x],1,3)==\"ts:\"){printf a[x]}};printf \"\\n\"}\' | grep \"ts:A:\" | sed \"s/ts:A://g\" > rna.strands
##BUG FIX HERE: bamToBed output misses commas at end of field 11 and 12!!! Needs Correction by awk script, otherwise bed12ToGTF.awk or bed12ToGTF_addscore-tacoFake.awk will miss terminal exons, which leads to missing UTRs!
bamToBed -split -bed12 -i rna.input.bam | awk 'BEGIN{OFS=\"\\t\";FS=OFS;}{\$11=\$11\",\";\$12=\$12\",\";print}' > rna.bed12
awk -v infile=rna.strands \'BEGIN{OFS=\"\\t\";FS=\"\\t\";while(getline l < infile){split(l,a,\"\\t\");h[a[1]]=a[2]}} {if(\$10==1){frame=\".\"} else{frame=\$6};counter[\$4]++;\$4=\$4\"_\"counter[\$4];if(h[\$4]==\"-\" && h[\$4]!=\"\"){if(\$6==\"+\"){\$6=\"-\"} else{\$6=\"+\"}};if(frame==\".\"){\$6=\".\"};print}' rna.bed12 | awk \'{n=split(\$11,a,\",\");for(x=1;x<n;x++){if(a[x]==0){t=1} else{t=0}};if(t==0){print}}\' > rna.stranded.bed12
";
		}

if($gtf ne "")  {
		$gtf2 = "input.gtf";
$COMMAND = "$COMMAND
cp ../$gtf $gtf2
";
		}
else {$gtf2 = "";}

if($rna eq ""){$files = "MODELS1.bed12";} else{$files = "MODELS1.bed12 rna.stranded.bed12";}

$COMMAND = "$COMMAND
##merge protein mappings + denovo transcritome using TACO and stringtie (add gtf here if provided)
##TACO
source $mamba/activate TACO
cat $files | awk -f $scr/bed12ToGTF_addscore-tacoFake.awk > evidence.stringtie-like.gtf
ls evidence.stringtie-like.gtf $gtf2 > TACO.gtf.list2
taco_run --max-isoforms 7 -p $cpu -o output2 TACO.gtf.list2
source $mamba/deactivate
##Stringtie
stringtie --merge -o evidence.merge.gtf TACO.gtf.list2
stringtie --merge -G output2/assembly.gtf -o MODELS2.gtf evidence.merge.gtf output2/assembly.gtf
sed \"s/MSTRG./MSTRG/g\" MODELS2.gtf | awk \'BEGIN{OFS=\"\\t\";FS=OFS;}{gsub(/\\./,\"_\",\$9);gsub(\"_p\",\".p\",\$9);print}\' > MODELS2.gtf.t 
mv MODELS2.gtf.t MODELS2.gtf 
gtfToGenePred MODELS2.gtf stdout | genePredToBed stdin MODELS2.bed12
##Transdecoder ORFs
gtf_to_alignment_gff3.pl MODELS2.gtf > MODELS2.gff
gtf_genome_to_cdna_fasta.pl MODELS2.gtf asm.fa > MODELS2.mrna.fa
seqtk seq -l 0 MODELS2.mrna.fa | split -l 4000 -
ls x?? | awk \'{print \"TransDecoder.LongOrfs -S -m 80 -t \"\$1}\' | parallel -j $cpu > parallel.ORFs.log 2>&1
cat x??.transdecoder_dir/longest_orfs.pep  > MODELS2.longest_orfs.pep
cat x??.transdecoder_dir/longest_orfs.gff3  > MODELS2.longest_orfs.gff3
rm -rf x?? x??.transdecoder_dir
";

#SCREEN ORFs and TRANSFER ORFs to genome
if($pdb eq "") {$pdb = $pro;} 
$COMMAND = "$COMMAND
seqtk seq ../$pdb > proteins.reference.faa
";

$COMMAND = "$COMMAND
##score ORFs by last alignment to references
lastdb -P $cpu -p PROTDB proteins.reference.faa
lastal -P $cpu -p BL80 -m100 -K 1 PROTDB MODELS2.longest_orfs.pep > MODELS2.maf
maf-convert blasttab MODELS2.maf| sort -k1,1V -k12,12rn| awk \'{n=split(\$1,a,\"\.\");name=a[1];for(x=2;x<n;x++){name=name\".\"a[x]};print name\"\\t\"\$0}\'| awk \'{if(o!=\$1){print};o=\$1}\' > MODELS2.bestORF
touch MODELS2.bestORF.remove MODELS2.bestORF.keeplist
awk \'{if(\$13>=200){print \$2 > \"MODELS2.bestORF.keeplist\";print \$1\".p\" > \"MODELS2.bestORF.remove\"}}\' MODELS2.bestORF
##Create list of best ORFs
awk \'{if(\$3==\"CDS\") {split(\$0,d,\"=\");print \$1\"\\t\"\$5-\$4+1\"\\t\"d[3]}}\' MODELS2.longest_orfs.gff3 | sort -k1,1 -k2,2rn | awk \'{if(oldid!=\$1) {split(\$3,a,\";\");print a[1]};oldid=\$1}\' > MODELS2.transdecoder.keeplist
grep -vFf MODELS2.bestORF.remove MODELS2.transdecoder.keeplist | cat - MODELS2.bestORF.keeplist > MODELS2.transdecoder.keeplist2
grep -w -F -f MODELS2.transdecoder.keeplist2 MODELS2.longest_orfs.gff3 > MODELS2.longest_orfs.gff3_BEST2
##shift CDS to ATG codons
grep -w CDS MODELS2.longest_orfs.gff3_BEST2 | awk \'{print \$1\"\\t\"\$4-1\"\\t\"\$5}' > MODELS2.start.bed
bedtools getfasta -split -s -bed MODELS2.start.bed -fi MODELS2.mrna.fa -fo MODELS2.start.fa
awk \'{if(substr(\$1,1,1)==\">\"){split(\$0,a,/[>:]/)} else{seq=toupper(\$1);for(x=1;x<=length(seq);x=x+3){if(substr(seq,x,3)==\"ATG\"){n=x;break}};print a[2]\"\\t\"n-1; }}\' MODELS2.start.fa > MODELS2.start.shift
awk -v infile=MODELS2.start.shift -v maxshiftpercent=20 \'BEGIN{OFS=\"\\t\";FS=\"\\t\";while(getline l < infile){split(l,a,\"\\t\");h[a[1]]=a[2]}} {if(\$3==\"CDS\" && 100*h[\$1]/(\$5-\$4-1)<=maxshiftpercent &&  h[\$1]>0){\$4=\$4+h[\$1];shift++};print} END{print \"shifted \"shift\" start-codons!\" > \"/dev/stderr\"}\' MODELS2.longest_orfs.gff3_BEST2 > MODELS2.longest_orfs.gff3_BEST2.corrATG
##transfer ORF coordinates to genome
cdna_alignment_orf_to_genome_orf.pl MODELS2.longest_orfs.gff3_BEST2.corrATG MODELS2.gff MODELS2.mrna.fa > MODELS2.CDS2.gff3
awk -f  $scr/TRANSDECODER-newversion-GFF3toGTF.awk MODELS2.CDS2.gff3 > MODELS2.CDS2.gtf
gtfToGenePred MODELS2.CDS2.gtf MODELS2.CDS2.gp
genePredToBed MODELS2.CDS2.gp MODELS2.CDS2.bed12
bash $scr/CDS_gtfToBed12 MODELS2.CDS2.gtf > MODELS2.CDS2only.bed12
bedtools getfasta -split -name -s -bed MODELS2.CDS2.bed12 -fi asm.fa -fo /dev/stdout > MODELS2.CDS2.mRNA.fa
bedtools getfasta -split -name -s -bed MODELS2.CDS2only.bed12 -fi asm.fa -fo /dev/stdout > MODELS2.CDS2.fa
$scr/TRANSLATE.sh MODELS2.CDS2.fa > MODELS2.CDS2.faa
##CLEAN-UP
ls MODELS2* | grep -v MODELS2.CDS2 | xargs rm
";		

##optional BUSCO2
if($odb ne "")  {
$COMMAND = "$COMMAND
##TEST BY BUSCO using $odb
run_BUSCO.py -i MODELS2.CDS2.faa -o BUSCO2 -l ../$odb -m proteins -c $cpu -f
mv run_BUSCO2/full_table_BUSCO2.tsv ALLMODELS.BUSCO.tsv
";
                }

$COMMAND = "$COMMAND
##functional annotation using custom PROTDB
lastal -P $cpu -p BL80 -m 100 -K1 PROTDB MODELS2.CDS2.faa > MODELS2.CDS2_vs_PROTDB.maf
grep \">\" proteins.reference.faa | awk \'{l=length(\$1);print substr(\$0,2,l-1)\"\\t\"substr(\$0,l+2)}\' > PROT-DB.desc.txt
maf-convert tab MODELS2.CDS2_vs_PROTDB.maf| grep -v \'^#\'| sed \"s/_frame.\\t/\\t\/g\" | sed \"s/EG2=//g\" | sort --buffer-size=128G --parallel=16 -k7,7V -k13,13g | awk \'{if(\$7!=o && \$13<=1){print};o=\$7}\' | awk -v infile=PROT-DB.desc.txt \'BEGIN{while(getline l<infile){split(l,a,\"\\t\");h[a[1]]=a[2]}} {print \$7\"\\t\"h[\$2]\" | \"\$2\"\\t\"\$1\" e-val: \"\$13}' > MODELS2.CDS2_vs_PROTDB.description.txt
";

if($egg eq 0)  {
$COMMAND = "$COMMAND
##functional annotation with EggNog
emapper.py --cpu $cpu --mp_start_method forkserver --data_dir $eggnog -o out --output_dir ./ --temp_dir ./ --override -m diamond --dmnd_ignore_warnings -i MODELS2.CDS2.faa --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 --itype proteins --tax_scope auto --target_orthologs all --go_evidence non-electronic --pfam_realign none --report_orthologs --decorate_gff yes --excel  > emapper.out  2>emapper.err
mv out.emapper.annotations ALLMODELS.eggnog.description.txt
";
		}

$COMMAND = "$COMMAND
##Final Clean-UP
mv MODELS2.CDS2.bed12 ALLMODELS.bed12
mv MODELS2.CDS2.faa ALLMODELS.faa
mv MODELS2.CDS2.fa ALLMODELS.cds.fa
mv MODELS2.CDS2.mRNA.fa ALLMODELS.mRNA.fa
mv MODELS2.CDS2_vs_PROTDB.description.txt ALLMODELS.lastp.description.txt
";

if($utr eq "true") {
$COMMAND = "$COMMAND
#POLISH UTRs
mawk -f $scr/clean-bed12-utrs.awk ALLMODELS.bed12 > ALLMODELS.bed12.polishUTR
mv ALLMODELS.bed12.polishUTR ALLMODELS.bed12
";
                }

$COMMAND = "$COMMAND
#CLUSTER MODELS BY CDS
$scr/CLUSTERBYCDS.sh ALLMODELS.bed12 0.1
ls | grep -Ev \'asm.fa|ALLMODELS.BUSCO.tsv|ALLMODELS.bed12|ALLMODELS.faa|ALLMODELS.mRNA.fa|ALLMODELS.cds.fa|ALLMODELS.lastp.description.txt|ALLMODELS.eggnog.description.txt\' | xargs rm -rf
";

$COMMAND = "$COMMAND
##create Final DB that contains all results
bash $scr/CREATE-DB.sh ALLMODELS.bed12 ALLMODELS.bed12.model_clusters.tsv ALLMODELS.lastp.description.txt ALLMODELS.eggnog.description.txt ALLMODELS.BUSCO.tsv > ALLMODELS-FINAL.bedDB
rm -f ALLMODELS.bed12 ALLMODELS.bed12.model_clusters.tsv ALLMODELS.cds.fa ALLMODELS.eggnog.description.txt ALLMODELS.faa ALLMODELS.lastp.description.txt ALLMODELS.mRNA.fa ALLMODELS.BUSCO.tsv orf.length

sort -k14,14n -k5,5rn ALLMODELS-FINAL.bedDB | awk \'{if(o!=\$14){print};o=\$14}\' > BESTMODELS-FINAL.bedDB
";

##summarize BUSCO from BESTMODELS-FINAL.bedDB
if($odb ne "")  {
$COMMAND = "$COMMAND
OUTPUT=\$(ls ../$odb/hmms/ | wc -l)
cut -f 39-42 BESTMODELS-FINAL.bedDB | grep -v '\\-' | grep  -E 'Complete|Duplicated|busco1' | cut -f 1 | sort | uniq    | awk -v max=\$OUTPUT 'BEGIN{i=0} {i++} END{print \"\\nAnalysis of BUSCOs ( n=\"max\" ) in BESTMODELS-FINAL.bedDB\\n\C: n=\"i-1\" / \"100*(i-1)/max\" percent\"}'
cut -f 39-42 BESTMODELS-FINAL.bedDB | grep -v '\\-' | grep  -E 'Complete|Duplicated|busco1' | cut -f 1 | sort | uniq -u | awk -v max=\$OUTPUT 'BEGIN{i=0} {i++} END{print \"S: n=\"i-1\" / \"100*(i-1)/max\" percent\"}'
cut -f 39-42 BESTMODELS-FINAL.bedDB | grep -v '\\-' | grep  -E 'Complete|Duplicated|busco1' | cut -f 1 | sort | uniq -d | awk -v max=\$OUTPUT 'BEGIN{i=0} {i++} END{print \"D: n=\"i\" / \"100*i/max\" percent\"}'
cut -f 39-42 BESTMODELS-FINAL.bedDB | grep -v '\\-' | grep -Ev 'Complete|Duplicated'        | cut -f 1 | sort | uniq    | awk -v max=\$OUTPUT 'BEGIN{i=0} {i++} END{print \"F: n=\"i-1\" / \"100*(i-1)/max\" percent\"}'
cut -f 39-42 BESTMODELS-FINAL.bedDB | grep -v '\\-' | grep  -E 'Complete|Duplicated|busco1|Fragmented' | cut -f 1 | sort | uniq | awk -v max=\$OUTPUT 'BEGIN{i=0} {i++} END{print \"M: n=\"max-(i-1)\" / \"100-100*(i-1)/max\" percent\\n\"}'
";
                }

$COMMAND="$COMMAND
##OUTPUT fasta files of mRNA, CDS and aminoacid, output GFF3 and GTF
python $scr/bedDB_to_gtf_gff.py BESTMODELS-FINAL.bedDB --fasta --genome asm.fa --gff -o BESTMODELS-FINAL.gff3
mv BESTMODELS-FINAL.gff3 BESTMODELS-FINAL-including-mRNA+CDS.gff3
python $scr/bedDB_to_gtf_gff.py BESTMODELS-FINAL.bedDB --gff -o BESTMODELS-FINAL.gff3
python $scr/bedDB_to_gtf_gff.py BESTMODELS-FINAL.bedDB -o BESTMODELS-FINAL.gtf
rm asm.fa asm.fa.fai

date
##END HANNO pipeline
";
##CREATE HANNO pipeline END


##OUTPUT HANNO PIPELINE SCRIPT
print "$COMMAND";
