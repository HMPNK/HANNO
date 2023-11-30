$SCRIPTS/CDS_gtfToBed12 $1 | sort -k1,1V -k2,2n > $1.bed12
bedtools intersect -sorted -e -f $2 -F $2 -s -split -wao -a $1.bed12 -b $1.bed12| cut -f 4,16 | sort -u | awk '{if(o!=$1){i++;list[i]=$1}h[$1]=h[$1]";"$2;o=$1;} END{for(x=1;x<=i;x++){print list[x]"\t"h[list[x]]}}' >  $1.intersection.tsv
mawk -f $SCRIPTS/get_clusters.awk $1.intersection.tsv | sort -k2,2n > $1.model_clusters.tsv
awk -v infile=$1.bed12 'BEGIN{while(getline l < infile){split(l,a,"\t");s[a[4]]=a[5]}} {print $0"\t"s[$1];}' $1.model_clusters.tsv | sort -k2,2n -k3,3rn| awk '{if(o!=$2){print $1};o=$2}' > $1.bestmodel.tsv
grep -wF -f $1.bestmodel.tsv $1 > $1.clustered.gtf


$SCRIPTS/CDS_gtfToBed12 $1.clustered.gtf > $1.clustered.cds.bed12

bedtools getfasta -s -split -name -bed  $1.clustered.cds.bed12  -fi $3 -fo /dev/stdout | cut -f1 -d '(' | $SCRIPTS/TRANSLATE.sh > $1.clustered.cds.faa

gtfToGenePred $1.clustered.gtf $1.clustered.gp
gtfToGenePred $1.clustered.gtf $1.clustered.gp
genePredToBed $1.clustered.gp $1.clustered.bed12
