awk -f $SCRIPTS/bed12ToGTF.awk $1 | $SCRIPTS/CDS_gtfToBed12 | sort -k1,1V -k2,2n > $1.cds.bed12
#USE bedtools version 2.27.1 !!! Had strange issues if using v2.31.1 ....
bedtools intersect -sorted -e -f $2 -F $2 -s -split -wao -a $1.cds.bed12 -b $1.cds.bed12 | cut -f 4,16 | sort -u | awk '{if(o!=$1){i++;list[i]=$1}h[$1]=h[$1]";"$2;o=$1;} END{for(x=1;x<=i;x++){print list[x]"\t"h[list[x]]}}' >  $1.intersection.tsv
mawk -f $SCRIPTS/get_clusters.awk $1.intersection.tsv | sort -k2,2n | sort -u > $1.model_clusters.tsv
rm -f $1.cds.bed12 $1.intersection.tsv
