#a better solution by H.Kuhl
BEGIN{OFS="\t";FS="\t";}
{

##BED12 to GENEPRED:
n=split($0,a,"\t");
if(n!=12){print "No BED12, expected 12 fields...Exit";exit};
$0="";
#switch fields
$1=a[4];$2=a[1];$3=a[6];$4=a[2];$5=a[3];$6=a[7];$7=a[8];$8=a[10];
score=a[5];
color=a[9];
#recalc exon boundaries
n1=split(a[11],b,",");

n2=split(a[12],c,",");
if(n1!=n2) {print "Problem with number of exons in BED12" ;exit}

#convert exon coords
$9=$4;
for(x=2;x<=n1;x++){es=$4+c[x];$9=$9","es};
$9=$9",";

#convert cds coords and set score to cd length
$10=$4+b[1];
for(x=2;x<=n2;x++){ee=$4+c[x]+b[x];$10=$10","ee;};
$10=$10",";

#do not need score here:
score=".";

delete a; delete b; delete c;
##GENEPREDTO GTF

i=0;j=1;frame=0;frameout=0;cdslen=0;cstart=0;cdscount=0

#read exons
n1=split($9,Estart,",");
n2=split($10,Eend,",");
if(n1!=n2){print "exon starts/ends count not matching";exit;}
merge=$9 $10;
n3=split(merge,allE,",");
asort(allE);

#create CDS from allE
i++;allC[i]=$6
for(x=1;x<n3;x++){if(allE[x]>$6 && allE[x]<$7){i++;;if(i==1){cstart=x};allC[i]=allE[x]}};
i++;allC[i]=$7

cend=cstart+i

if($3=="+"){
for(x=1;x<n3;x=x+2) 	{
			print $2"\tModel\texon\t"allE[x]+1"\t"allE[x+1]"\t"score"\t"$3"\t.\tgene_id \""$1"\"; transcript_id \""$1"\";"
			if(x>=cstart && x<=cstart+i)	{
							cdscount++;
									if(cdscount>1){frame=cdslen-int(cdslen/3)*3;if(frame==3){frame=0};if(frame==1){frameout=2} else if(frame==2){frameout=1} else{frameout=0}}
									cdslen=cdslen+allC[j+1]-allC[j];
									print $2"\tModel\tCDS\t"allC[j]+1"\t"allC[j+1]"\t"score"\t"$3"\t"frameout"\tgene_id \""$1"\"; transcript_id \""$1"\";";
							j=j+2
							}
			};
}

else if($3=="-"){
j=cend;
for(x=n3-1;x>0;x=x-2)     {
                        print $2"\tModel\texon\t"allE[x-1]+1"\t"allE[x]"\t"score"\t"$3"\t.\tgene_id \""$1"\"; transcript_id \""$1"\";"
                        if(x>=cstart && x<=cstart+i)	{
							cdscount++;
									if(cdscount>1){frame=cdslen-int(cdslen/3)*3;if(frame==3){frame=0};if(frame==1){frameout=2} else if(frame==2){frameout=1} else{frameout=0}}
									cdslen=cdslen+allC[j]-allC[j-1];
									print $2"\tModel\tCDS\t"allC[j-1]+1"\t"allC[j]"\t"score"\t"$3"\t"frameout"\tgene_id \""$1"\"; transcript_id \""$1"\";";
							j=j-2
							}
                        };
}
}
