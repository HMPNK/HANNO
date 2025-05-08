#simply remove UTRs not overlapping with CDS, will keep only one 5'UTR and one 3'UTR (helps cleaning annotation, if diverged mRNAs are used!)

BEGIN{
OFS="\t";
FS="\t";
#set max allowed distance of terminal UTR exons from CDS-start/end
if(cdsdist==""){cdsdist=300};
}

{
if(substr($1,1,1)=="#"){print;}
else{
ac=split($11,a,",");
bc=split($12,b,",");
cstart=$7;
cend=$8;
#analyse exons with cds overlap (maybe add distance parameter later)
for(x=1;x<=$10;x++)	{
				estart=$2+b[x];
				eend=$2+b[x]+a[x];
#			printf estart"\t"eend
			if((estart<=cstart-cdsdist && eend<=cstart-cdsdist) || (estart>cend+cdsdist && eend>cend+cdsdist) ){rem[x]=1;} else {}
			}

#change bed12
for(x=1;x<=$10;x++)     {
			if(rem[x]!=1) {exoncount++;field11=field11 a[x]",";field12=field12 b[x]",";}
			}
#correct transcript start/end
ac=split(field11,a,",");
bc=split(field12,b,",");
newstart=$2+b[1];
newend=$2+b[ac-1]+a[bc-1];
for(x=1;x<=bc-1;x++){field12new=field12new b[x]-b[1]",";}
#print "\n"newstart"\t"newend"\t"exoncount"\t"field11"\t"field12new;

#change bed12 and output
$2=newstart;
$3=newend;
$10=bc-1;
$11=field11;
$12=field12new;
print;

#printf "\n";
delete rem;exoncount="";field11="";field12="";field12new="";
}
}
