#!/usr/bin/awk -f
{
if($3=="mRNA") {i++;score=$6;;};
if($3=="exon" || $3=="CDS")
{
split($0,d,"\t");
x=split(d[9],n,/[=;,]/);
#print n[1];print n[2];print n[3];print n[4];print n[5];print n[6];print n[7];print n[8];
print d[1]"\t"d[2]"\t"d[3]"\t"d[4]"\t"d[5]"\t"score"\t"d[7]"\t"d[8]"\tgene_id \"helixG"i"\"; transcript_id \"helixT"i"\";";

}



}
