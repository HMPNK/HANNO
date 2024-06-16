#!/usr/bin/awk -f
{
if($3=="mRNA") {i++;score=$6;;};

if($3=="exon" || $3=="CDS")
  {
  split($0,d,"\t");
  x=split(d[9],n,/[=;,]/);
  print d[1]"\t"d[2]"\t"d[3]"\t"d[4]"\t"d[5]"\t"score"\t"d[7]"\t"d[8]"\tgene_id \"helixG"i"\"; transcript_id \"helixT"i"\";";
  }
}
