{
if($3=="exon" || $3=="CDS") 	{
				split($0,d,"\t");
				split(d[9],name,";");

				gsub("ID=cds.","",name[1]);
				gsub("ID=","",name[1]);
				count=split(name[1],name2,".");
				name2[3]=name2[1]"."name2[2];

				if(d[3]=="exon"){
				print d[1]"\t"d[2]"\t"d[3]"\t"d[4]"\t"d[5]"\t"d[6]"\t"d[7]"\t"d[8]"\tgene_id \""name2[3]"\"; transcript_id \""name2[3]"\";";}
				else if(d[3]=="CDS"){
				print d[1]"\t"d[2]"\t"d[3]"\t"d[4]"\t"d[5]"\t"d[6]"\t"d[7]"\t"d[8]"\tgene_id \""name2[3]"\"; transcript_id \""name2[3]"\";";}



				}





}
