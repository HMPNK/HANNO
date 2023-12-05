#create input data if not already existing (prevents from awk outputting error)
touch ALLMODELS.bed12.model_clusters.tsv; touch ALLMODELS.lastp.description.txt; touch ALLMODELS.eggnog.description.txt; touch ALLMODELS.BUSCO.tsv
#get cds length of each model
cut -f 1-12 $1 | awk -f $SCRIPTS/bed12ToGTF_addscore.awk | awk '{if($3=="CDS"){split($10,a,"\"");print a[2]"\t"$6}}' | uniq -c | awk '{print $2"\t"$3"\t"$1}' > orf.length
#sort by sequence id and coordinates
sort -k1,1V -k2,2n $1 | \
awk 'BEGIN{FS="\t";OFS=FS;}{$13=$4;print}' | \
#add orf.length
awk -v infile=orf.length 'BEGIN{FS="\t";OFS=FS;while(getline l < infile){n=split(l,a,FS);h[a[1]]=a[2]}} {x=$4;$5=h[x];print}' | \
#add CDS cluster info
awk -v infile=$2 'BEGIN{FS="\t";OFS=FS;while(getline l < infile){n=split(l,a,FS);h[a[1]]=a[2]}} {x=$13;if(h[x]!=""){if(t[h[x]]==""){i++;cluster[h[x]]=i;t[h[x]]++};$14=cluster[h[x]];print} else{$14="-";print}}' | \
#add lastp info
awk -v infile=$3 'BEGIN{FS="\t";OFS=FS;while(getline l < infile){gsub(" e-val: ","\t",l);gsub(/ \| /,"\t",l);n=split(l,a,FS);gsub(/\([\+\-]\)/,"",a[1]);h[a[1]]=a[2] FS a[3] FS a[4] FS a[5];}} {x=$13;if(h[x]!=""){$15=h[x];print} else{$15="-";for(x=1;x<=3;x++){$15=$15 FS "-"};print}}' | \
#add eggnog info
awk -v infile=$4 'BEGIN{FS="\t";OFS=FS; while(getline l < infile){split(l,a,FS); gsub(/\([\+\-]\)/,"",a[1]); h[a[1]]=l;}} {x=$13; if(h[x]!=""){$19=h[x]; print} else {$19="-";for(x=1;x<=21-1;x++){$19=$19 FS "-"};print} }' | \
#add busco info
awk -v infile=$5 'BEGIN{FS="\t";OFS=FS; while(getline l < infile){split(l,a,FS); gsub(/\([\+\-]\)/,"",a[3]); h[a[3]]=l;}} {x=$13; if(h[x]!=""){$40=h[x]; print}else {$40="-";for(x=1;x<=5-1;x++){$40=$40 FS "-"};print}}' | \
#rename bead name field with gene.transcript number
awk 'BEGIN{FS="\t";OFS=FS;}{tr[$14]++;$4="hanno.g"$14".t"tr[$14];print}' | \
#remove duplicate columns 19 and 42
cut -f 1-18,20-41,43- | \
#adjust scoring (orf length and add lastp,eggnog,busco scores; multiply by 2 if busco "Duplicated" or "Complete"; multiply by 1.5 if busco "Fragmented"); add original orf-length to column 43 and total exon length to 44, add cds exon count to 45 and total  exon count to 46
awk -v infile=orf.length 'BEGIN{FS="\t";OFS=FS;while(getline l < infile){n=split(l,a,FS);h[a[1]]=a[3]}}{$43=$5;n=split($11,a,",");for(x=1;x<=n;x++){$44=$44+a[x]};$45=h[$13];$46=$10;$5=$5+$17+$21+$41;if($40=="Duplicated" || $40=="Complete"){$5=2*$5} else if($40=="Fragmented"){$5=1.5*$5};print }' | \
#Chrom	mRNAstart	mRNAend	Name	score	strand	CDSstart	CDSend	itemRGB	blockCount	blockLen	blockstart	origName	CDSCLUSTER	lastP	description	besthitProt	score	evalue	seed_ortholog	evalue	score	eggNOG_OGs	max_annot_lvl	COG_category	Description	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction	PFAMs	busco1	busco2	busco3	orflen	mRNAlen	cds_exons_count	mRNA_exon_count
awk 'BEGIN{print "##Chrom\tmRNAstart\tmRNAend\tName\tscore\tstrand\tCDSstart\tCDSend\titemRGB\tblockCount\tblockLen\tblockstart\torigName\tCDSCLUSTER\tlastP\tdescription\tbesthitProt\tscore\tevalue\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\tbusco1\tbusco2\tbusco3\torflen\tmRNAlen\tcds_exons_count\tmRNA_exon_count"} {print}'
