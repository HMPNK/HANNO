#clusters bedtools intersect outputs to get gene loci
BEGIN{printf "\n"  > "/dev/stderr"}
{
n=split($2,a,";");
j++;
names[j]=$1;
#search clusters
for(x=2;x<=n;x++){if(c[a[x]]!="") {cluster=c[a[x]];break} else{cluster=""}};
#if match assign cluster
if(cluster!="") { for(x=2;x<=n;x++){c[a[x]]=cluster;} }
#if no match begin new cluster
else { i++; for(x=2;x<=n;x++){c[a[x]]=i;}}
}

END	{
	printf "CLUSTERS: "i"\n"  > "/dev/stderr"
	printf "\n"  > "/dev/stderr"
	for(x=1;x<=j;x++){print names[x]"\t"c[names[x]]};
	}
