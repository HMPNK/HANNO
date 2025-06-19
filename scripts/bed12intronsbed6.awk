#get intron boundaries from bed12 and write as bed6
{
n1=split($11,a,",");
n2=split($12,b,",");
for(x=2;x<n1;x++)	{
			print $1"\t"$2+a[x-1]+b[x-1]"\t"$2+b[x]"\t"$4"_intron_"x-1"\t.\t"$6
			}
}
