{
split($4,a,/[\._]/);
split($10,b,".");
if(a[1]a[2]==b[1]b[2] && a[1]a[2]a[3]!=$10) {
			e[$4]=e[$4]","$9;
			s[$4]=s[$4]","$8;
			if(t[$4]==""){i++;l[i]=$4};
			}
}

END{
for(x=1;x<=i;x++) {
			id=l[x];
			n1=split(s[id],a,",");
			n2=split(e[id],b,",");
			for(y=1;y<=n1;y++) {
				for(z=1;z<=n1;z++){if(b[y]<a[z]){if(b[y]!="" && a[z]!=""){print id"\t"b[y]"\t"a[z]}}}
				}
		}
}
