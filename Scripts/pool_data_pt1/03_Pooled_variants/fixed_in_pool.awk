#!/usr/bin/awk -f

NR==FNR{
aSwe[$1 "\t" $2]=$3; aCat[$1 "\t" $2]=$4 } 

NR!=FNR && FNR>2 && $4 != "0/0/0/0"{

Base[1]="A";
Base[2]="C";
Base[3]="G";
Base[4]="T";

split($4, a, "/"); 

for(i in a){
	if(a[i]==0){counter++} else{base=i; cov=a[i]}}; 

if(counter == 3){
	if(Base[base] == aSwe[$1 "\t" $2] && ($1 "\t" $2) in aSwe){print $1 "\t" $2 "\t" cov "\t" 1 "\t" aSwe[$1 "\t" $2] "\t" aCat[$1 "\t" $2]} 
	else if (Base[base] == aCat[$1 "\t" $2] && ($1 "\t" $2) in aCat) {print $1 "\t" $2 "\t" cov "\t" 0 "\t" aSwe[$1 "\t" $2] "\t" aCat[$1 "\t" $2]}  
	}
	
	
	counter = 0; 
}
