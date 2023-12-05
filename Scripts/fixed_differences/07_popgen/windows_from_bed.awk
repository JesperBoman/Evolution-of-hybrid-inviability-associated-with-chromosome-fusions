#u!/usr/bin/awk -f

#	WINDOWS FROM BED FILE	#
#	Made by Jesper Boman		#

#Creates windows from a bed file. Will make shorter windows then supplied for the scaffold-ending windows. 
#Outputs windows in 0-based coordinate system for start coordinate and 1-based for end coordinate as in bed-files

#Updated 23/8 - 2020
#Changelog: Minor fixes - changed where the window is written

#Updated 29/10 - 2020
#Changelog: Changed the header

#Updated 11/12 - 2022
#Modified from windows_from_fai4

#Usage: awk -f windows_from_bed.awk input.file.bed > output.bed

function ceil(x, y){y=int(x); return(x>y?y+1:y)}; #Ed Morton's ceil

{
#Change window size here

regionID = $1 "_" $2 "-" $3
window_size = 10000;
scaf=$1;
len=$3;
tot_windows=ceil((len-$2)/window_size);


while(i<tot_windows){
	if($2+window_size*(i+1)<len){
		print scaf "\t" $2+window_size*i "\t"  $2+window_size*(i+1) "\t" regionID "_" i+1;
		}		
	else { 
		print scaf "\t" $2+window_size*i "\t" len "\t" regionID "_" i+1;
		}
		i++;
		}
i=0;
}

