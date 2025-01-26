@links=("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4508928&format=file&file=GSM4508928%5Fadrenal%5Ffiltered%2Eseurat%2ERDS%2Egz","https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4508930&format=file&file=GSM4508930%5Fcerebrum%5Ffiltered%2Eseurat%2ERDS%2Egz","https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4508931&format=file&file=GSM4508931%5Feye%5Ffiltered%2Eseurat%2ERDS%2Egz","https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4508932&format=file&file=GSM4508932%5Fheart%5Ffiltered%2Eseurat%2ERDS%2Egz","https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4508933&format=file&file=GSM4508933%5Fintestine%5Ffiltered%2Eseurat%2ERDS%2Egz","https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4508934&format=file&file=GSM4508934%5Fkidney%5Ffiltered%2Eseurat%2ERDS%2Egz","https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4508935&format=file&file=GSM4508935%5Fliver%5Ffiltered%2Eseurat%2ERDS%2Egz","https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4508936&format=file&file=GSM4508936%5Flung%5Ffiltered%2Eseurat%2ERDS%2Egz","https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4508937&format=file&file=GSM4508937%5Fmuscle%5Ffiltered%2Eseurat%2ERDS%2Egz","https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4508938&format=file&file=GSM4508938%5Fpancreas%5Ffiltered%2Eseurat%2ERDS%2Egz","https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4508939&format=file&file=GSM4508939%5Fplacenta%5Ffiltered%2Eseurat%2ERDS%2Egz","https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4508940&format=file&file=GSM4508940%5Fspleen%5Ffiltered%2Eseurat%2ERDS%2Egz","https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4508941&format=file&file=GSM4508941%5Fstomach%5Ffiltered%2Eseurat%2ERDS%2Egz","https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4508942&format=file&file=GSM4508942%5Fthymus%5Ffiltered%2Eseurat%2ERDS%2Egz");
foreach $linkone (@links){
	if($linkone=~ /file\&file\=(.*)/){
	$tmpname=$1;
	$tmpname=~ s/\%5F/\_/g;
	$tmpname=~ s/\%2D/\-/g;
	$tmpname=~ s/\%2E/\./g;
	system("wget '".$linkone."' -q -o /dev/null -O ".$tmpname);
	system("gunzip ".$tmpname);
	}
}
