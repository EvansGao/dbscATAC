@links=("https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Data/scATAC-Healthy-Hematopoiesis-191120.rds","https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Disease-Data/scATAC-All-Hematopoiesis-MPAL-191120.rds");
foreach $linkone (@links){
	if($linkone=~ /\-Data\/(.*)/){
	$tmpname=$1;
	$tmpname=~ s/\%5F/\_/g;
	$tmpname=~ s/\%2D/\-/g;
	$tmpname=~ s/\%2E/\./g;
	system("wget '".$linkone."' -q -o /dev/null -O ".$tmpname);
##	system("gunzip ".$tmpname);
	}
}