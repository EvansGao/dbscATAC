open AA,">cellsID.txt";
my $firstdir="download/enhancer";
opendir(DIR,$firstdir) or "can't open the file";
@firstdir=readdir DIR;
@firstdir=grep{$_ ne "." && $_ ne ".."} @firstdir;
foreach $spe (@firstdir){
$speup=uc($spe);
	$seconddir=$firstdir."/".$spe;
	opendir(DIR,$seconddir) or "can't open the file";
	@seconddir=readdir DIR;
	@seconddir=grep{$_ ne "." && $_ ne ".."} @seconddir;
	$rankNum=1;
	foreach $cellfile (@seconddir){
	$cellname=$cellfile;
	$cellname=~ s/\.bed//g;
	print AA $cellname."\t".$spe."\t".$speup.sprintf("%03d",$rankNum)."\n";
	$rankNum++;
	}
}
close AA;