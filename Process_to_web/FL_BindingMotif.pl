$start = time;
#!note genome version
@spevers=("hg38","mm10","araTha1","dm6");
@spes=("hs","mm","at","dm");
$cutoff=3;
foreach $spe (@spes){
	my $dir="F:/Motif/".$spe."/sites";
	opendir(DIR,$dir) or "can't open the file";
	@dirs=readdir DIR;
	@dirs=grep{$_ ne "." && $_ ne ".."} @dirs;
	open AA,">Motifpre_".$spe.".bed";
	foreach $file (@dirs){
	$jasparID=$file;
	$jasparID=~ s/\.tsv$//g;
	open BB,$dir."/".$file;
		while(<BB>){
		s/\r|\n//g;
		@temp=split/\t/,$_;
		if(($temp[5]/100)>=$cutoff){
		print AA $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[3].":".$jasparID.":".($temp[5]/100).":".$temp[6]."\n";
		}

		
		}
	close BB;
	}
	close AA;
	system("bedtools sort -i Motifpre_".$spe.".bed>Motif_".$spe.".bed");
}
$duration = time - $start;
print "All are done: $duration s\n";

#