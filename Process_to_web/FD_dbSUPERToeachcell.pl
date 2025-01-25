$start = time;
mkdir("download/dbSUPER");
@spes=("hs","mm");
	foreach $spe (@spes){
	mkdir("download/dbSUPER/".$spe);
	my $dir="download/enhancer/".$spe;
	opendir(DIR,$dir) or "can't open the file";
	@dir=readdir DIR;
	@dir=grep{$_ ne "." && $_ ne ".."} @dir;
		foreach $cellfile (@dir){
		$cellname=$cellfile;
		$cellname=~ s/\.bed$//g;
		system("bedtools intersect -a ".$dir."/".$cellfile." -b dbSUPER_hg38_".$spe.".bed -f 0.1 -wa -wb>download/dbSUPER/".$spe."/".$cellfile);
		}
	}

$duration = time - $start;
print "All are done: $duration s\n";