$start = time;
mkdir("download/SEdbSUPER");
@spes=("hs","mm");
	foreach $spe (@spes){
	mkdir("download/SEdbSUPER/".$spe);
	my $dir="download/superenhancer/".$spe;
	opendir(DIR,$dir) or "can't open the file";
	@dir=readdir DIR;
	@dir=grep{$_ ne "." && $_ ne ".."} @dir;
		foreach $cellfile (@dir){
		$cellname=$cellfile;
		$cellname=~ s/\.bed$//g;
		system("bedtools intersect -a ".$dir."/".$cellfile." -b dbSUPER_hg38_".$spe.".bed -f 0.1 -wa -wb>download/SEdbSUPER/".$spe."/".$cellfile);
		}
	}

$duration = time - $start;
print "All are done: $duration s\n";