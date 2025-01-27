$start = time;
mkdir("download/SEGeneHancer");
@spes=("hs");
	foreach $spe (@spes){
	mkdir("download/SEGeneHancer/".$spe);
	my $dir="download/superenhancer/".$spe;
	opendir(DIR,$dir) or "can't open the file";
	@dir=readdir DIR;
	@dir=grep{$_ ne "." && $_ ne ".."} @dir;
		foreach $cellfile (@dir){
		$cellname=$cellfile;
		$cellname=~ s/\.bed$//g;
		system("bedtools intersect -a ".$dir."/".$cellfile." -b GH_hg38.bed -wa -wb>download/SEGeneHancer/".$spe."/".$cellfile);
		}
	}

$duration = time - $start;
print "All are done: $duration s\n";