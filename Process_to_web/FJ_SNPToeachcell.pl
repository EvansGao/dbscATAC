$start = time;
mkdir("download/SNP");
@spes=("hs");
	foreach $spe (@spes){
	mkdir("download/SNP/".$spe);
	my $dir="download/enhancer/".$spe;
	opendir(DIR,$dir) or "can't open the file";
	@dir=readdir DIR;
	@dir=grep{$_ ne "." && $_ ne ".."} @dir;
		foreach $cellfile (@dir){
		$cellname=$cellfile;
		$cellname=~ s/\.bed$//g;
		system("bedtools intersect -a ".$dir."/".$cellfile." -b SNPs_hg38.bed -wa -wb>download/SNP/".$spe."/".$cellfile);
		}
	}

$duration = time - $start;
print "All are done: $duration s\n";