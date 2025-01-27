$start = time;
mkdir("download/SEHEDD");
@spes=("hs");
	foreach $spe (@spes){
	mkdir("download/SEHEDD/".$spe);
	my $dir="download/superenhancer/".$spe;
	opendir(DIR,$dir) or "can't open the file";
	@dir=readdir DIR;
	@dir=grep{$_ ne "." && $_ ne ".."} @dir;
		foreach $cellfile (@dir){
		$cellname=$cellfile;
		$cellname=~ s/\.bed$//g;
		system("bedtools intersect -a ".$dir."/".$cellfile." -b HEDD_".$spe.".bed -F 1.0 -wa -wb>download/SEHEDD/".$spe."/".$cellfile);
		}
	}

$duration = time - $start;
print "All are done: $duration s\n";