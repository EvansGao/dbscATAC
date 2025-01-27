$start = time;
mkdir("download/SEMotif");
@spes=("hs","mm","dm");
	foreach $spe (@spes){
	mkdir("download/SEMotif/".$spe);
	my $dir="download/superenhancer/".$spe;
	opendir(DIR,$dir) or "can't open the file";
	@dir=readdir DIR;
	@dir=grep{$_ ne "." && $_ ne ".."} @dir;
		foreach $cellfile (@dir){
		$cellname=$cellfile;
		$cellname=~ s/\.bed$//g;
			if(!(-e "download/SEMotif/".$spe."/".$cellfile)){
			system("bedtools intersect -a ".$dir."/".$cellfile." -b Motif/Motif_".$spe.".bed -F 1.0 -wa -wb>download/SEMotif/".$spe."/".$cellfile);
			}
		}
	}

$duration = time - $start;
print "All are done: $duration s\n";