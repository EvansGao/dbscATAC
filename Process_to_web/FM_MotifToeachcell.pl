$start = time;
mkdir("download/Motif");
@spes=("hs","mm","at","dr","dm");
	foreach $spe (@spes){
	mkdir("download/Motif/".$spe);
	my $dir="download/enhancer/".$spe;
	opendir(DIR,$dir) or "can't open the file";
	@dir=readdir DIR;
	@dir=grep{$_ ne "." && $_ ne ".."} @dir;
		foreach $cellfile (@dir){
		$cellname=$cellfile;
		$cellname=~ s/\.bed$//g;
			if(!(-e "download/Motif/".$spe."/".$cellfile)){
			system("bedtools intersect -a ".$dir."/".$cellfile." -b Motif/Motif_".$spe.".bed -F 1.0 -wa -wb>download/Motif/".$spe."/".$cellfile);
			}
		}
	}

$duration = time - $start;
print "All are done: $duration s\n";