$start = time;
#!note genome version
#@spes=("human","mouse","fly");
#@spes=("gg","mm","at","dr");
@spes=("hs","mm","dm");
$cutoff=3;
foreach $spe (@spes){
	my $dir="AllSUEPs/".$spe;
	mkdir("download/SETargetGene");
	mkdir("download/SETargetGene/".$spe);
	opendir(DIR,$dir) or "can't open the file";
	@dirs=readdir DIR;
	@dirs=grep{$_ ne "." && $_ ne ".."} @dirs;
	foreach $file (@dirs){
	@enhancers=();
	%hashenhancerTogenes=();
	$cell=$file;
	$cell=~ s/\_EP\.txt//g;
	open BB,$dir."/".$file;
		while(<BB>){
		s/\r|\n//g;
		@temp=split/\_|\$/,$_;
		if(!exists $hashenhancerTogenes{$temp[0]}){
		$hashenhancerTogenes{$temp[0]}=$temp[1].":".$temp[2].":";
		push @enhancers,$temp[0];
		}elsif(exists $hashenhancerTogenes{$temp[0]} && !($hashenhancerTogenes{$temp[0]}=~ /\:$temp[2]\:/)){
		$hashenhancerTogenes{$temp[0]}.=";".$temp[1].":".$temp[2].":";
		}
		}
	close BB;
	open AA,">download/SETargetGene/".$spe."/".$cell.".bed";
	foreach $enh (@enhancers){
	print AA $enh."\t".$hashenhancerTogenes{$enh}."\n";
	}
	close AA;
	}
}
$duration = time - $start;
print "All are done: $duration s\n";

#