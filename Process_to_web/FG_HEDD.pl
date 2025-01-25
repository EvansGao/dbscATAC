	my $dir="HEDDraw/";
	opendir(DIR,$dir) or "can't open the file";
	@dirs=readdir DIR;
	@dirs=grep{$_ ne "." && $_ ne ".."} @dirs;
	@enhancers=();
	%hashenhToDisease=();
	%hashenhTotargetGene=();
	foreach $file (@dirs){
	$disease=$file;
	$disease=~ s/^DataDownload\-|\_EnhancerScore\.txt//g;
	open BB,$dir."/".$file;
	while(<BB>){
	s/\r|\n//g;
	@temp=split/\t/,$_;
	$enh="chr".$temp[1]."\t".$temp[2]."\t".$temp[3];
	if($temp[0] ne "EnhID"){
		if(!exists $hashenhToDisease{$enh}){
		$hashenhToDisease{$enh}=$disease;
		push @enhancers,$enh;
		@genes=split/\;\s+/,$temp[4];
			foreach $gene (@genes){
				if(!exists $hashenhTotargetGene{$enh}){
				$hashenhTotargetGene{$enh}=$gene;
				}elsif(exists $hashenhTotargetGene{$enh} && !($hashenhTotargetGene{$enh}=~ /$gene/)){
				$hashenhTotargetGene{$enh}.=",".$gene;
				}
			}
		}else{
		$hashenhToDisease{$enh}.=",".$disease;
		@genes=split/\;\s+/,$temp[4];
			foreach $gene (@genes){
				if(!exists $hashenhTotargetGene{$enh}){
				$hashenhTotargetGene{$enh}=$gene;
				}elsif(exists $hashenhTotargetGene{$enh} && !($hashenhTotargetGene{$enh}=~ /$gene/)){
				$hashenhTotargetGene{$enh}.=",".$gene;
				}
			}
		}
	}
	}
	close BB;
	}

	open AA,">HEDDpre_hs.bed";
	foreach $enh (@enhancers){
	print AA $enh."\t".$hashenhToDisease{$enh}."|".$hashenhTotargetGene{$enh}."\n";
	}
	close AA;
	system("bedtools sort -i HEDDpre_hs.bed>HEDD_hs.bed");
	unlink("./HEDDpre_hs.bed");


#