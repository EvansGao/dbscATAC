@spes=("hs","mm");
foreach $spe (@spes){
	my $dir="dbSUPERraw/".$spe;
	opendir(DIR,$dir) or "can't open the file";
	@dirs=readdir DIR;
	@dirs=grep{$_ ne "." && $_ ne ".."} @dirs;
	open AA,">dbSUPERpre_".$spe.".bed";
	foreach $file (@dirs){
	$cell=$file;
	$cell=~ s/^\d+\_|\.bed$//g;
	open BB,$dir."/".$file;
		while(<BB>){
		s/\r|\n//g;
		@temp=split/\t/,$_;
		if($temp[4]=~ /^SE\_\d+|^mSE\_\d+/){
		print AA $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[4].":".$cell."\n";
		}elsif($temp[3]=~ /^SE\_\d+|^mSE\_\d+/){
		print AA $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[3].":".$cell."\n";
		}

		
		}
	close BB;
	}
	close AA;
	system("bedtools sort -i dbSUPERpre_".$spe.".bed>dbSUPER_".$spe.".bed");
	system("liftOver dbSUPER_".$spe.".bed /data/gts/dbscATAC/project/liftover/hg19ToHg38.over.chain dbSUPER_hg38_".$spe.".bed unMapped");
	
}


#