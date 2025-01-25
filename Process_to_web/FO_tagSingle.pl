$start = time;
open IDID,"cellsID.txt";
%hashcellspeTOid=();
while(<IDID>){
s/\r|\n//g;
@tmp=split/\t/,$_;
$hashcellspeTOid{$tmp[0]."\t".$tmp[1]}=$tmp[2];
}
close IDID;

mkdir("single");
#@spes=("hs","mm","dm");
#@spes=("hs","mm","dm");
#@spes=("gg","mm","at","dr");
@spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
#@spes=("mm","at","dr","gg","hs","mam","pt");
foreach $spe (@spes){
	mkdir("single/".$spe);
	my $dir="download/enhancer/".$spe;
	opendir(DIR,$dir) or "can't open the file";
	@dirs=readdir DIR;
	@dirs=grep{$_ ne "." && $_ ne ".."} @dirs;
	foreach $file (@dirs){
		$cell=$file;
		$cell=~ s/\.bed//g;
#SNP
		%hashenhSNP=();
		if(-e "download/SNP/".$spe."/".$file){
		open SNPSNP,"download/SNP/".$spe."/".$file;
			while(<SNPSNP>){
			s/\r|\n//g;
			@tmp=split/\t/,$_;
			if(!exists $hashenhSNP{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$hashenhSNP{$tmp[0].":".$tmp[1]."-".$tmp[2]}=$tmp[7].":".$tmp[4].":".$tmp[5];
			}else{
			$hashenhSNP{$tmp[0].":".$tmp[1]."-".$tmp[2]}.=";".$tmp[7].":".$tmp[4].":".$tmp[5];
			}

			}
		close SNPSNP;
		}
#GeneHancer
		%hashenhGH=();
		if(-e "download/GeneHancer/".$spe."/".$file){
		open GH,"download/GeneHancer/".$spe."/".$file;
			while(<GH>){
			s/\r|\n//g;
			@tmp=split/\t/,$_;
			if(!exists $hashenhGH{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$hashenhGH{$tmp[0].":".$tmp[1]."-".$tmp[2]}=$tmp[7].":".$tmp[4].":".$tmp[5].":".$tmp[6];
			}else{
			$hashenhGH{$tmp[0].":".$tmp[1]."-".$tmp[2]}.=";".$tmp[7].":".$tmp[4].":".$tmp[5].":".$tmp[6];
			}

			}
		close GH;
		}
#dbSUPER
		%hashenhdbSUPER=();
		if(-e "download/dbSUPER/".$spe."/".$file){
		open DBSUPER,"download/dbSUPER/".$spe."/".$file;
			while(<DBSUPER>){
			s/\r|\n//g;
			@tmp=split/\t/,$_;
			if(!exists $hashenhdbSUPER{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$hashenhdbSUPER{$tmp[0].":".$tmp[1]."-".$tmp[2]}=$tmp[7].":".$tmp[4].":".$tmp[5].":".$tmp[6];
			}else{
			$hashenhdbSUPER{$tmp[0].":".$tmp[1]."-".$tmp[2]}.=";".$tmp[7].":".$tmp[4].":".$tmp[5].":".$tmp[6];
			}

			}
		close DBSUPER;
		}
#Motif
		%hashenhMotif=();
		if(-e "download/Motif/".$spe."/".$file){
		open MOTIF,"download/Motif/".$spe."/".$file;
			while(<MOTIF>){
			s/\r|\n//g;
			@tmp=split/\t/,$_;
			if(!exists $hashenhMotif{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$hashenhMotif{$tmp[0].":".$tmp[1]."-".$tmp[2]}=$tmp[7].":".$tmp[4].":".$tmp[5].":".$tmp[6];
			}else{
			$hashenhMotif{$tmp[0].":".$tmp[1]."-".$tmp[2]}.=";".$tmp[7].":".$tmp[4].":".$tmp[5].":".$tmp[6];
			}

			}
		close MOTIF;
		}		

#HEDD
		%hashenhHEDD=();
		if(-e "download/HEDD/".$spe."/".$file){
		open HEDD,"download/HEDD/".$spe."/".$file;
			while(<HEDD>){
			s/\r|\n//g;
			@tmp=split/\t/,$_;
			if(!exists $hashenhHEDD{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$hashenhHEDD{$tmp[0].":".$tmp[1]."-".$tmp[2]}=$tmp[7].":".$tmp[4].":".$tmp[5].":".$tmp[6];
			}else{
			$hashenhHEDD{$tmp[0].":".$tmp[1]."-".$tmp[2]}.=";".$tmp[7].":".$tmp[4].":".$tmp[5].":".$tmp[6];
			}

			}
		close HEDD;
		}
		
#TargetGene
		%hashenhTargetGene=();
		if(-e "download/TargetGene/".$spe."/".$file){
		open TargetGene,"download/TargetGene/".$spe."/".$file;
			while(<TargetGene>){
			s/\r|\n//g;
			@tmp=split/\t/,$_;
			if($tmp[0] ne ""){
			$hashenhTargetGene{$tmp[0]}=$tmp[1];
			}
			}
		close TargetGene;
		}
	
		open AA,"download/enhancer/".$spe."/".$file;
		open BB,">"."single/".$spe."/".$file;
		$num=1;
		while(<AA>){
		s/\r|\n//g;
		@tmp=split/\t/,$_;
		$ID=$hashcellspeTOid{$cell."\t".$spe}."-".sprintf("%05d",$num);
		$str=$ID."\t".$tmp[0].":".$tmp[1]."-".$tmp[2];
			if(exists $hashenhSNP{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$str.="\t".$hashenhSNP{$tmp[0].":".$tmp[1]."-".$tmp[2]};
			}else{
			$str.="\t"."NA";
			}
			if(exists $hashenhGH{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$str.="\t".$hashenhGH{$tmp[0].":".$tmp[1]."-".$tmp[2]};
			}else{
			$str.="\t"."NA";
			}
			if(exists $hashenhdbSUPER{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$str.="\t".$hashenhdbSUPER{$tmp[0].":".$tmp[1]."-".$tmp[2]};
			}else{
			$str.="\t"."NA";
			}
			if(exists $hashenhMotif{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$str.="\t".$hashenhMotif{$tmp[0].":".$tmp[1]."-".$tmp[2]};
			}else{
			$str.="\t"."NA";
			}
			if(exists $hashenhHEDD{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$str.="\t".$hashenhHEDD{$tmp[0].":".$tmp[1]."-".$tmp[2]};
			}else{
			$str.="\t"."NA";
			}
			if(exists $hashenhTargetGene{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$str.="\t".$hashenhTargetGene{$tmp[0].":".$tmp[1]."-".$tmp[2]};
			}else{
			$str.="\t"."NA";
			}
			
			
		print BB $str."\n";
		$num++;
		}
		close AA;
		close BB;
	
	}
}

$duration = time - $start;
print "All are done: $duration s\n";