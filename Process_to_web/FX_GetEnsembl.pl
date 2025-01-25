system("rm -R Ensembl");
#@spes=("hs","mm","at","dr");
#@spevers=("hg38","mm10","TAIR10","danRer10");
@spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
mkdir("Ensembl");
for($i=0;$i<scalar(@spes);$i++){
open AA,"../project/Refgene/ensembl_".$spevers[$i].".txt";
%chrgene=();
%hashgene=();
%hashgeneExonstart=();
%hashgeneName=();
%hashgeneStrand=();
%hashgeneStartEnd=();
@chrarray=();
while(<AA>){
chomp($_);
@tmp=split/\t/,$_;
if(!($tmp[6]=~ /Name/i)){
	if(!exists $chrgene{"chr".$tmp[6]}){
	$chrgene{"chr".$tmp[6]}=$tmp[0];
	push @chrarray,"chr".$tmp[6];
	}else{
		if(!($chrgene{"chr".$tmp[6]}=~ /$tmp[0]/)){
		$chrgene{"chr".$tmp[6]}.=",".$tmp[0];
		}
	}
	
}
$hashgeneExonstart{$tmp[0]."\t".$tmp[1]}=$tmp[2];

if(!exists $hashgene{$tmp[0]}){
$hashgene{$tmp[0]}=$tmp[1];
$hashgeneName{$tmp[0]}=$tmp[7];
		if($tmp[5] eq "1"){
		$hashgeneStrand{$tmp[0]}="+";
		}else{
		$hashgeneStrand{$tmp[0]}="-";
		}

$hashgeneStartEnd{$tmp[0]}=$tmp[3]."\t".$tmp[4];
}else{
$hashgene{$tmp[0]}.=",".$tmp[1];
}

}
close AA;
@chrarray=sort{$a<=>$b} @chrarray;

mkdir("Ensembl/".$spes[$i]);
foreach $chrom (@chrarray){
	if(length($chrom)<8){
	open $chrom,">Ensembl/".$spes[$i]."/".$chrom.".bed";
	@gene=split/\,/,$chrgene{$chrom};
	foreach $singlegene (@gene){
		@exonstarts=split/\,/,$hashgene{$singlegene};
		@exonstarts=sort{$a<=>$b} @exonstarts;
		$exonS="";
		$exonE="";
		for($j=0;$j<scalar(@exonstarts);$j++){
		$exonS.=",".$exonstarts[$j];
		$exonE.=",".$hashgeneExonstart{$singlegene."\t".$exonstarts[$j]};
		}
		$exonS=~ s/^\,//g;
		$exonE=~ s/^\,//g;
		print $chrom $chrom."\t".$hashgeneStartEnd{$singlegene}."\t".$hashgeneStrand{$singlegene}."\t".$hashgeneName{$singlegene}."\t".$exonS."\t".$exonE."\t".$singlegene."\n";
	}
	close $chrom;
	system("bedtools sort -i Ensembl/".$spes[$i]."/".$chrom.".bed >Ensembl/".$spes[$i]."/".$chrom);
	unlink("Ensembl/".$spes[$i]."/".$chrom.".bed");
	}
}

}