if(-d "./Genome"){
system("rm -R Genome");
}
mkdir("Genome");
#@spes=("ss","sc","gg","dr","ce","rn","dm");
#@spevers=("susScr3","sacCer3","galGal4","danRer10","ce10","rn5","dm3");
#@spevers=("susScr3","sc3","gg4","dr10","ce10","rn5","dm3");

@spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");

for($i=0;$i<scalar(@spes);$i++){
%hashchr=();
open AA,"../project/chromsizes/".$spevers[$i].".chrom.sizes";
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
	if($tmp[0] ne "" && !($tmp[0]=~ /random|chrUn/i)){
	$hashchr{$tmp[0]}="";
	}
}
close AA;

open CC,">".$spevers[$i].".bed";
if(-e "../project/Refgene/Refgene_".$spevers[$i].".txt"){
open BB,"../project/Refgene/Refgene_".$spevers[$i].".txt";
while(<BB>){
s/\r|\n//g;
@tmp=split/\t/,$_;
	if(exists $hashchr{$tmp[2]}){
	print CC $tmp[2]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[3]."\t".$tmp[12]."\t".$tmp[9]."\t".$tmp[10]."\n";
	}
}
close BB;
}elsif(-e "../project/Refgene/genes_exon_".$spevers[$i].".txt"){
@genepositions=();
%hashgenepositionTOstrand=();
%hashgenepositionTOname=();
%hashgenepositionTOexonStart=();
%hashgenepositionTOexonEnd=();
open BB,"../project/Refgene/genes_exon_".$spevers[$i].".txt";
while(<BB>){
s/\r|\n//g;
@tmp=split/\t/,$_;
$tmp[4]="chr".$tmp[4];
if($tmp[8] eq "1"){
$tmp[8]="+";
}elsif($tmp[8] eq "-1"){
$tmp[8]="-";
}
	if(!($tmp[0]=~ /Gene/i) && !exists $hashgenepositionTOstrand{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}){
	$hashgenepositionTOstrand{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}=$tmp[8];
	$hashgenepositionTOname{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}=$tmp[7];
	$hashgenepositionTOexonStart{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}=$tmp[2];
	$hashgenepositionTOexonEnd{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}=$tmp[3];
	push @genepositions,$tmp[4]."\t".$tmp[5]."\t".$tmp[6];
	}elsif(!($tmp[0]=~ /Gene/i) && exists $hashgenepositionTOstrand{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}){
	$hashgenepositionTOexonStart{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}.=",".$tmp[2];
	$hashgenepositionTOexonEnd{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}.=",".$tmp[3];
	}
}
close BB;
	foreach $geneposition (@genepositions){
	print CC $geneposition."\t".$hashgenepositionTOstrand{$geneposition}."\t".$hashgenepositionTOname{$geneposition}."\t".$hashgenepositionTOexonStart{$geneposition}."\t".$hashgenepositionTOexonEnd{$geneposition}."\n";
	}

}

close CC;
system("bedtools sort -i ".$spevers[$i].".bed>".$spevers[$i]."sort.bed");
mkdir("Genome/".$spes[$i]);
%hashchr=();
@chrs=();
open CHR,$spevers[$i]."sort.bed";
while(<CHR>){
s/\r|\n//g;
@tmp=split/\t/,$_;
	if($tmp[0] ne "" && !exists $hashchr{$tmp[0]}){
	$hashchr{$tmp[0]}="";
	push @chrs,$tmp[0];
	$tmpchr=$tmp[0];
	open $tmpchr,">"."Genome/".$spes[$i]."/".$tmpchr;
	print $tmpchr $_."\n";
	}elsif($tmp[0] ne "" && exists $hashchr{$tmp[0]}){
	$tmpchr=$tmp[0];
	print $tmpchr $_."\n";
	}
}
close CHR;
	foreach $chrom (@chrs){
	close $chrom;
	}
}
