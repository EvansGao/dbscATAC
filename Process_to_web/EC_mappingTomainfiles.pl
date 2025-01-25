$start = time;
use File::Copy;
%hashGene=();
#$spe="human";
#$speup="HS";
#$spelow="hs";

#@builds=("galGal5","danRer10","TAIR10");
#@spes=("chicken","zebrafish","arabidopsis");
#@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
@builds=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
@spes=("marmoset","mouse","cynomolgus","vervet","human","chicken","zebrafish","chimp","rhesus","fly","arabidopsis","rice","maize");
@speups=("CJ","MM","MF","CS","HS","GG","DR","PT","RM","DM","AT","OSJ","ZM");
@spelows=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
#@fullnames=("Callithrix_jacchus","Mus_musculus","Macaca_fascicularis","Chlorocebus_sabaeus","Homo_sapiens","Gallus_gallus","Danio_rerio","Pan_troglodytes","Macaca_mulatta","Drosophila_melanogaster","Arabidopsis_thaliana","Oryza_sativa_Japonica","Zea_mays");


#@builds=("gg","mm","at","dr");
#@spes=("chicken","mouse","arabidopsis","zebrafish");
#@speups=("GG","MM","AT","DR");
#@spelows=("gg","mm","at","dr");


mkdir("AllEPs");
for($i=0;$i<scalar(@builds);$i++){
mkdir("AllEPs/".$builds[$i]);
if(-d "download/interaction/".$builds[$i]){
%hashUniprotTOensembl=();
open EE,"webgenes/uniprot".$spes[$i].".txt";
while(<EE>){
s/\r|\n//g;
@temp=split/\t/,$_;
$temp[3]=~ s/^\s+|\s+$//g;
$temp[4]=~ s/^\s+|\s+$//g;
	if(!exists $hashUniprotTOensembl{$temp[3]} && $temp[3] ne ""){
	$hashUniprotTOensembl{$temp[3]}=$temp[1];
	}
	if(!exists $hashUniprotTOensembl{$temp[4]} && $temp[4] ne ""){
	$hashUniprotTOensembl{$temp[4]}=$temp[1];
	}
}
close EE;
my $dir="download/interaction/".$builds[$i];
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
@cells=();
%hashnameRank=();
$rankNum=1;
%hashgeneEnhPositions=();
%hashgenetoScore=();
%hashegimap=();
%hashegicell=();

foreach $file (@dir){
$cell=$file;
$cell=~ s/\_interaction\.txt$//g;
open EP, ">AllEPs/".$builds[$i]."/".$cell."_EP.txt";
open CELL,$dir."/".$file;
while(<CELL>){
s/\r|\n//g;
@tmp=split/\:|\-|\_|\||\t/,$_;
$tmpensembl="Unkown";
	if(exists $hashUniprotTOensembl{$tmp[8]}){
	$tmpensembl=$hashUniprotTOensembl{$tmp[8]};
	}elsif(!exists $hashUniprotTOensembl{$tmp[8]} && exists $hashUniprotTOensembl{$tmp[7]}){
	$tmpensembl=$hashUniprotTOensembl{$tmp[7]};
	}

print EP $tmp[0].":".$tmp[1]."-".$tmp[2]."_".$tmpensembl."\$".$tmp[7]."\$".$tmp[8]."\$".$tmp[4].":".$tmp[5]."-".$tmp[6]."\t".$tmp[10]."\n";
}
close CELL;
close EP;
}

%hashGene=();
open AA,"geneMapdata".$spes[$i].".txt";
%hash=();
while(<AA>){
chomp($_);
@temp=split/\t/,$_;
if($temp[0] ne ""){
$hashGene{$temp[0]}=$temp[1];
}
}
close AA;

mkdir("enhs");
mkdir("enhs/".$builds[$i]);
my $dir="AllEPs/".$builds[$i];
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
@cells=();
%hashnameRank=();
$rankNum=1;
%hashgeneEnhs=();
@genes=();
%hashgeneEnhPositions=();
%hashgenetoScore=();
open RANK,">cell_number".$builds[$i].".txt";
%hashegimap=();
%hashegicell=();
foreach $file (@dir){
$cell=$file;
$cell=~ s/\_EP\.txt//g;
$celluc=$cell;
$celluc=~s/\-|\_//g;
$celluc=uc($celluc);
$rankNumTwodigtal = sprintf("%02d",$rankNum);
$hashnameRank{$cell}=$speups[$i].$rankNumTwodigtal;
#print $hashnameRank{$cell}."\n";
print RANK $hashnameRank{$cell}."\t".$celluc."\t".$cell."\t"."Cell"."\n";
	open CELL,$dir."/".$file;
	open ENHPOS,">enhs/".$builds[$i]."/".$hashnameRank{$cell};
	$n=1;
	%hashexist=();
	while(<CELL>){
	s/\r|\n//g;
	@tmp=split/\t/,$_;
	@temp=split/\:|\-|\_|\$/,$tmp[0];
	if(!exists $hashGene{$temp[3]}){next;}
#chr3L:4467113-4467853_FBgn0035543$CG15020$chr3L$4421662$	20.445884
	if(!exists $hashexist{$temp[0].":".$temp[1]."-".$temp[2]}){
	$nn=sprintf("%05d",$n);
	$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}=$hashnameRank{$cell}."-".$nn;
	$noChr=$temp[0];
	$noChr=~ s/chr//g;
	print ENHPOS $hashnameRank{$cell}."-".sprintf("%05d",$n)."\t".$noChr.":".$temp[1]."-".$temp[2]."\n";
		if(!exists $hashgeneEnhs{$temp[3]}){
		$hashgeneEnhs{$temp[3]}=$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}.",";
		$hashgeneEnhPositions{$temp[3]}=$noChr.":".$temp[1]."-".$temp[2];
		$hashgenetoScore{$temp[3]}=$tmp[1];
		push @genes,$temp[3];
		}else{
			if(!($hashgeneEnhs{$temp[3]}=~ /$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}\,/)){
			$hashgeneEnhs{$temp[3]}.=$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}.",";
			$hashgeneEnhPositions{$temp[3]}.=",".$noChr.":".$temp[1]."-".$temp[2];
			$hashgenetoScore{$temp[3]}.=",".$tmp[1];
			}
		}
	$n++;
	}else{
	$nn=sprintf("%05d",$n);
	$noChr=$temp[0];
	$noChr=~ s/chr//g;
		if(!exists $hashgeneEnhs{$temp[3]}){
		$hashgeneEnhs{$temp[3]}=$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}.",";
		$hashgeneEnhPositions{$temp[3]}=$noChr.":".$temp[1]."-".$temp[2];
		$hashgenetoScore{$temp[3]}=$tmp[1];
		push @genes,$temp[3];
		}else{
			if(!($hashgeneEnhs{$temp[3]}=~ /$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}\,/)){
			$hashgeneEnhs{$temp[3]}.=$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}.",";
			$hashgeneEnhPositions{$temp[3]}.=",".$noChr.":".$temp[1]."-".$temp[2];
			$hashgenetoScore{$temp[3]}.=",".$tmp[1];
			}
		}
	
	
	
	
	}
	
	if(!exists $hashegimap{$temp[3]} && !exists $hashegicell{$temp[3]}){
	$hashegimap{$temp[3]}=$hashGene{$temp[3]};
	$hashegicell{$temp[3]}=$hashnameRank{$cell};
	}elsif(exists $hashegimap{$temp[3]} && !($hashegicell{$temp[3]}=~ /$hashnameRank{$cell}/)){
	$hashegicell{$temp[3]}.=",".$hashnameRank{$cell};
	}
}
close CELL;
close ENHPOS;
$rankNum++;
}
close RANK;

open DATA,">EPEP".$spes[$i].".txt";
open EGIMAP,">geneinfo".$spes[$i].".txt";
foreach $gene (@genes){
$hashgeneEnhs{$gene}=~ s/\,$//g;
#print DATA $transcript."\t".$hashtranscriptGene{$transcript}."\t".$hashtranscriptEnhs{$transcript}."\t".$hashtranscriptEnhPositions{$transcript}."\n";	
print DATA $gene."\t".$hashgenetoScore{$gene}."\t".$hashgeneEnhs{$gene}."\n";
print EGIMAP $gene."\t".$hashegicell{$gene}."\t".$hashGene{$gene}."\n";
}
close DATA;
close EGIMAP;

}
}
$duration = time - $start;
print "All are done: $duration s\n";