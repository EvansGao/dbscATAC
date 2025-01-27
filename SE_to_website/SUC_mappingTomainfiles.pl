$start = time;
use File::Copy;
%hashGene=();
#$spe="human";
#$speup="HS";
#$spelow="hs";

#@builds=("galGal5","danRer10","TAIR10");
#@spes=("chicken","zebrafish","arabidopsis");
#@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
@builds=("mm","hs","dm");
@spes=("mouse","human","fly");
@speups=("MM","HS","DM");
@spelows=("mm","hs","dm");
@spevers=("mm10","hg38","dm6");
#@builds=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
#@spes=("marmoset","mouse","cynomolgus","vervet","human","chicken","zebrafish","chimp","rhesus","fly","arabidopsis","rice","maize");
#@speups=("CJ","MM","MF","CS","HS","GG","DR","PT","RM","DM","AT","OSJ","ZM");
#@spelows=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");

#@fullnames=("Callithrix_jacchus","Mus_musculus","Macaca_fascicularis","Chlorocebus_sabaeus","Homo_sapiens","Gallus_gallus","Danio_rerio","Pan_troglodytes","Macaca_mulatta","Drosophila_melanogaster","Arabidopsis_thaliana","Oryza_sativa_Japonica","Zea_mays");


#@builds=("gg","mm","at","dr");
#@spes=("chicken","mouse","arabidopsis","zebrafish");
#@speups=("GG","MM","AT","DR");
#@spelows=("gg","mm","at","dr");
if(-d "./AllSUEPs"){
system("rm -R AllSUEPs");
}

mkdir("AllSUEPs");
for($i=0;$i<scalar(@builds);$i++){
mkdir("AllSUEPs/".$builds[$i]);
if(-d "download/SEgene/".$builds[$i]){
%hashnameTOensembl=();
open EE,"../project/Refgene/ensembl_".$spevers[$i].".txt";
$nn=0;
while(<EE>){
s/\r|\n//g;
$temp[3]=~ s/^\s+|\s+$//g;
@temp=split/\t/,$_;
	if(!exists $hashnameTOensembl{$temp[7]} && $temp[7] ne ""){
	$hashnameTOensembl{$temp[7]}=$temp[0];
	}
#if($nn==100){
#print "nice\n";
#print $hashnameTOensembl{"Gm24984"}."\n";
#}
#$nn++;
}
close EE;
%hashnameTOregion=();
%hashnameTOregionlen=();
open AA,"ROSErefseq/".$spevers[$i]."_refseq.ucsc";
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if($tmp[1] ne "" && !exists $hashnameTOregionlen{$tmp[12]}){
$hashnameTOregionlen{$tmp[12]}=$tmp[5]-$tmp[4];
$hashnameTOregion{$tmp[12]}=$tmp[2].":".$tmp[4]."-".$tmp[5];
}elsif($tmp[1] ne "" && exists $hashnameTOregionlen{$tmp[12]} && $hashnameTOregionlen{$tmp[12]}<($tmp[5]-$tmp[4])){
$hashnameTOregionlen{$tmp[12]}=$tmp[5]-$tmp[4];
$hashnameTOregion{$tmp[12]}=$tmp[2].":".$tmp[4]."-".$tmp[5];
}

}
close AA;

open UNIPROT,"../project/uniprot/uniprot_genes_".$spes[$i].".txt";
%hashnameTOuniprot=();
while(<UNIPROT>){
s/\r|\n//g;
@tmp=split/\t/,$_;
@tmpnames=split/\s+/,$tmp[4];
$hashnameTOuniprot{$tmpnames[0]}=$tmp[0];
}
close UNIPROT;

my $dir="download/SEgene/".$builds[$i];
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
$cell=~ s/\.bed$//g;
open EP, ">AllSUEPs/".$builds[$i]."/".$cell."_SUEP.txt";
open CELL,$dir."/".$file;
while(<CELL>){
s/\r|\n//g;
@tmp=split/\t/,$_;
$tmpensembl="Unkown";
	if(exists $hashnameTOensembl{$tmp[3]}){
	$tmpensembl=$hashnameTOensembl{$tmp[3]};
	}
$uniprot="Unkown";
	if(exists $hashnameTOuniprot{$tmp[3]}){
	$uniprot=$hashnameTOuniprot{$tmp[3]};
	}
print EP $tmp[0].":".$tmp[1]."-".$tmp[2]."_".$tmpensembl."\$".$tmp[3]."\$".$uniprot."\$".$hashnameTOregion{$tmp[3]}."\t"."1"."\n";
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

mkdir("SEenhs");
#enhs
mkdir("SEenhs/".$builds[$i]);
my $dir="AllSUEPs/".$builds[$i];
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
@cells=();
%hashnameRank=();
$rankNum=1;
%hashgeneSEs=();
@genes=();
%hashgeneEnhPositions=();
%hashgenetoScore=();
open RANK,">SEcell_number".$builds[$i].".txt";
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
	open ENHPOS,">SEenhs/".$builds[$i]."/".$hashnameRank{$cell};
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
		if(!exists $hashgeneSEs{$temp[3]}){
		$hashgeneSEs{$temp[3]}=$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}.",";
		$hashgeneEnhPositions{$temp[3]}=$noChr.":".$temp[1]."-".$temp[2];
		$hashgenetoScore{$temp[3]}=$tmp[1];
		push @genes,$temp[3];
		}else{
			if(!($hashgeneSEs{$temp[3]}=~ /$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}\,/)){
			$hashgeneSEs{$temp[3]}.=$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}.",";
			$hashgeneEnhPositions{$temp[3]}.=",".$noChr.":".$temp[1]."-".$temp[2];
			$hashgenetoScore{$temp[3]}.=",".$tmp[1];
			}
		}
	$n++;
	}else{
	$nn=sprintf("%05d",$n);
	$noChr=$temp[0];
	$noChr=~ s/chr//g;
		if(!exists $hashgeneSEs{$temp[3]}){
		$hashgeneSEs{$temp[3]}=$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}.",";
		$hashgeneEnhPositions{$temp[3]}=$noChr.":".$temp[1]."-".$temp[2];
		$hashgenetoScore{$temp[3]}=$tmp[1];
		push @genes,$temp[3];
		}else{
			if(!($hashgeneSEs{$temp[3]}=~ /$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}\,/)){
			$hashgeneSEs{$temp[3]}.=$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}.",";
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

open DATA,">SUEPSUEP".$spes[$i].".txt";
open EGIMAP,">SEgeneinfo".$spes[$i].".txt";
foreach $gene (@genes){
$hashgeneSEs{$gene}=~ s/\,$//g;
#print DATA $transcript."\t".$hashtranscriptGene{$transcript}."\t".$hashtranscriptSEs{$transcript}."\t".$hashtranscriptEnhPositions{$transcript}."\n";	
print DATA $gene."\t".$hashgenetoScore{$gene}."\t".$hashgeneSEs{$gene}."\n";
print EGIMAP $gene."\t".$hashegicell{$gene}."\t".$hashGene{$gene}."\n";
}
close DATA;
close EGIMAP;

}
}
$duration = time - $start;
print "All are done: $duration s\n";