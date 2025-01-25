$start = time;
use List::Util qw[min max sum];
#use Statistics::Basic qw<median>;
use List::MoreUtils qw(uniq);
use List::Util qw(shuffle);
@spes=("arabidopsis","marmoset","rice","maize","vervet","cynomolgus");
@builds=("TAIR10","calJac3","IRGSP1","B73v4","ChlSab1_1","macFas6");
@upstreams=("250","5000","800","3000","4000","5000");
#@downstreams=("25","500","500","500","200","280");
for($ii=0;$ii<scalar(@spes);$ii++){

$species=$spes[$ii];
$gbuild=$builds[$ii];
$upstream=$upstreams[$ii];
$downstream=$downstreams[$ii];
open CHROM,"/data/gts/dbscATAC/project/chromsizes/".$gbuild.".chrom.sizes";
%hashchr=();
while(<CHROM>){
s/\r|\n//g;
@tmp=split/\t/,$_;
$hashchr{$tmp[0]}="";
}
close CHROM;

open CODE,"/data/gts/dbscATAC/project/uniprot/uniprot_genes_".$species.".txt";
%hashname=();
%hashnameTOstandard=();
@uniprots=();
while(<CODE>){
s/\r|\n//g;
@tmp=split/\t/,$_;
@tmpnames=split/\s+/,$tmp[4];
if($tmp[0] ne ""){
push @uniprots,uc($tmp[0]);
}
foreach $tmpname (@tmpnames){
$tmpname=~ s/ZEAMMB73\_//g;
$hashname{uc($tmpname)}=uc($tmp[0]);
$hashnameTOstandard{uc($tmpname)}=$tmpnames[0];
}
}
close CODE;
print $hashnameTOstandard{"LARP4B"}."good";

open GENE,"/data/gts/dbscATAC/project/Refgene/genes_exon_".$gbuild.".txt";
%hashpro=();
open PRO,">standard_promotorpre.bed";
open PROEXON,">standard_promotor_exonpre.bed";
%hashrefgeneuniprot=();
$numgene=0;
%hashnumgene=();
while(<GENE>){
s/\r|\n//g;
@tmp=split/\t/,$_;
$tmp[2]=~ s/^\s+|\s+$//g;
$tmp[3]=~ s/^\s+|\s+$//g;
$tmp[4]=~ s/^\s+|\s+$//g;
$tmp[5]=~ s/^\s+|\s+$//g;
$tmp[6]=~ s/^\s+|\s+$//g;
$tmp[7]=~ s/^\s+|\s+$//g;
$tmp[0]=uc($tmp[0]);
$tmp[7]=uc($tmp[7]);
$uniprot="";
$standard="";
if($tmp[7] eq "" && exists $hashname{$tmp[0]}){
$uniprot=$hashname{$tmp[0]};
$standard=$hashnameTOstandard{$tmp[0]};
}elsif($tmp[7] ne "" && exists $hashname{$tmp[7]}){
$uniprot=$hashname{$tmp[7]};
$standard=$hashnameTOstandard{$tmp[7]};
}elsif($tmp[7] ne "" && !exists $hashname{$tmp[7]} && exists $hashname{$tmp[0]}){
$uniprot=$hashname{$tmp[0]};
$standard=$hashnameTOstandard{$tmp[0]};
}

$tmp[4]="chr".$tmp[4];
if(exists $hashchr{$tmp[4]} && $tmp[8] eq "1" && ($tmp[5]-$upstream)>0 && $uniprot ne ""){
	$hashrefgeneuniprot{$uniprot}="";
		if(!exists $hashpro{$tmp[4]."\t".($tmp[5]-$upstream)."\t".$tmp[5]."\t".$standard}){
		$hashpro{$tmp[4]."\t".($tmp[5]-$upstream)."\t".$tmp[5]."\t".$standard}="";
		print PRO $tmp[4]."\t".($tmp[5]-$upstream)."\t".$tmp[5]."\t".$standard.":".$uniprot."\n";
		}
	print PROEXON $tmp[4]."\t".($tmp[5]-$upstream)."\t".$tmp[5]."\t1"."\n";
	print PROEXON $tmp[4]."\t".$tmp[2]."\t".$tmp[3]."\t1"."\n";
}elsif(exists $hashchr{$tmp[4]} && $tmp[8] eq "-1" && $tmp[6]>0 && $uniprot ne ""){
	$hashrefgeneuniprot{$uniprot}="";
		if(!exists $hashpro{$tmp[4]."\t".$tmp[6]."\t".($tmp[6]+$upstream)."\t".$standard}){
		$hashpro{$tmp[4]."\t".$tmp[6]."\t".($tmp[6]+$upstream)."\t".$standard}="";
		print PRO $tmp[4]."\t".$tmp[6]."\t".($tmp[6]+$upstream)."\t".$standard.":".$uniprot."\n";
		}
	print PROEXON $tmp[4]."\t".$tmp[6]."\t".($tmp[6]+$upstream)."\t1"."\n";
	print PROEXON $tmp[4]."\t".$tmp[2]."\t".$tmp[3]."\t1"."\n";
}
if(!exists $hashnumgene{$tmp[0]}){
$hashnumgene{$tmp[0]}="";
$numgene++;
}
}
close GENE;
#print $numgene."good\n";

close PRO;
close PROEXON;
open LEFT,">uniprot_".$species."_genes_left.txt";
foreach $tmpuniprot (@uniprots){
	if(!exists $hashrefgeneuniprot{$tmpuniprot}){
		print LEFT $tmpuniprot."\n";
#		$leftgenes.=$hashname{$tmpname}.";";
#		if(index($leftgenes,$hashname{$tmpname})<0){
#		$leftgenes.=$hashname{$tmpname}.";";
#		}
#	if(exists $hashname{$tmpname} && !exists $hashrefgenename{$tmpname}){
#	print "OK";
	}
}
close LEFT;


system("bedtools sort -i standard_promotorpre.bed>standard_promotor_".$gbuild.".bed");
system("bedtools sort -i standard_promotor_exonpre.bed>standard_promotor_exonsort.bed");
system("bedtools merge -i standard_promotor_exonsort.bed -c 4 -o mean>standard_promotor_exon_".$gbuild.".bed");
unlink("standard_promotorpre.bed");
unlink("standard_promotor_exonpre.bed");
unlink("standard_promotor_exonsort.bed");
unlink("uniprot_".$species."_genes_left.txt");
}

