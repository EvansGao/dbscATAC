$start = time;
use List::Util qw[min max sum];
#use Statistics::Basic qw<median>;
use List::MoreUtils qw(uniq);
use List::Util qw(shuffle);
@spes=("human","human","mouse","mouse","chicken","zebrafish","chimp","rhesus","fly");
@builds=("hg19","hg38","mm9","mm10","galGal6","danRer10","panTro5","rheMac10","dm6");
@upstreams=("5000","5000","5000","5000","2000","2800","5000","5000","250");
@downstreams=("500","500","500","500","200","280","500","500","25");
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
@uniprots=();
while(<CODE>){
s/\r|\n//g;
@tmp=split/\t/,$_;
@tmpnames=split/\s+/,$tmp[4];
if($tmp[0] ne ""){
push @uniprots,$tmp[0];
}
foreach $tmpname (@tmpnames){
$hashname{$tmpname}=$tmp[0];
}
}
close CODE;


open GENE,"/data/gts/dbscATAC/project/Refgene/Refgene_".$gbuild.".txt";
%hashpro=();
open PRO,">standard_promotorpre.bed";
open PROEXON,">standard_promotor_exonpre.bed";
%hashrefgeneuniprot=();
while(<GENE>){
s/\r|\n//g;
@tmp=split/\t/,$_;
$tmp[9]=~ s/\,$//g;
$tmp[10]=~ s/\,$//g;
$tmp[12]=~ s/^\s+|\s+$//g;
if(exists $hashchr{$tmp[2]} && $tmp[3] eq "+" && ($tmp[4]-$upstream)>0 && exists $hashname{$tmp[12]}){
	$hashrefgeneuniprot{$hashname{$tmp[12]}}="";
		if(!exists $hashpro{$tmp[2]."\t".($tmp[4]-$upstream)."\t".($tmp[4]+$downstream)."\t".$tmp[12]}){
		$hashpro{$tmp[2]."\t".($tmp[4]-$upstream)."\t".($tmp[4]+$downstream)."\t".$tmp[12]}="";
		print PRO $tmp[2]."\t".($tmp[4]-$upstream)."\t".($tmp[4]+$downstream)."\t".$tmp[12].":".$hashname{$tmp[12]}."\n";
		}
	print PROEXON $tmp[2]."\t".($tmp[4]-$upstream)."\t".($tmp[4]+$downstream)."\t1"."\n";
	@exonstarts=split/\,/,$tmp[9];
	@exonends=split/\,/,$tmp[10];
	for($i=0;$i<scalar(@exonstarts);$i++){
	print PROEXON $tmp[2]."\t".$exonstarts[$i]."\t".$exonends[$i]."\t1"."\n";
	}
}elsif(exists $hashchr{$tmp[2]} && $tmp[3] eq "-" && ($tmp[5]-$downstream)>0 && exists $hashname{$tmp[12]}){
	$hashrefgeneuniprot{$hashname{$tmp[12]}}="";
		if(!exists $hashpro{$tmp[2]."\t".($tmp[5]-$downstream)."\t".($tmp[5]+$upstream)."\t".$tmp[12]}){
		$hashpro{$tmp[2]."\t".($tmp[5]-$downstream)."\t".($tmp[5]+$upstream)."\t".$tmp[12]}="";
		print PRO $tmp[2]."\t".($tmp[5]-$downstream)."\t".($tmp[5]+$upstream)."\t".$tmp[12].":".$hashname{$tmp[12]}."\n";
		}
	print PROEXON $tmp[2]."\t".($tmp[5]-$downstream)."\t".($tmp[5]+$upstream)."\t1"."\n";
	@exonstarts=split/\,/,$tmp[9];
	@exonends=split/\,/,$tmp[10];
	for($i=0;$i<scalar(@exonstarts);$i++){
	print PROEXON $tmp[2]."\t".$exonstarts[$i]."\t".$exonends[$i]."\t1"."\n";
	}
}
}
close GENE;
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

