$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
use 5.010;

@projects=("33164753","28825706","31375813","30078704","29539636","CRCSignac","PBMCSignac","31594952","33184180","31594952R","30166440","33208554","33427645","SignacmmBrain","31792411","33633365","29706549","29497140","34049878");
@projspe=("hs","hs","hs","mm","dm","hs","hs","mm","hs","hs","hs","mm","hs","mm","hs","hs","hs","mm","mm");
@ids=("PMID:33164753","PMID:28825706","PMID:31375813","PMID:30078704","PMID:29539636","Signac","Signac","PMID:31594952","PMID:33184180","PMID:31594952","PMID:30166440","PMID:33208554","PMID:33427645","Signac","PMID:31792411","PMID:33633365","PMID:29706549","PMID:29497140","PMID:34049878");
mkdir("download");
mkdir("download/enhancer");
mkdir("download/enhancer/hs");
mkdir("download/enhancer/mm");
mkdir("download/enhancer/dm");
mkdir("download/promoter");
mkdir("download/promoter/hs");
mkdir("download/promoter/mm");
mkdir("download/promoter/dm");
mkdir("download/interaction");
mkdir("download/interaction/hs");
mkdir("download/interaction/mm");
mkdir("download/interaction/dm");
mkdir("download/hundred");
mkdir("download/hundred/hs");
mkdir("download/hundred/mm");
mkdir("download/hundred/dm");
mkdir("download/peak");
mkdir("download/peak/hs");
mkdir("download/peak/mm");
mkdir("download/peak/dm");

open TRUEC,"SampleAnnotation0916.txt";
%hashcellTOtruecell=();
%hashcellTOdescribe=();
%hashcellTOtissue=();
%hashcellTOdisease=();
while(<TRUEC>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if($tmp[2] ne ""){
$hashcellTOtruecell{$tmp[0]}=$tmp[2];
$hashcellTOdescribe{$tmp[0]}=$tmp[3];
$hashcellTOtissue{$tmp[0]}=$tmp[4];
$hashcellTOdisease{$tmp[0]}=$tmp[5];
}
}
close TRUEC;

open AA,">SampleAnnotationpre.txt";
for($i=0;$i<scalar(@projects);$i++){
#@enhfiles=();
%hashcellTOnum=();
%hashcellTOfragnum=();
my $dir=$projects[$i]."/peak";
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
	foreach $file (@dir){
		if($file=~ /^(\d+)\_(\d+)\_(.*)\_$projects[$i]\.txt/i){
		$hashcellTOnum{$3}=$1;
		$hashcellTOfragnum{$3}=$2;
		}
	}
my $dir=$projects[$i]."/e";
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
	foreach $file (@dir){
		if($file=~ /(.*)\_$projects[$i]\_enh\.bed/i){
		$cellname=$1;
		$truename="";
		$describe=$cellname;
		$tissue="";
		$disease="";
		if(exists $hashcellTOtruecell{$cellname}){
		$truename=$hashcellTOtruecell{$cellname};
		$describe=$hashcellTOdescribe{$cellname};
		$tissue=$hashcellTOtissue{$cellname};
		$disease=$hashcellTOdisease{$cellname};
		}

		$enhnum=&LINES($dir."/".$file);
		$pronum=&PROLINES($projects[$i]."/p/".$cellname."_".$projects[$i]."_pro.bed");
		$internum=&LINES($projects[$i]."/i/".$cellname."_".$projects[$i]."_interaction.txt");
		
		print AA $cellname."\t".$projspe[$i]."\t".$truename."\t".$describe."\t".$tissue."\t".$disease."\t".$hashcellTOnum{$cellname}."\t".$ids[$i]."\t".$hashcellTOfragnum{$cellname}."\t".$enhnum."\t".$pronum."\t".$internum."\n";
#		copy($dir."/".$file, "download/enhancer/".$projspe[$i]."/".$hashcellTOtruecell{$cellname}.".bed");
#		copy($projects[$i]."/i/".$cellname."_".$projects[$i]."_interaction.txt", "download/interaction/".$projspe[$i]."/".$hashcellTOtruecell{$cellname}."_interaction.txt");
#		copy($projects[$i]."/hundred/hundred_".$cellname."_".$projects[$i].".txt", "download/hundred/".$projspe[$i]."/".$hashcellTOtruecell{$cellname}.".txt");
		}
	}

}
close AA;
$duration = time - $start;
print "All are done: $duration s\n";

sub  LINES()
{   
    my ($path)=@_;
$linenum=0;
open FILE,$path;
while (<FILE>) { $linenum++; }
close FILE;
return $linenum;
}


sub  PROLINES()
{   
    my ($path)=@_;
$linenum=0;
open FILE,$path;
%hashtmp=();
while (<FILE>) {
s/\r|\n//g;
s/\|$//g;
@tmp=split/\t/,$_;
@geneinfo=split/\|/,$tmp[3];
#@geneinfo=@geneinfo[1..$#geneinfo];
if($geneinfo[1] ne ""){
$hashtmp{$geneinfo[1]}="";
}
#foreach $gene (@geneinfo){
#$hashtmp{$gene}="";
#}
}
close FILE;

return scalar(keys %hashtmp);
}