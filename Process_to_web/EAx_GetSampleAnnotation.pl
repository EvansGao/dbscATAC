#@projects=("GSE179705","GSE173834","GSE151230","GSE131688");
#@projspe=("hs","hs","hs","mm","dm","hs","hs","mm","hs","hs","hs","mm","hs","mm","hs","hs","hs","mm","mm");
#@ids=("GEO:GSE179705","GEO:GSE173834","GEO:GSE151230","GEO:GSE131688");
#@builds=("galGal6","mm10","TAIR10","danRer10");
#@spes=("gg","mm","at","dr");

@projects=("GSE204761","GSE196794","GSE199739","GSE228270","GSE214082","GSE192772","GSE179705","GSE151230","GSE111586","GSE131688","GSE155178","GSE163697","GSE149683","GSE214132","10X_MouseBrain","GSE173834","GSE129785","GSE139369","GSE162690","10X_HumanPBMC");
#@projspe=("hs","hs","hs","mm","dm","hs","hs","mm","hs","hs","hs","mm","hs","mm","hs","hs","hs","mm","mm");
@ids=("GEO:GSE204761:35915179","GEO:GSE196794:35831300","GEO:GSE199739:35687698","GEO:GSE228270:37794584","GEO:GSE214082:37555319","GEO:GSE192772:37468639","GEO:GSE179705:NA","GEO:GSE151230:34170284","GEO:GSE111586:30078704","GEO:GSE131688:31639368","GEO:GSE155178:33964211","GEO:GSE163697:34987221","GEO:GSE149683:33184180","GEO:GSE214132:36482454","10X:10X_MouseBrain:10X_MouseBrain","GEO:GSE173834:34099698","GEO:GSE129785:31375813","GEO:GSE139369:37597510","GEO:GSE162690:33633365","10X:10X_HumanPBMC:10X_HumanPBMC");

@builds=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
@spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
%hashbuildTOspe=();
for($i=0;$i<scalar(@builds);$i++){
$hashbuildTOspe{$builds[$i]}=$spes[$i];
}

$commdir="/data/gts/dbscATAC/project";
@celltypes=();
%hashcelltypeTOspe=();
%hashcelltypeTOsource=();
%hashcelltypeTOcellnum=();
%hashcelltypeTOfragnum=();
%hashcelltypeTOinternum=();
%hashcelltypeTOenhnum=();
%hashcelltypeTOpronum=();

for($i=0;$i<scalar(@projects);$i++){

#get cell type information
my $interdir=$commdir."/".$projects[$i]."Yes/i";
opendir(DIR,$interdir) or "can't open the file";
@interdir=readdir DIR;
@interdir=grep{$_ ne "." && $_ ne ".."} @interdir;
foreach $inter (@interdir){
@interinfo=split/\./, $inter;
$tmpcelltype=$inter;
$tmpcelltype=~ s/\_interaction\.txt//g;
push @celltypes,$tmpcelltype;
$spe=$interinfo[2];
$spe=~ s/\_interaction//g;
$hashcelltypeTOspe{$tmpcelltype}=$hashbuildTOspe{$spe};
$hashcelltypeTOsource{$tmpcelltype}=$projects[$i];
$hashcelltypeTOinternum{$tmpcelltype}=&LINES($interdir."/".$inter);

}

my $enhdir=$commdir."/".$projects[$i]."Yes/e";
opendir(DIR,$enhdir) or "can't open the file";
@enhdir=readdir DIR;
@enhdir=grep{$_ ne "." && $_ ne ".."} @enhdir;
foreach $enh (@enhdir){
@enhinfo=split/\./, $enh;
$tmpcelltype=$enh;
$tmpcelltype=~ s/\_enh\.bed//g;
$hashcelltypeTOenhnum{$tmpcelltype}=&LINES($enhdir."/".$enh);
#print $hashcelltypeTOenhnum{$enhinfo[0]}."OK\n";
}

my $prodir=$commdir."/".$projects[$i]."Yes/p";
opendir(DIR,$prodir) or "can't open the file";
@prodir=readdir DIR;
@prodir=grep{$_ ne "." && $_ ne ".."} @prodir;
foreach $pro (@prodir){
@proinfo=split/\./, $pro;
$tmpcelltype=$pro;
$tmpcelltype=~ s/\_pro\.bed//g;
$hashcelltypeTOpronum{$tmpcelltype}=&PROLINES($prodir."/".$pro);
}


my $peakdir=$commdir."/".$projects[$i]."Yes/peak";
opendir(DIR,$peakdir) or "can't open the file";
@peakdir=readdir DIR;
@peakdir=grep{$_ ne "." && $_ ne ".."} @peakdir;
	foreach $file (@peakdir){
		if($file=~ /^(\d+)\_(\d+)\_(.*)\.txt/i){
		$hashcelltypeTOfragnum{$3}=$2;
		$hashcelltypeTOcellnum{$3}=$1;
#		print $3."\n";
		}
	}


#open AA, $commdir."/".$projects[$i]."Yes/sta.txt";
#while(<AA>){
#s/\r|\n//g;
#@tmp=split/\t/,$_;
#if($tmp[0] ne ""){
#$hashcelltypeTOcellnum{$tmp[0]}=$tmp[1];
#}
#
#}
#close AA;

}

open BB,">SampleAnnotationprenew.txt";
foreach $celltype (@celltypes){
print BB $celltype."\t".$hashcelltypeTOspe{$celltype}."\t".$celltype."\t".""."\t".""."\t".""."\t".$hashcelltypeTOcellnum{$celltype}."\t".$hashcelltypeTOsource{$celltype}."\t".$hashcelltypeTOfragnum{$celltype}."\t".$hashcelltypeTOenhnum{$celltype}."\t".$hashcelltypeTOpronum{$celltype}."\t".$hashcelltypeTOinternum{$celltype}."\n";
#print BB $celltype."\t".$hashcelltypeTOspe{$celltype}."\t".$celltype."\t".$describe."\t".$tissue."\t".$disease."\t".$hashcellTOnum{$cellname}."\t".$ids[$i]."\t".$hashcellTOfragnum{$cellname}."\t".$enhnum."\t".$pronum."\t".$internum."\n";


}
close BB;


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
}
close FILE;
return scalar(keys %hashtmp);
}
