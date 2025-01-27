$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
use 5.010;
#@spevers=("mm10","TAIR10","danRer10","galGal5","hg38","rheMac10","panTro5");
#@spes=("mm","at","dr","gg","hs","mam","pt");
@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
@spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");

%hashspeverTOspe=();
for($i=0;$i<scalar(@spevers);$i++){
$hashspeverTOspe{$spevers[$i]}=$spes[$i];
}


#@projects=("GSE179705","GSE173834","GSE151230","GSE131688");
#@projspe=("hs","hs","hs","mm","dm","hs","hs","mm","hs","hs","hs","mm","hs","mm","hs","hs","hs","mm","mm");
#@ids=("GEO:GSE179705","GEO:GSE173834","GEO:GSE151230","GEO:GSE131688");
##@projects=("GSE149683","GSE214082","GSE111586","GSE199739","10X_MouseBrain");
#@projspe=("hs","hs","hs","mm","dm","hs","hs","mm","hs","hs","hs","mm","hs","mm","hs","hs","hs","mm","mm");
##@ids=("GEO:GSE149683","GEO:GSE214082","GEO:GSE111586","10X:10X_MouseBrain");

@projects=("GSE111586","GSE199739","GSE149683","GSE214082","GSE163697","10X_MouseBrain");
#@projspe=("hs","hs","hs","mm","dm","hs","hs","mm","hs","hs","hs","mm","hs","mm","hs","hs","hs","mm","mm");
@ids=("GEO:GSE111586:30078704","GEO:GSE199739:35687698","GEO:GSE149683:33184180","GEO:GSE214082:37555319","GEO:GSE163697:34987221","10X:10X_MouseBrain:10X_MouseBrain");
@spevers=("mm10","mm10","hg38","mm10","dm6","mm10");

#system("rm -R download");
mkdir("download");
system("rm -R download/superenhancer");
system("rm -R download/SEgene");
mkdir("download/superenhancer");
mkdir("download/SEgene");
mkdir("download/deeptools");
#mkdir("download/enhancer/hs");
#mkdir("download/enhancer/mm");
#mkdir("download/enhancer/dm");
#mkdir("download/promoter");
#mkdir("download/promoter/hs");
#mkdir("download/promoter/mm");
#mkdir("download/promoter/dm");
#mkdir("download/interaction");
#mkdir("download/interaction/hs");
#mkdir("download/interaction/mm");
#mkdir("download/interaction/dm");
#mkdir("download/hundred");
#mkdir("download/hundred/hs");
#mkdir("download/hundred/mm");
#mkdir("download/hundred/dm");
#mkdir("download/peak");
#mkdir("download/peak/hs");
#mkdir("download/peak/mm");
#mkdir("download/peak/dm");
#mkdir("download/matrix");
#mkdir("download/matrix/hs");
#mkdir("download/matrix/mm");
#mkdir("download/matrix/dm");
open TRUEC,"SampleAnnotationSE.txt";
%hashcellTOtruecell=();
%hashspe=();
while(<TRUEC>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if($tmp[2] ne ""){
$hashcellTOtruecell{$tmp[0]}=$tmp[2];
}
}
close TRUEC;
#open AA,">SampleAnnotationpre.txt";

%hashspever=();
for($i=0;$i<scalar(@projects);$i++){
#@enhfiles=();
%hashcelltypeTOspe=();
%hashcelltypeSE=();
%hashcelltype=();
my $rdsdir="../project/".$projects[$i]."Yes/rds";
opendir(DIR,$rdsdir) or "can't open the file";
@rdsdir=readdir DIR;
@rdsdir=grep{$_ ne "." && $_ ne ".."} @rdsdir;
	foreach $rds (@rdsdir){
	@rdsinfo=split/\./, $rds;
	$hashcelltypeTOspe{$rdsinfo[0]}=$rdsinfo[2];
	}

%hashNMidTOgenename=();
open AA,"ROSErefseq/".$spevers[$i]."_refseq.ucsc";
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if($tmp[1] ne ""){
$hashNMidTOgenename{$tmp[1]}=$tmp[12];
}

}
close AA;



my $SEtruedir="../project/".$projects[$i]."Yes/SEtrue";
opendir(DIR,$SEtruedir) or "can't open the file";
@SEtruedir=readdir DIR;
@SEtruedir=grep{$_ ne "." && $_ ne ".."} @SEtruedir;
foreach $SEtrue (@SEtruedir){
$celltype=$SEtrue;
$celltype=~ s/\_Super\.bed//g;
	open AA,$SEtruedir."/".$SEtrue;
	while(<AA>){
	s/\r|\n//g;
	@tmp=split/\t/,$_;
	if($tmp[3] ne ""){
	$hashcelltypeSE{$celltype."\t".$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}=$tmp[3];
	}
	}
	close AA;
	if(!exists $hashcelltype{$celltype}){
	$hashcelltype{$celltype}="";
	}
$spe=$hashcelltypeTOspe{$celltype};
#print $spe."\n";
	if(!exists $hashspe{$spe}){
	$hashspe{$spe}="";
	mkdir("download/superenhancer/".$hashspeverTOspe{$spe});
	mkdir("download/SEgene/".$hashspeverTOspe{$spe});
	mkdir("download/deeptools/".$hashspeverTOspe{$spe});
	}
$truename=$hashcellTOtruecell{$celltype.".".$projects[$i].".".$spe};
#print $hashspeverTOspe{$spe}."\n";
copy("../project/".$projects[$i]."Yes/SEtrue/".$celltype."_Super.bed", "download/superenhancer/".$hashspeverTOspe{$spe}."/".$truename.".bed");

if(-e "../project/".$projects[$i]."Yes/deeptools/".$celltype.".".$projects[$i].".".$spe."_HeatmapSE.png"){
copy("../project/".$projects[$i]."Yes/deeptools/".$celltype.".".$projects[$i].".".$spe."_HeatmapSE.png", "download/deeptools/".$hashspeverTOspe{$spe}."/".$truename."_HeatmapSE.png");
}
if(-e "../project/".$projects[$i]."Yes/deeptools/".$celltype.".".$projects[$i].".".$spe."_HeatmapE.png"){
copy("../project/".$projects[$i]."Yes/deeptools/".$celltype.".".$projects[$i].".".$spe."_HeatmapE.png", "download/deeptools/".$hashspeverTOspe{$spe}."/".$truename."_HeatmapE.png");
}
if(-e "../project/".$projects[$i]."Yes/tmp/".$celltype.".".$projects[$i].".".$spe."/".$celltype."_Plot_points.png"){
copy("../project/".$projects[$i]."Yes/tmp/".$celltype.".".$projects[$i].".".$spe."/".$celltype."_Plot_points.png", "download/deeptools/".$hashspeverTOspe{$spe}."/".$truename."_Plotpoints.png");
}


}

my $SEgenedir="../project/".$projects[$i]."Yes/SEgene";
opendir(DIR,$SEgenedir) or "can't open the file";
@SEgenedir=readdir DIR;
@SEgenedir=grep{$_ ne "." && $_ ne ".."} @SEgenedir;
foreach $SEgene (@SEgenedir){
$celltype=$SEgene;
$celltype=~ s/\_SuperStitched\_REGION\_TO\_GENE\.txt//g;
$spe=$hashcelltypeTOspe{$celltype};
$truename=$hashcellTOtruecell{$celltype.".".$projects[$i].".".$spe};
	if(exists $hashcelltype{$celltype}){
	open AA,$SEgenedir."/".$SEgene;
	open BB,">download/SEgene/".$hashspeverTOspe{$spe}."/".$truename.".bed";
#	print "download/SEgene/".$hashspeverTOspe{$spe}."/".$truename.".bed\n";
	while(<AA>){
	s/\r|\n//g;
	@tmp=split/\t/,$_;
	@NMs=split/\,/,$tmp[6];
	@genenames=map {$hashNMidTOgenename{$_}} @NMs;
	@genenames=uniq(@genenames);
	if(exists $hashcelltypeSE{$celltype."\t".$tmp[1]."\t".$tmp[2]."\t".$tmp[3]}){
	print BB $tmp[1]."\t".$tmp[2]."\t".$tmp[3]."\t".$hashNMidTOgenename{$tmp[8]}."\t".join(",",@genenames)."\n";
#		if($truename eq "Vip"){
#		print $hashNMidTOgenename{"NR_046233"}."good\n";
#		}
	
	}
	}
	close AA;
	close BB;
	}


}


#print $enhinfo[0].$hashcellTOtruecell{$enhinfo[0]}."OK\n";
#$truename=$hashcellTOtruecell{$enhinfo[0].".".$enhinfo[1].".".$spe};
#copy("../project/".$projects[$i]."Yes/rds/".$enhinfo[0].".".$projects[$i].".".$spe.".rds", "download/matrix/".$hashspeverTOspe{$spe}."/".$truename.".rds");
#copy("../project/".$projects[$i]."Yes/p/".$enhinfo[0].".".$projects[$i].".".$spe."_pro.bed", "download/promoter/".$hashspeverTOspe{$spe}."/".$truename."_pro.bed");
#copy("../project/".$projects[$i]."Yes/e/".$enhinfo[0].".".$projects[$i].".".$spe."_enh.bed", "download/enhancer/".$hashspeverTOspe{$spe}."/".$truename.".bed");
#copy("../project/".$projects[$i]."Yes/i/".$enhinfo[0].".".$projects[$i].".".$spe."_interaction.txt", "download/interaction/".$hashspeverTOspe{$spe}."/".$truename."_interaction.txt");
#copy("../project/".$projects[$i]."Yes/hundred/".$enhinfo[0].".".$projects[$i].".".$spe.".txt", "download/hundred/".$hashspeverTOspe{$spe}."/".$truename.".txt");
#copy("../project/".$projects[$i]."Yes/peak/".$hashcellTOnum{$enhinfo[0]}."_".$hashcellTOfragnum{$enhinfo[0]}."_".$enhinfo[0].".".$projects[$i].".".$spe.".txt", "download/peak/".$hashspeverTOspe{$spe}."/".$truename.".txt");


}
#system("tar -cjf matrix.tbz download/matrix");
#close AA;
$duration = time - $start;
print "All are done: $duration s\n";