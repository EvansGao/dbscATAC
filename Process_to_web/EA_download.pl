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
@projects=("GSE204761","GSE196794","GSE199739","GSE228270","GSE214082","GSE192772","GSE179705","GSE151230","GSE111586","GSE131688","GSE155178","GSE163697","GSE149683","GSE214132","10X_MouseBrain","GSE173834","GSE129785","GSE139369","GSE162690","10X_HumanPBMC");
#@projspe=("hs","hs","hs","mm","dm","hs","hs","mm","hs","hs","hs","mm","hs","mm","hs","hs","hs","mm","mm");
#@ids=("GEO:GSE204761","GEO:GSE196794","GEO:GSE199739","GEO:GSE228270","GEO:GSE214082","GEO:GSE192772","GEO:GSE179705","GEO:GSE151230","GEO:GSE111586","GEO:GSE131688","GEO:GSE155178","GEO:GSE163697","GEO:GSE149683","GEO:GSE214132","10X:10X_MouseBrain","GEO:GSE173834");
@ids=("GEO:GSE204761:35915179","GEO:GSE196794:35831300","GEO:GSE199739:35687698","GEO:GSE228270:37794584","GEO:GSE214082:37555319","GEO:GSE192772:37468639","GEO:GSE179705:NA","GEO:GSE151230:34170284","GEO:GSE111586:30078704","GEO:GSE131688:31639368","GEO:GSE155178:33964211","GEO:GSE163697:34987221","GEO:GSE149683:33184180","GEO:GSE214132:36482454","10X:10X_MouseBrain:10X_MouseBrain","GEO:GSE173834:34099698","GEO:GSE129785:31375813","GEO:GSE139369:37597510","GEO:GSE162690:33633365","10X:10X_HumanPBMC:10X_HumanPBMC");
#system("rm -R download");
mkdir("download");
mkdir("download/enhancer");
#mkdir("download/enhancer/hs");
#mkdir("download/enhancer/mm");
#mkdir("download/enhancer/dm");
mkdir("download/promoter");
#mkdir("download/promoter/hs");
#mkdir("download/promoter/mm");
#mkdir("download/promoter/dm");
mkdir("download/interaction");
#mkdir("download/interaction/hs");
#mkdir("download/interaction/mm");
#mkdir("download/interaction/dm");
mkdir("download/hundred");
#mkdir("download/hundred/hs");
#mkdir("download/hundred/mm");
#mkdir("download/hundred/dm");
mkdir("download/peak");
#mkdir("download/peak/hs");
#mkdir("download/peak/mm");
#mkdir("download/peak/dm");
mkdir("download/matrix");
#mkdir("download/matrix/hs");
#mkdir("download/matrix/mm");
#mkdir("download/matrix/dm");
open TRUEC,"SampleAnnotation.txt";
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

for($i=0;$i<scalar(@projects);$i++){
#@enhfiles=();
%hashcellTOnum=();
%hashcellTOfragnum=();
my $dir="../project/".$projects[$i]."Yes/peak";
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
	foreach $file (@dir){
		if($file=~ /^(\d+)\_(\d+)\_(.*)\_$projects[$i]/i){
		$hashcellTOnum{$3}=$1;
		$hashcellTOfragnum{$3}=$2;
		}
	}

my $enhdir="../project/".$projects[$i]."Yes/e";
opendir(DIR,$enhdir) or "can't open the file";
@enhdir=readdir DIR;
@enhdir=grep{$_ ne "." && $_ ne ".."} @enhdir;
foreach $enh (@enhdir){
@enhinfo=split/\./, $enh;
$spe=$enhinfo[2];
$spe=~ s/\_enh//g;
if(!exists $hashspe{$spe}){
$hashspe{$spe}="";
mkdir("download/enhancer/".$hashspeverTOspe{$spe});
mkdir("download/promoter/".$hashspeverTOspe{$spe});
mkdir("download/interaction/".$hashspeverTOspe{$spe});
mkdir("download/hundred/".$hashspeverTOspe{$spe});
mkdir("download/peak/".$hashspeverTOspe{$spe});
mkdir("download/matrix/".$hashspeverTOspe{$spe});
}

#print $enhinfo[0].$hashcellTOtruecell{$enhinfo[0]}."OK\n";
$truename=$hashcellTOtruecell{$enhinfo[0].".".$enhinfo[1].".".$spe};
copy("../project/".$projects[$i]."Yes/rds/".$enhinfo[0].".".$projects[$i].".".$spe.".rds", "download/matrix/".$hashspeverTOspe{$spe}."/".$truename.".rds");
copy("../project/".$projects[$i]."Yes/p/".$enhinfo[0].".".$projects[$i].".".$spe."_pro.bed", "download/promoter/".$hashspeverTOspe{$spe}."/".$truename."_pro.bed");
copy("../project/".$projects[$i]."Yes/e/".$enhinfo[0].".".$projects[$i].".".$spe."_enh.bed", "download/enhancer/".$hashspeverTOspe{$spe}."/".$truename.".bed");
copy("../project/".$projects[$i]."Yes/i/".$enhinfo[0].".".$projects[$i].".".$spe."_interaction.txt", "download/interaction/".$hashspeverTOspe{$spe}."/".$truename."_interaction.txt");
copy("../project/".$projects[$i]."Yes/hundred/hundred_".$enhinfo[0].".".$projects[$i].".".$spe.".txt", "download/hundred/".$hashspeverTOspe{$spe}."/".$truename.".txt");
copy("../project/".$projects[$i]."Yes/peak/".$hashcellTOnum{$enhinfo[0]}."_".$hashcellTOfragnum{$enhinfo[0]}."_".$enhinfo[0].".".$projects[$i].".".$spe.".txt", "download/peak/".$hashspeverTOspe{$spe}."/".$truename.".txt");

}

}
#system("tar -cjf matrix.tbz download/matrix");
#close AA;
$duration = time - $start;
print "All are done: $duration s\n";