$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
use 5.010;
#@spevers=("mm10","TAIR10","danRer10","galGal5","hg38","rheMac10","panTro5");
#@spes=("mm","at","dr","gg","hs","mam","pt");
@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
@spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
%hashspeverTOspe=();
system("rm -R download/markers");
system("rm -R rnanewsub");
mkdir("download/markers");
mkdir("rnanewsub");

for($i=0;$i<scalar(@spevers);$i++){
$hashspeverTOspe{$spevers[$i]}=$spes[$i];
}



#@projects=("GSE179705","GSE173834","GSE151230","GSE131688");
#@projspe=("hs","hs","hs","mm","dm","hs","hs","mm","hs","hs","hs","mm","hs","mm","hs","hs","hs","mm","mm");
#@ids=("GEO:GSE179705","GEO:GSE173834","GEO:GSE151230","GEO:GSE131688");
###@projects=("GSE204761","GSE196794","GSE199739","GSE228270","GSE214082","GSE192772","GSE179705","GSE151230","GSE111586","GSE131688","GSE155178","GSE163697","GSE149683","GSE214132","10X_MouseBrain","GSE173834");
#@projspe=("hs","hs","hs","mm","dm","hs","hs","mm","hs","hs","hs","mm","hs","mm","hs","hs","hs","mm","mm");
###@ids=("GEO:GSE204761","GEO:GSE196794","GEO:GSE199739","GEO:GSE228270","GEO:GSE214082","GEO:GSE192772","GEO:GSE179705","GEO:GSE151230","GEO:GSE111586","GEO:GSE131688","GEO:GSE155178","GEO:GSE163697","GEO:GSE149683","GEO:GSE214132","10X:10X_MouseBrain","GEO:GSE173834");


@projects=("GSE204761","GSE196794","GSE228270","GSE214082","GSE192772","GSE179705","GSE151230","GSE111586","GSE155178","GSE163697","GSE149683","GSE214132","10X_MouseBrain","GSE173834","GSE129785","GSE139369","GSE162690","10X_HumanPBMC");
#@projspe=("hs","hs","hs","mm","dm","hs","hs","mm","hs","hs","hs","mm","hs","mm","hs","hs","hs","mm","mm");
@ids=("GEO:GSE204761","GEO:GSE196794","GEO:GSE228270","GEO:GSE214082","GEO:GSE192772","GEO:GSE179705","GEO:GSE151230","GEO:GSE111586","GEO:GSE155178","GEO:GSE163697","GEO:GSE149683","GEO:GSE214132","10X:10X_MouseBrain","GEO:GSE173834","GEO:GSE129785","GEO:GSE139369","GEO:GSE162690","10X:10X_HumanPBMC");

#mkdir("download/enhancer/hs");
#mkdir("download/enhancer/mm");
#mkdir("download/enhancer/dm");

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
%hashspe=();
open FLAG,">markercellflag.txt";
for($i=0;$i<scalar(@projects);$i++){
%hashcells=();
my $enhdir="../project/".$projects[$i]."Yes/e";
opendir(DIR,$enhdir) or "can't open the file";
@enhdir=readdir DIR;
@enhdir=grep{$_ ne "." && $_ ne ".."} @enhdir;
foreach $enhfile (@enhdir){
$celloriname=$enhfile;
$celloriname=~ s/\_enh\.bed//g;
$hashcells{$celloriname}=$enhfile;
@enhinfo=split/\./, $enhfile;
$spe=$enhinfo[2];
$spe=~ s/\_enh//g;
if(!exists $hashspe{$spe}){
$hashspe{$spe}="";
mkdir("download/markers/".$hashspeverTOspe{$spe});
}
}

my $dir="../project/".$projects[$i]."Yes";
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
	foreach $file (@dir){
		if($file=~ /^rnanewmarkers.*$projects[$i].*\.txt$/i){
		$rnanewsub="NA";
		$rnanewsubtmp=$file;
		$rnanewsubtmp=~ s/txt$/rds/g;
		$rnanewsubtmp=~ s/^rnanewmarkers/rnanewsub/g;
		if(-e $dir."/".$rnanewsubtmp){
		$rnanewsub=$rnanewsubtmp;
		copy($dir."/".$rnanewsub, "rnanewsub/".$rnanewsub);
		}
		@info=split/\./,$file;
		$info[0]=~ s/rnanewmarkers\_//g;
		%hashfileexist=();
		@filecells=();
		$filenum=10;
		open AA, $dir."/".$file;
		while(<AA>){
		s/\r|\n//g;
		@tmp=split/\t/,$_;
		$celltmpname=$info[0]."_".$tmp[6].".".$projects[$i].".".$info[2];
		if($projects[$i] eq "GSE192772" || $projects[$i] eq "GSE151230" || $projects[$i] eq "GSE173834"){
		$celltmpname=$tmp[6].".".$projects[$i].".".$info[2];
		}elsif($projects[$i] eq "GSE179705"){
		$celltmpname=$projects[$i]."_".$tmp[6].".".$projects[$i].".".$info[2];
#		print $celltmpname."aa\n";
		}
			if(exists $hashcells{$celltmpname} && !exists $hashfileexist{$celltmpname}){
			$hashfileexist{$celltmpname}=$filenum;

			print FLAG $tmp[6]."\t".$hashspeverTOspe{$info[2]}."\t".$hashcellTOtruecell{$celltmpname}."\t".$celltmpname."\t".$file."\t".$rnanewsub."\n";
			open $filenum,">download/markers/".$hashspeverTOspe{$info[2]}."/".$hashcellTOtruecell{$celltmpname}.".txt";
			print $filenum $_."\n";
			push @filecells,$celltmpname;
			$filenum++;
			}elsif(exists $hashcells{$celltmpname} && exists $hashfileexist{$celltmpname}){
			$tmpfilenum=$hashfileexist{$celltmpname};
			print $tmpfilenum $_."\n";
			}
		}
		close AA;
		foreach $filecell (@filecells){
		$tmpfilenum=$hashfileexist{$filecell};
		close $tmpfilenum;
		}
		}
	}
}
close FLAG;

$duration = time - $start;
print "All are done: $duration s\n";