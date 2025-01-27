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

open TRUEC,"SampleAnnotation.txt";
%hashspecellTOall=();
while(<TRUEC>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if($tmp[2] ne ""){
$hashspecellTOall{$tmp[1]."\t".$tmp[2]}=$_;
}
}
close TRUEC;

@celltypes=();
%hashcelltypeTOmarker=();
my $commdir="/data/gts/dbscATAC/DataProcess/download/markers";
opendir(DIR,$commdir) or "can't open the file";
@firstdir=readdir DIR;
@firstdir=grep{$_ ne "." && $_ ne ".."} @firstdir;
for($i=0;$i<scalar(@firstdir);$i++){
	my $seconddir=$commdir."/".$firstdir[$i];
	opendir(DIR,$seconddir) or "can't open the file";
	@seconddir=readdir DIR;
	@seconddir=grep{$_ ne "." && $_ ne ".."} @seconddir;
	for($j=0;$j<scalar(@seconddir);$j++){
	$tmpcelltype=$seconddir[$j];
	$tmpcelltype=~ s/\.txt$//g;
	$hashcelltypeTOmarker{$firstdir[$i]."\t".$tmpcelltype}=&LINES($commdir."/".$firstdir[$i]."/".$seconddir[$j]);
	push @celltypes,$firstdir[$i]."\t".$tmpcelltype;
	}
}

open BB,">SampleAnnotationM.txt";
foreach $celltype (@celltypes){
print BB $hashspecellTOall{$celltype}."\t".$hashcelltypeTOmarker{$celltype}."\n";
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


$duration = time - $start;
print "All are done: $duration s\n";
