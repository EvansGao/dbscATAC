#@projects=("GSE179705","GSE173834","GSE151230","GSE131688");
#@projspe=("hs","hs","hs","mm","dm","hs","hs","mm","hs","hs","hs","mm","hs","mm","hs","hs","hs","mm","mm");
#@ids=("GEO:GSE179705","GEO:GSE173834","GEO:GSE151230","GEO:GSE131688");
#@builds=("galGal6","mm10","TAIR10","danRer10");
#@spes=("gg","mm","at","dr");

@projects=("GSE111586","GSE199739","GSE149683","GSE214082","GSE163697","10X_MouseBrain");
#@projspe=("hs","hs","hs","mm","dm","hs","hs","mm","hs","hs","hs","mm","hs","mm","hs","hs","hs","mm","mm");
@ids=("GEO:GSE111586:30078704","GEO:GSE199739:35687698","GEO:GSE149683:33184180","GEO:GSE214082:37555319","GEO:GSE163697:34987221","10X:10X_MouseBrain:10X_MouseBrain");


@builds=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
@spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
%hashbuildTOspe=();
for($i=0;$i<scalar(@builds);$i++){
$hashbuildTOspe{$builds[$i]}=$spes[$i];
}

open TRUEC,"SampleAnnotation.txt";
%hashcellTOall=();
while(<TRUEC>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if($tmp[2] ne ""){
$hashcellTOall{$tmp[0]}=$_;
}
}
close TRUEC;


$commdir="/data/gts/dbscATAC/project";
@celltyperaws=();
%hashcelltyperawTOse=();
for($i=0;$i<scalar(@projects);$i++){
%hashcelltypeTOspe=();
my $rdsdir="../project/".$projects[$i]."Yes/rds";
opendir(DIR,$rdsdir) or "can't open the file";
@rdsdir=readdir DIR;
@rdsdir=grep{$_ ne "." && $_ ne ".."} @rdsdir;
	foreach $rds (@rdsdir){
	@rdsinfo=split/\./, $rds;
	$hashcelltypeTOspe{$rdsinfo[0]}=$rdsinfo[2];
	}
#get cell type information
my $SEtruedir=$commdir."/".$projects[$i]."Yes/SEtrue";
opendir(DIR,$SEtruedir) or "can't open the file";
@SEtruedir=readdir DIR;
@SEtruedir=grep{$_ ne "." && $_ ne ".."} @SEtruedir;
foreach $SEtrue (@SEtruedir){
$celltype=$SEtrue;
$celltype=~ s/\_Super\.bed//g;
$celltyperaw=$celltype.".".$projects[$i].".".$hashcelltypeTOspe{$celltype};
print $celltyperaw."\n";
$senum=&LINES($SEtruedir."/".$SEtrue);
push @celltyperaws,$celltyperaw;
$hashcelltyperawTOse{$celltyperaw}=$senum;
}
}

open BB,">SampleAnnotationSE.txt";
foreach $celltyperaw (@celltyperaws){
print BB $hashcellTOall{$celltyperaw}."\t".$hashcelltyperawTOse{$celltyperaw}."\n";
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
