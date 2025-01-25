$start = time;
open AA,"SampleAnnotation.txt";
@groups=();
%hashgroupTocells=();
%hashgroupTocellenhs=();
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if(!exists $hashgroupTocells{$tmp[1]."\t".$tmp[4]}){
$hashgroupTocells{$tmp[1]."\t".$tmp[4]}=$tmp[2];
$hashgroupTocellenhs{$tmp[1]."\t".$tmp[4]}=&LINES($tmp[1],$tmp[2].".bed");
push @groups,$tmp[1]."\t".$tmp[4];
}else{
$hashgroupTocells{$tmp[1]."\t".$tmp[4]}.=";".$tmp[2];
$hashgroupTocellenhs{$tmp[1]."\t".$tmp[4]}.=";".&LINES($tmp[1],$tmp[2].".bed");
}
}

open SPETISSUE,">speciestissueinfo.txt";
foreach $group (@groups){
print SPETISSUE $group."\t".$hashgroupTocells{$group}."\t".$hashgroupTocellenhs{$group}."\n";
}
close SPETISSUE;


$duration = time - $start;
print "All are done: $duration s\n";

sub  LINES()
{   
    my ($spe,$filename)=@_;
$linenum=0;
open FILE,"download/enhancer/".$spe."/".$filename;
while (<FILE>) { $linenum++; }
close FILE;
return $linenum;
}