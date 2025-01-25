$start = time;

open CC,"counts.mm";
%hashpeakindexTOcellindex=();
while(<CC>){
s/\r|\n//g;
@tmp=split/\s+/,$_;
if(!/\%/ && $tmp[2] < 10000){
#print DD $hashindexTOpeak{$tmp[0]}."\t".$hashindexTOcell{$tmp[1]}."\t"."1"."\n";
if(!exists $hashpeakindexTOcellindex{$tmp[0]}){
$hashpeakindexTOcellindex{$tmp[0]}=$tmp[1];
}else{
$hashpeakindexTOcellindex{$tmp[0]}.=";".$tmp[1];
}

}
}
close CC;


open BB,"barcodes.txt";
$celln=1;
%hashindexTOcell=();
while(<BB>){
s/\r|\n//g;
if($_ ne ""){
$hashindexTOcell{$celln}=$_;
$celln++;
}
}
close BB;



open AA,"peaks.txt";
open BB,">peakstmp.bed";
$peakn=1;
%hashindexTOpeak=();
while(<AA>){
s/\r|\n//g;
s/\_/\t/g;
if($_ ne ""){
$hashindexTOpeak{$peakn}=$_;
print BB $_."\t".$peakn."\n";
$peakn++;
}
}
close AA;
close BB;
system("bedtools sort -i peakstmp.bed>peaks.bed");
unlink("peakstmp.bed");

open DD,"peaks.bed";
open EE,">fragments.tsv";
while(<DD>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if($_ ne ""){
@cellindexs=split/\;/,$hashpeakindexTOcellindex{$tmp[3]};
foreach $cellindex (@cellindexs){
print EE $tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\t".$hashindexTOcell{$cellindex}."\t"."1"."\n";
}

}
}
close DD;
close EE;

#bgzip -c fragments.tsv> integrated.fragments.tsv.gz
#system("bedtools sort -i fragments.tsv>integrated.fragments.tsv");
#unlink("fragments.tsv");
#system("bgzip -c integrated.fragments.tsv > integrated.fragments.tsv.gz");
system("bgzip -c fragments.tsv > integrated.fragments.tsv.gz");
system("tabix -p vcf integrated.fragments.tsv.gz");
#system("sort -k4,4 integrated.fragments.tsv > integratedsort.fragments.tsv");
#system("gzip integratedsort.fragments.tsv");
unlink("integrated.fragments.tsv");
unlink("peaks.bed");
#bgzip -c integrated.fragments.tsv > integrated.fragments.tsv.gz
##tabix -p vcf integrated.fragments.tsv.gz

#gunzip integratedsort.fragments.tsv.gz
#bgzip -c integratedsort.fragments.tsv > integratedsort.fragments.tsv.gz
##tabix -p vcf integratedsort.fragments.tsv.gz

$duration = time - $start;
print "All are done: $duration s\n";