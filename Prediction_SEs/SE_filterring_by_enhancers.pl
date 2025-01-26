$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
mkdir("SEtrue");
mkdir("SEgene");
my $dir="e";
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
foreach $file (@dir){
@names=split/\./,$file;
$standardname=$file;
$standardname=~ s/\_enh\.bed$//g;
if(-e "tmp/".$standardname."/".$names[0]."_SuperStitched.table.txt"){
system("cp tmp/".$standardname."/".$names[0]."_SuperStitched_REGION_TO_GENE.txt SEgene");
open AA,"tmp/".$standardname."/".$names[0]."_SuperStitched.table.txt";
open BB,">tmp/".$standardname."/".$names[0]."_PutativeSuper.bed";
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if($tmp[1]=~ /^chr/ && $tmp[8] eq "1"){
print BB $tmp[1]."\t".$tmp[2]."\t".$tmp[3]."\t".$tmp[6]."\n";
}

}
close AA;
close BB;
system("bedtools sort -i tmp/".$standardname."/".$names[0]."_PutativeSuper.bed>SEtrue/".$names[0]."_PutativeSupersort.bed");
system("bedtools intersect -a SEtrue/".$names[0]."_PutativeSupersort.bed -b e/".$file." -wa -wb>SEtrue/".$names[0]."_intersect.bed");
#system("bedtools intersect -a SEtrue/".$names[0]."_PutativeSupersort.bed -b e/".$file." -wa -wb -f 0.1>SEtrue/".$names[0]."_intersect.bed");
open AA,"SEtrue/".$names[0]."_intersect.bed";
open BB,">SEtrue/".$names[0]."_Super.bed";
%hashSuper=();
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if(!exists $hashSuper{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\t".$tmp[3]}){
$hashSuper{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\t".$tmp[3]}="";
print BB $tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\t".$tmp[3]."\n";
}

}
close AA;
close BB;
unlink("tmp/".$standardname."/".$names[0]."_PutativeSuper.bed");
unlink("SEtrue/".$names[0]."_PutativeSupersort.bed");
unlink("SEtrue/".$names[0]."_intersect.bed");
}
}
