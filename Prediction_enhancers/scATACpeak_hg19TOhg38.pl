$origfile= $ARGV[0];
open AA,$origfile;
open BB,">".$origfile."_hg19.bed";
while(<AA>){
s/\r|\n//g;
if(/^chr\d+\_\d+|^chrX\_\d+|^chrY\_\d+/i){
@tmp=split/\-|\_/,$_;
print BB $tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\t".$_."\n";
}
}
close AA;
close BB;
system("liftOver ".$origfile."_hg19.bed ../liftover/hg19ToHg38.over.chain ".$origfile."_hg38.bed unMapped");
open AA,$origfile."_hg38.bed";
open BB,">".$origfile."_hg38";
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
print BB $tmp[3]."\t".$tmp[0]."_".$tmp[1]."_".$tmp[2]."\n";

}
close AA;
close BB;
unlink($origfile."_hg19.bed");
unlink($origfile."_hg38.bed");