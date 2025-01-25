$start = time;
open AA,"scCells.txt";
%hashmarkTOproject=();
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if($tmp[0] ne ""){
$hashmarkTOproject{$tmp[0]}=$tmp[1];
}

}
close AA;


open AA,"cellsID.txt";
open BB,">SampleAnnotationQuestion.txt";
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
@info=split/\_/,$tmp[0];
$mark=@info[$#info];
$tmp[0]=~ s/^\s+|\s+$//g;
$tmp[1]=~ s/^\s+|\s+$//g;
$tmp[2]=~ s/^\s+|\s+$//g;
if($tmp[0] ne ""){
print BB $tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\t".$hashmarkTOproject{$mark}."\t".$hashmarkTOproject{$mark}."\t".$hashmarkTOproject{$mark}."\t".$hashmarkTOproject{$mark}."\t".$hashmarkTOproject{$mark}."\n";
}
}
close AA;
close BB;
$duration = time - $start;
print "All are done: $duration s\n";