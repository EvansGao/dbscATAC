#@spes=("chicken","human","mouse","arabidopsis");
#@spetags=("ENSDARG","ENSG","ENSMUSG","AT");
@spes=("marmoset","mouse","cynomolgus","vervet","human","chicken","zebrafish","chimp","rhesus","fly","arabidopsis","rice","maize");
#@spetags=("ENSGALG","ENSDARG","AT","ENSMUSG");
#$spe="fly";
#$spetag="FBgn";
for($i=0;$i<scalar(@spes);$i++){
open AA,"geneMapdata".$spes[$i].".txt";
%hashIDtoInfo=();
@genes=();
while(<AA>){
s/\r|\n//g;
@temp=split/\t/,$_;
#if(($temp[0]=~ /$spetags[$i]/i) && $temp[0] ne ""){
if($temp[0] ne ""){
$hashIDtoInfo{$temp[0]}=$temp[1];
push @genes,$temp[0];
}
}
close AA;


%hashGeneIDtoCell=();
open CELL,"geneinfo".$spes[$i].".txt";
while(<CELL>){
s/\r|\n//g;
@temp=split/\t/,$_;
if(!exists $hashGeneIDtoCell{$temp[0]}){
$hashGeneIDtoCell{$temp[0]}=$temp[1];
}
}
close CELL;



open II,">geneforsearch".$spes[$i].".txt";
foreach $gene (@genes){
if(!exists $hashGeneIDtoCell{$gene}){
print II $gene."\t".$hashIDtoInfo{$gene}."\t"."NO"."\n";
}else{
print II $gene."\t".$hashIDtoInfo{$gene}."\t".$hashGeneIDtoCell{$gene}."\n";
}
}
close II;

}

