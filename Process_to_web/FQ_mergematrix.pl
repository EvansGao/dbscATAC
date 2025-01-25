@mtxnames=("Liver_Endothelial_II_cells_MM1","Lung_Endothelial_II_cells_MM1");
@barc=("AAA","BBB","CCC","DDD","EEE","FFF","GGG","HHH","JJJ","KKK");
%hashcellpeak=();
%hashcellTOlength=();
open LAB,">label.txt";
open BB,">allpeakspre.bed";
open CC,">mtx.txt";
for($i=0;$i<scalar(@mtxnames);$i++){
open AA,"download/hundred/mm/".$mtxnames[$i].".txt";
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
@peakinfo=split/\_/,$tmp[0];
if($tmp[0] ne ""){
print BB join("\t",@peakinfo)."\t".$tmp[0]."|".$mtxnames[$i]."\n";
$hashcellpeak{$tmp[0]."|".$mtxnames[$i]}=join("\t",@tmp[1..$#tmp]);
	if(!exists $hashcellTOlength{$mtxnames[$i]}){
	$hashcellTOlength{$mtxnames[$i]}=scalar(@tmp)-1;
	}
}
}
close AA;
for($j=0;$j<$hashcellTOlength{$mtxnames[$i]};$j++){
print LAB $barc[$i].$j."\t".$mtxnames[$i]."\n";
print CC "\t".$barc[$i].$j;
}

}
print CC "\n";
close BB;
close LAB;

system("bedtools sort -i allpeakspre.bed>allpeakspresort.bed");
system("bedtools merge -i allpeakspresort.bed -c 4 -o collapse>allpeaks.bed");
unlink("allpeakspresort.bed");
unlink("allpeakspre.bed");
open AA,"allpeaks.bed";
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
%hashtmpcell=();
	@tmpcellpeaks=split/\,/,$tmp[3];
	foreach $tmpcellpeak (@tmpcellpeaks){
		$tmpcell=$tmpcellpeak;
		$tmpcell=~ s/^.*\|//g;
		if(exists $hashcellpeak{$tmpcellpeak}){
		$hashtmpcell{$tmpcell}=$hashcellpeak{$tmpcellpeak};
		}
	}
	print CC $tmp[0]."_".$tmp[1]."_".$tmp[2];
	foreach $mtxname (@mtxnames){
		if(exists $hashtmpcell{$mtxname}){
		print CC "\t".$hashtmpcell{$mtxname};
		}else{
		$tmpzero="";
			for($j=0;$j<$hashcellTOlength{$mtxname};$j++){
			$tmpzero.="0"."\t";
			}
		$tmpzero=~ s/\t$//g;
		print CC "\t".$tmpzero;
		}
	
	}
	print CC "\n";
}
close AA;
close CC;

