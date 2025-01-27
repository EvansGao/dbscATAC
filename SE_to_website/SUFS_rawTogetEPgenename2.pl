use List::MoreUtils qw(uniq);
#@spes=("human","mouse","fly");
#@spes=("mm","dr","at","gg");
@spes=("hs","mm","dm");
mkdir("SEEPgenename2");
for($i=0;$i<scalar(@spes);$i++){
#print $spes[$i]."\n";
mkdir("SEEPgenename2/".$spes[$i]);
my $dir="AllSUEPs/".$spes[$i];
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
foreach $filedir (@dir){
if($filedir ne "." && $filedir ne ".."){
	$cell=$filedir;
	$cell=~ s/\_EP\.txt//g;
	open CELL,$dir."/".$filedir;
	@enhs=();
	%hashenh=();
	while(<CELL>){
	s/\r|\n//g;
	@tmp=split/\t/,$_;
	@temp=split/\:|\-|\_|\$/,$tmp[0];
		if(!exists $hashenh{$temp[0]."\t".$temp[1]."\t".$temp[2]}){
		$hashenh{$temp[0]."\t".$temp[1]."\t".$temp[2]}=$temp[4];
		push @enhs,$temp[0]."\t".$temp[1]."\t".$temp[2];
		}else{
		$hashenh{$temp[0]."\t".$temp[1]."\t".$temp[2]}.=",".$temp[4];
		}
	
	}
	close CELL;
	open AA,">".$cell.".bed";
	foreach $enh (@enhs){
	@genenames=split/\,/,$hashenh{$enh};
	@genenames=uniq(@genenames);
	print AA $enh."\t".join(";",@genenames)."\n";
	}
	close AA;
	system("bedtools sort -i ".$cell.".bed>SEEPgenename2/".$spes[$i]."/".$cell.".bed");
	unlink($cell.".bed");
}
}


}
