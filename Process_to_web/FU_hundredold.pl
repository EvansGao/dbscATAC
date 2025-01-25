#@tissues=("hs");
@tissues=("hs","mm","dm");
foreach $tissue (@tissues){
my $dir="/data/gts/SingleCell/Cicero/download/hundredbed/".$tissue;
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
@celltypes=();
%hashcelltypeTOnum=();
foreach $file (@dir){
if($file=~ /\.txt/){
$tmpcelltype=$file;
$tmpcelltype=~ s/\.txt//g;
push @celltypes,$tmpcelltype;
$hashcelltypeTOnum{$tmpcelltype}=&GetNum("/data/gts/SingleCell/Cicero/download/hundredbed/".$tissue."/".$tmpcelltype.".txt");
}
}
#@celltypes=("Arteries_C","B_cell_A","Macrophages_A","Capillary_1_A","MyoFibrob_A","Ciliated_cells_A","Pericytes_A","Pericytes_B","T_cell_A","T_cell_B");
open AA,">correlation".$tissue.".txt";
for($i=0;$i<scalar(@celltypes);$i++){
	for($j=$i;$j<scalar(@celltypes);$j++){
	print AA $celltypes[$i]."\t".$celltypes[$j]."\t".&Similarity($celltypes[$i],$celltypes[$j],$tissue)."\n";
	}
}
close AA;
}


sub  Similarity()
{   
    my ($fileA,$fileB,$tissue)=@_;
    $dirA="/data/gts/SingleCell/Cicero/download/hundredbed/".$tissue."/".$fileA.".txt";
    $dirB="/data/gts/SingleCell/Cicero/download/hundredbed/".$tissue."/".$fileB.".txt";
	$numdirA=$hashcelltypeTOnum{$fileA};
	$numdirB=$hashcelltypeTOnum{$fileB};
	system("bedtools intersect -a ".$dirA." -b ".$dirB.">tmpoverlap.bed");
	$numoverlap=&GetNum("tmpoverlap.bed");
#	unlink("tmpoverlap.bed");
#	system("cat ".$dirA." ".$dirB.">tmpcat.bed");
#	system("bedtools sort -i tmpcat.bed>tmpcatsort.bed");
#	system("bedtools merge -i tmpcatsort.bed>tmpunion.bed");
#	unlink("tmpcat.bed");
#	unlink("tmpcatsort.bed");
#	$lengthunion=&GetNum("tmpunion.bed");
#	unlink("tmpunion.bed");
	return $numoverlap/($numdirA+$numdirB-$numoverlap);
}

sub   GetNum()
{   
    my ($file)=@_;
    $num=0;
    open FILE,$file;
	while(<FILE>){
	chomp($_);
	if($_ ne ""){
	$num ++;
	}
	}
    close FILE;
    return $num;
    
    }