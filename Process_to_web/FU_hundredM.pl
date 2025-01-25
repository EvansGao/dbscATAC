@tissues=("hs","mm","dm");
foreach $tissue (@tissues){
my $dir="/data/gts/SingleCell/Cicero/download/hundredbed/".$tissue;
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
@celltypes=();
foreach $file (@dir){
if($file=~ /\.txt/){
$tmpcelltype=$file;
$tmpcelltype=~ s/\.txt//g;
push @celltypes,$tmpcelltype;
}
}
@groups=();
for($i=0;$i<scalar(@celltypes);$i++){
	for($j=$i;$j<scalar(@celltypes);$j++){
	push @groups,$celltypes[$i]."\t".$celltypes[$j];
	}
}

open AA,">correlation".$tissue.".txt";
$piece=30;
my @processes=();
my $pronum=0;
for($ii=0;$ii<scalar(@groups);$ii++){
	@cells=split/\t/,$groups[$ii];
   $processes[$pronum]=fork();
   if($processes[$pronum]){   
   print "experiment: ".$groups[$ii]."\n";
   }else{ 
    Similarity($cells[0],$cells[1],$tissue);         # child handles 
    exit 0; 
  }
  if(($pronum+1)%$piece==0){
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}
  }elsif(($pronum+1)%$piece!=0 && ($pronum+1)==scalar(@groups)){
  	for($k=$piece*int(($pronum+1)/$piece);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*int(($pronum+1)/$piece);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}
  }
	$pronum++;
}

close AA;
}


sub  Similarity()
{   
    my ($fileA,$fileB,$tissue)=@_;
    $dirA="/data/gts/SingleCell/Cicero/download/hundredbed/".$tissue."/".$fileA.".txt";
    $dirB="/data/gts/SingleCell/Cicero/download/hundredbed/".$tissue."/".$fileB.".txt";
	system("bedtools intersect -a ".$dirA." -b ".$dirB.">".$fileA.$fileB."tmpoverlap.bed");
	$lengthoverlap=&GetNum($fileA.$fileB."tmpoverlap.bed");
	unlink($fileA.$fileB."tmpoverlap.bed");
	system("cat ".$fileA." ".$fileB.">".$fileA.$fileB."tmpcat.bed");
	system("bedtools sort -i ".$fileA.$fileB."tmpcat.bed>".$fileA.$fileB."tmpcatsort.bed");
	system("bedtools merge -i ".$fileA.$fileB."tmpcatsort.bed>".$fileA.$fileB."tmpunion.bed");
	unlink($fileA.$fileB."tmpcat.bed");
	unlink($fileA.$fileB."tmpcatsort.bed");
	$lengthunion=&GetNum($fileA.$fileB."tmpunion.bed");
	unlink($fileA.$fileB."tmpunion.bed");
	print AA $fileA."\t".$fileB."\t".($lengthoverlap/$lengthunion)."\n";
#	return $lengthoverlap/$lengthunion;
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