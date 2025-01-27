$start = time;
#@species=("Gallus_gallus","Mus_musculus","Arabidopsis_thaliana","Danio_rerio");
#@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
#@spenames=("marmoset","mouse","cynomolgus","vervet","human","chicken","zebrafish","chimp","rhesus","fly","arabidopsis","rice","maize");
#@spes=("gg","mm","at","dr");
#@tracks=("combined","P300","POLR2A","Histone","TF-binding","DHS","FAIRE","MNase-seq","GRO-seq","CAGE","MPRA","STARR-seq","CHIA-PET");
#@tracks=("Consensus");
#@firstchrs=("chr1","chr1","chr1","chr1");

@species=("marmoset","mouse","cynomolgus","vervet","human","chicken","zebrafish","chimp","rhesus","fly","arabidopsis","rice","maize");
@spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
@tracks=("Consensus");
@firstchrs=("chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr2R","chr1","chr1","chr1");


%hashspetrackTocells=();
%hashspetrackTocellenhs=();
for($i=0;$i<scalar(@species);$i++){
if(-d "SEcells/".$spes[$i]){
	my $dir="SEcells/".$spes[$i];
	opendir(DIR,$dir) or "can't open the file";
	@dir=readdir DIR;
	@dir=grep{$_ ne "." && $_ ne ".."} @dir;
	foreach $file (@dir){
		print $file."\n";
		my $seconddir="SEcells/".$spes[$i]."/".$file."/".$firstchrs[$i];
		opendir(DIR,$seconddir) or "can't open the file";
		@seconddir=readdir DIR;
		@seconddir=grep{$_ ne "." && $_ ne ".."} @seconddir;
		foreach $track (@seconddir){
			if($track eq "combined"){
			$track="Consensus";
			}
			
			if(!exists $hashspetrackTocells{$spes[$i]."\t".$track}){
			$hashspetrackTocells{$spes[$i]."\t".$track}=$file;
			$hashspetrackTocellenhs{$spes[$i]."\t".$track}=&LINES($spes[$i],$file.".bed");
			}else{
			$hashspetrackTocells{$spes[$i]."\t".$track}.=";".$file;
			$hashspetrackTocellenhs{$spes[$i]."\t".$track}.=";".&LINES($spes[$i],$file.".bed");
			}
		}
	}
}
}

open SPETRACK,">SEspeciestrackinfo.txt";
for($i=0;$i<scalar(@species);$i++){
if(-d "SEcells/".$spes[$i]){
	for($j=0;$j<scalar(@tracks);$j++){
		if(exists $hashspetrackTocells{$spes[$i]."\t".$tracks[$j]}){
		print SPETRACK $species[$i]."\t".$tracks[$j]."\t".$hashspetrackTocells{$spes[$i]."\t".$tracks[$j]}."\t".$hashspetrackTocellenhs{$spes[$i]."\t".$tracks[$j]}."\n";
		
		}
	}
}
}
close SPETRACK;

$duration = time - $start;
print "All are done: $duration s\n";

sub  LINES()
{   
    my ($spe,$filename)=@_;
$linenum=0;
open FILE,"download/superenhancer/".$spe."/".$filename;
while (<FILE>) { $linenum++; }
close FILE;
return $linenum;
}
