$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );

if(-d "./SEcells"){
system("rm -R SEcells");
}

#$spename="fly";
#$spestr="dm";
#@chrarray=("chr2L","chr2LHet","chr2R","chr2RHet","chr3L","chr3LHet","chr3R","chr3RHet","chr4","chrU","chrUextra","chrX","chrXHet","chrYHet","chrM");

#$spename="mouse";
#$spestr="mm";
#@chrarray=("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY");
#@spenames=("chicken","mouse","arabidopsis","zebrafish");
@spenames=("mouse","human","fly");
@spestrs=("mm","hs","dm");
#@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
@spechrs=("chr1:chr2:chrX:chr3:chr4:chr5:chr6:chr7:chr10:chr8:chr14:chr9:chr11:chr13:chr12:chr15:chr16:chr17:chrY:chr18:chr19","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chrX:chr8:chr9:chr11:chr10:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr20:chr19:chrY:chr22:chr21:chrM","chr3R:chr3L:chr2R:chrX:chr2L:chrY:chr4:chrM");
#$spename="human";
#$spestr="hs";



#@spestrs=("gg","mm","at","dr");

#@spechrs=("chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chr20:chr21:chr22:chr23:chr24:chr25:chr26:chr27:chr28:chr30:chrZ:chr31:chr32:chr33:chrW:chrM","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chrX:chrY","chr1:chr2:chr3:chr4:chr5:chrMt:chrPt","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chr20:chr21:chr22:chr23:chr24:chr25:chrM");
#@chrarray=("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chrX","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr20","chrY","chr19","chr22","chr21");

for($i=0;$i<scalar(@spenames);$i++){
@celllines=();
#@chrarray=("chr2L","chr2LHet","chr2R","chr2RHet","chr3L","chr3LHet","chr3R","chr3RHet","chr4","chrU","chrUextra","chrX","chrXHet","chrYHet","chrM");
#my $dir="/zp1/data/tgao/EnhancerAtlas2/standarddmfalseWeb/".$spestr."/".$spestr."allcells";
my $dir="download/superenhancer/".$spestrs[$i];
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
foreach $filedir (@dir){
push @celllines,$filedir;
}
#print scalar(@celllines);
mkdir("SEcells");
mkdir("SEcells/".$spestrs[$i]);
@chrarray=split/\:/,$spechrs[$i];
#@celllines=("416B","BAT","Bone_marrow","CD4+CD8+","CD19+","CMP","EpiSC","ESC_NPC","Forelimb_E13","GMP","HFSC","Striatum","WAT","Brain","C3H10Thalf","CD43-","Cerebellum","CH12","Cortex","ESC_Bruce4","G1E-ER4","G1E","Heart","Kidney","3T3-L1","BMCE","Large_intestine_epithelial","Limb_E14.5","Liver_E14.5","NIH-3T3","NPC");
#@celllines=("L3_eye_discs","embryo_2-4");
$piece=16;
my @processes=();
my $pronum=0;
for($ii=0;$ii<scalar(@celllines);$ii++){
   $processes[$pronum]=fork();
   if($processes[$pronum]){   
#   print $ii."\t".$celllines[$ii]."\n";
   }else{ 
   webformat($celllines[$ii],$dir);         # child handles 
#print $celllines[$ii]."\t".$dir."\n";
    exit 0; 
  }
  if(($pronum+1)%$piece==0){
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}
  }elsif(($pronum+1)%$piece!=0 && ($pronum+1)==scalar(@celllines)){
  	for($k=$piece*int(($pronum+1)/$piece);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*int(($pronum+1)/$piece);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}
  }
	$pronum++;
}

$duration = time - $start;
print "All are done: $duration s\n";

}

sub  webformat(){
my ($name, $dir)=@_;
%hashpositionToID=();
%hashpositionToName=();

$name=~ s/\.bed//g;
$filepath="AllSUEPs/".$spestrs[$i]."/".$name."_SUEP.txt";
print $filepath."\n";
	if (-e $filepath) {
	print "Yes";
	open PRE,$filepath;
		while(<PRE>){
		chomp($_);
		@temp=split/\:|\-|\_|\$/,$_;
		$position=$temp[0].":".$temp[1]."-".$temp[2];
			if(!exists $hashpositionToID{$position}){
			$hashpositionToID{$position}=$temp[3];
			$hashpositionToName{$position}=$temp[4];
			}else{
			$hashpositionToID{$position}.=",".$temp[3];
			$hashpositionToName{$position}.=",".$temp[4];
			}
		}
	close PRE;
	}

mkdir("SEcells/".$spestrs[$i]."/".$name);
%hashchrTOall=();
#print $dir."/".$name.".bed";
open DD,$dir."/".$name.".bed";
while(<DD>){
s/\r|\n//g;
@tmp=split/\t/,$_;
#print "good";
$tmpID="";
$tmpName="";
	if(exists $hashpositionToID{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
	$tmpID=$hashpositionToID{$tmp[0].":".$tmp[1]."-".$tmp[2]};
	$tmpName=$hashpositionToName{$tmp[0].":".$tmp[1]."-".$tmp[2]};
	}
	if($tmp[0] ne "" && !exists $hashchrTOall{$tmp[0]}){
	$hashchrTOall{$tmp[0]}=$_."\t".$tmp[1]."\t".$tmp[2]."\t".$tmp[3]."\t".$tmpID."\t".$tmpName."\n";
	}elsif($tmp[0] ne "" && exists $hashchrTOall{$tmp[0]}){
	$hashchrTOall{$tmp[0]}.=$_."\t".$tmp[1]."\t".$tmp[2]."\t".$tmp[3]."\t".$tmpID."\t".$tmpName."\n";
	}
}
close DD;


foreach $chromo (@chrarray){
	if(exists $hashchrTOall{$chromo}){
	mkdir("SEcells/".$spestrs[$i]."/".$name."/".$chromo);
	open $chromo,">"."SEcells/".$spestrs[$i]."/".$name."/".$chromo."/combined";
	print $chromo $hashchrTOall{$chromo};
	close $chromo;
	}

}

}

