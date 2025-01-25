$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );

if(-d "./cells"){
system("rm -R cells");
}

#$spename="fly";
#$spestr="dm";
#@chrarray=("chr2L","chr2LHet","chr2R","chr2RHet","chr3L","chr3LHet","chr3R","chr3RHet","chr4","chrU","chrUextra","chrX","chrXHet","chrYHet","chrM");

#$spename="mouse";
#$spestr="mm";
#@chrarray=("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY");
#@spenames=("chicken","mouse","arabidopsis","zebrafish");
@spenames=("marmoset","mouse","cynomolgus","vervet","human","chicken","zebrafish","chimp","rhesus","fly","arabidopsis","rice","maize");
@spestrs=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
#@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
@spechrs=("chr1:chr2:chr3:chr4:chr5:chr6:chr7:chrX:chr10:chr11:chr8:chr9:chr12:chr13:chr14:chr15:chr16:chr17:chr21:chr19:chr22:chr18:chr20:chrY","chr1:chr2:chrX:chr3:chr4:chr5:chr6:chr7:chr10:chr8:chr14:chr9:chr11:chr13:chr12:chr15:chr16:chr17:chrY:chr18:chr19","chrMT:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chr1:chr20:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chrX","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chr20:chr21:chr22:chr23:chr24:chr25:chr26:chr27:chr28:chr29:chrMT:chrX:chrY","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chrX:chr8:chr9:chr11:chr10:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr20:chr19:chrY:chr22:chr21:chrM","chr1:chr2:chr3:chr4:chrZ:chr5:chr7:chr6:chr8:chr9:chr10:chr12:chr11:chr13:chr14:chr20:chr15:chr18:chr17:chr19:chr27:chr33:chr21:chrW:chr24:chr31:chr23:chr26:chr22:chr28:chr25:chr16:chr30:chr32:chrM","chr4:chr7:chr5:chr3:chr6:chr2:chr1:chr9:chr16:chr20:chr8:chr17:chr14:chr13:chr18:chr12:chr19:chr15:chr23:chr21:chr10:chr11:chr24:chr22:chr25:chrM","chr1:chr3:chr4:chr5:chr6:chr7:chrX:chr8:chr12:chr10:chr11:chr2B:chr9:chr2A:chr13:chr14:chr15:chr17:chr16:chr18:chr20:chr19:chr22:chr21:chrY","chr1:chr2:chr5:chr3:chr6:chr4:chr7:chrX:chr8:chr9:chr11:chr12:chr14:chr15:chr13:chr10:chr17:chr16:chr20:chr18:chr19:chrY","chr3R:chr3L:chr2R:chrX:chr2L:chrY:chr4:chrM","chr1:chr2:chr3:chr4:chr5:chrMt:chrPt","chr10:chr11:chr12:chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chrMt:chrPt","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chrMt:chrPt");
#$spename="human";
#$spestr="hs";



#@spestrs=("gg","mm","at","dr");

#@spechrs=("chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chr20:chr21:chr22:chr23:chr24:chr25:chr26:chr27:chr28:chr30:chrZ:chr31:chr32:chr33:chrW:chrM","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chrX:chrY","chr1:chr2:chr3:chr4:chr5:chrMt:chrPt","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chr20:chr21:chr22:chr23:chr24:chr25:chrM");
#@chrarray=("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chrX","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr20","chrY","chr19","chr22","chr21");

for($i=0;$i<scalar(@spenames);$i++){
@celllines=();
#@chrarray=("chr2L","chr2LHet","chr2R","chr2RHet","chr3L","chr3LHet","chr3R","chr3RHet","chr4","chrU","chrUextra","chrX","chrXHet","chrYHet","chrM");
#my $dir="/zp1/data/tgao/EnhancerAtlas2/standarddmfalseWeb/".$spestr."/".$spestr."allcells";
my $dir="./download/enhancer/".$spestrs[$i];
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
foreach $filedir (@dir){
push @celllines,$filedir;
}
#print scalar(@celllines);
mkdir("cells");
mkdir("cells/".$spestrs[$i]);
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
$filepath="./AllEPs/".$spestrs[$i]."/".$name."_EP.txt";
	if (-e $filepath) {
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

mkdir("cells/".$spestrs[$i]."/".$name);
%hashchrTOall=();
print $dir."/".$name.".bed";
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
	mkdir("cells/".$spestrs[$i]."/".$name."/".$chromo);
	open $chromo,">"."cells/".$spestrs[$i]."/".$name."/".$chromo."/combined";
	print $chromo $hashchrTOall{$chromo};
	close $chromo;
	}

}

}

