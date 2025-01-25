$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
my $dir="R";
if(-d "j"){
system("rm -R j");
system("rm -R c");
system("rm -R e");
system("rm -R p");
system("rm -R i");
}
mkdir("j");
mkdir("c");
mkdir("e");
mkdir("p");
mkdir("i");
#Macaca fascicularis
@builds=("hg38","mm10","galGal6","danRer10","panTro5","rheMac10","calJac3","dm6","ChlSab1_1","macFas6");
@spes=("hs","mm","gg","dr","pt","rm","cj","dm","cs","mf");
%hashbuildTOspe=();
for($kk=0;$kk<scalar(@builds);$kk++){
$hashbuildTOspe{$builds[$kk]}=$spes[$kk];
}


opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
@celltypes=();
%hashcelltype=();
foreach $file (@dir){
$tmpcelltype=$file;
	if($tmpcelltype=~ /\_cicero\.txt/i){
	$tmpcelltype=~ s/\_cicero\.txt|\_jaccard\.txt//g;
		if(!exists $hashcelltype{$tmpcelltype} && !(-e "e/".$tmpcelltype."_enh.bed")){
		$hashcelltype{$tmpcelltype}="";
		push @celltypes,$tmpcelltype;
		}
	}
}

#@celltypes=("Lung_Alveolar_macrophages_30078704");
$nameA=int(rand(100000));
$nameB=int(rand(100000));
foreach $celltype (@celltypes){
open $nameA,"./R/".$celltype."_jaccard.txt";
@celltypeinfo=split/\./,$celltype;
@tmpRegions=();
%hashtmpRegionSig=();
$lengsum=0;
$sigsum=0;
while(<$nameA>){
s/\r|\n|\"//g;
@tmp=split/\_|\s+/,$_;
if($tmp[0] ne "" && $tmp[0] ne "x"){
push @tmpRegions,$tmp[0]."\t".$tmp[1]."\t".$tmp[2];
$hashtmpRegionSig{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}=$tmp[3];
$lengsum+=$tmp[2]-$tmp[1];
$sigsum+=($tmp[2]-$tmp[1])*$tmp[3];
}
}
close $nameA;

if($lengsum==0){
open $nameA,">".$celltype.".error";
close $nameA;
}else{
$meansig=$sigsum/$lengsum;
open $nameA,">j/".$celltype."_jaccardATACpre.bed";
	foreach $tmpRegion (@tmpRegions){
	print $nameA $tmpRegion."\t".(10*$hashtmpRegionSig{$tmpRegion}/$meansig)."\n";
	}
close $nameA;
}

system("bedtools sort -i j/".$celltype."_jaccardATACpre.bed>j/".$celltype."_jaccardATAC.bed");
system("bedtools subtract -a j/".$celltype."_jaccardATAC.bed -b ../perl/standard_promotor_exon_".$celltypeinfo[2].".bed -A>j/".$celltype."_jaccardenhpre.bed");
if(exists $hashbuildTOspe{$celltypeinfo[2]} && $hashbuildTOspe{$celltypeinfo[2]} =~ /hs|mm|dm/){
system("bedtools subtract -a j/".$celltype."_jaccardenhpre.bed -b /data/gts/dbscATAC/project/silencer/silencer".$hashbuildTOspe{$celltypeinfo[2]}.".bed -N -f 0.1>j/".$celltype."_jaccardenhfilter.bed");
system("bedtools sort -i j/".$celltype."_jaccardenhfilter.bed>j/".$celltype."_jaccardenh.bed");
unlink("j/".$celltype."_jaccardenhfilter.bed");
}else{
system("bedtools sort -i j/".$celltype."_jaccardenhpre.bed>j/".$celltype."_jaccardenh.bed");
}

unlink("j/".$celltype."_jaccardenhpre.bed");
system("bedtools intersect -a j/".$celltype."_jaccardATAC.bed -b ../perl/standard_promotor_".$celltypeinfo[2].".bed -wa -wb -f 0.1>j/".$celltype."_jaccardprointersect.bed");
open $nameB,"j/".$celltype."_jaccardprointersect.bed";
@pros=();
%hashproTOscore=();
%hashproTOgene=();
while(<$nameB>){
s/\r|\n|\"//g;
@tmp=split/\t/,$_;
if(!exists $hashproTOgene{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}){
$hashproTOgene{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}=$tmp[7]."|";
push @pros,$tmp[0]."\t".$tmp[1]."\t".$tmp[2];
$hashproTOscore{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}=sprintf("%.4f",$tmp[3]);
}elsif(exists $hashproTOgene{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]} && !($hashproTOgene{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}=~ /$tmp[7]/)){
$hashproTOgene{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}.=$tmp[7]."|";
}
}
close $nameB;
open $nameB,">j/".$celltype."_jaccardpropre.bed";
foreach $pro (@pros){
print $nameB $pro."\t".$hashproTOscore{$pro}."|".$hashproTOgene{$pro}."\n";
}
close $nameB;
system("bedtools sort -i j/".$celltype."_jaccardpropre.bed>j/".$celltype."_jaccardpro.bed");

open $nameA,"./R/".$celltype."_cicero.txt";
print $celltype."\n";
@pairs=();
%hashpairTOscore=();
open $nameB,">c/".$celltype."_ciceroATACpre.bed";
%hashcATAC=();
while(<$nameA>){
s/\r|\n|\"//g;
@tmp=split/\_|\,/,$_;
if($tmp[0]){
push @pairs,$tmp[1]."\t".$tmp[2]."\t".$tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[6];
$hashpairTOscore{$tmp[1]."\t".$tmp[2]."\t".$tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}=$tmp[7];
}
if($tmp[0] ne "" && !exists $hashcATAC{$tmp[1]."\t".$tmp[2]."\t".$tmp[3]}){
$hashcATAC{$tmp[1]."\t".$tmp[2]."\t".$tmp[3]}="";
print $nameB $tmp[1]."\t".$tmp[2]."\t".$tmp[3]."\t1"."\n";
}

if($tmp[0] ne "" && !exists $hashcATAC{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}){
$hashcATAC{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}="";
print $nameB $tmp[4]."\t".$tmp[5]."\t".$tmp[6]."\t1"."\n";
}
}
close $nameA;
close $nameB;
system("bedtools sort -i c/".$celltype."_ciceroATACpre.bed>c/".$celltype."_ciceroATAC.bed");
system("bedtools intersect -a j/".$celltype."_jaccardenh.bed"." -b c/".$celltype."_ciceroATAC.bed"." -u -wa >e/".$celltype."_enh.bed");
system("bedtools intersect -a j/".$celltype."_jaccardpro.bed"." -b c/".$celltype."_ciceroATAC.bed"." -u -wa >p/".$celltype."_propre.bed");
open $nameA,"e/".$celltype."_enh.bed";
%hashenh=();
while(<$nameA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if(!exists $hashenh{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}){
$hashenh{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}=sprintf("%.4f",$tmp[3]);
}

}
close $nameA;
open $nameA,"p/".$celltype."_propre.bed";
%hashpro=();
while(<$nameA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if(!exists $hashpro{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}){
$hashpro{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}=$tmp[3];
}

}
close $nameA;
%hashnorepeat=();
open $nameB,">i/".$celltype."_interaction.txt";
open $nameA,">p/".$celltype."_pro.bed";
%hashtruepro=();
foreach $pair (@pairs){
#	print "OK";
@tmp=split/\t/,$pair;
	if(exists $hashenh{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]} && exists $hashpro{$tmp[3]."\t".$tmp[4]."\t".$tmp[5]}){
		print "Good";
		if(!exists $hashtruepro{$tmp[3]."\t".$tmp[4]."\t".$tmp[5]}){
		$hashtruepro{$tmp[3]."\t".$tmp[4]."\t".$tmp[5]}="";
		print $nameA $tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$hashpro{$tmp[3]."\t".$tmp[4]."\t".$tmp[5]}."\n";
		}
		@allpros=split/\|/,$hashproTOgene{$tmp[3]."\t".$tmp[4]."\t".$tmp[5]};
		@allpros=grep{$_ ne "" } @allpros;
		foreach $pro (@allpros){
		$allstr=$tmp[0].":".$tmp[1]."-".$tmp[2]."|".$hashenh{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}."\t".$tmp[3].":".$tmp[4]."-".$tmp[5]."|".$pro."|"."\t".$hashpairTOscore{$pair};
			if(!exists $hashnorepeat{$allstr}){
			$hashnorepeat{$allstr}="";
			print $nameB $allstr."\n";
			}
		}
	}
	if(exists $hashenh{$tmp[3]."\t".$tmp[4]."\t".$tmp[5]} && exists $hashpro{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}){
		if(!exists $hashtruepro{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}){
		$hashtruepro{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}="";
		print $nameA $tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\t".$hashpro{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}."\n";
		}
		@allpros=split/\|/,$hashproTOgene{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]};
		@allpros=grep{$_ ne "" } @allpros;
		foreach $pro (@allpros){
		$allstr=$tmp[3].":".$tmp[4]."-".$tmp[5]."|".$hashenh{$tmp[3]."\t".$tmp[4]."\t".$tmp[5]}."\t".$tmp[0].":".$tmp[1]."-".$tmp[2]."|".$pro."|"."\t".$hashpairTOscore{$pair};
			if(!exists $hashnorepeat{$allstr}){
			$hashnorepeat{$allstr}="";
			print $nameB $allstr."\n";
			}
		
		}

	}
}
close $nameB;
close $nameA;
unlink("j/".$celltype."_jaccardATACpre.bed");
unlink("j/".$celltype."_jaccardATAC.bed");
unlink("j/".$celltype."_jaccardprointersect.bed");
unlink("j/".$celltype."_jaccardpropre.bed");
unlink("c/".$celltype."_ciceroATACpre.bed");
unlink("p/".$celltype."_propre.bed");
}
$duration = time - $start;
print "All are done: $duration s\n";
#sub   NORM()
#{   
#	my ($GSM,$numtype)=@_;
#	@tmpRegions=();
#	%hashtmpRegionSig=();
#	open $GSM,$hashtmpGSMTotargetfile{$GSM};
#	$lengsum=0;
#	$sigsum=0;
#	while(<$GSM>){
#	s/\r|\n//g;
#	@tmp=split/\t|\s+/,$_;
#		if($tmp[0]=~ /^chr/ && ($tmp[2]-$tmp[1])<=2500 && ($tmp[2]-$tmp[1])>1 && $tmp[$numtype]>0 && $tmp[1]>=0){
#			push @tmpRegions,$tmp[0]."\t".$tmp[1]."\t".$tmp[2];
#			$hashtmpRegionSig{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}=$tmp[$numtype];
#			$lengsum+=$tmp[2]-$tmp[1];
#			$sigsum+=($tmp[2]-$tmp[1])*$tmp[$numtype];
#		}
#	
#	}
#	close $GSM;
#	if($lengsum==0){
#	open $GSM,">".$GSM.".wrongzero";
#	close $GSM;
#	system("mv ".$GSM.".wrongzero Res".$spestr);
#	next;
#	}
#	$meansig=$sigsum/$lengsum;
#	open $GSM,">ENH-".$spestrup."-".$GSM.".bed";
#	foreach $tmpRegion (@tmpRegions){
#	print $GSM $tmpRegion."\t".(10*$hashtmpRegionSig{$tmpRegion}/$meansig)."\n";
#	}
#	close $GSM;
#}