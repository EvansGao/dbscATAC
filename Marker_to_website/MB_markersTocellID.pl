open AA,">cellsIDMarker.txt";
my $firstdir="download/markers";
opendir(DIR,$firstdir) or "can't open the file";
@firstdir=readdir DIR;
@firstdir=grep{$_ ne "." && $_ ne ".."} @firstdir;
foreach $spe (@firstdir){
$speup=uc($spe);
	$seconddir=$firstdir."/".$spe;
	opendir(DIR,$seconddir) or "can't open the file";
	@seconddir=readdir DIR;
	@seconddir=grep{$_ ne "." && $_ ne ".."} @seconddir;
	$rankNum=1;
	foreach $cellfile (@seconddir){
	$cellname=$cellfile;
	$cellname=~ s/\.bed//g;
	print AA $cellname."\t".$spe."\t".$speup.sprintf("%03d",$rankNum)."\n";
	$rankNum++;
	}
}
close AA;


##my $dir="/data/gts/dbscATAC/DataProcess/download/markers";
###my $dir="E:/dbSpatial/perlweb/markers";
##opendir(DIR,$dir) or "can't open the file";
##@dir=readdir DIR;
##
##@allelements=();
##%hashInfoToGene=();
##%hashInfoToCelltype=();
##%hashInfoToTissue=();
##%hashInfoToSpecies=();
##%hashInfoToL2FC=();
##%hashInfoToPvalue=();
##%hashInfoToPCT1=();
##%hashInfoToPCT2=();
##%hashInfoToMethod=();
##%hashInfoToGSM=();
##%hashInfoToPMID=();
##%hashInfoToSample=();
##
###Top	p.value	FDR	logFC	feats	cluster	ranking
###p_val	avg_log2FC	pct.1	pct.2	p_val_adj	cluster	gene
###1.79070029562646e-178	0.630245068880859	1	1	2.72025281908615e-174	1	MT-ND4
##
##foreach $file (@dir){
##@info=split/\_/,$file;
##$species=$info[0];
##$tissue=$info[1];
##$pmid=$info[2];
##$gsm=$info[3];
##$sample=$file;
##$sample=~ s/^.*GSM\d+\_|_giotto_markers\.txt|_seurat_markers\.txt//g;
##	open MARK,$dir."/".$file;
##	while(<MARK>){
##	s/\r|\n//g;
##	@tmp=split/\t/,$_;
##	$tmp[5]="C".$tmp[5];
##		if($file=~ /Human.*seurat/i && $tmp[0] ne "p_val"){
##		$element=$tmp[6]."\t".$tmp[5]."\t".$gsm;
##			if($tmp[4]=~ /e/){
##			@pvalues=split/e/,$tmp[4];
##			$tmp[4]=sprintf("%.3f",$pvalues[0])."e".$pvalues[1];
##			}elsif(!($tmp[4]=~ /e/) && $tmp[4]>0){
##			$tmp[4]=sprintf("%.3e",$tmp[4]);
##			}
##			if(!exists $hashInfoToGene{$element}){
##			push @allelements,$element;
##			$hashInfoToGene{$element}=$tmp[6];
##			$hashInfoToCelltype{$element}=$sample."_".$tmp[5];
##			$hashInfoToTissue{$element}=$tissue;
##			$hashInfoToSpecies{$element}=$species;
##			$hashInfoToL2FC{$element}=sprintf("%.3f",$tmp[1]);
##			$hashInfoToPvalue{$element}=$tmp[4];
##			$hashInfoToPCT1{$element}=$tmp[2];
##			$hashInfoToPCT2{$element}=$tmp[3];
##			$hashInfoToMethod{$element}="Seurat";
##			$hashInfoToGSM{$element}=$gsm;
##			$hashInfoToPMID{$element}=$pmid;
##			$hashInfoToSample{$element}=$sample;
##			}else{
##			$hashInfoToL2FC{$element}.=";".sprintf("%.3f",$tmp[1]);
##			$hashInfoToPvalue{$element}.=";".$tmp[4];
##			$hashInfoToPCT1{$element}=$tmp[2];
##			$hashInfoToPCT2{$element}=$tmp[3];
##			$hashInfoToMethod{$element}.=";"."Seurat";
##			}
##		
##		}elsif($file=~ /Human.*giotto/i && $tmp[0] ne "Top"){
##		$element=$tmp[4]."\t".$tmp[5]."\t".$gsm;
##			if($tmp[1]=~ /e/){
##			@pvalues=split/e/,$tmp[1];
##			$tmp[1]=sprintf("%.3f",$pvalues[0])."e".$pvalues[1];
##			}elsif(!($tmp[1]=~ /e/) && $tmp[1]>0){
##			$tmp[1]=sprintf("%.3e",$tmp[1]);
###			print $tmp[1]."OK";
##			}
##		
##			if(!exists $hashInfoToGene{$element}){
##			push @allelements,$element;
##			$hashInfoToGene{$element}=$tmp[4];
##			$hashInfoToCelltype{$element}=$sample."_".$tmp[5];
##			$hashInfoToTissue{$element}=$tissue;
##			$hashInfoToSpecies{$element}=$species;
##			$hashInfoToL2FC{$element}=sprintf("%.3f",$tmp[3]);
##			$hashInfoToPvalue{$element}=$tmp[1];
##			$hashInfoToPCT1{$element}="NA";
##			$hashInfoToPCT2{$element}="NA";
##			$hashInfoToMethod{$element}="Giotto";
##			$hashInfoToGSM{$element}=$gsm;
##			$hashInfoToPMID{$element}=$pmid;
##			$hashInfoToSample{$element}=$sample;
##			}else{
##			$hashInfoToL2FC{$element}.=";".sprintf("%.3f",$tmp[3]);
##			$hashInfoToPvalue{$element}.=";".$tmp[1];
###			$hashInfoToPCT1{$element}=$tmp[2];
###			$hashInfoToPCT2{$element}=$tmp[3];
##			$hashInfoToMethod{$element}.=";"."Giotto";
##			}
##		}elsif($file=~ /Mouse.*seurat/i && $tmp[0] ne "p_val"){
##		
##		}elsif($file=~ /Mouse.*giotto/i && $tmp[0] ne "p_val"){
##		
##		}
##	
##	
##	}
##	close MARK;
##
##}
##
##%hashInfoToDBid=();
##open AA,">human_markers.txt";
##open BB,">mouse_markers.txt";
##$DBid=1;
##foreach $element (@allelements){
##@basic=split/\t/,$element;
##print AA "dbSp-".sprintf("%06d",$DBid)."\t".$hashInfoToGene{$element}."\t".$hashInfoToCelltype{$element}."\t".$hashInfoToTissue{$element}."\t".$hashInfoToSpecies{$element}."\t".$hashInfoToL2FC{$element}."\t".$hashInfoToPvalue{$element}."\t".$hashInfoToPCT1{$element}."\t".$hashInfoToPCT2{$element}."\t".$hashInfoToMethod{$element}."\t".$hashInfoToGSM{$element}."\t".$hashInfoToPMID{$element}."\t".$hashInfoToSample{$element}."\n";
##$DBid++;
##}
##close AA;
##close BB;
