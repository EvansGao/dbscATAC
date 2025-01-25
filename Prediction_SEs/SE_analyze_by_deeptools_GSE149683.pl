$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
#sleep(10000);
my $dir="fragmentnew";
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
%hashsize=();
#%hashspecies=();
foreach $file (@dir){
#@fileinfo=split/\./,$file;
$size = -s "fragmentnew/".$file;
@SEfiles=split/\./,$file;
#$cicerofile=~ s/\.rds$//g;
if($file=~ /\.bed\.tar\.gz$/){
#if($file=~ /\.bed/ && !(-e "SE/".$SEfiles[0]."_SuperStitched.table.txt")){
#if($file=~ /\.bed/){
$tmpcelltype=$file;
$tmpcelltype=~ s/\.bed\.tar\.gz//g;
	if(-d "tmp/".$tmpcelltype && !(-e "deeptools/".$tmpcelltype."_HeatmapE.png" && -e "deeptools/".$tmpcelltype."_HeatmapSE.png")){
#	if(-d "tmp/".$tmpcelltype){
	$hashsize{$tmpcelltype}=$size;
	}
#push @celltypes,$tmpcelltype;
}
}
@celltypes=keys %hashsize;
@celltypes=sort {$hashsize{$a} <=> $hashsize{$b}} @celltypes;
#@celltypes=("30078704_celltype_Heart_Cardiomyocytes.GSE111586.mm10");
print join("Good\n",@celltypes);
#print $file."\n";
#@celltypes=@celltypes[0..2];

#30078704_celltype_Cerebellum_Purkinje_cells_SuperStitched.table.txt
#30078704_celltype_Cerebellum_Purkinje_cells.GSE149683.mm10
#mkdir("R");
#mkdir("hundred");
#mkdir("peak");
$piece=5;
my @processes=();
my $pronum=0;

if(-d "fragmentnew"){
print "OK\n";
mkdir("deeptools");
#mkdir("SE");
#mkdir("SEfig");
#mkdir("SEgene");

for($ii=0;$ii<scalar(@celltypes);$ii++){
   $processes[$pronum]=fork();
   if($processes[$pronum]){
   print "experiment: ".$celltypes[$ii]."\n";
   }else{ 
    BIGWIG($celltypes[$ii]); # child handles
    exit 0; 
  }
  if(($pronum+1)%$piece==0){
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}
  }elsif(($pronum+1)%$piece!=0 && ($pronum+1)==scalar(@tissues)){
  	for($k=$piece*int(($pronum+1)/$piece);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*int(($pronum+1)/$piece);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}
  }
	$pronum++;
}


}
#@celltypes=("Acinar_cells_pancreas_33184180","Antigen_presenting_cells_thymus_33184180","Astrocytes_cerebrum_2_33184180");
#@celltypes=@celltypes[20..$#celltypes];
sub   BIGWIG()
{   
	my ($cluster)=@_;
@cellinfo=split/\./,$cluster;
$spe=$cellinfo[2];
$filename=$cluster;
$filename=~ s/\.|\+|\_|\-//g;
$SE=$filename."RES";

system("tar -zxvf fragmentnew/".$cluster.".bed.tar.gz");
open $filename, "fragmentnew/".$cluster.".bed";
open $SE,">deeptools/".$cluster."tmp.bed";
while(<$filename>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if($tmp[0]=~ /chr\d+$|chrX$|chrY$/){
print $SE $tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\t".$tmp[5]."\n";
}
}
close $filename;
close $SE;

#system("awk -F '\t' '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$6}' fragment/".$cluster.".bed>tmp/".$cluster."tmp.bed");
system("bedtools sort -i deeptools/".$cluster."tmp.bed>deeptools/".$cluster."tmpsort.bed");
system("bedtools merge -i deeptools/".$cluster."tmpsort.bed -c 4 -o sum>deeptools/".$cluster.".bedgraph");
system("bedGraphToBigWig deeptools/".$cluster.".bedgraph /data/gts/dbscATAC/project/chromsizes/".$spe.".chrom.sizes deeptools/".$cluster.".bigwig");
unlink("fragment/".$cluster.".bed");
unlink("deeptools/".$cluster."tmp.bed");
unlink("deeptools/".$cluster."tmpsort.bed");
unlink("deeptools/".$cluster.".bedgraph");
#print "computeMatrix scale-regions -S tmp/".$cluster.".bigwig -R e/".$cluster."_enh.bed -a 10000 -b 10000 --regionBodyLength 600 --maxThreshold 66 --missingDataAsZero -out tmp/".$cluster."_SE.tab.gz\n";
#computeMatrix scale-regions -S tmp/30078704_celltype_WholeBrain_Purkinje_cells.GSE111586.mm10.bigwig -R e/30078704_celltype_WholeBrain_Purkinje_cells.GSE111586.mm10_enh.bed -a 10000 -b 10000 --regionBodyLength 600 --maxThreshold 66 --missingDataAsZero -out tmp/30078704_celltype_WholeBrain_Purkinje_cells.GSE111586.mm10_SE.tab.gz
system("computeMatrix scale-regions -S deeptools/".$cluster.".bigwig -R e/".$cluster."_enh.bed -a 10000 -b 10000 --regionBodyLength 1000 --maxThreshold 500 --missingDataAsZero -out deeptools/".$cluster."_E.tab.gz");
system("plotHeatmap -m deeptools/".$cluster."_E.tab.gz -out deeptools/".$cluster."_HeatmapE.png --colorMap summer --startLabel \"start         \" --endLabel \"        end\"  --regionsLabel \"\" --heatmapHeight 12 --heatmapWidth 8 --samplesLabel  \"Typical Enhancer\"");
#plotHeatmap -m tmp/30078704_celltype_WholeBrain_SOM+_Interneurons.GSE111586.mm10_E.tab.gz -out deeptools/30078704_celltype_WholeBrain_SOM+_Interneurons.GSE111586.mm10_HeatmapE.png --startLabel "start         " --endLabel "        end"  --regionsLabel "" --heatmapHeight 10 --samplesLabel  "Typical Enhancer"


open $filename, "tmp/".$cluster."/".$cellinfo[0]."_SuperStitched.table.txt";
open $SE,">deeptools/".$cluster."_SEtmp.bed";
while(<$filename>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if($tmp[1]=~ /chr/){
print $SE $tmp[1]."\t".$tmp[2]."\t".$tmp[3]."\t".$tmp[6]."\t".$tmp[4]."\n";
}
}
close $filename;
close $SE;
system("bedtools sort -i deeptools/".$cluster."_SEtmp.bed>deeptools/".$cluster."_SE.bed");
#bedtools sort -i tmp/30078704_celltype_WholeBrain_Purkinje_cells.GSE111586.mm10_SEtmp.bed>tmp/30078704_celltype_WholeBrain_Purkinje_cells.GSE111586.mm10_SE.bed
system("computeMatrix scale-regions -S deeptools/".$cluster.".bigwig -R deeptools/".$cluster."_SE.bed -a 10000 -b 10000 --binSize 200 --regionBodyLength 18000 --maxThreshold 500 --missingDataAsZero -out deeptools/".$cluster."_SE.tab.gz");

#system("computeMatrix scale-regions -S bw/".$cluster.".bigwig -R SE/".$cluster."_hg19_Super.bed -a 10000 -b 10000 --binSize 200 --regionBodyLength 18000 --missingDataAsZero -out tmp/".$cluster."_SUPER.tab.gz");
#computeMatrix scale-regions -S tmp/30078704_celltype_WholeBrain_Purkinje_cells.GSE111586.mm10.bigwig -R tmp/30078704_celltype_WholeBrain_Purkinje_cells.GSE111586.mm10_SE.bed -a 10000 -b 10000 --binSize 200 --regionBodyLength 18000 --maxThreshold 500 --missingDataAsZero -out tmp/30078704_celltype_WholeBrain_Purkinje_cells.GSE111586.mm10_SE.tab.gz

system("plotHeatmap -m deeptools/".$cluster."_SE.tab.gz -out deeptools/".$cluster."_HeatmapSE.png --colorMap summer --startLabel \"        start\" --endLabel \"end        \"  --regionsLabel \"\" --heatmapHeight 12 --heatmapWidth 8 --samplesLabel  \"Super-Enhancer\"");
unlink("deeptools/".$cluster."_SEtmp.bed");
unlink("deeptools/".$cluster."_SE.bed");
unlink("deeptools/".$cluster."_E.tab.gz");
unlink("deeptools/".$cluster."_SE.tab.gz");
unlink("deeptools/".$cluster.".bigwig");
unlink("deeptools/".$cluster."_enh.bed");
#	system("plotHeatmap -m ".$cluster."_SE.tab.gz -out heatmap/".$cluster."_Heatmap.png --yMin 0 --yMax 12 --startLabel \"start         \" --endLabel \"        end\"  --regionsLabel \"\" --heatmapHeight 15 --samplesLabel  \"Typical Enhancer\"");
#	plotHeatmap -m CD14+_Monocytes_10X_v1_pbmc_10k_SE.tab.gz -out heatmap/CD14+_Monocytes_10X_v1_pbmc_10k_Heatmap.png --yMin 0 --startLabel "start         " --endLabel "        end"  --regionsLabel "" --heatmapHeight 15 --samplesLabel  "Typical Enhancer"
	
#	unlink("tmp/".$cluster."tmp.bed");
#	unlink("tmp/".$cluster."bedgraph");
}



$duration = time - $start;
print "All are done: $duration s\n";



#foreach $celltype (@celltypes){
#system("Rscript ./JCicero.R -m ".$celltype);
#}
