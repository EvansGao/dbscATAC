$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
#sleep(10000);
my $dir="fragment";
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
%hashsize=();
#%hashspecies=();
foreach $file (@dir){
#@fileinfo=split/\./,$file;
$size = -s "fragment/".$file;
@SEfiles=split/\./,$file;
#$cicerofile=~ s/\.rds$//g;
if($file=~ /\.bed$/){
#if($file=~ /\.bed/ && !(-e "SE/".$SEfiles[0]."_SuperStitched.table.txt")){
#if($file=~ /\.bed/){
$tmpcelltype=$file;
$tmpcelltype=~ s/\.bed//g;
#0	if(-d "tmp/".$tmpcelltype){
	if(-d "tmp/".$tmpcelltype && !(-e "deeptools/".$tmpcelltype."_HeatmapE.png" && -e "deeptools/".$tmpcelltype."_HeatmapSE.png")){
	$hashsize{$tmpcelltype}=$size;
	}
#push @celltypes,$tmpcelltype;
}
}
@celltypes=keys %hashsize;
@celltypes=sort {$hashsize{$a} <=> $hashsize{$b}} @celltypes;
#@celltypes=("30078704_celltype_WholeBrain_Purkinje_cells.GSE111586.mm10","30078704_celltype_WholeBrain_SOM+_Interneurons.GSE111586.mm10");
print join("Good\n",@tmpcelltype);
#print $file."\n";
#@celltypes=@celltypes[0..2];

#30078704_celltype_Cerebellum_Purkinje_cells_SuperStitched.table.txt
#30078704_celltype_Cerebellum_Purkinje_cells.GSE149683.mm10
#mkdir("R");
#mkdir("hundred");
#mkdir("peak");
$piece=10;
my @processes=();
my $pronum=0;

if(-d "fragment"){
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
#system("tar -zxvf fragment/".$cluster.".bed.tar.gz");
system("awk -F '\t' '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$6}' fragment/".$cluster.".bed>tmp/".$cluster."tmp.bed");
system("bedtools sort -i tmp/".$cluster."tmp.bed>tmp/".$cluster."tmpsort.bed");
system("bedtools merge -i tmp/".$cluster."tmpsort.bed -c 4 -o sum>tmp/".$cluster.".bedgraph");
system("bedGraphToBigWig tmp/".$cluster.".bedgraph /data/gts/dbscATAC/project/chromsizes/".$spe.".chrom.sizes tmp/".$cluster.".bigwig");
#unlink("fragment/".$cluster.".bed");
unlink("tmp/".$cluster."tmp.bed");
unlink("tmp/".$cluster."tmpsort.bed");
unlink("tmp/".$cluster.".bedgraph");
#print "computeMatrix scale-regions -S tmp/".$cluster.".bigwig -R e/".$cluster."_enh.bed -a 10000 -b 10000 --regionBodyLength 600 --maxThreshold 66 --missingDataAsZero -out tmp/".$cluster."_SE.tab.gz\n";
#computeMatrix scale-regions -S tmp/30078704_celltype_WholeBrain_Purkinje_cells.GSE111586.mm10.bigwig -R e/30078704_celltype_WholeBrain_Purkinje_cells.GSE111586.mm10_enh.bed -a 10000 -b 10000 --regionBodyLength 600 --maxThreshold 66 --missingDataAsZero -out tmp/30078704_celltype_WholeBrain_Purkinje_cells.GSE111586.mm10_SE.tab.gz
system("computeMatrix scale-regions -S tmp/".$cluster.".bigwig -R e/".$cluster."_enh.bed -a 10000 -b 10000 --regionBodyLength 600 --maxThreshold 500 --missingDataAsZero -out tmp/".$cluster."_E.tab.gz");
system("plotHeatmap -m tmp/".$cluster."_E.tab.gz -out deeptools/".$cluster."_HeatmapE.png --colorMap summer --startLabel \"start         \" --endLabel \"        end\"  --regionsLabel \"\" --heatmapHeight 12 --heatmapWidth 8 --samplesLabel  \"Typical Enhancer\"");
#plotHeatmap -m tmp/30078704_celltype_WholeBrain_SOM+_Interneurons.GSE111586.mm10_E.tab.gz -out deeptools/30078704_celltype_WholeBrain_SOM+_Interneurons.GSE111586.mm10_HeatmapE.png --startLabel "start         " --endLabel "        end"  --regionsLabel "" --heatmapHeight 10 --samplesLabel  "Typical Enhancer"
$filename=$cluster;
$filename=~ s/\.|\+|\_|\-//g;
$SE=$filename."RES";
open $filename, "tmp/".$cluster."/".$cellinfo[0]."_AllStitched.table.txt";
open $SE,">tmp/".$cluster."_SEtmp.bed";
while(<$filename>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if($tmp[1]=~ /^chr/ && $tmp[8] eq "1"){
print $SE $tmp[1]."\t".$tmp[2]."\t".$tmp[3]."\t".$tmp[6]."\t".$tmp[4]."\n";
}
}
close $filename;
close $SE;
system("bedtools sort -i tmp/".$cluster."_SEtmp.bed>tmp/".$cluster."_SE.bed");
#bedtools sort -i tmp/30078704_celltype_WholeBrain_Purkinje_cells.GSE111586.mm10_SEtmp.bed>tmp/30078704_celltype_WholeBrain_Purkinje_cells.GSE111586.mm10_SE.bed
system("computeMatrix scale-regions -S tmp/".$cluster.".bigwig -R tmp/".$cluster."_SE.bed -a 10000 -b 10000 --binSize 200 --regionBodyLength 18000 --maxThreshold 500 --missingDataAsZero -out tmp/".$cluster."_SE.tab.gz");

#system("computeMatrix scale-regions -S bw/".$cluster.".bigwig -R SE/".$cluster."_hg19_Super.bed -a 10000 -b 10000 --binSize 200 --regionBodyLength 18000 --missingDataAsZero -out tmp/".$cluster."_SUPER.tab.gz");
#computeMatrix scale-regions -S tmp/30078704_celltype_WholeBrain_Purkinje_cells.GSE111586.mm10.bigwig -R tmp/30078704_celltype_WholeBrain_Purkinje_cells.GSE111586.mm10_SE.bed -a 10000 -b 10000 --binSize 200 --regionBodyLength 18000 --maxThreshold 500 --missingDataAsZero -out tmp/30078704_celltype_WholeBrain_Purkinje_cells.GSE111586.mm10_SE.tab.gz

system("plotHeatmap -m tmp/".$cluster."_SE.tab.gz -out deeptools/".$cluster."_HeatmapSE.png --colorMap summer --startLabel \"        start\" --endLabel \"end        \"  --regionsLabel \"\" --heatmapHeight 12 --heatmapWidth 8 --samplesLabel  \"Super-Enhancer\"");
unlink("tmp/".$cluster."_SEtmp.bed");
unlink("tmp/".$cluster."_SE.bed");
unlink("tmp/".$cluster."_E.tab.gz");
unlink("tmp/".$cluster."_SE.tab.gz");
unlink("tmp/".$cluster.".bigwig");
unlink("tmp/".$cluster."_enh.bed");
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
