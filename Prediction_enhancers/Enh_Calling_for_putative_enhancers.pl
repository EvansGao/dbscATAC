$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
my $dir="rds";
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
%hashsize=();
#%hashspecies=();
foreach $file (@dir){
#@fileinfo=split/\./,$file;
$size = -s "rds/".$file;
$cicerofile=$file;
$cicerofile=~ s/\.rds$//g;
if($file=~ /\.rds/ && !(-e "R/".$cicerofile."_cicero.txt")){
$tmpcelltype=$file;
$tmpcelltype=~ s/\.rds//g;
$hashsize{$tmpcelltype}=$size;
#push @celltypes,$tmpcelltype;
}
}
@celltypes=keys %hashsize;
@celltypes=sort {$hashsize{$a} <=> $hashsize{$b}} @celltypes;
print join("Good\n",@celltypes);
#print $file."\n";
#@celltypes=@celltypes[0..1];


mkdir("R");
mkdir("hundred");
mkdir("peak");
$piece=2;
my @processes=();
my $pronum=0;

if(-d "fragment"){
mkdir("tmp");
mkdir("SE");
mkdir("SEfig");
mkdir("SEgene");
for($ii=0;$ii<scalar(@celltypes);$ii++){
   $processes[$pronum]=fork();
   if($processes[$pronum]){
   print "experiment: ".$celltypes[$ii]."\n";
   }else{ 
    ALL($celltypes[$ii]); # child handles
    exit 0; 
  }
  if(($pronum+1)%$piece==0){
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}
  }elsif(($pronum+1)%$piece!=0 && ($pronum+1)==scalar(@GSMforSRRctrTRUE)){
  	for($k=$piece*int(($pronum+1)/$piece);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*int(($pronum+1)/$piece);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}
  }
	$pronum++;
}


}else{
for($ii=0;$ii<scalar(@celltypes);$ii++){
   $processes[$pronum]=fork();
   if($processes[$pronum]){
   print "experiment: ".$celltypes[$ii]."\n";
   }else{ 
    CICERO($celltypes[$ii]); # child handles
    exit 0; 
  }
  if(($pronum+1)%$piece==0){
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}
  }elsif(($pronum+1)%$piece!=0 && ($pronum+1)==scalar(@GSMforSRRctrTRUE)){
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



sub   CICERO()
{   
	my ($cell)=@_;
	system("Rscript ./JCicero_all.R -m ".$cell);
}

sub   ALL()
{   
	my ($cluster)=@_;
	system("Rscript ./JCicero_all.R -m ".$cluster);
	system("bedtools sort -i fragment/".$cluster.".bed>tmp/".$cluster."_sorted.bed");
	@cellinfo=split/\./,$cluster;
	$spe=$cellinfo[2];
	system("bedToBam -i tmp/".$cluster."_sorted.bed -g /data/gts/dbscATAC/project/chromsizes/".$spe.".chrom.sizes > tmp/".$cluster.".bam");
##	#samtools view tmp/cluster_10_sorted.bam | head -5
	system("samtools sort -o tmp/".$cluster."_sort.bam tmp/".$cluster.".bam");
	system("samtools index tmp/".$cluster."_sort.bam");
	unlink("tmp/".$cluster.".bam");
	unlink("tmp/".$cluster."_sorted.bed");
	$size=2.91e9;
	$distance=12500;
	if($spe eq "mm10"){
	$size=2.7e9;
	$distance=11250;
	}
	system("macs2 callpeak -t tmp/".$cluster."_sort.bam -f BAM -g ".$size." --nomodel --nolambda -n ".$cluster." -q 0.01");
#	macs2 callpeak -t tmp/GSM6598797_CD_C3.GSE214082.mm10_sort.bam -f BAM -g 2.7e9 --nomodel --nolambda -n GSM6598797_CD_C3.GSE214082.mm10 -q 0.01
	rename $cluster."_peaks.narrowPeak", $cluster.".bed" or warn "Rename _peaks to bed failed: $!\n";
	unlink($cluster."_summits.bed");
	unlink($cluster."_peaks.xls");
	if (-e $cluster.".bed") { 
	system("awk -F '\t' '{print \$1\"\t\"\$4\"\t\"\"\t\"\$2\"\t\"\$3\"\t\"\"\t\"\".\"\"\t\"\"\t\"\$4}' ".$cluster.".bed>".$cluster.".gff");
#	system("rm -rf ".$cluster);	
	mkdir("tmp/".$cluster);
	system("ROSE_main.py -g ".$spe." -i ".$cluster.".gff -r tmp/".$cluster."_sort.bam --custom='/data/gts/dbSpatial/superenhancer/ROSE/annotation/".$spe."_refseq.ucsc' -o tmp/".$cluster." -s ".$distance." -t 2500");
#ROSE_main.py -g hg38 -i A549_ip_1_peaks.narrowPeak.gff -r tmp/A549_ip_1.sorted.bam  --custom='/data/gts/dbSpatial/superenhancer/ROSE/annotation/hg38_refseq.ucsc' -o tmp/A549_ip_1_peaks.narrowPeak  -s 11500  -t 2500
#	unlink($cluster.".bed");
#/data/gts/dbscATAC/superenhancer/annotation/hg38_refseq.ucsc
	wait();
	sleep(60);
	unlink("tmp/".$cluster."_sort.bam");
	unlink("tmp/".$cluster."_sort.bam.bai");
	unlink($cluster.".gff");
	system("mv tmp/".$cluster."/".$cellinfo[0]."_SuperStitched.table.txt SE");
	system("mv tmp/".$cluster."/".$cellinfo[0]."_SuperStitched_REGION_TO_GENE.txt SEgene");
	system("mv tmp/".$cluster."/".$cellinfo[0]."_Plot_points.png SEfig");
	system("rm -rf tmp/".$cluster);
#GSM6598797_CD_C5

#GSM6598797_CD_C5.GSE214082.mm10
#ROSE_main.py -g ".$spe." -i ".$cluster.".gff -r tmp/".$cluster."_sort.bam --custom='/data/gts/dbSpatial/superenhancer/ROSE/annotation/".$spe."_refseq.ucsc' -o tmp/".$cluster." -s 12500
#ROSE_main.py -g mm10 -i GSM6598797_CD_C5.GSE214082.mm10.gff -r tmp/GSM6598797_CD_C5.GSE214082.mm10_sort.bam --custom='/data/gts/dbSpatial/superenhancer/ROSE/annotation/mm10_refseq.ucsc' -o tmp/GSM6598797_CD_C5.GSE214082.mm10 -s 12500
#mv tmp/GSM6598797_CD_C5.GSE214082.mm10_sort.bam ./
#mv tmp/GSM6598797_CD_C5.GSE214082.mm10_sort.bam.bai ./
	#python ROSE_main.py -g hg19 -i SRR029348_hg19.gff -r SRR029348_sorted.bam -o SRR029348 -s 12500
	#python ROSE_main.py -g hg19 -i SRR029348_hg19.gff -r SRR029348_sorted.bam -o SRR029348 -s 12500
	}

}

#foreach $celltype (@celltypes){
#system("Rscript ./JCicero.R -m ".$celltype);
#}
$duration = time - $start;
print "All are done: $duration s\n";



#foreach $celltype (@celltypes){
#system("Rscript ./JCicero.R -m ".$celltype);
#}
