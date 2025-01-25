$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );

system("rm -R seq");
mkdir("seq");
@builds=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
@spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
@spechrs=("chr1:chr2:chr3:chr4:chr5:chr6:chr7:chrX:chr10:chr11:chr8:chr9:chr12:chr13:chr14:chr15:chr16:chr17:chr21:chr19:chr22:chr18:chr20:chrY","chr1:chr2:chrX:chr3:chr4:chr5:chr6:chr7:chr10:chr8:chr14:chr9:chr11:chr13:chr12:chr15:chr16:chr17:chrY:chr18:chr19","chrMT:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chr1:chr20:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chrX","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chr20:chr21:chr22:chr23:chr24:chr25:chr26:chr27:chr28:chr29:chrMT:chrX:chrY","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chrX:chr8:chr9:chr11:chr10:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr20:chr19:chrY:chr22:chr21:chrM","chr1:chr2:chr3:chr4:chrZ:chr5:chr7:chr6:chr8:chr9:chr10:chr12:chr11:chr13:chr14:chr20:chr15:chr18:chr17:chr19:chr27:chr33:chr21:chrW:chr24:chr31:chr23:chr26:chr22:chr28:chr25:chr16:chr30:chr32:chrM","chr4:chr7:chr5:chr3:chr6:chr2:chr1:chr9:chr16:chr20:chr8:chr17:chr14:chr13:chr18:chr12:chr19:chr15:chr23:chr21:chr10:chr11:chr24:chr22:chr25:chrM","chr1:chr3:chr4:chr5:chr6:chr7:chrX:chr8:chr12:chr10:chr11:chr2B:chr9:chr2A:chr13:chr14:chr15:chr17:chr16:chr18:chr20:chr19:chr22:chr21:chrY","chr1:chr2:chr5:chr3:chr6:chr4:chr7:chrX:chr8:chr9:chr11:chr12:chr14:chr15:chr13:chr10:chr17:chr16:chr20:chr18:chr19:chrY","chr3R:chr3L:chr2R:chrX:chr2L:chrY:chr4:chrM","chr1:chr2:chr3:chr4:chr5:chrMt:chrPt","chr10:chr11:chr12:chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chrMt:chrPt","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chrMt:chrPt");
#@spenames=("marmoset","mouse","cynomolgus","vervet","human","chicken","zebrafish","chimp","rhesus","fly","arabidopsis","rice","maize");
@linkkeys=("https://ftp.ensembl.org/pub/release-80/fasta/callithrix_jacchus/dna/Callithrix_jacchus.C_jacchus3.2.1.dna.chromosome","https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome","https://ftp.ensembl.org/pub/release-112/fasta/macaca_fascicularis/dna/Macaca_fascicularis.Macaca_fascicularis_6.0.dna.primary_assembly","https://ftp.ensembl.org/pub/release-112/fasta/chlorocebus_sabaeus/dna/Chlorocebus_sabaeus.ChlSab1.1.dna.chromosome","https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome","https://ftp.ensembl.org/pub/release-106/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.chromosome","https://ftp.ensembl.org/pub/release-80/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome","https://ftp.ensembl.org/pub/release-112/fasta/pan_troglodytes/dna/Pan_troglodytes.Pan_tro_3.0.dna.chromosome","https://ftp.ensembl.org/pub/release-112/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.primary_assembly","https://ftp.ensembl.org/pub/release-111/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly","https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome","https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.chromosome","https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-49/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.chromosome");
#https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-49/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.chromosome.1.fa.gz
#https://ftp.ensembl.org/pub/release-102/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.1.fa.gz
#https://ftp.ensembl.org/pub/release-111/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.2L.fa.gz
#https://ftp.ensembl.org/pub/release-112/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.primary_assembly.1.fa.gz
#chimphttp://useast.ensembl.org/biomart/martview/54c798e06daa89b37378af92404d7adf 112 Pan_tro_3.0
#zebrafishhttp://may2015.archive.ensembl.org/biomart/martview/84ea87f76305b6cae6ab8ca3f1c69826 GRCz10/danRer10
#chickenhttps://apr2022.archive.ensembl.org/biomart/martview/58809e781a3afb67a2db9bf606be1aad  GRCg6a/galGal6
#cynomolgushttp://useast.ensembl.org/biomart/martview/54c798e06daa89b37378af92404d7adf 112
#maizehttps://nov2020-plants.ensembl.org/biomart/martview/906720c0d1dfd493210d4f6750047291
#flyhttp://jan2024.archive.ensembl.org/biomart/martview/087e12dede9a9f9f681702b754afcae5
#marmosethttp://may2015.archive.ensembl.org/biomart/martview/84ea87f76305b6cae6ab8ca3f1c69826
#rhesushttps://useast.ensembl.org/biomart/martview/45a1148f4a7082997e275fa16019ece2  112   Mmul_10
#ricehttps://plants.ensembl.org/biomart/martview/8cfcdb4ff2f653e19c58e6ddbb8f97b4   59
#vervethttps://useast.ensembl.org/biomart/martview/d908531fed14c15ccf599d64c3dbd042 112
#@spes=("human","chicken","zebrafish","arabidopsis","mouse");
#@fullnames=("Callithrix_jacchus","Mus_musculus","Macaca_fascicularis","Chlorocebus_sabaeus","Homo_sapiens","Gallus_gallus","Danio_rerio","Pan_troglodytes","Macaca_mulatta","Drosophila_melanogaster","Arabidopsis_thaliana","Oryza_sativa_Japonica","Zea_mays");

#@builds=("TAIR10");
#@spes=("at");
#@spechrs=("chr1:chr2:chr3:chr4:chr5:chrMt:chrPt");
#@linkkeys=("https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome");


for($i=0;$i<scalar(@builds);$i++){
mkdir("seq/".$spes[$i]);
@chrarray=split/\:/,$spechrs[$i];
$piece=8;
my @processes=();
my $pronum=0;
for($ii=0;$ii<scalar(@chrarray);$ii++){
   $processes[$pronum]=fork();
   if($processes[$pronum]){   
#   print $ii."\t".$celllines[$ii]."\n";
   }else{ 
   Seqindex($chrarray[$ii],$spes[$i],$linkkeys[$i]);         # child handles 
    exit 0; 
  }
  if(($pronum+1)%$piece==0){
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}
  }elsif(($pronum+1)%$piece!=0 && ($pronum+1)==scalar(@chrarray)){
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


$duration = time - $start;
print "All are done: $duration s\n";


sub  Seqindex(){
my ($chrom, $spe, $linkkey)=@_;
$chromNum=$chrom;
$chromNum=~s/^chr//g;
#$wholelink=$linkkey.$chrom.".fa.gz";
#$tmpname="";
#if($wholelink=~ /dna\/(.*)$/){
#$tmpname=$1;
#}
#https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
#ftp://ftp.ensembl.org/pub/release-80/fasta/callithrix_jacchus/dna/Callithrix_jacchus.C_jacchus3.2.1.dna.chromosome.1.fa.gz
#print "wget '".$linkkey.".".$chromNum.".fa.gz' -q -o /dev/null -O seq/".$spe."/".$chrom.".fa.gz\n";
if(!(-e "seq/".$spe."/".$chrom.".fa")){
system("wget '".$linkkey.".".$chromNum.".fa.gz' -q -o /dev/null -O seq/".$spe."/".$chrom.".fa.gz");
system("gunzip seq/".$spe."/".$chrom.".fa.gz");
system("sed -i '1d' seq/".$spe."/".$chrom.".fa");
system("sed -i '1"."s\/^\/>".$chrom."\\n\/' seq/".$spe."/".$chrom.".fa");
if(-e "seq/".$spe."/".$chrom.".fa"){
system("samtools faidx seq/".$spe."/".$chrom.".fa");
}
}

}

