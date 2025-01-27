$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
use 5.010;
open IDID,"cellsIDMarker.txt";
%hashcellspeTOid=();
while(<IDID>){
s/\r|\n//g;
@tmp=split/\t/,$_;
$hashcellspeTOid{$tmp[0]."\t".$tmp[1]}=$tmp[2];
}
close IDID;

open TRUEC,"SampleAnnotation.txt";
%hashcellTOtissue=();
%hashcellTOresource=();
%hashspe=();
while(<TRUEC>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if($tmp[2] ne ""){
$hashcellTOtissue{$tmp[2]}=$tmp[4];
$hashcellTOresource{$tmp[2]}=$tmp[7];
}
}
close TRUEC;



@spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
@spefullnames=("Callithrix_jacchus (marmoset)","Mus_musculus (mouse)","Macaca_fascicularis (cynomolgus)","Chlorocebus_sabaeus (vervet)","Homo_sapiens (human)","Gallus_gallus (chicken)","Danio_rerio (zebrafish)","Pan_troglodytes (chimp)","Macaca_mulatta (rhesus)","Drosophila_melanogaster (fly)","Arabidopsis_thaliana (arabidopsis)","Oryza_sativa_Japonica (rice)","Zea_mays (maize)");
@shorts=("marmoset","mouse","cynomolgus","vervet","human","chicken","zebrafish","chimp","rhesus","fly","arabidopsis","rice","maize");
%hashspeTOfullname=();
%hashspeTOver=();
%hashspeTOshort=();
for($i=0;$i<scalar(@spes);$i++){
$hashspeTOfullname{$spes[$i]}=$spefullnames[$i];
$hashspeTOver{$spes[$i]}=$spevers[$i];
$hashspeTOshort{$spes[$i]}=$shorts[$i];
}

my $dir="/data/gts/dbscATAC/DataProcess/download/markers";
#my $dir="E:/dbSpatial/perlweb/markers";
opendir(DIR,$dir) or "can't open the file";
@dirs=readdir DIR;
@dirs=grep{$_ ne "." && $_ ne ".."} @dirs;
foreach $filedir (@dirs){
open GENE,"/data/gts/dbscATAC/project/Refgene/ensembl_".$hashspeTOver{$filedir}.".txt";
%hashgeneTOensembl=();
%hashgeneTOensemblname=();
while(<GENE>){
s/\r|\n//g;
@tmp=split/\t/,$_;

if($tmp[7] ne "" && !exists $hashgeneTOensembl{$tmp[7]}){
$hashgeneTOensembl{$tmp[7]}=$tmp[0].":chr".$tmp[6].":".$tmp[3].":".$tmp[4].":".$tmp[5];
$hashgeneTOensemblname{$tmp[7]}=$tmp[0];
}elsif($tmp[7] eq "" && $tmp[0] ne ""){
$hashgeneTOensembl{$tmp[0]}=$tmp[0].":chr".$tmp[6].":".$tmp[3].":".$tmp[4].":".$tmp[5];
$hashgeneTOensemblname{$tmp[0]}=$tmp[0];
}
}
close GENE;


open SPRCIES,">markers".$filedir.".txt";
my $spedir="/data/gts/dbscATAC/DataProcess/download/markers/".$filedir;
opendir(DIR,$spedir) or "can't open the file";
@cells=readdir DIR;
@cells=grep{$_ ne "." && $_ ne ".."} @cells;
	foreach $file (@cells){
	$num=1;
	$egintfile=$file;
	$egintfile=~ s/(\.txt)/\_interaction$1/g;
	%hashenhinfo=();
		if(-e "/data/gts/dbscATAC/DataProcess/download/interaction/".$filedir."/".$egintfile){
			open INT,"/data/gts/dbscATAC/DataProcess/download/interaction/".$filedir."/".$egintfile;
			while(<INT>){
			s/\r|\n//g;
			@tmp=split/\t|\:|\-|\|/,$_;
			if(!exists $hashenhinfo{$tmp[7]}){
			$hashenhinfo{$tmp[7]}=$tmp[0].":".$tmp[1].":".$tmp[2].":".$tmp[3].":".$tmp[10];
			}elsif(exists $hashenhinfo{$tmp[7]} && index($hashenhinfo{$tmp[7]}, $tmp[0].":".$tmp[1].":".$tmp[2].":".$tmp[3].":".$tmp[10]) == -1){
			$hashenhinfo{$tmp[7]} .= ";".$tmp[0].":".$tmp[1].":".$tmp[2].":".$tmp[3].":".$tmp[10];
			}
			}
			close INT;
		}
		
	open UNIPROT,"/data/gts/dbscATAC/project/uniprot/uniprot_genes_".$hashspeTOshort{$filedir}.".txt";
	%hashaliasTOname=();
	while(<UNIPROT>){
	s/\r|\n//g;
	@tmp=split/\t/,$_;
	@genenames=split/\s+/,$tmp[4];
		for($i=0;$i<scalar(@genenames);$i++){
		$genenames[$i]=~ s/ZEAMMB73\_//g;
		$hashaliasTOname{$genenames[$i]}=$genenames[0];
		}
	}
	close UNIPROT;
	$segintfile=$file;
	$segintfile=~ s/txt$/bed/g;
	%hashseinfo=();
		if(-e "/data/gts/dbscATAC/DataProcess/download/SEgene/".$filedir."/".$segintfile){
			open INT,"/data/gts/dbscATAC/DataProcess/download/SEgene/".$filedir."/".$segintfile;
			while(<INT>){
			s/\r|\n//g;
			@tmp=split/\t|\:|\-|\|/,$_;
			$tmpgenestr=$tmp[3].",".$tmp[4];
			@tmpgenes=split/\,/,$tmpgenestr;
			@tmpgenes=uniq(@tmpgenes);
			foreach $tmpgene (@tmpgenes){
				if($tmpgene eq "Mtss1l"){
				print $tmpgene."\t".$hashaliasTOname{$tmpgene}."\t".$filedir."\n";
				}
				if(exists $hashaliasTOname{$tmpgene}){
				$tmpgene=$hashaliasTOname{$tmpgene};
				}
				
				if(!exists $hashseinfo{$tmpgene}){
				$hashseinfo{$tmpgene}=$tmp[0].":".$tmp[1].":".$tmp[2];
				}elsif(exists $hashseinfo{$tmpgene} && index($hashseinfo{$tmpgene}, $tmp[0].":".$tmp[1].":".$tmp[2]) == -1){
				$hashseinfo{$tmpgene} .= ";".$tmp[0].":".$tmp[1].":".$tmp[2];
				}
			}

			}
			close INT;
		}
	
	
	$cell=$file;
	$cell=~ s/\.txt//g;
		open MARK,$spedir."/".$file;
		while(<MARK>){
		s/\r|\n//g;
		@tmp=split/\t/,$_;
		$ID=$hashcellspeTOid{$file."\t".$filedir}."-M".sprintf("%04d",$num);
		$enhinfo="NA";
		if(exists $hashenhinfo{$tmp[7]}){
		$enhinfo=$hashenhinfo{$tmp[7]};
		}
		$seinfo="NA";
		if(exists $hashseinfo{$tmp[7]}){
		$seinfo=$hashseinfo{$tmp[7]};
		}

		if($tmp[5]=~ /e/){
		@pvalues=split/e/,$tmp[5];
		$tmp[5]=sprintf("%.3f",$pvalues[0])."e".$pvalues[1];
		}elsif(!($tmp[5]=~ /e/) && $tmp[5]>0){
		$tmp[5]=sprintf("%.3e",$tmp[5]);
		}
		$tmpensembl="NA";
		if(exists $hashgeneTOensemblname{$tmp[7]}){
		$tmpensembl=$hashgeneTOensemblname{$tmp[7]};
		}
		print SPRCIES $ID."\t".$tmp[7]."\t".$tmpensembl."\t".$filedir."\t".$tmp[6]."\t".$cell."\t".$hashcellTOtissue{$cell}."\t".$hashspeTOfullname{$filedir}."\t".sprintf("%.3f",$tmp[2])."\t".$tmp[5]."\t".sprintf("%.3f",$tmp[3])."\t".sprintf("%.3f",$tmp[4])."\t"."Signac"."\t".$hashcellTOresource{$cell}."\t".$hashgeneTOensembl{$tmp[7]}."\t".$enhinfo."\t".$seinfo."\n";
		$num++;
#print AA "dbSp-".sprintf("%06d",$DBid)."\t".$hashInfoToGene{$element}."\t".$hashInfoToCelltype{$element}."\t".$hashInfoToTissue{$element}."\t".$hashInfoToSpecies{$element}."\t".$hashInfoToL2FC{$element}."\t".$hashInfoToPvalue{$element}."\t".$hashInfoToPCT1{$element}."\t".$hashInfoToPCT2{$element}."\t".$hashInfoToMethod{$element}."\t".$hashInfoToGSM{$element}."\t".$hashInfoToPMID{$element}."\t".$hashInfoToSample{$element}."\n";

		}
	}
close SPRCIES;
}


$duration = time - $start;
print "All are done: $duration s\n";






