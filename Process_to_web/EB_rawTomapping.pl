use List::MoreUtils qw(uniq);
#@spes=("human","chicken","zebrafish","arabidopsis");
#@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
#@builds=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
#marmoset (calJac3),mouse (mm10),cynomolgus (macFas6),vervet (ChlSab1_1),human (hg38),chicken (galGal6), zebrafish (danRer10), chimp (panTro5), rhesus (rheMac10), fly (dm6), arabidopsis (TAIR10), rice (IRGSP1), maize (B73v4)
#human (hg38), mouse (mm10), marmoset (calJac3), cynomolgus (macFas6),vervet (ChlSab1_1), chicken (galGal6), zebrafish (danRer10), chimp (panTro5), rhesus (rheMac10), fly (dm6), arabidopsis (TAIR10), rice (IRGSP1), maize (B73v4)
@spes=("marmoset","mouse","cynomolgus","vervet","human","chicken","zebrafish","chimp","rhesus","fly","arabidopsis","rice","maize");
#chimphttp://useast.ensembl.org/biomart/martview/54c798e06daa89b37378af92404d7adf 112 Pan_tro_3.0
#zebrafishhttp://may2015.archive.ensembl.org/biomart/martview/84ea87f76305b6cae6ab8ca3f1c69826 GRCz10/danRer10
#chickenhttps://apr2022.archive.ensembl.org/biomart/martview/58809e781a3afb67a2db9bf606be1aad  GRCg6a/galGal6
#cynomolgushttp://useast.ensembl.org/biomart/martview/54c798e06daa89b37378af92404d7adf 112
#maizehttps://nov2020-plants.ensembl.org/biomart/martview/906720c0d1dfd493210d4f6750047291
#flyhttp://jan2024.archive.ensembl.org/biomart/martview/087e12dede9a9f9f681702b754afcae5
#marmosethttp://may2015.archive.ensembl.org/biomart/martview/84ea87f76305b6cae6ab8ca3f1c69826
#rhesushttps://useast.ensembl.org/biomart/martview/45a1148f4a7082997e275fa16019ece2  112   Mmul_10
#ricehttps://plants.ensembl.org/biomart/martview/8cfcdb4ff2f653e19c58e6ddbb8f97b4   59
#arabidopsishttps://plants.ensembl.org/biomart/martview/8cfcdb4ff2f653e19c58e6ddbb8f97b4 59
#vervethttps://useast.ensembl.org/biomart/martview/d908531fed14c15ccf599d64c3dbd042 112
#@spes=("human","chicken","zebrafish","arabidopsis","mouse");
#@fullnames=("Callithrix_jacchus","Mus_musculus","Macaca_fascicularis","Chlorocebus_sabaeus","Homo_sapiens","Gallus_gallus","Danio_rerio","Pan_troglodytes","Macaca_mulatta","Drosophila_melanogaster","Arabidopsis_thaliana","Oryza_sativa_Japonica","Zea_mays");
for($i=0;$i<scalar(@spes);$i++){
open AA,"webgenes/genes".$spes[$i].".txt";
%hashIDs=();
@IDs=();
%hashtogene=();
%hashtoprotein= ();
%hashtoPosition= ();
%hashtoTranscriptName= ();
%hashtoGeneName= ();
#%hashto
while(<AA>){
s/\r|\n//g;
@temp=split/\t/,$_;
if($spes[7] eq "maize"){
$temp[7]=$temp[1];
}
#print $temp[6]."\n";
if(!($temp[0]=~ /Transcript/i) && $temp[0] ne ""){
#if($temp[0]=~ /^EN|^FB|^AT|^Zm/i){	
	if(!exists $hashIDs{$temp[0]}){
	$hashIDs{$temp[0]}=$temp[0];
	push @IDs,$temp[0];
	}
if(!exists $hashtogene{$temp[0]}){
$hashtogene{$temp[0]}=$temp[1];
}elsif(exists $hashtogene{$temp[0]} && !($hashtogene{$temp[0]}=~ /$temp[1]/)){
$hashtogene{$temp[0]}.=";".$temp[1];
}

if(!exists $hashtoprotein{$temp[0]}){
$hashtoprotein{$temp[0]}=$temp[2];
}elsif(exists $hashtoprotein{$temp[0]} && !($hashtoprotein{$temp[0]}=~ /$temp[2]/)){
$hashtoprotein{$temp[0]}.=";".$temp[2];
}


$hashtoPosition{$temp[0]}="chr".$temp[8].":".$temp[3]."-".$temp[4].":".$temp[5];
$hashtoTranscriptName{$temp[0]}=$temp[7];
$hashtoGeneName{$temp[0]}=$temp[6];

#print $hashtoGeneName[$temp[0]]."\n";
}
}
close AA;
%hashtoPDB= ();
%hashtoEntrezGene= ();
%hashtoHGNC= ();

open BB,"webgenes/pdb".$spes[$i].".txt";
while(<BB>){
s/\r|\n//g;
@temp=split/\t/,$_;
if(!($temp[0]=~ /Transcript/i) && $temp[0] ne ""){
	if(!exists $hashIDs{$temp[0]}){
	$hashIDs{$temp[0]}=$temp[0];
	push @IDs,$temp[0];
	}
	
if(!exists $hashtoHGNC{$temp[0]}){
$hashtoHGNC{$temp[0]}=$temp[2];
}elsif(exists $hashtoHGNC{$temp[0]} && !($hashtoHGNC{$temp[0]}=~ /$temp[2]/)){
$hashtoHGNC{$temp[0]}.=";".$temp[2];
}

if(!exists $hashtoEntrezGene{$temp[0]}){
$hashtoEntrezGene{$temp[0]}=$temp[3];
}elsif(exists $hashtoEntrezGene{$temp[0]} && !($hashtoEntrezGene{$temp[0]}=~ /$temp[3]/)){
$hashtoEntrezGene{$temp[0]}.=";".$temp[3];
}

if(!exists $hashtoPDB{$temp[0]}){
$hashtoPDB{$temp[0]}=$temp[4];
}elsif(exists $hashtoPDB{$temp[0]} && !($hashtoPDB{$temp[0]}=~ /$temp[4]/)){
$hashtoPDB{$temp[0]}.=";".$temp[4];
}

}
}
close BB;

%hashRefSeqmRNA=();
%hashRefSeqncRNA=();
%hashChEMBL=();
open CC,"webgenes/ref".$spes[$i].".txt";
while(<CC>){
s/\r|\n//g;
@temp=split/\t/,$_;
if(!($temp[0]=~ /Transcript/i) && $temp[0] ne ""){
	if(!exists $hashIDs{$temp[0]}){
	$hashIDs{$temp[0]}=$temp[0];
	push @IDs,$temp[0];
	}
	
if(!exists $hashChEMBL{$temp[0]}){
$hashChEMBL{$temp[0]}=$temp[4];
}elsif(exists $hashChEMBL{$temp[0]} && !($hashChEMBL{$temp[0]}=~ /$temp[4]/)){
$hashChEMBL{$temp[0]}.=";".$temp[4];
}

if(!exists $hashRefSeqmRNA{$temp[0]}){
$hashRefSeqmRNA{$temp[0]}=$temp[2];
}elsif(exists $hashRefSeqmRNA{$temp[0]} && !($hashRefSeqmRNA{$temp[0]}=~ /$temp[2]/)){
$hashRefSeqmRNA{$temp[0]}.=";".$temp[2];
}

if(!exists $hashRefSeqncRNA{$temp[0]}){
$hashRefSeqncRNA{$temp[0]}=$temp[3];
}elsif(exists $hashRefSeqncRNA{$temp[0]} && !($hashRefSeqncRNA{$temp[0]}=~ /$temp[3]/)){
$hashRefSeqncRNA{$temp[0]}.=";".$temp[3];
}
}
}
close CC;

%hashRefSeqProteinID= ();
%hashUCSCID=();
%hashWikiGeneID=();
open DD,"webgenes/UCSC".$spes[$i].".txt";
while(<DD>){
s/\r|\n//g;
@temp=split/\t/,$_;
if(!($temp[0]=~ /Transcript/i) && $temp[0] ne ""){
	if(!exists $hashIDs{$temp[0]}){
	$hashIDs{$temp[0]}=$temp[0];
	push @IDs,$temp[0];
	}
	
if(!exists $hashRefSeqProteinID{$temp[0]}){
$hashRefSeqProteinID{$temp[0]}=$temp[2];
}elsif(exists $hashRefSeqProteinID{$temp[0]} && !($hashRefSeqProteinID{$temp[0]}=~ /$temp[2]/)){
$hashRefSeqProteinID{$temp[0]}.=";".$temp[2];
}

if(!exists $UCSCID{$temp[0]}){
$UCSCID{$temp[0]}=$temp[3];
}elsif(exists $UCSCID{$temp[0]} && !($UCSCID{$temp[0]}=~ /$temp[3]/)){
$UCSCID{$temp[0]}.=";".$temp[3];
}

if(!exists $hashWikiGeneID{$temp[0]}){
$hashWikiGeneID{$temp[0]}=$temp[4];
}elsif(exists $hashWikiGeneID{$temp[0]} && !($hashWikiGeneID{$temp[0]}=~ /$temp[4]/)){
$hashWikiGeneID{$temp[0]}.=";".$temp[4];
}
}
}
close DD;

open UNIPROT,"/data/gts/dbscATAC/project/uniprot/uniprot_genes_".$spes[$i].".txt";
%hashuniprotTOnames=();
while(<UNIPROT>){
s/\r|\n//g;
@tmp=split/\t/,$_;
@genenames=split/\s+/,$tmp[4];
for($m=0;$m < scalar(@genenames);$m++){
$genenames[$m]=~ s/ZEAMMB73\_//g;
}
	if($tmp[0] ne ""){
#	$hashuniprotTOnames{$tmp[0]}=join(";",@genenames[1..$#genenames]);
	$hashuniprotTOnames{$tmp[0]}=join(";",@genenames);
	}
}
close UNIPROT;



%hashSwissProtID= ();
%hashTrEMBLAccession=();
%hashSwissProtAccession=();
%hashAlias=();
open EE,"webgenes/uniprot".$spes[$i].".txt";
while(<EE>){
s/\r|\n//g;
@temp=split/\t/,$_;
#print $temp[2];

if(!($temp[0]=~ /Transcript/i) && $temp[0] ne ""){
	if(!exists $hashIDs{$temp[0]}){
	$hashIDs{$temp[0]}=$temp[1];
	push @IDs,$temp[0];
	}

	if(!exists $hashAlias{$temp[0]}){
	$hashAlias{$temp[0]}="";
		if(exists $hashuniprotTOnames{$temp[2]}){
		$hashAlias{$temp[0]}=$hashuniprotTOnames{$temp[2]};
		}
		if(exists $hashuniprotTOnames{$temp[3]}){
		$hashAlias{$temp[0]}.=";".$hashuniprotTOnames{$temp[3]};
		}
		if(exists $hashuniprotTOnames{$temp[4]}){
		$hashAlias{$temp[0]}.=";".$hashuniprotTOnames{$temp[4]};
		}
	}elsif(exists $hashAlias{$temp[0]}){
		if(exists $hashuniprotTOnames{$temp[2]}){
		$hashAlias{$temp[0]}.=";".$hashuniprotTOnames{$temp[2]};
		}
		if(exists $hashuniprotTOnames{$temp[3]}){
		$hashAlias{$temp[0]}.=";".$hashuniprotTOnames{$temp[3]};
		}
		if(exists $hashuniprotTOnames{$temp[4]}){
		$hashAlias{$temp[0]}.=";".$hashuniprotTOnames{$temp[4]};
		}
	}

if(!exists $hashSwissProtID{$temp[0]}){
$hashSwissProtID{$temp[0]}=$temp[2];
}elsif(exists $hashSwissProtID{$temp[0]} && !($hashSwissProtID{$temp[0]}=~ /$temp[2]/)){
$hashSwissProtID{$temp[0]}.=";".$temp[2];
}

if(!exists $hashTrEMBLAccession{$temp[0]}){
$hashTrEMBLAccession{$temp[0]}=$temp[3];
}elsif(exists $hashTrEMBLAccession{$temp[0]} && !($hashTrEMBLAccession{$temp[0]}=~ /$temp[3]/)){
$hashTrEMBLAccession{$temp[0]}.=";".$temp[3];
}

if(!exists $hashSwissProtAccession{$temp[0]}){
$hashSwissProtAccession{$temp[0]}=$temp[4];
}elsif(exists $hashSwissProtAccession{$temp[0]} && !($hashSwissProtAccession{$temp[0]}=~ /$temp[4]/)){
$hashSwissProtAccession{$temp[0]}.=";".$temp[4];
}

}
}
close EE;

@genes=();
%hashgeneinfo=();
foreach $ID (@IDs){
if($ID ne ""){
$hashtogene{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashtoprotein{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashtoPosition{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashtoTranscriptName{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashtoGeneName{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashtoHGNC{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashChEMBL{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashRefSeqmRNA{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashRefSeqncRNA{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashtoEntrezGene{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashtoPDB{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashRefSeqProteinID{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$UCSCID{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashFlyBasenameGene{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashSwissProtID{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashAlias{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashTrEMBLAccession{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashSwissProtAccession{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$gene=$hashtogene{$ID};

	if(!exists $hashgeneinfo{$gene}){
	$hashgeneinfo{$gene}=$ID.";".$hashtoprotein{$ID}.";".$hashtoTranscriptName{$ID}.";".$hashtoGeneName{$ID}.";".$hashtoHGNC{$ID}.";".$hashChEMBL{$ID}.";".$hashRefSeqmRNA{$ID}.";".$hashRefSeqncRNA{$ID}.";".$hashtoEntrezGene{$ID}.";".$hashtoPDB{$ID}.";".$hashRefSeqProteinID{$ID}.";".$UCSCID{$ID}.";".$hashWikiGeneID{$ID}.";".$hashSwissProtID{$ID}.";".$hashAlias{$ID}.";".$hashTrEMBLAccession{$ID}.";".$hashSwissProtAccession{$ID};
	push @genes,$gene;
	}else{
	$hashgeneinfo{$gene}.=";".$ID.";".$hashtoprotein{$ID}.";".$hashtoTranscriptName{$ID}.";".$hashtoGeneName{$ID}.";".$hashtoHGNC{$ID}.";".$hashChEMBL{$ID}.";".$hashRefSeqmRNA{$ID}.";".$hashRefSeqncRNA{$ID}.";".$hashtoEntrezGene{$ID}.";".$hashtoPDB{$ID}.";".$hashRefSeqProteinID{$ID}.";".$UCSCID{$ID}.";".$hashWikiGeneID{$ID}.";".$hashSwissProtID{$ID}.";".$hashAlias{$ID}.";".$hashTrEMBLAccession{$ID}.";".$hashSwissProtAccession{$ID};
	}
#print $hashUniProtID{$ID}."\n";

}
}


open FF,">geneMapdata".$spes[$i].".txt";
foreach $gene (@genes){
@geneinfo=split/\;/,$hashgeneinfo{$gene};
@geneinfo=uniq(@geneinfo);
@geneinfo=grep {$_ ne ""} @geneinfo;
	if($gene ne ""){
	print FF $gene."\t".join(",",@geneinfo)."\n";
	}
}
close FF;

}