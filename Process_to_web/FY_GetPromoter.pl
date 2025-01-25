use List::Util qw[min max];

@spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
@upstreams=("5000","5000","5000","4000","5000","2000","2800","5000","5000","250","250","800","3000");
@downstreams=("500","500","500","400","500","200","280","500","500","25","25","80","300");
#@spes=("mm","at","dr","gg","hs","mam","pt");
#@spevers=("mm10","TAIR10","danRer10","galGal6","hg38","rheMac10","panTro5");
#@upstreams=("5000","250","2800","2000","5000","5000","5000");
#@downstreams=("500","25","280","200","500","500","500");
for($i=0;$i<scalar(@spes);$i++){
open CHROM,"../project/chromsizes/".$spevers[$i].".chrom.sizes";
%hashchr=();
while(<CHROM>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if(length($tmp[0])<8){
$hashchr{$tmp[0]}="";
}
}
close CHROM;

if(-e "../project/Refgene/Refgene_".$spevers[$i].".txt"){
open AA,"../project/Refgene/Refgene_".$spevers[$i].".txt";
%hashpro=();
open PRO,">promoterpre".$spes[$i].".bed";
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
$tmp[9]=~ s/\,$//g;
$tmp[10]=~ s/\,$//g;
$tmp[12]=~ s/^\s+|\s+$//g;
if(exists $hashchr{$tmp[2]} && $tmp[3] eq "+" && ($tmp[4]-$upstreams[$i])>0 && $tmp[12] ne ""){
		if(!exists $hashpro{$tmp[2]."\t".($tmp[4]-$upstreams[$i])."\t".($tmp[4]+$downstreams[$i])."\t".$tmp[12]}){
		$hashpro{$tmp[2]."\t".($tmp[4]-$upstreams[$i])."\t".($tmp[4]+$downstreams[$i])."\t".$tmp[12]}="";
		print PRO $tmp[2]."\t".($tmp[4]-$upstreams[$i])."\t".($tmp[4]+$downstreams[$i])."\t".$tmp[12]."\n";
		}
}elsif(exists $hashchr{$tmp[2]} && $tmp[3] eq "-" && ($tmp[5]-$downstreams[$i])>0 && $tmp[12] ne ""){
		if(!exists $hashpro{$tmp[2]."\t".($tmp[5]-$downstreams[$i])."\t".($tmp[5]+$upstreams[$i])."\t".$tmp[12]}){
		$hashpro{$tmp[2]."\t".($tmp[5]-$downstreams[$i])."\t".($tmp[5]+$upstreams[$i])."\t".$tmp[12]}="";
		print PRO $tmp[2]."\t".($tmp[5]-$downstreams[$i])."\t".($tmp[5]+$upstreams[$i])."\t".$tmp[12]."\n";
		}
}
}
close AA;
close PRO;
}elsif(-e "../project/Refgene/genes_exon_".$spevers[$i].".txt"){
open AA,"../project/Refgene/genes_exon_".$spevers[$i].".txt";
%hashpro=();
open PRO,">promoterpre".$spes[$i].".bed";
while(<AA>){
s/\r|\n//g;
s/\r|\n//g;
@tmp=split/\t/,$_;
$tmp[2]=~ s/^\s+|\s+$//g;
$tmp[3]=~ s/^\s+|\s+$//g;
$tmp[4]=~ s/^\s+|\s+$//g;
$tmp[4]="chr".$tmp[4];
$tmp[5]=~ s/^\s+|\s+$//g;
$tmp[6]=~ s/^\s+|\s+$//g;
$tmp[7]=~ s/^\s+|\s+$//g;
if(exists $hashchr{$tmp[4]} && $tmp[8] eq "1" && ($tmp[5]-$upstreams[$i])>0  && $tmp[7] ne ""){
		if(!exists $hashpro{$tmp[4]."\t".($tmp[5]-$upstreams[$i])."\t".($tmp[5]+$downstreams[$i])."\t".$standard}){
		$hashpro{$tmp[4]."\t".($tmp[5]-$upstreams[$i])."\t".($tmp[5]+$downstreams[$i])."\t".$tmp[7]}="";
		print PRO $tmp[4]."\t".($tmp[5]-$upstreams[$i])."\t".($tmp[5]+$downstreams[$i])."\t".$tmp[7]."\n";
		}
}elsif(exists $hashchr{$tmp[4]} && $tmp[8] eq "-1" && ($tmp[6]-$downstreams[$i])>0 && $tmp[7] ne ""){
		if(!exists $hashpro{$tmp[4]."\t".($tmp[6]-$downstreams[$i])."\t".($tmp[6]+$upstreams[$i])."\t".$standard}){
		$hashpro{$tmp[4]."\t".($tmp[6]-$downstreams[$i])."\t".($tmp[6]+$upstreams[$i])."\t".$tmp[7]}="";
		print PRO $tmp[4]."\t".($tmp[6]-$downstreams[$i])."\t".($tmp[6]+$upstreams[$i])."\t".$tmp[7]."\n";
		}
}
}
close AA;
close PRO;
}

system("bedtools sort -i promoterpre".$spes[$i].".bed>promotersort".$spes[$i].".bed");
open CC,"promotersort".$spes[$i].".bed";
open DD,">promoter".$spes[$i].".bed";
@peaks=();
%hashpeak=();
while(<CC>){
chomp($_);
@temp=split/\t/,$_;
$temp[3]=~ s/\(.*//g;
if($temp[0] ne "" && !exists $hashpeak{$temp[0]."\t".$temp[1]."\t".$temp[2]}){
push @peaks,$temp[0]."\t".$temp[1]."\t".$temp[2];
$hashpeak{$temp[0]."\t".$temp[1]."\t".$temp[2]}=$temp[3];
}elsif(exists $hashpeak{$temp[0]."\t".$temp[1]."\t".$temp[2]} && !($hashpeak{$temp[0]."\t".$temp[1]."\t".$temp[2]}=~ /$temp[3]/)){
$hashpeak{$temp[0]."\t".$temp[1]."\t".$temp[2]}.=",".$temp[3];
}

}
close CC;
foreach $peak (@peaks){
@names=split/\,/,$hashpeak{$peak};
@names=sort { lc($a) cmp lc($b) } @names;
$hashpeak{$peak}=join(",",@names);
print DD $peak."\t".$hashpeak{$peak}."\n";
}
close DD;
unlink("promoterpre".$spes[$i].".bed");
unlink("promotersort".$spes[$i].".bed");

}


