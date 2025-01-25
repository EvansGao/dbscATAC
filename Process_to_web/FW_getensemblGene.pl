#@spes=("hs","mm","at","dr");
#@spechrs=("chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chr20:chr21:chr22:chrX:chrY","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chrX:chrY","chr1:chr2:chr3:chr4:chr5:chrMt:chrPt","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chr20:chr21:chr22:chr23:chr24:chr25:chrM");


@spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
#@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
@spechrs=("chr1:chr2:chr3:chr4:chr5:chr6:chr7:chrX:chr10:chr11:chr8:chr9:chr12:chr13:chr14:chr15:chr16:chr17:chr21:chr19:chr22:chr18:chr20:chrY","chr1:chr2:chrX:chr3:chr4:chr5:chr6:chr7:chr10:chr8:chr14:chr9:chr11:chr13:chr12:chr15:chr16:chr17:chrY:chr18:chr19","chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chr1:chr20:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chrX","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chr20:chr21:chr22:chr23:chr24:chr25:chr26:chr27:chr28:chr29:chrMT:chrX:chrY","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chrX:chr8:chr9:chr11:chr10:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr20:chr19:chrY:chr22:chr21:chrM","chr1:chr2:chr3:chr4:chrZ:chr5:chr7:chr6:chr8:chr9:chr10:chr12:chr11:chr13:chr14:chr20:chr15:chr18:chr17:chr19:chr27:chr33:chr21:chrW:chr24:chr31:chr23:chr26:chr22:chr28:chr25:chr16:chr30:chr32:chrM","chr4:chr7:chr5:chr3:chr6:chr2:chr1:chr9:chr16:chr20:chr8:chr17:chr14:chr13:chr18:chr12:chr19:chr15:chr23:chr21:chr10:chr11:chr24:chr22:chr25:chrM","chr1:chr3:chr4:chr5:chr6:chr7:chrX:chr8:chr12:chr10:chr11:chr2B:chr9:chr2A:chr13:chr14:chr15:chr17:chr16:chr18:chr20:chr19:chr22:chr21:chrY","chr1:chr2:chr5:chr3:chr6:chr4:chr7:chrX:chr8:chr9:chr11:chr12:chr14:chr15:chr13:chr10:chr17:chr16:chr20:chr18:chr19:chrY","chr3R:chr3L:chr2R:chrX:chr2L:chrY:chr4:chrM","chr1:chr2:chr3:chr4:chr5:chrMt:chrPt","chr10:chr11:chr12:chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chrMt:chrPt","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chrMt:chrPt");

for($i=0;$i<scalar(@spes);$i++){
@tmpchrs=split/\:/,$spechrs[$i];
%hashchr=();
foreach $tmpchr (@tmpchrs){
$hashchr{$tmpchr}="";
}
open CELL,"./ensemblgenes/ensemblgenesori".$spes[$i].".txt";
open BB,">./ensemblgenes/ensemblgenes".$spes[$i].".txt";
while(<CELL>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if(exists $hashchr{"chr".$tmp[4]}){
print BB $_."\n";
}
}
close CELL;
close BB;
}
