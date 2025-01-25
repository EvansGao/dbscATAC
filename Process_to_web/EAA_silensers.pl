@species=("hs","mm","dm");
foreach $spe (@species){
open AA,$spe.".bed";
open REDUN,">silencerredun".$spe.".txt";
#open BB,">silencer".$spe.".txt";
#open CC,">silencer".$spe."one.txt";
#open DD,">silencer".$spe."more.txt";
#@tissues=();
%hashtissueTOsilencer=();
@silencers=();
%hashsilencerTOnum=();
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if($tmp[9]=~ /throughput/ && $tmp[0]=~ /^chr\d+$|^chrX$|^chrY$|^chr3r$/){
print REDUN $_."\n";
#push @silencers,$_;
#	if(!exists ){}
if(!exists $hashsilencerTOnum{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}){
$hashsilencerTOnum{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}=1;
push @silencers,$tmp[0]."\t".$tmp[1]."\t".$tmp[2];
}else{
$hashsilencerTOnum{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}+=1;
}
}
}
close AA;
close REDUN;
open PRE,">presilencer".$spe.".bed";
foreach $silencer (@silencers){
	print PRE $silencer."\t".$hashsilencerTOnum{$silencer}."\n";
#	if($hashsilencerTOnum{$silencer}==1){
#	print CC $silencer."\t".$hashsilencerTOnum{$silencer}."\n";
#	}else{
#	print DD $silencer."\t".$hashsilencerTOnum{$silencer}."\n";
#	}
}
#close BB;
#close CC;
#close DD;
close PRE;
system("bedtools sort -i presilencer".$spe.".bed>sortsilencer".$spe.".bed");
#bedtools sort -i presilencerhs.bed>sortsilencerhs.bed
system("bedtools merge -i sortsilencer".$spe.".bed -c 4 -o mean>silencer".$spe.".bed");
#unlink("presilencer".$spe.".bed");
#unlink("sortsilencer".$spe.".bed");
}
#open AA,"Mus_musculus.bed";
#open BB,">silencermm.txt";
#@silencers=();
#while(<AA>){
#s/\r|\n//g;
#@tmp=split/\t/,$_;
#if($tmp[9]=~ /throughput/){
#print BB $_."\n";
#push @silencers,$_;
#}
#}
#close AA;
#close BB;