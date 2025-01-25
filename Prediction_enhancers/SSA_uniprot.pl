##open AA,"../uniprot/uniprotkb_organism_id_3702_2023_12_20.tsv";
##open BB,">uniprot_genes_arabidopsis.txt";
##while(<AA>){
##s/\r|\n//g;
##@tmp=split/\t/,$_;
###if($tmp[1] eq "reviewed"){
##if($tmp[5] =~ /Arabidopsis\sthaliana/){
##print BB $tmp[0]."\t".$tmp[2]."\t".$tmp[1]."\t".$tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[6]."\n";
##}
##}
##close AA;
##close BB;
##
##open AA,"../uniprot/uniprotkb_organism_id_7227_2023_10_12.tsv";
##open BB,">uniprot_genes_fly.txt";
##while(<AA>){
##s/\r|\n//g;
##@tmp=split/\t/,$_;
###if($tmp[1] eq "reviewed"){
##if($tmp[5] =~ /Drosophila\smelanogaster/){
##print BB $tmp[0]."\t".$tmp[2]."\t".$tmp[1]."\t".$tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[6]."\n";
##}
##}
##close AA;
##close BB;
##
##open AA,"../uniprot/uniprotkb_organism_id_7955_2023_12_14.tsv";
##open BB,">uniprot_genes_zebrafish.txt";
##while(<AA>){
##s/\r|\n//g;
##@tmp=split/\t/,$_;
###if($tmp[1] eq "reviewed"){
##if($tmp[5] =~ /Danio\srerio/){
##print BB $tmp[0]."\t".$tmp[2]."\t".$tmp[1]."\t".$tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[6]."\n";
##}
##}
##close AA;
##close BB;
##
##open AA,"../uniprot/uniprotkb_organism_id_9606_2023_10_11.tsv";
##open BB,">uniprot_genes_human.txt";
##while(<AA>){
##s/\r|\n//g;
##@tmp=split/\t/,$_;
###if($tmp[1] eq "reviewed"){
##if($tmp[5] =~ /Homo\ssapiens/){
##print BB $tmp[0]."\t".$tmp[2]."\t".$tmp[1]."\t".$tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[6]."\n";
##}
##}
##close AA;
##close BB;
##
##open AA,"../uniprot/uniprotkb_organism_id_10090_2023_10_11.tsv";
##open BB,">uniprot_genes_mouse.txt";
##while(<AA>){
##s/\r|\n//g;
##@tmp=split/\t/,$_;
###if($tmp[1] eq "reviewed"){
##if($tmp[5] =~ /Mus\smusculus/){
##print BB $tmp[0]."\t".$tmp[2]."\t".$tmp[1]."\t".$tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[6]."\n";
##}
##}
##close AA;
##close BB;
##
##open AA,"../uniprot/uniprotkb_organism_id_9598_2023_10_12.tsv";
##open BB,">uniprot_genes_chimpanzee.txt";
##while(<AA>){
##s/\r|\n//g;
##@tmp=split/\t/,$_;
###if($tmp[1] eq "reviewed"){
##if($tmp[5] =~ /Pan\stroglodytes/){
##print BB $tmp[0]."\t".$tmp[2]."\t".$tmp[1]."\t".$tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[6]."\n";
##}
##}
##close AA;
##close BB;
##
##open AA,"../uniprot/uniprotkb_organism_id_9544_2023_10_12.tsv";
##open BB,">uniprot_genes_rhesusMacaque.txt";
##while(<AA>){
##s/\r|\n//g;
##@tmp=split/\t/,$_;
###if($tmp[1] eq "reviewed"){
##if($tmp[5] =~ /Macaca\smulatta/){
##print BB $tmp[0]."\t".$tmp[2]."\t".$tmp[1]."\t".$tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[6]."\n";
##}
##}
##close AA;
##close BB;
##
##open AA,"../uniprot/uniprotkb_organism_id_9031_2023_11_27.tsv";
##open BB,">uniprot_genes_chicken.txt";
##while(<AA>){
##s/\r|\n//g;
##@tmp=split/\t/,$_;
###if($tmp[1] eq "reviewed"){
##if($tmp[5] =~ /Gallus\sgallus/){
##print BB $tmp[0]."\t".$tmp[2]."\t".$tmp[1]."\t".$tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[6]."\n";
##}
##}
##close AA;
##close BB;
##
##open AA,"../uniprot/uniprotkb_organism_id_9483_2024_04_05.tsv";
##open BB,">uniprot_genes_marmoset.txt";
##while(<AA>){
##s/\r|\n//g;
##@tmp=split/\t/,$_;
###if($tmp[1] eq "reviewed"){
##if($tmp[5] =~ /Callithrix\sjacchus/){
##print BB $tmp[0]."\t".$tmp[2]."\t".$tmp[1]."\t".$tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[6]."\n";
##}
##}
##close AA;
##close BB;

##open AA,"../uniprot/uniprotkb_organism_id_39947_2024_04_09.tsv";
##open BB,">uniprot_genes_rice.txt";
##while(<AA>){
##s/\r|\n//g;
##@tmp=split/\t/,$_;
###if($tmp[1] eq "reviewed"){
##if($tmp[5] =~ /Oryza\ssativa/){
##print BB $tmp[0]."\t".$tmp[2]."\t".$tmp[1]."\t".$tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[6]."\n";
##}
##}
##close AA;
##close BB;

#open AA,"../uniprot/uniprotkb_organism_id_4577_2024_06_21.tsv";
#open BB,">uniprot_genes_maize.txt";
#while(<AA>){
#s/\r|\n//g;
#@tmp=split/\t/,$_;
##if($tmp[1] eq "reviewed"){
#if($tmp[5] =~ /Zea\smays.*Maize/){
#print BB $tmp[0]."\t".$tmp[2]."\t".$tmp[1]."\t".$tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[6]."\n";
#}
#}
#close AA;
#close BB;


#open AA,"../uniprot/uniprotkb_organism_id_60711_2024_06_27.tsv";
#open BB,">uniprot_genes_vervet.txt";
#while(<AA>){
#s/\r|\n//g;
#@tmp=split/\t/,$_;
##if($tmp[1] eq "reviewed"){
#if($tmp[5] =~ /Chlorocebus\ssabaeus.*Green\smonkey/i){
#print BB $tmp[0]."\t".$tmp[2]."\t".$tmp[1]."\t".$tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[6]."\n";
#}
#}
#close AA;
#close BB;

open AA,"../uniprot/uniprotkb_organism_id_9541_2024_06_28.tsv";
open BB,">uniprot_genes_cynomolgus.txt";
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
#if($tmp[1] eq "reviewed"){
if($tmp[5] =~ /Macaca\sfascicularis.*Crab/i){
print BB $tmp[0]."\t".$tmp[2]."\t".$tmp[1]."\t".$tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[6]."\n";
}
}
close AA;
close BB;



open AA,"uniprot_genes_human.txt";
%hashgene=();
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if($tmp[0] ne ""){
$hashgene{$tmp[0]}="";
}
}
close AA;

open AA,"../uniprot/uniprot_human_genes.txt";
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if(!exists $hashgene{$tmp[0]}){
print $tmp[0]."\n";
}
}
close AA;

#marmoset

