##https://atlas.gs.washington.edu/mouse-atac/data/

#@links=("http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/BoneMarrow_62216.bam","http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/BoneMarrow_62216.bam.bai");

#@links=("http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/BoneMarrow_62016.bam","http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/BoneMarrow_62016.bam.bai");
#@links=("http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/Cerebellum_62216.bam","http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/Cerebellum_62216.bam.bai");

#@tissues=("BoneMarrow","","","","","","","");
#@tissues=("BoneMarrow","Heart","LargeIntestine");
#@downloads=("BoneMarrow_62016;BoneMarrow_62216","HeartA_62816","LargeIntestineA_62816;LargeIntestineB_62816");

#@tissues=("LargeIntestine","Kidney","Liver","Lung","PreFrontalCortex","Spleen","Testes","Thymus","WholeBrain");
#@downloads=("LargeIntestineA_62816;LargeIntestineB_62816","Kidney_62016","Liver_62016","Lung1_62216;Lung2_62216","PreFrontalCortex_62216","Spleen_62016","Testes_62016","Thymus_62016","WholeBrainA_62216;WholeBrainA_62816");

#all tissues
#@tissues=("Spleen","Testes","Thymus","WholeBrain");
#@tissues=("Kidney","Liver","PreFrontalCortex","BoneMarrow","Cerebellum","Heart","LargeIntestine","Lung","SmallIntestine");
@tissues=("Spleen","Testes","Thymus","WholeBrain","Kidney","Liver","PreFrontalCortex","BoneMarrow","Cerebellum","Heart","LargeIntestine","Lung","SmallIntestine");
@downloads=("Spleen_62016","Testes_62016","Thymus_62016","WholeBrainA_62216;WholeBrainA_62816","Kidney_62016","Liver_62016","PreFrontalCortex_62216","BoneMarrow_62016;BoneMarrow_62216","Cerebellum_62216","HeartA_62816","LargeIntestineA_62816;LargeIntestineB_62816","Lung1_62216;Lung2_62216","SmallIntestine_62816");
%hashsample=();
for($i=0;$i<scalar(@tissues);$i++){
$hashsample{$tissues[$i]}=$downloads[$i];
}


$piece=4;
my @processes=();
my $pronum=0;

for($ii=0;$ii<scalar(@tissues);$ii++){
   $processes[$pronum]=fork();
   if($processes[$pronum]){
   print "experiment: ".$tissues[$ii]."\n";
   }else{ 
    ALL($tissues[$ii]); # child handles
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


#@tissues=("Cerebellum");
#@downloads=("Cerebellum_62216");
##1170806
##1170584
#@tissues=("Heart","SmallIntestine");
#@downloads=("HeartA_62816","SmallIntestine_62816");
#@tissues=("Heart","SmallIntestine","LargeIntestine");
#@downloads=("HeartA_62816","SmallIntestine_62816","LargeIntestineA_62816;LargeIntestineB_62816");
#wget "http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/HeartA_62816.bam">& out &
sub   ALL()
{   
	my ($tissue)=@_;
@samples=split/\;/,$hashsample{$tissue};
$bamstr="";
	if(scalar(@samples)>1){
		for($j=0;$j<scalar(@samples);$j++){
		system("wget 'http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/".$samples[$j].".bam' -q -o /dev/null -O ".$samples[$j].".bam");
		system("wget 'http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/".$samples[$j].".bam.bai' -q -o /dev/null -O ".$samples[$j].".bam.bai");
		$bamstr.=" ".$samples[$j].".bam";
		}
		system("samtools merge ".$tissue.".bam".$bamstr);
		system("samtools index ".$tissue.".bam");
		for($j=0;$j<scalar(@samples);$j++){
		system("rm -R ".$samples[$j].".bam");
		}
	}else{
		system("wget 'http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/".$samples[0].".bam' -q -o /dev/null -O ".$tissue.".bam");
		system("wget 'http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/".$samples[0].".bam.bai' -q -o /dev/null -O ".$tissue.".bam.bai");
	}

#samtools index BoneMarrow.bam
system("sync; echo 3 > /proc/sys/vm/drop_caches");
system("sinto fragments -b ".$tissue.".bam -p 8 -f ".$tissue."_fragments.bed --barcode_regex \"[^:]\*\"");
#system("sinto fragments -b ".$tissue.".bam -p 4 -f ".$tissue."_fragments.bed --barcode_regex \"[^:]\*\"");
#sinto fragments -b BoneMarrow.bam -p 8 -f BoneMarrow_fragments.bed --barcode_regex "[^:]*"
#sinto fragments -b Heart.bam -p 8 -f Heart_fragments.bed --barcode_regex "[^:]*"
system("rm -R ".$tissue.".bam");
system("rm -R ".$tissue.".bam.bai");
#print "sinto fragments -b ".$tissue.".bam -p 8 -f ".$tissue."_fragments.bed --barcode_regex \"[^:]\*\"";
system("liftOver ".$tissue."_fragments.bed mm9Tomm10.over.chain ".$tissue."_fragments.mm10.bed unMapped");
system("sort -k1,1 -k2,2n ".$tissue."_fragments.mm10.bed >".$tissue."_fragments.sort.bed");
system("rm -R ".$tissue."_fragments.bed");
system("bgzip -@ 4 ".$tissue."_fragments.sort.bed");
#bgzip -@ 8 BoneMarrow_fragments.sort.bed
system("tabix -p bed ".$tissue."_fragments.sort.bed.gz");
#tabix -p bed BoneMarrow_fragments.sort.bed.gz
system("rm -R ".$tissue."_fragments.sort.bed");
system("rm -R ".$tissue."_fragments.mm10.bed");
}

#liftOver Cerebellum_fragments.sort.bed mm9Tomm10.over.chain Cerebellum_fragments.sort.mm10.bed unMapped
#sort -k1,1 -k2,2n Cerebellum_fragments.sort.mm10.bed >Cerebellum_fragments.sort.mm10.sort.bed
#samtools merge BoneMarrow.bam BoneMarrow_62016.bam BoneMarrow_62216.bam
#samtools index BoneMarrow.bam
#