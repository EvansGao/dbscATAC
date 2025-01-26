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
@tissues=("Kidney","Liver","PreFrontalCortex","Spleen");
@downloads=("Kidney_62016","Liver_62016","PreFrontalCortex_62216","Spleen_62016");

#@tissues=("Cerebellum");
#@downloads=("Cerebellum_62216");
##1170806
##1170584
#@tissues=("Heart","SmallIntestine");
#@downloads=("HeartA_62816","SmallIntestine_62816");
#@tissues=("Heart","SmallIntestine","LargeIntestine");
#@downloads=("HeartA_62816","SmallIntestine_62816","LargeIntestineA_62816;LargeIntestineB_62816");
#wget "http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/HeartA_62816.bam">& out &
for($i=0;$i<scalar(@tissues);$i++){
@samples=split/\;/,$downloads[$i];
$bamstr="";
	if(scalar(@samples)>1){
		for($j=0;$j<scalar(@samples);$j++){
		system("wget 'http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/".$samples[$j].".bam' -q -o /dev/null -O ".$samples[$j].".bam");
		system("wget 'http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/".$samples[$j].".bam.bai' -q -o /dev/null -O ".$samples[$j].".bam.bai");
		$bamstr.=" ".$samples[$j].".bam";
		}
		system("samtools merge ".$tissues[$i].".bam".$bamstr);
		system("samtools index ".$tissues[$i].".bam");
		for($j=0;$j<scalar(@samples);$j++){
		system("rm -R ".$samples[$j].".bam");
		}
	}else{
		system("wget 'http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/".$samples[0].".bam' -q -o /dev/null -O ".$tissues[$i].".bam");
		system("wget 'http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/".$samples[0].".bam.bai' -q -o /dev/null -O ".$tissues[$i].".bam.bai");
	}

#samtools index BoneMarrow.bam
system("sync; echo 3 > /proc/sys/vm/drop_caches");
system("sinto fragments -b ".$tissues[$i].".bam -p 8 -f ".$tissues[$i]."_fragments.bed --barcode_regex \"[^:]\*\"");
#system("sinto fragments -b ".$tissues[$i].".bam -p 4 -f ".$tissues[$i]."_fragments.bed --barcode_regex \"[^:]\*\"");
#sinto fragments -b Cerebellum.bam -p 8 -f Cerebellum_fragments.bed --barcode_regex "[^:]*"
#sinto fragments -b Heart.bam -p 8 -f Heart_fragments.bed --barcode_regex "[^:]*"
system("rm -R ".$tissues[$i].".bam");
system("rm -R ".$tissues[$i].".bam.bai");
#print "sinto fragments -b ".$tissues[$i].".bam -p 8 -f ".$tissues[$i]."_fragments.bed --barcode_regex \"[^:]\*\"";
system("liftOver ".$tissues[$i]."_fragments.bed mm9Tomm10.over.chain ".$tissues[$i]."_fragments.mm10.bed unMapped");
system("sort -k1,1 -k2,2n ".$tissues[$i]."_fragments.mm10.bed >".$tissues[$i]."_fragments.sort.bed");
system("rm -R ".$tissues[$i]."_fragments.bed");
system("bgzip -@ 4 ".$tissues[$i]."_fragments.sort.bed");
#bgzip -@ 8 BoneMarrow_fragments.sort.bed
system("tabix -p bed ".$tissues[$i]."_fragments.sort.bed.gz");
#tabix -p bed BoneMarrow_fragments.sort.bed.gz
system("rm -R ".$tissues[$i]."_fragments.sort.bed");

}

#liftOver Cerebellum_fragments.sort.bed mm9Tomm10.over.chain Cerebellum_fragments.sort.mm10.bed unMapped
#sort -k1,1 -k2,2n Cerebellum_fragments.sort.mm10.bed >Cerebellum_fragments.sort.mm10.sort.bed
#samtools merge BoneMarrow.bam BoneMarrow_62016.bam BoneMarrow_62216.bam
#samtools index BoneMarrow.bam
#