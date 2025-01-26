@links=("http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/BoneMarrow_62216.bam","http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/BoneMarrow_62216.bam.bai");

#@links=("http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/BoneMarrow_62016.bam","http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/BoneMarrow_62016.bam.bai");
#@links=("http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/Cerebellum_62216.bam","http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/Cerebellum_62216.bam.bai");
foreach $linkone (@links){
	if($linkone=~ /bams\/(.*)/){
	$tmpname=$1;
	system("wget '".$linkone."' -q -o /dev/null -O ".$tmpname);
#	system("gunzip ".$tmpname);
	}
}