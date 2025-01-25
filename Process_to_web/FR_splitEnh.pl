#@spes=("hs","mm","dm");
#@speups=("HS","MM","DM");

#@spes=("gg","mm","at","dr");
@spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
#@speups=("GG","MM","AT","DR");
@speups=("CJ","MM","MF","CS","HS","GG","DR","PT","RM","DM","AT","OSJ","ZM");
mkdir("enhdetail");
for($i=0;$i<scalar(@spes);$i++){
mkdir("enhdetail/".$spes[$i]);
my $dir="enhs/".$spes[$i];
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
#@dir=grep{$_ ne "." && $_ ne ".."} @dir;
	foreach $file (@dir){
	if($file ne "." && $file ne ".."){
	%hashnum=();
	@arraynum=();
		open AA,$dir."/".$file;
#		print $file;
		mkdir("enhdetail/".$spes[$i]."/".$file);
		while(<AA>){
		chomp($_);
		@temp=split/\t/,$_;
		$temp[0]=~ s/$speups[$i]\d\d-//g;
		$start=int(($temp[0]-1)/100);
		$end=$start+1;
			if(!exists $hashnum{$start}){
			push @arraynum,$start;
			$hashnum{$start}=$start;
			open $start,">enhdetail/".$spes[$i]."/".$file."/".($start*100+1)."-".($start*100+100);
			print $start $_."\n";
			}else{
			print $start $_."\n";
			}
		}
		close AA;
		foreach $num (@arraynum){
		close $num;
		}
	}
}

}

