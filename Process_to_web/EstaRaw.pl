open AA,"SampleAnnotationpre.txt";
open BB,">allsta.txt";
@spes=("hs","mm","dm");
@para=("celltype","Enhancers","scDatasets");
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
}
close AA;