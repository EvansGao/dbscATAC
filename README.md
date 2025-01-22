# dbscATAC: a resource of single-cell super-enhancers/enhancers and gene markers for scATAC-seq data
![Fig1](https://github.com/user-attachments/assets/184c25ed-988c-4671-9fdc-3a6be4d8c92b)
dbscATAC is a specialized single-cell database annotating 213,835 super-enhancers, 347,484 marker genes, 13,470,526 enhancers, and 10,402,346 enhancer-gene interactions from 1,668,076 single cells across 1,028 tissue/cell types in 13 species. . The predicted datasets are available in an interactive website http://www.enhanceratlas.org.

# Dependencies
To conduct identification of super-enhancers, typical enhancers, and gene markers, following softwares are required:
✯ R version > 4.2<br />
✯ Perl >5.16.3<br />
✯ deeptools > 2<br />

# Usage
perl EAGLE.pl -E <Enhancer> -G <Expression> -S <Species><br />
Example for prediction of EG interaction: perl EAGLE.pl -E inputexample/cell_enh.bed -G inputexample/cell_gene.txt -S human<br />
Essential Options<br />
-E: A tab-delineate file indicate the enhancer positions and signals with format as "\<chr\>\t\<start\>\t\<end\>\t\<signal\>"<br />
-G: A tab-delineate file displayed the expression value of ensembl genes with format as "\<ensembl\>\t\<value\>"<br />
-S: Speceis "human" or "mouse". The default is "human"<br />
  
# Softwares
To run EAGLE, following softwares are required:<br />
Perl v5.16.3<br />
Matlab R2017b<br />
