map=ENSG00000000460_snps.txt
genotypeFile=WHLBLD.snps_1.txt
gene=ENSG00000000460
outdir=ase_out
file=ENSG00000000460_simulation.txt

expressions=ENSG00000000460_ase.txt
perm=1
file2=ENSG00000000460_mapase.txt
java -jar ../ase.jar simulation -a $genotypeFile -m $map -g $gene -o $outdir -f $file -n 50

java -jar ../ase.jar mapase -a $genotypeFile -m $map -g $gene -o $outdir -f $file2 -b $expressions -p $perm
