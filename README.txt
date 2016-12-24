This program was compiled using Java 1.8.0_60. 

Usage: ase <function> [<args>] 

Functions:
<simulation>	This function performs a simulation for one gene and a randomly chosen SNP with all possible numbers of errors.
<mapase>	This function maps variants to ASE.

Options:
-m	Map from gene to SNPs
-a	Genotype file
-b	Expression file
-g	Gene name
-p	Number of permutations
-n 	Number of samples
-o	Output directory
-f	Output file name
-h	help statement

Sample usage is included in test_data/runASE.sh.

#################

INPUT 

The map file contains a list of SNPs in the first column.  For the simulations, a SNP is chosen randomly from this list to call ASE.  For the ASE mapping in real data, a p-value based on our algorithm is calculated for all SNPs in this list.  SNPs must have same identifier as in genotype file.

The genotype file contains the phrase "ID" followed by sample IDs separated with whitespace in the first line.  In subsequent lines, the first column contains a SNP id, and the following columns contains a genotype, which is assumed to be a value in [0,2] representing the number of reference alleles.  The genotype values may be floating point numbers if the genotypes are imputed.

The expression file contains the number of reads mapping to the reference and alternative allele for each heterozygous SNP in each individual.  The 7th column contains the subject ID.  The 9th column contains the number of reads mapping to the alternative allele.  The 11th column contains the total number of reads.  All reads must correspond to the gene of interest.

The gene name is a string identifying the gene tested.  This is used for defaut output naming.

#################

OUTPUT

Simulation output file columns:
1. Gene name
2. SNP ID
3. Minor allele-frequency of SNP
4. Number of individuals with ASE
5. Number of individuals heterozygous for SNP
6. Total number of individuals
7. Number of mismatches
8. P-value

ASE mapping algorithm output file:
1. Gene name
2. SNP ID
3. Number of individuals with ASE
4. Number of individuals heterozygous for SNP
5. Total number of individuals
6. Number of mismatches
7. P-value
