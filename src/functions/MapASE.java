package functions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;

import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;



import genome.Gene;
import genome.SNP;
import parse.ParseExpressions;
import parse.ParseGenotypes;
import parse.ParseMap;
import sample.ExpMatrix;
import sample.ExpSample;
import sample.GenoMatrix;
import sample.GenoSample;

public class MapASE {
	int[][] genotypes;
	int numSNPs;
	String[] snpids;
	
	int[] expressions;
	String geneName;
	
	String[] sampleids;
	int numSamples;
	
	Map<String, double[]> pmap;
	double[] pointwise;
	
	int perm;
	File outdir;
	String filename;
	
	
	public MapASE(InputStream map, InputStream genotypesInput,
			InputStream expressionsInput, String g, int p, File o, String f) throws IOException{
		perm=p;
		outdir = o;
		filename=f;
		
		ParseMap parsemap = new ParseMap();
		parsemap.parseMap(map, g);
		Gene gene = parsemap.getGene();
		geneName = gene.getId();
		List<SNP> snps = parsemap.getSNPs();
		numSNPs = snps.size();
		Map<String, SNP> snpMap = parsemap.getSnpMap();	

		//parse genotype and ASE data
		List<String> genosampleIDs = ParseGenotypes.parseGenotypes(genotypesInput,snps, snpMap);
		//Gene gene = new Gene(g);
		List<String> expsampleIDs = ParseExpressions.parseExpressions(expressionsInput, gene, outdir);
	
		genosampleIDs.retainAll(expsampleIDs);
		numSamples = genosampleIDs.size();
		
		//get subset of genotype and ASE data that have matching samples
		GenoMatrix genomat = new GenoMatrix(snps, genosampleIDs);
		ExpMatrix expmat = new ExpMatrix(gene.getExpsamples(), genosampleIDs);
		genotypes = genomat.getGenotypes();
		expressions = expmat.getExpressions();
		sampleids = expmat.getSampleids();
		snpids = genomat.getSnpids();
		
		if(genotypes[0].length!=expressions.length){
			System.out.println("Sample subset not working");
			System.exit(1);
		}

	}

	
	public void mapase() throws IOException {
		pmap = new HashMap<String, double[]>(); //possible p-values for all SNPs
		
		pointwise = new double[numSNPs]; //array to store pointwise p-values
		String[] line = pointwisePValue(); //populate pointwise array and return array of output for each snp
		
		double[] permPValues = new double[perm]; //list of minimum p-value for each permutation
		for(int p=1; p<=perm; p++){
			permPValues[p-1]= permutationPValue();
		}

		double[] permutationPValue = calcPValues(permPValues);
		
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+filename)));
		for(int i=0; i<numSNPs; i++){
			outfile.write(line[i]+"\t"+permutationPValue[i]+"\n");
		}
		outfile.close();

	}



	//Calculates pointwise p-values and populates $realPValues and $lines
	public String[] pointwisePValue(){
		String[] line=new String[numSNPs];
		for(int i=0; i<numSNPs; i++){
			//number of ones in ase samples
			int m=0;
			//number of ones in genotype samples
			int k=0;
			//number of mismatches
			int incorrect=0;
			
			//number of ASE and heterozygous
			int a = 0;
			//number of not ASE and heterozygous
			int b = 0;
			
			for (int j=0; j<numSamples;j++) {
				int hasASE = expressions[j];
				int isHetero = genotypes[i][j];
				
				if(hasASE==1){
					m=m+1; //increment count of number of individuals with ASE
				}
				if(isHetero==1){
					k=k+1; //increment count of number of individuals that are heterozygous
				}

				if(isHetero != hasASE){
					incorrect++; //increment counter for number of mismatches
				}	

				if(hasASE==1 && isHetero==1){
					a=a+1;
				}
				
				if(hasASE==0 && isHetero==1){
					b=b+1;
				}
				
			}

			//calculate p-value based on hypergeometric distribution
			ComputeSig sig = new ComputeSig(numSamples, m, k, incorrect);
			double[] poss = sig.significance();
			pmap.put(snpids[i], poss); //store possible p-values for snp i
			
			int j = (incorrect - Math.abs(m-k))/2; //index of pointwise p-value in array of possible p-values
			double p = poss[j]; //pointwise p-value for snp i
			pointwise[i] = p; // store pointwise p-value
			
			//proportion of heterozygous individuals w/i ASE subset
			double p1 = a*1.0/m;
			//proportion of heterozygous individuals w/i balanced subset
			double p2 = b*1.0/(numSamples - m);
			
			//calculate statistic
			double p3 = (p1+p2)/2;
			double s = (p1 - p2) / Math.sqrt(2.0*p3*(1-p3)/numSamples);
			
			line[i] = geneName+"\t"+snpids[i]+"\t"+a+"\t"+b+"\t"+p1+"\t"+p2+"\t"+s+"\t"+m+"\t"+k+"\t"+numSamples+"\t"+incorrect+"\t"+p; //create line of output for SNP
		}
		
		return line;

	}
	
	//returns minimum pointwise p-value from permuted data
	public double permutationPValue(){
		double min=1;
		shuffleArray(expressions);
		for(int i=0; i<numSNPs; i++){
		
			//number of ones in ase samples
			int m=0;
			//number of ones in genotype samples
			int k=0;

			int incorrect=0;
			for (int j=0; j<numSamples;j++) {
				
				int hasASE = expressions[j];
				int isHetero = genotypes[i][j];
				if(hasASE==1){
					m=m+1;
				}
				if(isHetero==1){
					k=k+1;
				}

				if(isHetero != hasASE){
					incorrect++;
				}	
			}
			int j = (incorrect - Math.abs(m-k))/2;
			if(pmap.get(snpids[i]).length-1<j){
				System.out.println("permutation p value wrong");
				System.exit(1);
			}
			double p = pmap.get(snpids[i])[j];
			if(p<min){
				min=p;
			}
		}
		return min;
	}
	



	static void shuffleArray(int[] ar)
	  {
	    // If running on Java 6 or older, use `new Random()` on RHS here
	    Random rnd = ThreadLocalRandom.current();
	    for (int i = ar.length - 1; i > 0; i--)
	    {
	      int index = rnd.nextInt(i + 1);
	      // Simple swap
	      int a = ar[index];
	      ar[index] = ar[i];
	      ar[i] = a;
	    }
	  }
/**
//For e-gene style calculation
	private double calcPValues(double minPointwise, double[] permPValues) {
		if(permPValues.length!=perm){
			System.out.println("Not all permutations completed");
			System.exit(1);
		}
		
		int count = 0; //number of permuted p-values that are equal to or lower than $p
		for(int j=0; j<perm; j++){
			if(permPValues[j]<=minPointwise){
				count++;
			}
		}
		
		return(count*1.0/perm);
	}
	

**/
	private double[] calcPValues(double[] permPValues) {
		int perm = permPValues.length;
		double[] permutationPValues = new double[numSNPs];
		for(int i=0; i<numSNPs; i++){
			double point = pointwise[i];
			int sig = 0;
			for(int j=0; j<perm; j++){
				if(permPValues[j]<=point){
					sig++;
				}
			}
			permutationPValues[i] = sig*1.0/perm;
		}
		return permutationPValues;
	}


}
