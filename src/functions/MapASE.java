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
	/**
	public MapASE(InputStream map, InputStream genotypesInput, String g, int p, File o, String f) throws IOException{
		perm=p;
		outdir = o;
		filename=f;
		
		//get gene and cis-snp info
		ParseMap parsemap = new ParseMap();
		parsemap.parseMap(map, g);
		Gene gene = parsemap.getGene();
		geneName = gene.getId();
		List<SNP> snps = parsemap.getSNPs();
		numSNPs = snps.size();
		Map<String, SNP> snpMap = parsemap.getSnpMap();	
		
		//parse genotype and ASE data
		List<String> genosampleIDs = ParseGenotypes.parseGenotypes(genotypesInput,snps, snpMap);
		numSamples = genosampleIDs.size();
		
		//get subset of genotype and ASE data that have matching samples
		GenoMatrix genomat = new GenoMatrix(snps, genosampleIDs);
		genotypes = genomat.getGenotypes();
		snpids = genomat.getSnpids();
		
		
		ExpMatrix expmat = new ExpMatrix(snps.get(0).getGenosamples());
		expressions = expmat.getExpressions();
		sampleids = expmat.getSampleids();
		
		if(genotypes[0].length!=expressions.length){
			System.out.println("Sample subset not working");
			System.exit(1);
		}
	}

**/
	public MapASE(GenoMatrix genotypes2, ExpMatrix expmat, String geneid, int perm2, File outfile, String filename2) {
		perm = perm2;
		outdir = outfile;
		filename = filename2;
		Gene g = new Gene(geneid);
		
	}

	public MapASE(int[][] geno, int[] is, String[] sampleids2, String[] snpids2, String geneid, int perm2, File outfile,
			String filename2) {
		genotypes = geno;
		expressions = is;
		sampleids = sampleids2;
		snpids = snpids2;
		perm = perm2;
		outdir = outfile;
		filename = filename2;
		Gene g = new Gene(geneid);
		numSNPs = snpids.length;
		numSamples = sampleids.length;
	}
	public void mapase() throws IOException {
		pmap = new HashMap<String, double[]>(); //possible p-values for all SNPs
		
		String line = pointwisePValue();
		String[] tokens = line.split("\\s+");
		System.out.println(line);
		int i = tokens.length-1;
		double minPointwise = Double.parseDouble(tokens[i]);
		
		double[] permPValues = new double[perm]; //list of minimum p-value for each permutation
		for(int p=1; p<=perm; p++){
			permPValues[p-1]= permutationPValue();
		}


		double genePValue = calcPValues(minPointwise, permPValues);
		
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+filename)));
		outfile.write(line+"\t"+genePValue+"\n");
		outfile.close();

	}



	//Calculates pointwise p-values and populates $realPValues and $lines
	public String pointwisePValue(){
		String line="";
		double min=1;
		for(int i=0; i<numSNPs; i++){
			//number of ones in ase samples
			int m=0;
			//number of ones in genotype samples
			int k=0;
			//number of mismatches
			int incorrect=0;
			
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

			}

			//calculate p-value based on hypergeometric distribution
			ComputeSig sig = new ComputeSig(numSamples, m, k, incorrect);
			double[] poss = sig.significance();
			int j = (incorrect - Math.abs(m-k))/2;
			double p = poss[j];
			if(p<min){
				min = p;
				line = geneName+"\t"+snpids[i]+"\t"+m+"\t"+k+"\t"+numSamples+"\t"+incorrect+"\t"+p;
			}
			
			pmap.put(snpids[i], poss);
			
		}
		
		return line;

	}
	
	//returns minimum pointwise p-value from permuted data
	public double permutationPValue(){
		double min=1;
		for(int i=0; i<numSNPs; i++){
			shuffleArray(expressions);
		
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
	



}
