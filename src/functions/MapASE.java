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
		pmap = new HashMap<String, double[]>();
		double[] pointwisePValues = new double[numSNPs];
		String[] lines = new String[numSNPs];
		double minPointwise = pointwisePValue(pointwisePValues, lines);
		
		
		int notSig=0;
		List<Double> permPValues = new ArrayList<Double>(); //list of minimum p-value for each permutation
		boolean adaptive=false;
		for(int p=1; p<=perm; p++){
			double min = permutationPValue();
			permPValues.add(min);
			
			if(min>=minPointwise){
				notSig++;
			}
			//New adaptive permutations
			BinomialTest bt = new BinomialTest();
			//System.out.println("bt param: "+p+"\t"+notSig);
			if(bt.binomialTest(p, notSig, .0000025, AlternativeHypothesis.GREATER_THAN, .0000025)){
				System.out.println("Permutations: "+ p);
				adaptive=true;
				break;
			}
		}
		
		//nothing is significant
		if(adaptive){
			BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+filename)));
			for(int i=0; i<numSNPs; i++){
				String line = lines[i]+"\t"+1.0+"\t"+0+"\n";
				outfile.write(line);
			}
			outfile.close();
		}
		else{
			System.out.println("Permutations: "+ perm);
			
			List<Double> pValues = calcPValues(pointwisePValues, permPValues); //p-values from permutation test
			
			BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+filename)));
			for(int i=0; i<numSNPs; i++){
				double pv = pValues.get(i);
				String line = lines[i]+"\t"+pv;

				if(isSignificant(pv)){
					line = line+"\t"+1;
				}
				else{
					line = line+"\t"+0;
				}

				outfile.write(line+"\n");
			}
			outfile.close();

		}

		
	}


	private void writePointwise(String[] lines) throws IOException {
		System.out.println("Writing pointwise p-values to file");
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+geneName+"_pointwise.txt")));
		for(int i=0; i<lines.length; i++){
			outfile.write(lines[i]+"\n");
		}
		outfile.close();
	}

	//Calculates pointwise p-values and populates $realPValues and $lines
	public double pointwisePValue(double[] realPValues, String[] lines){
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
			double p = sig.significance();
			realPValues[i]=p;
			if(p<min){
				min = p;
			}
			
			double[] poss = sig.possiblePValues();
			pmap.put(snpids[i], poss);

			String line = geneName+"\t"+snpids[i]+"\t"+m+"\t"+k+"\t"+numSamples+"\t"+incorrect+"\t"+p;
			lines[i] = line;
		}
		
		return min;

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


	private List<Double> calcPValues(double[] pointwisePValues, List<Double> permPValues) {
		if(permPValues.size()!=perm){
			System.out.println("Not all permutations completed");
			System.exit(1);
		}
		
		List<Double> pValues = new ArrayList<Double>();
		
		//iterate over pointwise p-values for SNPs
		for(int i=0; i<pointwisePValues.length; i++){
			double p = pointwisePValues[i]; //pointwise p-value for ith SNP in $snps
			int count = 0; //number of permuted p-values that are equal to or lower than $p
			for(int j=0; j<perm; j++){
				if(permPValues.get(j)<=p){
					count++;
				}
			}
			pValues.add(count*1.0/perm); //add permutation p-value for SNP i to $pValues
		}
		return pValues;
	}
	
	//test whether pvalue passes genome-wide significance threshold
	public boolean isSignificant(double z){
		if(z<=.0000025){
			return true;
		}
		return false;
	}



}
