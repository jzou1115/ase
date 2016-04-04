package functions;

import java.io.BufferedInputStream;
import java.util.Arrays;
import java.util.Collection;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;
import org.apache.commons.math3.util.CombinatoricsUtils;

import parse.ParseExpressions;
import parse.ParseGene;
import parse.ParseGenotypes;
import parse.ParseMap;
import parse.ParseSNP;
import parse.ParseSamples;
import sample.ExpSample;
import sample.GenoSample;
import genome.DiSNP;
import genome.Gene;
import genome.SNP;

public class Combinations {
	//initialized in getCombindations
	int combs;
	int[][] genotypes;
	String[] sampleids;
	String[] combids;
	
	//initialized in mapase
	HashMap<String, double[]> pmap;
	
	public Combinations(InputStream map, String gene, InputStream genotypes, InputStream expression, int p, File out, String f) throws IOException {
		
		ParseMap parsemap = new ParseMap();
		parsemap.parseMap(map, gene);
		Gene g = parsemap.getGene();
		List<SNP> snps = parsemap.getSNPs();
		Map<String, SNP> snpMap = parsemap.getSnpMap();
		
		ParseGenotypes.parseGenotypes(genotypes, snps, snpMap);
		ParseExpressions.parseExpressions(expression, g, out);
		
		getCombinations(g, snps, p, out, f);
		
		//final long startTime = System.currentTimeMillis();
		//List<SNP> combs = ParseGenotypes.parseGenotypes(new BufferedInputStream( new FileInputStream(new File(out+File.separator+g.getId()+"_combinations.txt"))));
		//final long endTime = System.currentTimeMillis();
		//System.out.println("Total time to create snp objects: " + (endTime - startTime) );
		
		mapase(g,p, out, f);
		
	}
	
	public void getCombinations(Gene g, List<SNP> snps, int perm, File outdir, String filename) throws IOException{
		System.out.println("Generating combinations of SNPs");
		System.out.println("Number of snps: "+snps.size());
		
		List<GenoSample> samples = snps.get(0).getGenosamples();
		int total = samples.size();
		sampleids = new String[total];
		for(int i=0; i<total; i++){
			sampleids[i]= samples.get(i).getID();
		}
		
		int[][] temp = new int[snps.size()][total];
		String[] snpids = new String[snps.size()];
		for(int i=0; i<snps.size();i++){
			SNP s = snps.get(i);
			snpids[i]=s.getId();
			List<GenoSample> genos = s.getGenosamples();
			for(int j=0; j<genos.size(); j++){
				temp[i][j] = genos.get(j).getHetero();
			}
		}
		
		final long startTime = System.currentTimeMillis();
	
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+g.getId()+"_combinations.txt")));
		String header = "SNP\t";
		for(String s:sampleids){
			header = header+s+"\t";
		}
		outfile.write(header.trim()+"\n");
	
		combs=(int) CombinatoricsUtils.binomialCoefficientDouble(snps.size(), 2);
		combids = new String[combs];
		genotypes = new int[combs][total];
		int snpNum = 0;
		for(int i=0; i<temp.length; i++){
			for(int j=i+1; j<temp.length; j++){
				int[] snp1 = temp[i];
				int[] snp2 = temp[j];
				String line= snpids[i]+"-"+snpids[j]+"\t";
				combids[snpNum] = snpids[i]+"-"+snpids[j]+"\t";
				combs++;
				
				for(int k=0; k<snp1.length; k++){
					if(snp1[k]+snp2[k]>0){
						line = line+1+"\t";
						genotypes[snpNum][k] = 1;
					}
					else{
						line = line+0+"\t";
						genotypes[snpNum][k] = 0;
					}
				}
				snpNum++;
				if(snpNum>=combs){
					break;
				}
				outfile.write(line.trim()+"\n");
			}
		}
		outfile.close();
		final long endTime = System.currentTimeMillis();
		System.out.println("Total time to create combinations of snps: " + (endTime - startTime) );
		
		
	}
	

	
	public void mapase(Gene g, int perm, File out, String filename) throws IOException{

		pmap = new HashMap<String, double[]>();
		getSubset(g);
		List<ExpSample> ase = g.getExpsamples();
		
		double[] pointwisePValues = new double[combs];
		String[] lines = new String[combs];
		pointwisePValue(g, pointwisePValues, lines);
		
		/**
		System.out.println("Writing pointwise p-values to file");
		BufferedWriter outfile0 = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(out+File.separator+g.getId()+"_pointwise.txt")));
		for(String l:lines){
			if(l!=null){
				outfile0.write(l+"\n");	
			}
		}
		outfile0.close();
		**/
		
		double minPointwise = minPValue(pointwisePValues);
		
		System.out.println("minpointwise: "+minPointwise);
		int notSig =0;
		List<Double> permPValues = new ArrayList<Double>();
		System.out.println("Starting permuations");
		boolean adaptive=false;
		for(int p=1; p<=perm; p++){
			double min = permutationPValue(ase);
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
			BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(out+File.separator+filename)));
			for(int i=0; i<genotypes.length; i++){
				String line = lines[i]+"\t"+1.0+"\t"+0+"\n";
				outfile.write(line);
			}
			outfile.close();
		}
		else{
			System.out.println("Permutations: "+ perm);
			
			List<Double> pValues = calcPValues(pointwisePValues, permPValues, perm); //p-values from permutation test
			
			BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(out+File.separator+filename)));
			for(int i=0; i<genotypes.length; i++){
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
	
	public double minPValue(double[] a){
		double min = 1.0;
		for(double d:a){
			if(d<min){
				min = d;
			}
		}
		return min;
	}
		
	//Calculates pointwise p-values and populates $realPValues and $lines
	public void pointwisePValue(Gene g, double[] pointwisePValues, String[] lines){
		String genename = g.getId();
		List<ExpSample> ase = g.getExpsamples();
		int subsetSize = sampleids.length;
		//System.out.println(genotypes.length+"\t"+genotypes[0].length);
		//System.out.println(subsetSize+"\t"+ase.size());
		//System.out.println(ase.size());
		
		for(int i=0; i<genotypes.length; i++){
			//number of ones in ase samples
			int m=0;
			//number of ones in genotype samples
			int k=0;
			//number of mismatches
			int incorrect=0;
			for (int j=0; j<ase.size();j++) {
				int hasASE = ase.get(j).getASE();
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
			//System.out.println("incorrect: "+incorrect+" m: "+m+" k: "+k);
			//calculate p-value based on hypergeometric distribution
			ComputeSig sig = new ComputeSig(subsetSize, m, k, incorrect);
			double p = sig.significance();
			pointwisePValues[i]=p;
			
			double[] poss = sig.possiblePValues();
			pmap.put(combids[i], poss);
			/**
			int j = (incorrect - Math.abs(m-k))/2;
			if(poss[j]!=p){
				System.out.println("Bug in cache hypergeometric");
				System.exit(1);
			}
			**/
			String line = genename+"\t"+combids[i]+"\t"+m+"\t"+k+"\t"+subsetSize+"\t"+incorrect+"\t"+p;
			lines[i]=line;
		}

	}
	
	//returns minimum pointwise p-value from permuted data
	public double permutationPValue(List<ExpSample> ase){
		List<Double> pvals = new ArrayList<Double>();
		for(int i=0; i<genotypes.length; i++){
			/**
			List<Integer> subsetGeno = new ArrayList<Integer>();
			for(int j=0; j<genotypes[i].length; j++){
				subsetGeno.add(genotypes[i][j]);
			}
			**/
			//Collections.shuffle(ase);
			
			//System.out.println(subsetGeno.size()+"\t"+ase.size());
			
			//number of ones in ase samples
			int m=0;
			//number of ones in genotype samples
			int k=0;

			int incorrect=0;
			for (int j=0; j<ase.size();j++) {
				int hasASE = ase.get(j).getASE();
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
			//System.out.println("incorrect: "+incorrect+" m: "+m+" k: "+k);
			int j = (incorrect - Math.abs(m-k))/2;
			//System.out.println("i: "+i+" j: "+j);
			double p = pmap.get(combids[i])[j];
		
			/**
			//System.out.println(s.getId()+"\t"+subsetGeno.size()+"\t"+m+"\t"+k+"\t"+incorrect);
			ComputeSig sig = new ComputeSig(subsetGeno.size(), m, k, incorrect);
			double p = sig.significance();
		**/	
			pvals.add(p);
		}
		return Collections.min(pvals);
	}
	

	//return subset of genosamples that also have expsamples
	public void getSubset(Gene g){
		List<ExpSample> exp = g.getExpsamples();
		
		//sampleids of expsamples
		String[] ids = new String[exp.size()];
		for(int i=0; i<exp.size(); i++){
			ids[i] = exp.get(i).getID();
		}

		//find subset of genosamples that have matching expsamples
		List<ExpSample> newexpsamples = new ArrayList<ExpSample>();
		List<String> newsampleids = new ArrayList<String>();
		List<Integer> subset = new ArrayList<Integer>();
		for(int j=0; j<ids.length; j++){
			String s = ids[j];
			int i = indexSample(s, sampleids);
			if(i>=0){
				subset.add(i);
				newsampleids.add(ids[j]);
				newexpsamples.add(exp.get(j));
			}
		}
		g.replaceSamples(newexpsamples);
		
		//remove genosamples that do not have expsample
		int[][] pruned = new int[genotypes.length][subset.size()];
		for(int i=0; i<genotypes.length; i++){
			for(int j=0; j<subset.size(); j++){
				pruned[i][j] = genotypes[i][subset.get(j)];
			}
		}
		
		genotypes = pruned;
		
		String[] ret = new String[newsampleids.size()];
		for(int i=0; i<newsampleids.size(); i++){
			ret[i] = newsampleids.get(i);
		}
		sampleids= ret;
	}
	
	public int indexSample(String id, String[] genoids){
		for(int i=0; i<genoids.length; i++){
			if(genoids[i].equals(id)){
				return i;
			}
		}
		return -1;
	}
	


	private List<Double> calcPValues(double[] pointwisePValues, List<Double> permPValues, int perm) {
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
