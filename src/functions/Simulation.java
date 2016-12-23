package functions;

import genome.Gene;
import genome.SNP;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import parse.ParseGenotypes;
import parse.ParseMap;
import sample.ExpSample;
import sample.GenoSample;

public class Simulation {
	Gene gene;
	List<SNP> snps;
	Map<String,SNP> snpMap;
	List<String> sampleNames;
	
	private File outdir;
	
	/**
	 * Reads input files and sets state for simulation class
	 * @param map Map file for genes and SNPs
	 * @param g Gene id
	 * @param genotypes Genotype file
	 * @param out Output directory
	 * @throws IOException
	 */
	public void setTestGene(InputStream map, String g, InputStream genotypes, File out) throws IOException{
		outdir = out;
				
		ParseMap parsemap = new ParseMap();
		parsemap.parseMap(map, g);
		gene = parsemap.getGene();
		snps = parsemap.getSNPs();
		snpMap = parsemap.getSnpMap();

		ParseGenotypes.parseGenotypes(genotypes,snps, snpMap);

	}
	

	
	public List<Integer> getSubset(int total, int num){
		Random rand = new Random();
		List<Integer> ret = new ArrayList<Integer>();
		int i=0;
		int j;
		while(i<num){
			j= rand.nextInt(total);
			if(!ret.contains(j)){
				ret.add(j);
				i++;
			}
		}

		return ret;
	}
/*	
	//return subset of genosamples that also have expsamples
	public List<GenoSample> getSubset(int total, SNP s, List<ExpSample> exp){
		Set<String> sampleids = new HashSet<String>();
		for(ExpSample e:exp){
			sampleids.add(e.getID());
		}
		
		List<GenoSample> geno = new ArrayList<GenoSample>(s.getGenosamples());
		Set<String> sampleids2 = new HashSet<String>();
		for(int j=0; j< geno.size(); j++){
			sampleids2.add(geno.get(j).getID());
		}
		
		sampleids.retainAll(sampleids2);
		
		List<GenoSample> subset = new ArrayList<GenoSample>();
		for(GenoSample g:geno){
			if(sampleids.contains(g.getID())){
				subset.add(g);
			}
		}
		return subset;
	}*/


	/**
	 * Calls ASE based on randomly chosen SNP for sample size n.  Runs ASE mapping algorithm on simulated data. Outputs p-value for all possible numbers of errors using closed form solution for all possible p-values.
	 * @param filename Output file name for simulation p-values
	 * @throws IOException
	 */
	public void mapASE(String filename, int n) throws IOException{

		//choose random SNP for simulation
		Random rand = new Random(17);
		int randInt = rand.nextInt(snps.size());
		SNP s = snps.get(randInt);
		
		//minor allele frequency of SNP s
		double f = calculateMAF(s);
		
		List<GenoSample> allGenotypes = s.getGenosamples();
		
		//Put random sample of size n in subset
		List<Integer> indices = getSubset(s.getGenosamples().size(), n);
		List<GenoSample> subset = new ArrayList<GenoSample>();
		for(int i:indices){
			subset.add(allGenotypes.get(i));
		}
		//ASE calls for subset of individuals
		List<ExpSample> ase = aseCall(subset);
		if(ase.size()!=subset.size()){
			System.out.println("Subset function malfunction");
			System.exit(1);
		}
		
		System.out.println("Number of individuals: "+ase.size());
		
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+filename)));

		//number of ones in ase samples
		int m=0;
		//number of ones in genotype samples
		int k=0;
		//number of mismatches between ASE status and genotype
		int incorrect=0;
		for (int i=0; i<subset.size();i++) {
			GenoSample g = subset.get(i);
			int hasASE = ase.get(i).getASE();
			
			if(hasASE==1){
				m=m+1;
			}
			int isHetero = g.getHetero();
			if(isHetero==1){
				k=k+1;
			}

			if(isHetero != hasASE){
				incorrect++;
			}	

		}
		System.out.println("Incorrect: "+incorrect);
		//set state for ComputeSig state
		System.out.println(subset.size()+"\t"+m+"\t"+k+"\t"+incorrect);
		
		ComputeSig sig = new ComputeSig(subset.size(), m, k, incorrect);
		//p-values for all possible numbers of errors
		double[] poss = sig.significance();
		for(int j=0 ; j< poss.length ; j ++){
			//number of mismatches
			int inc = 2*j + Math.abs(m-k);
			//p-value for inc mismatches
			double p = poss[j];
			String line = gene.toString()+"\t"+s.getId()+"\t"+f+"\t"+m+"\t"+k+"\t"+ase.size()+"\t"+inc+"\t"+p;
			outfile.write(line+"\n");
		}
		outfile.close();

	}

	/**
	 * Use Hardy Weinburg to find minor allele frequency of SNP
	 * @param s SNP
	 * @return minor allele frequency of SNP
	 */
	public double calculateMAF(SNP s){
		int hetero=0;
		
		List<GenoSample> genos = s.getGenosamples();
		for(int i=0; i<genos.size();i++){
			if(genos.get(i).getHetero()==1){
				hetero++;
			}
		}

		double pq = 0.5*hetero/genos.size();
		double discriminant = Math.sqrt(1-4*pq);
		
		double root = (1+discriminant)/2;
		double root2 = (1-discriminant)/2;

		if(root>=0 && root<=1){
			if(root<1-root){
				return root;
			}
			return 1-root;
		}
		
		if(root2>=0 && root2<=1){
			if(root2<1-root2){
				return root2;
			}
			return 1-root2;
		}
		return -1;
	}
	
	/**
	 * Use genotypes of a SNP to call ASE for individuals with zero errors
	 * @param genos List of genotypes
	 * @return List of ASE calls
	 */
	public List<ExpSample> aseCall(List<GenoSample> genos){
		List<ExpSample> ret = new ArrayList<ExpSample>();
		for(int i=0; i<genos.size();i++){
			GenoSample g = genos.get(i);
			ret.add(new ExpSample(g.getID(), g.getHetero()));
		}
		return ret;
	}
	
	public List<SNP> getSnps(){
		return snps;
	}
	public Gene getGene(){
		return gene;
	}


}
