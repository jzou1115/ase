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
	private Map<Integer,Double> sig;
	
	public void setTestGene(InputStream map2, String gene2, InputStream genotypes, File out) throws IOException{
		outdir = out;
				
		ParseMap parsemap = new ParseMap();
		parsemap.parseMap(map2, gene2);
		gene = parsemap.getGene();
		snps = parsemap.getSNPs();
		snpMap = parsemap.getSnpMap();

		ParseGenotypes.parseGenotypes(genotypes,snps, snpMap);

	}
	
	
	public List<SNP> getSnps(){
		return snps;
	}
	public Gene getGene(){
		return gene;
	}


	
	public Object[] getSubset(int total, int num){
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

		return ret.toArray();
	}
	
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
	}


	//mapase without permutations
	public void mapASE(String filename) throws IOException{

		Random rand = new Random(17);
		int randInt = rand.nextInt(snps.size());
		SNP s = snps.get(randInt);
		double f = calculateMAF(s);
		
		List<ExpSample> ase = aseCall(s);
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+filename)));


		//subset of genosamples that have expression samples
		List<GenoSample> subset = getSubset(ase.size(), s, ase);
		if(subset==null){
			System.out.println("Subset size 0");
			System.exit(1);
		}

		//m and k are used to approximate significance
		//number of ones in ase samples
		int m=0;
		//number of ones in genotype samples
		int k=0;

		int correct=0;
		int incorrect=0;
		for (int i=0; i<subset.size();i++) {
			GenoSample g = subset.get(i);
			String sampleID = g.getID();

			int hasASE = ase.get(i).getASE();
			if(hasASE==1){
				m=m+1;
			}
			int isHetero = g.getHetero();
			if(isHetero==1){
				k=k+1;
			}

			if(isHetero == hasASE){
				correct++;
			}
			else if(isHetero != hasASE){
				incorrect++;
			}	
			else{
				System.out.println("Missing data: "+sampleID );
			}
		}

		ComputeSig sig = new ComputeSig(subset.size(), m, k, incorrect);
		double[] poss = sig.significance();
		for(int j=0 ; j< poss.length ; j ++){
			int inc = 2*j + Math.abs(m-k);
			double p = poss[j];
			String line = gene.toString()+"\t"+s.toString()+"\t"+f+"\t"+inc+"\t"+p;
			outfile.write(line+"\n");
		}
		outfile.close();

	}

	public boolean isSignificant(int e, int incorrect, double z){
		if(incorrect<=e){
			if(z<=.0000025){
				return true;
			}
		}
		return false;
	}


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
	
	public List<ExpSample> aseCall(SNP s){
		List<GenoSample> genos = s.getGenosamples();
		List<ExpSample> ret = new ArrayList<ExpSample>();
		for(int i=0; i<genos.size();i++){
			GenoSample g = genos.get(i);
			ret.add(new ExpSample(g.getID(), g.getHetero()));
		}
		return ret;
	}

}
