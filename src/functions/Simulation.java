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
import parse.ParseSNP;
import parse.ParseSamples;
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
		for(String k:snpMap.keySet()){
			System.out.println(k);
		}

		ParseGenotypes.parseGenotypes(genotypes,snps, snpMap);
		System.out.println(gene.toString()+"\t"+snps.size());
	}
	
	
	public List<SNP> getSnps(){
		return snps;
	}
	public Gene getGene(){
		return gene;
	}


	
	public Object[] getSubset(int total, int num){
		//System.out.println(total);
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
	
/*	//simulation
	public int mapASE(int[] ase, int st, double f, String snpid) throws IOException{
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+"simulation"+st+"_"+snpid+"_"+f+".txt")));
		
		int variants=0;
		
		if(ase.length<sampleSize){
			System.out.println("Not enough samples");
			System.exit(1);
		}
		
		Object[] subset = getSubset(ase.length, sampleSize);

		for(SNP s:snps){
			List<GenoSample> gsamples = s.getGenosamples();
			if(gsamples.size()<sampleSize){
				continue;
			}
			
			int correct=0;
			int incorrect=0;
			for (int i=0; i<subset.length;i++) {
				int ind= (int) subset[i];
				String sampleID = gsamples.get(ind).getID();
				int isHetero = gsamples.get(ind).getHetero();
				if(isHetero == ase[ind]){
					correct++;
				}
				else if(isHetero != ase[ind]){
					incorrect++;
				}
				else{
					System.out.println("Missing data: "+sampleID );
				}
			}
			
			String line = s.getId()+"\t"+f+"\t"+correct+"\t"+incorrect;
			for(int e=0; e<=errors; e++){
				if(incorrect<=e){
					line = line+"\t"+1;
				}
				else{
					line = line+"\t"+0;
				}
			}

			outfile.write(line+"\n");

		}
		
		outfile.close();
		return variants;
	}
*/


	//mapase without permutations
	public void mapASE(String filename) throws IOException{
		//List<ExpSample> ase = gene.getExpsamples();

/**
		if(ase.size()<sampleSize){
			System.out.println("Not enough expsamples");
			System.exit(1);
		}
**/
		Random rand = new Random(17);
		int randInt = rand.nextInt(snps.size());
		SNP s = snps.get(randInt);
		double f = calculateMAF(s);
		
		List<ExpSample> ase = aseCall(s);
		System.out.println("ASE");
		for(int i=0; i<ase.size(); i++){
			System.out.println(ase.get(i).toString());
		}
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
/*	
	public int[] permute(int[] ase, int n){
		int a[] = new int[n];
		int ind[] = new int[n];
		
		for(int i=0;i<n;i++){
			ind[i] =0;
		}
		int index;
		Random rand = new Random();
		for(int i=0; i<n;i++){
			do{
				index = rand.nextInt(n);
			} while(ind[index] != 0);
			ind[index] = 1;
			a[i]= ase[index];
		}
		return a;
	}
	*/

	public double calculateMAF(SNP s){
		int hetero=0;
		System.out.println(s.getId());
		
		List<GenoSample> genos = s.getGenosamples();
		for(int i=0; i<genos.size();i++){
			System.out.println(genos.get(i).getHetero());
			if(genos.get(i).getHetero()==1){
				hetero++;
			}
		}
		

		System.out.println("hetero: "+hetero);
		double pq = 0.5*hetero/genos.size();
		System.out.println("pq: "+pq);
		double discriminant = Math.sqrt(1-4*pq);
		
		double root = (1+discriminant)/2;
		double root2 = (1-discriminant)/2;
		System.out.println("root1: "+root);
		System.out.println("root2: "+root2);
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
	
/**
 	public int[] getASE(Gene g){
		List<ExpSample> ase = g.getExpsamples();
		int[] ret = new int[ase.size()];
		for(int i=0; i<ase.size();i++){
			ret[i] = ase.get(i).getASE();
		}
		return ret;
	}
**/
/**
	public void randomSimulation() throws IOException{
		Random rand = new Random(17);
		int randInt = rand.nextInt(snps.size());
		SNP s = snps.get(randInt);
		double f = calculateMAF(s);
		
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+s.getId()+"_"+f+"_simulation.txt")));

		int[] ase = aseCall(s);	

		Map<Integer,Double> ret = new HashMap<Integer,Double>();
		for(int i=0; i<=errors;i++){
			ret.put(i, 0.0);
		}
		
		//adaptive permutation
		int adaptive=0; //number of error rates that have at least five snps that match permuted
		for(int i=0; i<perm;i++){
			int[] p = permute(ase, ase.length);
			//key: number of errors ; value: whether or not a snp was found
			Map<Integer,Integer> result = mapASE(p);
			for(int key:result.keySet()){
				if(ret.containsKey(key)){
					double value = ret.get(key);
					value = value+result.get(key);
					ret.put(key, value);
					//if 5th snp found for error rate, increment adaptive
					if(result.get(key)>0 && value==5){
						adaptive++;
					}

				}
			}
			//if all error rates have 5 snps found, stop permutations
			if(adaptive==ret.keySet().size()){
				//outfile.write(adaptive+" stopping permutations");
				System.out.println(adaptive+" stopping permutations");
				//make p-value 1.0 for all error rates
				for(int key:ret.keySet()){
					ret.put(key, 1.0);
					outfile.write(key+"\t"+1.0+"\n");
				}
				outfile.close();
			}	
		}
		
		System.out.println("Significance key set size: "+ret.keySet().size());
		for(int key:ret.keySet()){
			double val = ret.get(key);
			val = val/perm;
			ret.put(key, val);
			outfile.write(key+"\t"+val+"\n");
		}
		
		outfile.close();
		//sig = ret;


	}
	**/
/*	public double testSignificance(int[] geno, int[] ase, int mismatch) throws IOException{	
		
		int total =0; //number of permutations that have equal or less mismatches
		for(int i=0; i<perm;i++){

			int[] p = permute(geno, geno.length);
			int h = hamming(ase, p);
			
			if(h<=mismatch){
				total++;
			}
			
			if(total>25){
				break;
			}
		}
		
		if(total>25){
			return 1.0;
		}

		return total*1.0/perm;
	}*/
	

	/*private int hamming(int[] ase, int[] p){
		int errors=0;
		if(ase.length != p.length){
			return -1;
		}
		else{
			for(int i=1; i< ase.length; i++){
				if(ase[i]!=p[i]){
					errors=errors+1;
				}
			}
		}
		return errors;
	}
	*/
/*	
	private double calculatePValue(double[] perms) {
		int total = 0;
		for(double p:perms){
			if(p>1){
				total++;
			}
		}
		return 1.0*total/perms.length;
	}*/

/*	public boolean passThreshold(int incorrect, int correct){
		if (correct == 0) {
			return false;
		}
		else if (incorrect>errors) {
			return false;
		}
		return true;
	}*/
}
