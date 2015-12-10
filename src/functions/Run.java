package functions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import sample.*;
import genome.*;

public class Run {

	private Gene gene;
	private List<SNP> snps;	
	private int errors;
	private int perm;
	private int sampleSize;
	private File outdir;
	private Map<Integer,Double> sig;

	public Run(Gene g, List<SNP> s, int e, int p, int n, File out){
		snps=s;
		gene=g;
		errors= e;
		perm=p;
		sampleSize=n;
		outdir=out;
	}

	
	public List<SNP> getSnps(){
		return snps;
	}
	public Gene getGene(){
		return gene;
	}

	public int getErrors(){
		return errors;
	}
	public int getPerm(){
		return perm;
	}
	public int sampleSize(){
		return sampleSize;
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
	
	public Object[] getSubset(int total, int num, SNP s, List<ExpSample> exp){
		List<Integer> tried = new ArrayList<Integer>();
		Random rand = new Random();
		List<Integer> ret = new ArrayList<Integer>();
		int i=0;
		int j;
		while(i<num){
			if(tried.size()==total){
				return ret.toArray();
			}
			j= rand.nextInt(total);
			if(!tried.contains(j)){
				tried.add(j);
				String gtexID = exp.get(j).getID();
				if(s.getSample(gtexID)!=null){
					ret.add(j);
					i++;
				}
			}
		}

		return ret.toArray();
	}
	
	//simulation
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
	
	//testsignificance
	public Map<Integer,Integer> mapASE(int[] ase) throws IOException{
		Map<Integer,Integer> ret = new HashMap<Integer, Integer>();
		for(int i=0; i<=errors;i++){
			ret.put(i, 0);
		}
		
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
			if(passThreshold(incorrect, correct)){
				for(int j: ret.keySet()){
					if(j>=incorrect && ret.get(j)==0){
						ret.put(j, 1);
						variants++;
					}
				}
			}
			
			if(variants==(errors+1)){

				return ret;
			}

		}

		return ret;
	}

	/**
	
	//temp debug perm method
	public void mapASE(String st) throws IOException{

		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+st+"_permASE.txt")));

		List<ExpSample> ase = gene.getExpsamples();
		
		if(ase.size()<sampleSize){
			System.out.println("Not enough samples");
			System.exit(1);
		}
		
		for(SNP s:snps){
			Object[] subset = getSubset(ase.size(), sampleSize, s, gene.getExpsamples());
			if(subset!=null){
				int correct=0;
				int incorrect=0;
				int[] geno = new int[sampleSize];
				int[] aseSub = new int[sampleSize];
				int m=0;
				int k=0;
				for (int i=0; i<subset.length;i++) {
					int ind= (int) subset[i];
					ExpSample e = ase.get(ind);
					String sampleID = e.getID();
					GenoSample g = s.getSample(sampleID);

					int hasASE = e.getASE();
					if(hasASE==1){
						m++;
					}
					int isHetero = g.getHetero();
					if(isHetero==1){
						k++;
					}
					geno[i] = isHetero;
					aseSub[i] = hasASE;
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
				ComputeSig si = new ComputeSig(sampleSize, m, k, incorrect);
				double p = si.significance();
				
				double[] sig = testSignificance(geno,aseSub, incorrect);
				String line = s.getId()+"\t"+correct+"\t"+incorrect+"\t"+p;
				outfile.write(line+"\n");
				for(int a=0; a<sig.length;a++){
					outfile.write(sig[a]+"\n");
				}
				break;
			}

		}
		//outfile.write("variants: "+variants+"\n");
		outfile.close();
	}
	**/
	
	//mapase with permuations
	public void mapASE(String st) throws IOException{

		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+st+"_permASE.txt")));

		List<ExpSample> ase = gene.getExpsamples();
		
		if(ase.size()<sampleSize){
			System.out.println("Not enough samples");
			System.exit(1);
		}

		for(SNP s:snps){
			Object[] subset = getSubset(ase.size(), sampleSize, s, gene.getExpsamples());
			if(subset==null){
				continue;
			}
			int correct=0;
			int incorrect=0;
			int[] geno = new int[sampleSize];
			int[] aseSub = new int[sampleSize];
			for (int i=0; i<subset.length;i++) {
				int ind= (int) subset[i];
				ExpSample e = ase.get(ind);
				String sampleID = e.getID();
				GenoSample g = s.getSample(sampleID);
				
				int hasASE = e.getASE();
				int isHetero = g.getHetero();
				geno[i] = isHetero;
				aseSub[i] = hasASE;
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

			double sig = testSignificance(geno,aseSub, incorrect);
			String line = s.getId()+"\t"+correct+"\t"+incorrect+"\t"+sig;
	
			if(isSignificant(errors, incorrect, sig)){
				line = line+"\t"+1;
			}
			else{
				line = line+"\t"+0;
			}
		
			if(subset.length!=0){
				outfile.write(line+"\n");
			}
		}
		//outfile.write("variants: "+variants+"\n");
		outfile.close();
	}


	//mapase without permutations
	public void mapASE(String st, String filename) throws IOException{
		List<ExpSample> ase = gene.getExpsamples();
		if(ase.size()<sampleSize){
			System.out.println("Not enough expsamples");
			System.exit(1);
		}
	

		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+st+filename)));

		
		for(SNP s:snps){
			Object[] subset = getSubset(ase.size(), sampleSize, s, ase);
			if(subset==null || subset.length==0){
				System.out.println("Not enough genosamples in "+s.getId());
				continue;
			}
			
			//m and k are used to approximate significance
			//number of ones in ase samples
			int m=0;
			//number of ones in genotype samples
			int k=0;
			
			int correct=0;
			int incorrect=0;
			for (int i=0; i<subset.length;i++) {
				int ind= (int) subset[i];
				ExpSample e = ase.get(ind);
				String sampleID = e.getID();
				GenoSample g = s.getSample(sampleID);

				int hasASE = e.getASE();
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

			ComputeSig sig = new ComputeSig(sampleSize, m, k, incorrect);
			double p = sig.significance();
			String line = s.getId()+"\t"+correct+"\t"+incorrect+"\t"+p;
			if(isSignificant(errors, incorrect, p)){
				line = line+"\t"+1;
			}
			else{
				line = line+"\t"+0;
			}

			outfile.write(line+"\n");
		}


	
		//outfile.write("variants: "+variants+"\n");
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
	
	public int[] aseCall(SNP s){
		List<GenoSample> genos = s.getGenosamples();
		int[] ret = new int[genos.size()];
		for(int i=0; i<genos.size();i++){
			ret[i] = genos.get(i).getHetero();
		}
		return ret;
	}
	
	public int[] getASE(Gene g){
		List<ExpSample> ase = g.getExpsamples();
		int[] ret = new int[ase.size()];
		for(int i=0; i<ase.size();i++){
			ret[i] = ase.get(i).getASE();
		}
		return ret;
	}
	
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
	
	public double testSignificance(int[] geno, int[] ase, int mismatch) throws IOException{	
		
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
	}
	/**
	//temp debug perm method
	public double[] testSignificance(int[] geno, int[] ase, int mismatch) throws IOException{	
		double[] pvals = new double[perm];
		int total=0; //number of permutations that have equal or less mismatches
		for(int i=0; i<perm;i++){

			int[] p = permute(geno, geno.length);
			int h = hamming(ase, p);
			
			if(h<=mismatch){
				total++;
			}
			pvals[i] = total*1.0/(i+1);
			
		}
		
		return pvals;
	}
	**/
	private int hamming(int[] ase, int[] p){
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
	
	
	private double calculatePValue(double[] perms) {
		int total = 0;
		for(double p:perms){
			if(p>1){
				total++;
			}
		}
		return 1.0*total/perms.length;
	}

	public boolean passThreshold(int incorrect, int correct){
		if (correct == 0) {
			return false;
		}
		else if (incorrect>errors) {
			return false;
		}
		return true;
	}
}