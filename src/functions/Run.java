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

	/**
	public Run(Gene g, List<SNP> s, int e, int n, File out){
		snps=s;
		gene=g;
		errors= e;
		sampleSize=n;
		outdir=out;
	}
	**/
	
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
			if(passThreshold(incorrect, correct)){
				outfile.write(s.getId()+"\t"+correct+"\t"+incorrect+"\t"+1+"\n");
			}
			else{
				outfile.write(s.getId()+"\t"+correct+"\t"+incorrect+"\t"+0+"\n");	
			}
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
				if(ret.get(incorrect)==0){
					ret.put(incorrect, 1);
					variants++;
				}
			}
			
			if(variants==(errors+1)){
				return ret;
			}

		}
		
		return ret;
	}
	/**
	public int mapASE(int[] ase, String st) throws IOException{
		int variants=0;
		Object[] subset = getSubset(ase.length, sampleSize);
		
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+st+"_"+subset.length+"_simulation.txt")));
		for(SNP s:snps){
			List<GenoSample> gsamples = s.getGenosamples();
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
				outfile.write(s.getId()+"\t"+correct+"\t"+incorrect+"\t"+1+"\n");
			}
			else{
				outfile.write(s.getId()+"\t"+correct+"\t"+incorrect+"\t"+0+"\n");	
			}
		}
		outfile.close();
		return variants;
	}
	**/
	
	//mapase
	public int mapASE(String st) throws IOException{
		
		Map<Integer,Double> pval = testSignificance();
		
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+st+"_mapase.txt")));

		List<ExpSample> ase = gene.getExpsamples();
		
		if(ase.size()<sampleSize){
			System.out.println("Not enough samples");
			System.exit(1);
		}

		int variants=0;
		for(SNP s:snps){
			Object[] subset = getSubset(ase.size(), sampleSize, s, gene.getExpsamples());
			if(subset==null){
				continue;
			}
			int correct=0;
			int incorrect=0;
			for (int i=0; i<subset.length;i++) {
				int ind= (int) subset[i];
				ExpSample e = ase.get(ind);
				String sampleID = e.getID();
				GenoSample g = s.getSample(sampleID);
				
				int hasASE = e.getASE();
				int isHetero = g.getHetero();
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
			if(passThreshold(incorrect, correct)){
				variants++;
			}
			String line = s.getId()+"\t"+correct+"\t"+incorrect;
			for(int e=0; e<=errors; e++){
				if(isSignificant(e, incorrect, sig.get(e))){
					line = line+"\t"+1;
				}
				else{
					line = line+"\t"+0;
				}
			}
			outfile.write(line+"\n");
		}
		//outfile.write("variants: "+variants+"\n");
		outfile.close();
		return variants;
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
	
	public int randomSimulation() throws IOException{
		int pass=0;
		//outfile.write("Number of permutations: "+perm+"\n");
		//outfile.write("Number of samples: "+sampleSize+"\n");
		//outfile.write("Maximum number of errors: "+errors+"\n");
		
		Random rand = new Random(13);
		int randInt = rand.nextInt(snps.size());
		SNP s = snps.get(randInt);
	//	BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+gene.getId()+"_simulation_"+perm+"_"+sampleSize+"_"+errors+"_"+threshold+"_"+s.getId()+"_"+itr+".txt")));
		//outfile.write("SNP used: "+s.toString()+"\n");

		double f = calculateMAF(s);
		int[] ase = aseCall(s);	
		//int x= mapASE(ase, outdir+File.separator+s.getId());

		int total=0;
		double[] perms = new double[perm];
		for(int i=0; i<perm;i++){
			if(total>=5){
				//outfile.write("Not Significant");
				break;
			}
			int[] p = permute(ase, ase.length);
			int a = mapASE(p, i, f, s.getId());
			perms[i]=(double) a;
			if(a>1){
				total=total+1;
				//outfile.write("PermNum"+i+"\t"+a+"\t"+1+"\n");
				System.out.println("PermNum"+i+"\t"+a+"\t"+1);
			}
			else{
				//outfile.write("PermNum"+i+"\t"+a+"\t"+0+"\n");
				System.out.println("PermNum"+i+"\t"+a+"\t"+0);
			}
		}

		/**
	//	double p = calculatePValue(perms);

		//double mean= 1.0*total/perm;

		//int numSNPs = snps.size();

		//if(p<threshold){
			pass++;
		}
		outfile.write("GeneID\tMAF\tSimulationMean\tP-value\tTotalSNPs\n");
		outfile.write(s.getId()+"\t"+f+"\t"+mean+"\t"+p+"\t"+numSNPs+"\n");

		//outfile.write("power="+pass+"/"+perm+"="+1.0*pass/perm+"\n");
		 **/
		//outfile.close();
		return pass;
	}
	
	public Map<Integer,Double> testSignificance() throws IOException{
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+gene.getId()+"_significance.txt")));

		Random rand = new Random(13);
		int randInt = rand.nextInt(snps.size());
		SNP s = snps.get(randInt);

		int[] ase = aseCall(s);	

		Map<Integer,Double> ret = new HashMap<Integer,Double>();
		for(int i=0; i<=errors;i++){
			ret.put(i, 0.0);
		}
		
		for(int i=0; i<perm;i++){
			//adaptive permutation
			int adaptive=0;
			
			int[] p = permute(ase, ase.length);
			Map<Integer,Integer> result = mapASE(p);
			for(int key:result.keySet()){
				if(ret.containsKey(key)){
					double value = ret.get(key);
					value = value+result.get(key);
					ret.put(key, value);
				}
				if(ret.get(key)>5){
					adaptive++;
				}
			}
			
			if(adaptive==ret.keySet().size()){
				outfile.write(adaptive+" errors detected, stopping permutations");
				System.out.println(adaptive+" errors detected, stopping permutations");
				for(int key:ret.keySet()){
					ret.put(key, 1.0);
				}
				sig = ret;
				return ret;
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
		sig = ret;
		return ret;
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