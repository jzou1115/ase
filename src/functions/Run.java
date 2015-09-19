package functions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import sample.*;
import genome.*;

public class Run {

	private Gene gene;
	private List<SNP> snps;	
	private double threshold;
	private int errors;
	private int perm;
	private int sampleSize;
	private File outdir;

	public Run(Gene g, List<SNP> s, double t, int e, int p, int n, File out){
		snps=s;
		gene=g;
		threshold=t;
		errors= e;
		perm=p;
		sampleSize=n;
		outdir=out;
	}

	public Run(Gene g, List<SNP> s, double t, int e, int n, File out){
		snps=s;
		gene=g;
		threshold=t;
		errors= e;
		sampleSize=n;
		outdir=out;
	}
	
	public List<SNP> getSnps(){
		return snps;
	}
	public Gene getGene(){
		return gene;
	}
	public double getThreshold(){
		return threshold;
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
		System.out.println(total);
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
		System.out.println(total);
		Random rand = new Random();
		List<Integer> ret = new ArrayList<Integer>();
		int i=0;
		int j;
		while(i<num){
			if(tried.size()==total){
				return null;
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
	
	public int mapASE(int[] ase){
		int variants=0;
		Object[] subset = getSubset(ase.length, sampleSize);
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
				variants++;
			}
		}
		return variants;
	}
	
	public int mapASE(int[] ase, String st) throws IOException{
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+st+"_simulation.txt")));

		int variants=0;
		Object[] subset = getSubset(ase.length, sampleSize);
		outfile.write("subset len: "+subset.length+"\n");
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
				variants++;
				outfile.write(s.getId()+"\t"+incorrect+"\n");
			}
		}
		outfile.write("variants: "+variants+"\n");
		outfile.close();
		return variants;
	}
	
	public int mapASE(String st) throws IOException{
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+st+"_mapase.txt")));

		List<ExpSample> ase = gene.getExpsamples();
		//System.out.println("Num ExpSamples: "+ase.size());
		
		int samples;
		if(ase.size()<sampleSize){
			samples = ase.size();
		}
		else{
			samples = sampleSize;
		}
		
		
		outfile.write("subset len: "+samples+"\n");
		int variants=0;
		int total =0;
		for(SNP s:snps){
			total++;
			//System.out.println(s.getId()+"\t"+total);
			Object[] subset = getSubset(ase.size(), samples, s, gene.getExpsamples());
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
				outfile.write(s.getId()+"\t"+incorrect+"\n");
			}
		}
		outfile.write("variants: "+variants+"\n");
		outfile.close();
		return variants;
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
	
	
	public int allSimulations() throws IOException{
		int pass=0;
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+gene.getId()+"_simulation_"+perm+"_"+sampleSize+"_"+errors+"_"+threshold+".txt")));
		outfile.write("Number of permutations: "+perm+"\n");
		outfile.write("Number of samples: "+sampleSize+"\n");
		outfile.write("Maximum number of errors: "+errors+"\n");
		//outfile.write("Significance: "+threshold+"\n");
		outfile.write("GeneID\tMAF\tNumVariants\tSimulationMean\tP-value\tTotalSNPs\tReduction\n");
		for(SNP s: snps){
			double f = calculateMAF(s);
			//System.out.println(s.getId()+"\t"+f);
			int[] ase = aseCall(s);	
			int x= mapASE(ase, outdir+File.separator+s.getId());
			
			int total=0;
			double[] perms = new double[perm];
			for(int i=0; i<perm;i++){
				int[] p = permute(ase, ase.length);
				int a = mapASE(p);
				total=total+a;
				perms[i]=(double) a;
			}
			
			double p = calculatePValue(perms);
			
			double mean= 1.0*total/perm;
			
			int numSNPs = snps.size();
			double reduction = 1.0-1.0*x/numSNPs;
			
			if(p<threshold){
				pass++;
			}
			outfile.write(s.getId()+"\t"+f+"\t"+x+"\t"+mean+"\t"+p+"\t"+numSNPs+"\t"+reduction+"\n");
		}
		//outfile.write("power="+pass+"/"+perm+"="+1.0*pass/perm+"\n");
		outfile.close();
		return pass;
	}
	
	
	public int randomSimulation() throws IOException{
		int pass=0;
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+gene.getId()+"_simulation_"+perm+"_"+sampleSize+"_"+errors+"_"+threshold+".txt")));
		outfile.write("Number of permutations: "+perm+"\n");
		outfile.write("Number of samples: "+sampleSize+"\n");
		outfile.write("Maximum number of errors: "+errors+"\n");
		
		Random rand = new Random();
		int randInt = rand.nextInt(snps.size());
		SNP s = snps.get(randInt);
		outfile.write("SNP used: "+s.toString()+"\n");

		double f = calculateMAF(s);
		int[] ase = aseCall(s);	
		//int x= mapASE(ase, outdir+File.separator+s.getId());

		int total=0;
		double[] perms = new double[perm];
		for(int i=0; i<perm;i++){
			int[] p = permute(ase, ase.length);
			int a = mapASE(p);
			total=total+a;
			perms[i]=(double) a;
			System.out.println("PermNum"+i+"\t"+a);
		}

		double p = calculatePValue(perms);

		double mean= 1.0*total/perm;

		int numSNPs = snps.size();

		if(p<threshold){
			pass++;
		}
		outfile.write("GeneID\tMAF\tSimulationMean\tP-value\tTotalSNPs\n");
		outfile.write(s.getId()+"\t"+f+"\t"+mean+"\t"+p+"\t"+numSNPs+"\n");

		//outfile.write("power="+pass+"/"+perm+"="+1.0*pass/perm+"\n");
		outfile.close();
		return pass;
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