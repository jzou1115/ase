package functions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
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
	//private Map<SNP,List<GenoSample>> genoMap;
	
	//public Run(Gene g, List<SNP> s, double t, int e, int p, int n, Map<SNP, List<GenoSample>> geno){
	public Run(Gene g, List<SNP> s, double t, int e, int p, int n, File out){
		snps=s;
		gene=g;
		threshold=t;
		errors= e;
		perm=p;
		sampleSize=n;
		outdir=out;
		//genoMap = geno;
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
//	public Map<SNP,List<GenoSample>> getGenoMap(){
//		return genoMap;
//	}
	
	public Object[] getSubset(int total){
		System.out.println(total);
		Random rand = new Random();
		List<Integer> ret = new ArrayList<Integer>();
		int i=0;
		int j;
		while(i<sampleSize){
			j= rand.nextInt(total);
			if(!ret.contains(j)){
				ret.add(j);
				i++;
			}
		}

		return ret.toArray();
	}
	
	public int mapASE(int[] ase){
		int variants=0;
		Object[] subset = getSubset(ase.length);
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
				//System.out.println(gene.getId()+"\t"+s.getId());
				variants++;
			}
		}
		return variants;
	}
	
	public int mapASE(int[] ase, String st) throws IOException{
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(st+".txt")));

		int variants=0;
		Object[] subset = getSubset(ase.length);
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
				//System.out.println(gene.getId()+"\t"+s.getId());
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
		//System.out.println("calc map for "+s.getId());
		//if(genoMap.containsKey(s)){
			//List<GenoSample> genos = genoMap.get(s);
		List<GenoSample> genos = s.getGenosamples();
			for(int i=0; i<genos.size();i++){
				if(genos.get(i).getHetero()==1){
					hetero++;
				}
			}
			if(hetero>1.0*genos.size()/2){
				//System.out.println(1.0*(genos.size()-hetero)/genos.size());
				return 1.0*(genos.size()-hetero)/genos.size();
			}
			//System.out.println(1.0*hetero/genos.size());
			return 1.0*hetero/genos.size();
		//}
		//System.out.println("no key in genomap");
		//return 0.0;
	}
	
	public int[] aseCall(SNP s){
		//List<GenoSample> genos = genoMap.get(s);
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
			
			double p = calculatePValue(x, perms);
			
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
	
	
	private double calculatePValue(int x, double[] perms) {
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