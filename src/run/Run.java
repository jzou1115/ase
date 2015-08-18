package run;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import sample.*;
import genome.*;

public class Run {

	Gene gene;
	List<SNP> snps;	
	int threshold;
	int errors;
	int perm;
	int sampleSize;
	
	public Run(Gene g, List<SNP> s, int t, int e, int p, int n){
		snps=s;
		gene=g;
		threshold=t;
		errors= e;
		perm=p;
		sampleSize=n;
	}

	
	
	public int mapASE(int[] ase){
		int variants=0;
		for(SNP s:snps){
			List<GenoSample> gsamples = s.getGenosamples();
			int correct=0;
			int incorrect=0;
			for (int i=0; i<gsamples.size();i++) {
				String sampleID = gsamples.get(i).getID();
				int isHetero = gsamples.get(i).getHetero();
				if(isHetero == ase[i]){
					correct++;
				}
				else if(isHetero != ase[i]){
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
	
	public double simulate(SNP s, int[] ase) throws FileNotFoundException{
		int total=0;
		for(int i=0; i<perm;i++){
			int[] p = permute(ase, ase.length);
			int a = mapASE(p);
			total=total+a;
		}
		return 1.0*total/perm;
	}
	
	public double calculateMAF(SNP s){
		int hetero=0;
		
		List<GenoSample> genos = s.getGenosamples();
		for(int i=0; i<genos.size();i++){
			if(genos.get(i).getHetero()==1){
				hetero++;
			}
		}
		if(hetero>=genos.size()){
			return 1.0*(genos.size()-hetero)/genos.size();
		}
		return 1.0*hetero/genos.size();
	}
	
	public int[] aseCall(SNP s){
		List<GenoSample> genos = s.getGenosamples();
		int[] ret = new int[genos.size()];
		for(int i=0; i<genos.size();i++){
			ret[i] = genos.get(i).getHetero();
		}
		return ret;
	}
	
	public void allSimulations() throws IOException{
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(gene.getId()+"_simulation.txt")));
		outfile.write("GeneID\tMAF\tNumVariants\tSimulationMean\n");
		for(SNP s: snps){
			double f = calculateMAF(s);
			int[] ase = aseCall(s);
			double mean = simulate(s, ase);		
			int x= mapASE(ase);
			outfile.write(s.getId()+"\t"+f+"\t"+x+"\t"+mean+"\n");
		}
	}
	
	public boolean passThreshold(int incorrect, int correct){
		if (correct == 0) {
			return false;
		}
		else if (incorrect>threshold) {
			return false;
		}
		return true;
	}
}