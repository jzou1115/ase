package run;

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
			ArrayList<GenoSample> gsamples = s.getGenosamples();
			int correct=0;
			int incorrect=0;
			for (GenoSample g : gsamples) {
				String sampleID = g.getSampleID();
				int isHetero = g.getHetero();
				randInt = rand.nextInt(2);
				if(isHetero == randInt){
					correct++;
				}
				else if(isHetero != randInt){
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
	
	public int[] permute(int[] ase){
		
	}
	
	public double simulate(){
		
	}
	
	
	//TODO think about when correct AND incorrect can be 0, because this happened before
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