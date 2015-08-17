package run;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import sample.ExpSample;
import sample.GenoSample;
import gene.*;
import snp.*;

public class Run {
	int threshold;
	List<SNP> snps;
	Gene gene;
	
	int variants;
	
	public Run(Gene g, List<SNP> s, int t){
		snps=s;
		gene=g;
		threshold=t;
		variants = 0;
	}

	
	public int runSim(){
		Random rand = new Random();
		int randInt = rand.nextInt(2);

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
			if(passThreshold(incorrect)){
				//System.out.println(gene.getId()+"\t"+s.getId());
				variants++;
			}
		}
		return variants;
	}
	
	
	public int run(){
		Map<String, ExpSample> emap= gene.getExpsamples();
		for(SNP s:snps){
			ArrayList<GenoSample> gsamples= s.getGenosamples();
			int correct=0;
			int incorrect=0;
			for (GenoSample g : gsamples) {
				String sampleID = g.getSampleID();
				int isHetero = g.getHetero();
				//TODO: use a separate case for when emap doesn't exist, so both correct + incorrect = 0
				if(emap.containsKey(sampleID) && (isHetero == emap.get(sampleID).getASE())){
					correct++;
				}
				else if(emap.containsKey(sampleID) && (isHetero != emap.get(sampleID).getASE())){
					incorrect++;
				}
				else{
					System.out.println("Missing data: "+sampleID + "SNP id:" + s.getId());
				}
			}
			if(passThreshold(incorrect)){
				//System.out.println(gene.getId()+"\t"+s.getId());
				variants++;
			}
		}
		return variants;
	}
	
	
	public boolean passThreshold(int incorrect){
		if(incorrect>threshold){
			return false;
		}
		return true;
	}
	
	public boolean passThreshold(int correct, int incorrect){
		if(1.0*correct/incorrect >= threshold){
			return true;
		}
		return false;
	}
	

}