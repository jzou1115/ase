package run;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import sample.*;
import genome.*;

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
			if(passThreshold(incorrect, correct)){
				//System.out.println(gene.getId()+"\t"+s.getId());
				variants++;
			}
		}
		return variants;
	}
	
	
	public int run() throws Exception{
		Map<String, ExpSample> emap= gene.getExpsamples();
		for(SNP s:snps){
			ArrayList<GenoSample> gsamples= s.getGenosamples();
			int correct=0;
			int incorrect=0;
			for (GenoSample g : gsamples) {
				String sampleID = g.getSampleID();
				int isHetero = g.getHetero();
				if (emap.containsKey(sampleID)){
					if(isHetero == emap.get(sampleID).getASE()){
						correct++;
					}
					else{
						incorrect++;
					}
				}				 
				else {
					throw new Exception("Missing data: "+sampleID + "SNP id:" + s.getId());
				}
			}
			if(passThreshold(incorrect, correct)){
				//System.out.println(gene.getId()+"\t"+s.getId());
				variants++;
			}
		}
		return variants;
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