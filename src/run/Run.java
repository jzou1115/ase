package run;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import sample.ExpSample;
import sample.GenoSample;
import sample.Sample;
import gene.*;
import snp.*;
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
		Map<String, ExpSample> emap= gene.getExpsamples();
		Random rand = new Random();
		int randInt = rand.nextInt(2);

		for(SNP s:snps){
			Map<String, GenoSample> gmap= s.getGenosamples();
			int correct=0;
			int incorrect=0;
			for(String gtexId:gmap.keySet()){
				randInt = rand.nextInt(2);
				
				if(emap.containsKey(gtexId) & (gmap.get(gtexId).getHetero() == randInt)){
					correct++;
				}
				else if(emap.containsKey(gtexId) & (gmap.get(gtexId).getHetero() != randInt)){
					incorrect++;
				}
				else{
					System.out.println("Missing data: "+gtexId);
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
			Map<String, GenoSample> gmap= s.getGenosamples();
			int correct=0;
			int incorrect=0;
			for(String gtexId:gmap.keySet()){
				if(emap.containsKey(gtexId) & (gmap.get(gtexId).getHetero() == emap.get(gtexId).getASE())){
					correct++;
				}
				else if(emap.containsKey(gtexId) & (gmap.get(gtexId).getHetero() != emap.get(gtexId).getASE())){
					incorrect++;
				}
				else{
					System.out.println("Missing data: "+gtexId);
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