package run;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import sample.ExpSample;
import sample.GenoSample;
import sample.Sample;
import gene.*;
import snp.*;
import genome.*;

public class Run {
	double threshold;
	SNPgroup snps;
	GeneGroup genes;
	Map<Gene,List<SNP>> map;
	
	Map<String, Result> results;
	
	int sim;
	
	Map<String, Integer> counts;
	
	public Run(SNPgroup s, GeneGroup g, double t){
		snps=s;
		genes=g;
		threshold=t;
		
		results= new HashMap<String,Result>();

		snpsToGenes();
	}
	
	public void snpsToGenes(){
		map = new HashMap<Gene, List<SNP>>();
		
		List<SNP> snpList = snps.getSnps();
		//Collections.sort(snpList);
		
		List<Gene> geneList = genes.getGenes();
		//Collections.sort(geneList);
		
		for(SNP s:snpList){
			for(Gene g:geneList){
				if(match(s,g)){
					if(map.get(g)==null){
						List<SNP> val = new ArrayList<SNP>();
						val.add(s);
						map.put(g, val);
					}
					else{
						map.get(g).add(s);
					}
					break;
				}
			}
		}
		/**
		for(Gene g:map.keySet()){
			List<SNP> s = map.get(g);
			for(SNP x:s){
				System.out.println(x.getId()+"\t"+g.getId());
			}
		}
		**/
	
	}
	
	public boolean match(SNP s, Gene g){
		if(g.region.expand(100).contains(s.getLocation())){
			return true;
		}
		return false;
	}
	
	public void run(){
		for(Gene g:genes.getGenes()){
			List<SNP> snpList = map.get(g);
			Map<String, ExpSample> emap= g.getExpsamples();
			for(SNP s:snpList){
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
				Result r = new Result(correct, incorrect, s.getId());
				results.put(s.getId(), r);
				//System.out.println(r.toString());
			}
			
		}
	}
	
	public void runAll(int n){
		counts = new HashMap<String, Integer>();
		
		for(int i=0; i<n; i++){
			runSim();
		}
		
		for(String s:counts.keySet()){
			System.out.println(s+"\t"+counts.get(s));
		}
	}
	
	
	public boolean passThreshold(int correct, int incorrect){
		if(1.0*correct/incorrect >= threshold){
			return true;
		}
		return false;
	}
	
	public void runSim(){
		
		for(Gene g:genes.getGenes()){
			List<SNP> snpList = map.get(g);
			Map<String, ExpSample> emap= g.getExpsamples();
			int randInt= (int) Math.round(Math.random());
			for(SNP s:snpList){
				System.out.println(s.toString());
				Map<String, GenoSample> gmap= s.getGenosamples();
				int correct=0;
				int incorrect=0;
				for(String gtexId:gmap.keySet()){
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
				System.out.println(correct+"\t"+incorrect);
				if(counts.get(s.getId())==null){
					//System.out.println("Creating key in counts for: "+s.getId());
					counts.put(s.getId(), 0);
				}
				
				if(passThreshold(correct,incorrect)){
					int temp= counts.get(s.getId());
					temp++;
					counts.put(s.getId(), temp);
					//System.out.println("updating key in counts for: "+s.getId() + '\t' + temp);
				}
				else{
					;
				}

				//System.out.println(r.toString());
			}
			
		}
	}

}