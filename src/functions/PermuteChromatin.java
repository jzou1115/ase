package functions;
import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import genome.ChromState;
import genome.Gene;
import genome.GenomicCoordinate;
import genome.GenomicRegion;
import genome.SNP;
import parse.ParseChromState;

import org.apache.commons.math3.distribution.BinomialDistribution; //binomial distribution cdf
import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;

public class PermuteChromatin {
	
	File permFile;
	List<ChromState> chromatin;
	Map<Gene, File> chromatinFiles;
	List<SNP> snps;
	List<SNP> variants;
	List<Gene> genes;
	
	GenomicRegion chromosome;
	int chromNum;
	long min;
	long max;

	File outdir;
	String filename;
	
	/**
	public PermuteChromatin(File chrom, List<SNP> snp, List<SNP> variant, Gene gene, File out, String f) throws FileNotFoundException{
		snps = snp;
		variants = variant;
		
		chromNum = snps.get(0).getLocation().getChromosome();
		chromFile = new File(chrom,"chr"+chromNum);
		chromatin = ParseChromState.parseChromState(new BufferedInputStream( new FileInputStream(chrom+File.separator+"chr"+chromNum+File.separator+gene.getId()+"_chromatin.txt")));
		
		min = chromatin.get(0).getRegion().getStart().getCoord();
		max = chromatin.get(chromatin.size()-1).getRegion().getEnd().getCoord();
		
		chromosome = new GenomicRegion(new GenomicCoordinate(chromNum,min), new GenomicCoordinate(chromNum,max));
	
		outdir=out;
		filename=f;
	}
	**/
	
	//allows list of genes
	public PermuteChromatin(File chrom, List<SNP> snp, List<SNP> variant, List<Gene> gene, File out, String f) throws FileNotFoundException{
		snps = snp;
		variants = variant;
		genes = gene;
		
	//	chromNum = snps.get(0).getLocation().getChromosome();
	//	chromFile = new File(chrom,"chr"+chromNum);
		chromatinFiles = mapFiles(chrom);
		permFile=new File(chrom+File.separator+"perm");
//		chromatin = ParseChromState.parseChromState(new BufferedInputStream( new FileInputStream(chrom+File.separator+"chr"+chromNum+File.separator+gene.getId()+"_chromatin.txt")));
		
		//min = chromatin.get(0).getRegion().getStart().getCoord();
		//max = chromatin.get(chromatin.size()-1).getRegion().getEnd().getCoord();
		
//		chromosome = new GenomicRegion(new GenomicCoordinate(chromNum,min), new GenomicCoordinate(chromNum,max));
	
		outdir=out;
		filename=f;
	}
	
	public HashMap<Gene,File> mapFiles(File chrom){
		HashMap<Gene,File> ret = new HashMap<Gene,File>();
		for(Gene g:genes){
			String id = g.getId();
			String chr = g.getSNPs().get(0).getId().split("_")[0];
			File file = new File(chrom+File.separator+"chr"+chr+File.separator+id+"_chromatin.txt");
			if(file.exists()){
				ret.put(g, file);
			}
			else{
				System.out.println(id+" chromatin file missing");
				System.exit(1);
			}
		}
		return ret;
	}
	/**
	public void filtersnps(Gene g){
		GenomicRegion proximal = g.copyRegion().expand(250000);
		List<SNP> filtered = new ArrayList<SNP>();
		for(SNP s:snps){
			if(proximal.contains(s.getLocation())){
				filtered.add(s);
			}
		}
		snps = filtered;
	}

	public void filterchrom(Gene g){
		GenomicRegion proximal = g.copyRegion().expand(250000);
		List<ChromState> filtered = new ArrayList<ChromState>();
		for(ChromState c: chromatin){
			if(proximal.overlaps(c.getRegion())){
				filtered.add(c);
			}
		}
		System.out.println("Chromatin size: "+chromatin.size());
		chromatin = filtered;
		System.out.println("Filtered chromatin size: "+chromatin.size());
		System.out.println(filtered.size());
	}
	**/
	
	public void permute(int n) throws IOException{
	//	BufferedWriter file = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+"permutation"+n)));
		List<ChromState> chromatin2 = new ArrayList<ChromState>();
		for(int i=0; i<chromatin.size();i++){
			int ind;
			if((i+n)<chromatin.size()){
				ind = i+n;
			}
			else{
				int offset = (i+n)-chromatin.size();
				ind = offset;
			}
			
			ChromState c = chromatin.get(i);
			ChromState c2 = chromatin.get(ind);
			chromatin2.add(new ChromState(c2.getState(), c.copyRegion()));
			
		}
		
		if(chromatin2.size()!=chromatin.size()){
			System.out.print("error in chrom perm");
		}
		/**
		for(int i=0; i<chromatin.size(); i++){
			file.write(chromatin.get(i).toString()+"\n");
			file.write(chromatin2.get(i).toString()+"\n");
		}
		file.close();
		**/
		chromatin = chromatin2;
		Collections.sort(chromatin);
		
	}
	

	
	public Map<String, Double> mapVar(){
		Map<String, Double> ret = new HashMap<String, Double>();
		Map<SNP,ChromState> map = AssignChromState.assignStateSNP(variants, chromatin);

		for(SNP s:variants){
			if(map.containsKey(s)){
				String state = map.get(s).getState();
				if(!ret.containsKey(state)){
					ret.put(state, 0.0);
				}

				Double val = ret.get(state);
				val = val+1;
				ret.put(state, val);

			}
		}
		
		
		int numVar = variants.size();
		for(String chromState:ret.keySet()){
			Double prop = ret.get(chromState)*1.0/numVar;
			ret.put(chromState, prop);
		}

		return ret;
		
	}
	
	public Map<String, Double> mapVar(List<ChromState> chrom){
		Map<String, Double> ret = new HashMap<String, Double>();
		Map<SNP,ChromState> map = AssignChromState.assignStateSNP(variants, chrom);

		for(SNP s:variants){
			if(map.containsKey(s)){
				String state = map.get(s).getState();
				if(!ret.containsKey(state)){
					ret.put(state, 0.0);
				}

				Double val = ret.get(state);
				val = val+1;
				ret.put(state, val);

			}
		}
		

		int numVar = variants.size();
		for(String chromState:ret.keySet()){
			Double prop = ret.get(chromState)*1.0/numVar;
			ret.put(chromState, prop);
		}

		return ret;
		
	}
	
	public Map<String, Double> mapVar(List<SNP> var, List<ChromState> chrom){
		Map<String, Double> ret = new HashMap<String, Double>();
		Map<SNP,ChromState> map = AssignChromState.assignStateSNP(var, chrom);

		for(SNP s:var){
			if(map.containsKey(s)){
				String state = map.get(s).getState();
				if(!ret.containsKey(state)){
					ret.put(state, 0.0);
				}

				Double val = ret.get(state);
				val = val+1;
				ret.put(state, val);

			}
		}

		return ret;
		
	}
	
	public Map<String, Integer> countVar(){
		Map<String, Integer> ret = new HashMap<String, Integer>();
		Map<SNP,ChromState> map = AssignChromState.assignStateSNP(variants, chromatin);

		for(SNP s:variants){
			if(map.containsKey(s)){
				String state = map.get(s).getState();
				if(!ret.containsKey(state)){
					ret.put(state, 0);
				}

				int val = ret.get(state);
				val = val+1;
				ret.put(state, val);

			}
		}
		return ret;
		
	}
	
	public Map<String, Integer> countVar(Gene g) throws FileNotFoundException{
		List<SNP> var = g.getSNPs();
		List<ChromState> chrom = ParseChromState.parseChromState(new BufferedInputStream( new FileInputStream(chromatinFiles.get(g))));
		
		Map<String, Integer> ret = new HashMap<String, Integer>();
		Map<SNP,ChromState> map = AssignChromState.assignStateSNP(var, chrom);

		for(SNP s:var){
			if(map.containsKey(s)){
				String state = map.get(s).getState();
				if(!ret.containsKey(state)){
					ret.put(state, 0);
				}

				int val = ret.get(state);
				val = val+1;
				ret.put(state, val);

			}
		}
		return ret;
		
	}

	public double calculatePValue(int trials, int variants, double p){
		BinomialTest bin = new BinomialTest();
		double pval = bin.binomialTest(trials,variants,p, AlternativeHypothesis.GREATER_THAN);
		return pval;
	}

	public int getVarSize() {
		return variants.size();
	}

	public int getSnpsSize() {
		return snps.size();
	}

	public void testEnrichment(int p) throws IOException{

		Map<String,Integer> realvar = countVar();
		
		Map<String,Double> allPerm = new HashMap<String, Double>();
		
		//circular permutations to get null mean
		Random rand = new Random();
		for(int i=0; i<p; i++){
			if(i%1000==0){
				System.out.println(i);
			}
			permute(rand.nextInt(chromatin.size()/2));
			Map<String, Double> perm = mapVar();
			for(String chromState:perm.keySet()){
				if(! allPerm.containsKey(chromState)){
					allPerm.put(chromState, 0.0);
				}
				double numVar = perm.get(chromState) + allPerm.get(chromState);
				allPerm.put(chromState, numVar);
			}
		}
		//Calculate averages of each chromatin state
		for(String chromState:allPerm.keySet()){
			double avg = allPerm.get(chromState)*1.0/p;
			allPerm.put(chromState, avg);
		}
		
		//calculate p-values of chromatin state enrichment and write to file
		BufferedWriter file = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+filename+"_enrichment.txt")));
		for(String chromState:allPerm.keySet()){
			if(realvar.containsKey(chromState)){
				double pval = calculatePValue(variants.size(), realvar.get(chromState), allPerm.get(chromState));
				file.write(chromState+"\t"+ allPerm.get(chromState)*variants.size()+"\t"+realvar.get(chromState)+"\t"+pval+"\n");
			}
			else{
				double pval = calculatePValue(variants.size(), 0, allPerm.get(chromState));
				file.write(chromState+"\t"+ allPerm.get(chromState)*variants.size()+"\t"+0.0+"\t"+pval+"\n");
			}
		}
		file.close();
	}
	
	public void testEnrichmentGene(int p) throws IOException{
		
		//Real data
		Map<String,Integer> realvar = new HashMap<String,Integer>();
		for(Gene g:genes){
			Map<String,Integer> gvar = countVar(g);
			for(String key: gvar.keySet()){
				if(!realvar.keySet().contains(key)){
					realvar.put(key, 0);
				}
				int num = gvar.get(key) + realvar.get(key);
				realvar.put(key, num);
			}
		}
		System.out.println(realvar.toString());
		
		//Permuted data
		Map<String,Double> allPerm = new HashMap<String, Double>();
		
		File[] listOfFiles = permFile.listFiles();
		if(listOfFiles.length<p){
			System.out.println("Not enough gene chromatin state files");
			System.exit(1);
		}
		
		Random rand = new Random(13);
		int totalSNPs = 0;
		for(int i=0; i<p; i++){
			int fileNum = rand.nextInt(listOfFiles.length);
			File chromF = listOfFiles[fileNum];
			List<ChromState> chrom = ParseChromState.parseChromState(new BufferedInputStream( new FileInputStream(chromF)));
			int chr1 = chrom.get(0).getRegion().getStart().getChromosome();
			long start1 = chrom.get(0).getRegion().getStart().getCoord();
			
			for(Gene g:genes){
				totalSNPs = totalSNPs + g.getSNPs().size();
				File chromF2 = chromatinFiles.get(g);
				List<ChromState> chrom2 = ParseChromState.parseChromState(new BufferedInputStream( new FileInputStream(chromF2)));
				long start2 = chrom2.get(0).getRegion().getStart().getCoord();
				List<SNP> newvar = alignSNPs(g.getSNPs(), start2, chr1, start1);
				Map<String, Double> perm = mapVar(newvar, chrom);
				
				for(String chromState:perm.keySet()){
					if(realvar.containsKey(chromState)){
						if(realvar.get(chromState) <= perm.get(chromState)){
							if(! allPerm.containsKey(chromState)){
								allPerm.put(chromState, 0.0);
							}
							double temp = allPerm.get(chromState) +1;
							allPerm.put(chromState, temp);
						}
					}
					else{
						realvar.put(chromState, 0);
						if(realvar.get(chromState) <= perm.get(chromState)){
							if(! allPerm.containsKey(chromState)){
								allPerm.put(chromState, 0.0);
							}
							double temp = allPerm.get(chromState) +1;
							allPerm.put(chromState, temp);
						}
						
					}
				}
			}
		

		}
		/**
		
		//Calculate averages of each chromatin state
		for(String chromState:allPerm.keySet()){
			double avg = allPerm.get(chromState)*1.0/(p*totalSNPs);
			allPerm.put(chromState, avg);
		}
		System.out.println("Avg perm: "+ allPerm.toString());
		**/
		for(String chromState:allPerm.keySet()){
			double avg = allPerm.get(chromState)*1.0/(p);
			allPerm.put(chromState, avg);
		}
		
		//calculate p-values of chromatin state enrichment and write to file
		BufferedWriter file = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+filename+"_enrichment.txt")));
		for(String chromState:allPerm.keySet()){
			if(realvar.containsKey(chromState)){
				double pval = allPerm.get(chromState);
				//double pval = calculatePValue(variants.size(), realvar.get(chromState), allPerm.get(chromState));
				file.write(chromState+"\t"+realvar.get(chromState)+"\t"+pval+"\n");
			}
			else{
				double pval = allPerm.get(chromState);
				//double pval = calculatePValue(variants.size(), 0, allPerm.get(chromState));
				file.write(chromState+"\t"+0.0+"\t"+pval+"\n");
			}
		}
		file.close();
	}

	private List<SNP> alignSNPs(List<SNP> list, long start2, int chr, long start) {
		List<SNP> ret = new ArrayList<SNP>();
		long diff = start - start2;
		for(SNP s: list){
			SNP newSnp = new SNP(s.getId(), chr, s.getLocation().getCoord()+diff);
			ret.add(newSnp);
		}
		return ret;
	}

	private List<ChromState> alignChromStates(List<ChromState> permchrom, long start) {
		List<ChromState> ret = new ArrayList<ChromState>();
		long diff = start - permchrom.get(0).getRegion().getStart().getCoord();
		for(ChromState c: permchrom){
			GenomicCoordinate newstart = new GenomicCoordinate(c.getRegion().getChromosome(), c.getRegion().getStart().getCoord()+diff);
			GenomicCoordinate newend = new GenomicCoordinate(c.getRegion().getChromosome(), c.getRegion().getEnd().getCoord()+diff);
			ret.add(new ChromState(c.getState(), new GenomicRegion(newstart, newend)));
		}
		return ret;
	}
	
}
