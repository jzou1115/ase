package run;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import sample.ExpSample;
import sample.GenoSample;
import gene.*;
import snp.*;

public class Run {
	List<SNP> snps;
	Gene gene;

	/* RunType = 1 - use given ASE, 2 - assign ASE based on random snp, 3 - shuffle ASE, 4 - assign based on SNP */
	int runType;
	int threshold;
	
	/* Final counts of variants that work */
	public List<String> variantIds;
	
	/*If runtype = 4, then this is the SNP that the simulation will be run with */
	String snpId;
	
	int numberPermutations = 100;
	int significanceThreshold = 1;
	
	public Run(Gene g, List<SNP> s, int t, String id){
		snps=s;
		gene=g;
		threshold=t;
		runType = 4;
		variantIds = new ArrayList<String>();
		snpId = id;
	}
	
	public Run(Gene g, List<SNP> s, int t, int rType){
		snps=s;
		gene=g;
		threshold=t;
		runType = rType;
		variantIds = new ArrayList<String>();
	}
	
	
	public ExpSample find(ArrayList<ExpSample> samples, String id) {
		for (ExpSample s : samples){
			if (s.id.equalsIgnoreCase(id)){
				return s;
			}
		}
		return null;
	}
	
	ArrayList<ExpSample> getExpSamples( ){
		ArrayList<ExpSample> ASE = gene.esamples;
		switch (runType) {
			case 1:
				return ASE;
			case 2:
				gene.esamples = assignASE(getRandomSNP());
				return gene.esamples;
			case 3:
				return shuffle();
			case 4:
				gene.esamples = assignASE(getSpecificSNP(snpId));
				return gene.esamples;
			default:
				return ASE;
		}
	}

	private ArrayList<ExpSample> assignASE(SNP snp) {
		
		ArrayList<ExpSample> assignedASE = new ArrayList<ExpSample>();
		for (GenoSample s : snp.getGenosamples()){
			ExpSample esample = new ExpSample (s.id, s.getHetero());
			assignedASE.add(esample);
		}
		return assignedASE;
	}
	
	private SNP getRandomSNP(){
		Random random = new Random();
		return snps.get(random.nextInt(snps.size()));
	}
	
	private SNP getSpecificSNP(String snpId){
		for (SNP s : snps){
			if(s.getId().equalsIgnoreCase(snpId)){
				return s;
			}
		}
		System.out.println("In Run.java -- cannot find SNP associated with snpId");
		return null;
	}
	
	
	private ArrayList<ExpSample> shuffle() {
		ArrayList<ExpSample> ASE = gene.getExpsamples();
		
		ArrayList<ExpSample> copyArray = new ArrayList<ExpSample>();
		for (int i=0; i<ASE.size(); i++){
			copyArray.add(new ExpSample(ASE.get(i).id, ASE.get(i).hasASE));
		}
		
		Collections.shuffle(copyArray);
		
		for (int i = 0; i < copyArray.size(); i++){
			copyArray.get(i).id = ASE.get(i).id;
		}
		return copyArray;
	}

	public int run(){
		ArrayList<ExpSample> expSamples = getExpSamples();
		for(SNP s:snps){
			ArrayList<GenoSample> gsamples= s.getGenosamples();
			int correct=0;
			int incorrect=0;
			for (GenoSample g : gsamples) {
				String sampleID = g.getSampleID();
				int isHetero = g.getHetero();
				ExpSample expSample = this.find(expSamples, sampleID);
				if (expSample == null)
					System.out.println("ERROR: Missing data: "+sampleID + "SNP id:" + s.getId());
				else if(isHetero == expSample.getASE())
					correct++;
				else
					incorrect++;
			}
			if(passThreshold(incorrect, correct)){
				variantIds.add(s.getId());
			}
		}
		return variantIds.size();
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
	
	/* Runs Permutation Test */
	public boolean isStatisticallySignificant(){
		int numberBetterPermutations = 0;
		for (int i=0; i<numberPermutations; i++){
			Run r = new Run(gene, snps, threshold, 3);
			r.run();
			if (r.variantIds.size() != 0){
				numberBetterPermutations++;
			}
		}
		double pValue = 1.0 * numberBetterPermutations / numberPermutations;
		System.out.println(pValue);

		if (numberBetterPermutations >= significanceThreshold){
			return false;
		}
		return true;
	}
}