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
	
	/*this is the SNP that the simulation will be run with */
	public SNP snp;
	
	int numberPermutations = 1000;
	int significanceThreshold = 1;
	
	public Run(Gene g, List<SNP> s, int t, SNP snp){
		snps=s;
		gene=g;
		Collections.sort(g.esamples);
		threshold=t;
		runType = 2;
		variantIds = new ArrayList<String>();
		this.snp = snp;
	}
	
	public Run(Gene g, List<SNP> s, int t, int rType){
		snps=s;
		gene=g;
		threshold=t;
		runType = rType;
		variantIds = new ArrayList<String>();
	}
	
	ArrayList<ExpSample> getExpSamples(){
		switch (runType) {
			case 1:
				setRandomSNP();
			case 2:
				gene.esamples = assignASE();
			case 3:
				Collections.sort(gene.esamples);
				return gene.esamples;
			case 4:
				Collections.shuffle(gene.esamples);
				return gene.esamples;
			default:
				return gene.esamples;
		}
	}
	private ArrayList<ExpSample> assignASE() {
		ArrayList<ExpSample> assignedASE = new ArrayList<ExpSample>();
		for (GenoSample s : snp.getGenosamples()){
			ExpSample esample = new ExpSample (s.id, s.getHetero());
			assignedASE.add(esample);
		}
		return assignedASE;
	}
	
	private void setRandomSNP(){
		Random random = new Random();
		snp = snps.get(random.nextInt(snps.size()));
	}

	public void run(){
		ArrayList<ExpSample> esamples = getExpSamples();
		
		for(SNP s:snps){
			ArrayList<GenoSample> gsamples= s.getGenosamples();
			int correct=0;
			int incorrect=0;
			for (int i =0; i<gsamples.size(); i++){
				if (gsamples.get(i).isHetero == esamples.get(i).hasASE)
					correct ++;
				else
					incorrect ++;
			}
			if(passThreshold(incorrect, correct)){
				variantIds.add(s.getId());
			}
		}
	}

	private boolean passThreshold(int incorrect, int correct){
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
			Run r = new Run(gene, snps, threshold, 4);
			r.run();
			if (r.variantIds.size() != 0){
				numberBetterPermutations++;
			}
		}
		double pValue = 1.0 * numberBetterPermutations / numberPermutations;
		System.out.print(pValue + ",");

		if (numberBetterPermutations >= significanceThreshold){
			return false;
		}
		return true;
	}
}