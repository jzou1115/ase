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

	/* RunType = 1 - use given ASE, 2 - assign ASE based on random snp, 3 - used shuffle ASE */
	int runType;
	int threshold;
	/* Final counts of variants that work */
	public List<String> variantIds;

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
				gene.esamples = assign();
				return gene.esamples;
			case 3:
				return shuffle();
			default:
				return ASE;
		}
	}
	
	/*
	 * Chooses SNP randomly. Changes gene.esamples.
	 */
	private ArrayList<ExpSample> assign() {
		Random random = new Random();
		SNP snp = snps.get(random.nextInt(snps.size()));
		
		System.out.println("Assigning ASE based on SNP: " + snp.getId());
		
		ArrayList<ExpSample> assignedASE = new ArrayList<ExpSample>();
		for (GenoSample s : snp.getGenosamples()){
			ExpSample esample = new ExpSample (s.id, s.getHetero());
			assignedASE.add(esample);
		}
		return assignedASE;
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
}