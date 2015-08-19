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
	int threshold;
	List<SNP> snps;
	Gene gene;
	int variants;
	int runType;
	
	public Run(Gene g, List<SNP> s, int t, int rType){
		snps=s;
		gene=g;
		threshold=t;
		variants = 0;
		runType = rType;
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
		ArrayList<ExpSample> ASE = gene.getExpsamples();
		switch (runType) {
			case 1:
				return ASE;
			case 2:
				return shuffle(ASE);
			case 3:
				return assign();
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
		
		//TODO think about assigning this
		gene.esamples = assignedASE;
		
		return assignedASE;
	}

	private ArrayList<ExpSample> shuffle(ArrayList<ExpSample> ASE) {
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
				if ( expSample != null ){
					if(isHetero == expSample.getASE()){
						correct++;
					}
					else{
						incorrect++;
					}
				}
				else {
					System.out.println("ERROR: Missing data: "+sampleID + "SNP id:" + s.getId());
				}
			}
			if(passThreshold(incorrect, correct)){
				variants++;
			}
		}
		//System.out.println("There are " + variants + " variants in this permutation.");
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