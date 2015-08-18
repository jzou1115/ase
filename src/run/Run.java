package run;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import sample.*;
import genome.*;
import statistics.*;

public class Run {
	
	Gene gene;
	List<SNP> snps;
	int num;
	int errors;
	List<ExpSample> expressionData;
	List<GenoSample> genotypeData;
	
	public Run(Gene g, 	List<SNP> sn, int n, int e, List<ExpSample> exp, List<GenoSample> geno){
		gene=g;
		snps=sn;
		num=n;
		errors=e;
		expressionData = exp;
		genotypeData = geno;
	}

	
/**
	public List<SNP> runSim() throws FileNotFoundException{
		Random rand = new Random();
		
		int numSamples = 50;
		List<Integer> randSamples = new ArrayList<Integer>();
		while(randSamples.size()<num){
			int randInt = rand.nextInt(numSamples);
			if(!randSamples.contains(randInt)){
				randSamples.add(randInt);
			}
		}

		//BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("out.txt")));
		
		List<SNP> ret = new ArrayList<SNP>();
		for(SNP s:snps){
			int[] gmap= samples.get(s);
			int correct=0;
			int incorrect=0;
			
			for(int i: randSamples){
				int randInt = rand.nextInt(2);
				
				if(gmap[i] == randInt){
					correct++;
				}
				else{
					incorrect++;
				}
			}
			if(passThreshold(incorrect)){
				//System.out.println(gene.getId()+"\t"+s.getId());
				ret.add(s);
			}
		}
		return ret;
	}
	
	

	public boolean passThreshold(int incorrect){
		if(incorrect>errors){
			return false;
		}
		return true;
	}
	
**/

}