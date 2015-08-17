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

import sample.ExpSample;
import sample.GenoSample;
import sample.Sample;
import gene.*;
import snp.*;
import genome.*;

public class Run {
	int errors;
	List<SNP> snps;
	Gene gene;
	int num;
	Map<SNP, int[]> samples;
	
	public Run(Gene g, 	List<SNP> sn, Map<SNP, int[]> s, int e, int n){
		samples=s;
		gene=g;
		errors=e;
		num=n;
		snps=sn;
	}

	
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
	


}