package run;

import java.util.ArrayList;
import java.util.Map;

import sample.Sample;
import gene.*;
import snp.*;

public class Run {
	float threshold;
	SNPgroup snps;
	GeneGroup genes;
	Map<SNP,Gene> map;
	
	HashMap<String, Result> results;
	
	Run(SNPgroup s, GeneGroup g, Map<SNP,Gene> m, float t){
		snps=s;
		genes=g;
		map=m;
		threshold=t;
	}
	
	public void run(){

	}
}

