package parse;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;

import genome.Gene;
import sample.ExpSample;

public class ParseExpressions {
	
	public static List<String> parseExpressions(InputStream expressions, Gene g, File outdir) throws IOException{
		System.out.println("Reading expressions");
		BufferedReader br = new BufferedReader(new InputStreamReader(expressions));
		String line = br.readLine(); //skip header
		
		Map<String, Integer> reads = new HashMap<String,Integer>(); //Map GTEx id of individual to total number of reads for individual
		Map<String, Integer> left = new HashMap<String,Integer>(); //Map GTEx id of individual to ase call based on binomial test
		List<String> ret = new ArrayList<String>();
		try {
			//iterate over lines in file (SNPs used to call ASE)
			while((line = br.readLine()) != null){
				String[] tokens = line.split("\\s+");
				if(tokens.length==24){
					String gtexid = tokens[6].trim();
					
					if(!tokens[15].equals("NA")){
						double p = Double.parseDouble(tokens[15]);
						int totalReads = Integer.parseInt(tokens[10]);
						int refReads = Integer.parseInt(tokens[8]);
						int altReads = Integer.parseInt(tokens[9]);
						
						//initialize total reads and hap 1 reads to zero
						if(!reads.containsKey(gtexid)){
							reads.put(gtexid, 0);
							left.put(gtexid, 0);
						}
						
						//add total read count
						int newTotalReads = reads.get(gtexid) + totalReads;
						reads.put(gtexid, newTotalReads);
						
						//add reads to left (hap 1)
						String genotype = tokens[18];
						if(isLeft(genotype)){
							int newLeftReads = left.get(gtexid) + altReads;
							left.put(gtexid, newLeftReads);
						}
					
					}
					

				} 
			}

			br.close();


			//determine whether each individual has ASE or not
			//BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+g.getId()+"_ase.txt")));		
		
			for(String sample: reads.keySet()){
				int totalReads = reads.get(sample);
				int leftReads = left.get(sample);
				
				//binomial test for ASE call
				//binomialTest(int numberOfTrials, int numberOfSuccesses, double probability, AlternativeHypothesis alternativeHypothesis)
				BinomialTest bt = new BinomialTest();
				double p = bt.binomialTest(totalReads, leftReads, 0.05, AlternativeHypothesis.TWO_SIDED);
				
				ExpSample expsamp;
				if(p < .05){
					expsamp = new ExpSample(sample, 1);
				}
				else{
					expsamp = new ExpSample(sample, 0);
				}
				
				g.addSample(expsamp);
				ret.add(sample);
				//outfile.write(sample+"\t"+call+"\n");
		
			}
			//outfile.close();
		

		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		return ret;
		
	}

	private static boolean isLeft(String string) {
		String[] tokens = string.split(";");
		String[] tokens2 = tokens[1].split("|");
		if(!isInteger(tokens2[0])){
			System.out.println("Error: Cannot parse genotype "+string);
			System.exit(1);
		}
		if(Integer.parseInt(tokens2[0])==1){
			return true;
		}
		return false;
	}

}
