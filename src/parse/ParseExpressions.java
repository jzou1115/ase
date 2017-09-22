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
		Map<String, Integer> aseCall = new HashMap<String,Integer>(); //Map GTEx id of individual to ase call based on binomial test
		List<String> ret = new ArrayList<String>();
		try {
			//iterate over lines in file (SNPs used to call ASE)
			while((line = br.readLine()) != null){
				try{
					String[] tokens = line.split("\\s+");
					String gtexid = tokens[6].trim();
					
					double p = Double.parseDouble(tokens[15]);
					int totalReads = Integer.parseInt(tokens[10]);

					if(!reads.containsKey(gtexid)){
						reads.put(gtexid, 0);
					}
					
					if(totalReads > reads.get(gtexid)){
						reads.put(gtexid, totalReads);
						if(p<.05){
							aseCall.put(gtexid, 1);
						}
						else{
							aseCall.put(gtexid, 0);
						}
					}
					

				} catch (Exception e){
					e.printStackTrace();
					System.exit(1);
				}
			}

			br.close();


			//determine whether each individual has ASE or not
			BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+g.getId()+"_ase.txt")));		
			
			for(String sample: aseCall.keySet()){
				ExpSample expsamp;
				int call = aseCall.get(sample);
				if(call==0){
					expsamp = new ExpSample(sample, 0);
				}
				else{
					expsamp = new ExpSample(sample, 1);
				}
				g.addSample(expsamp);
				
				outfile.write(sample+"\t"+call+"\n");
		
			}
			outfile.close();
		

		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		return ret;
		
	}


}
