package parse;

import genome.Gene;
import genome.GenomicCoordinate;
import genome.SNP;
import sample.ExpSample;
import sample.GenoSample;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class ParseGene {

	public static Gene parseGene(String line){
		String[] tokens = line.split("\t");
		String id = tokens[0];
		int chr = Integer.parseInt(tokens[1]);
		long s = Long.parseLong(tokens[2]);
		long e = Long.parseLong(tokens[3]);
		GenomicCoordinate start = new GenomicCoordinate(chr, s);
		GenomicCoordinate end = new GenomicCoordinate(chr, e);
		return new Gene(id, start, end);
	}
	public static List<Gene> readGeneGroup(InputStream in){
		List<Gene> genes = new ArrayList<Gene>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(in));
		String line;
		try {
			while((line = reader.readLine()) != null){
				try{
					genes.add(parseGene(line));
				} catch (Exception e){
					//do nothing
					//System.out.println("cannot add Gene");
				}
			}
			reader.close();
		} catch (IOException e) {
			System.out.println("no lines");
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return genes;
	}
	
	public static List<ExpSample> parseExpressions(InputStream expressions, Gene g, File outdir) throws IOException{
		System.out.println("Reading expressions");
		BufferedReader br = new BufferedReader(new InputStreamReader(expressions));
		String line;
		
		Map<String, Integer> reads = new HashMap<String,Integer>();
		Map<String, Integer> ref = new HashMap<String,Integer>();
		try {
			while((line = br.readLine()) != null){
				try{
					String[] tokens = line.split("\\s+");
					String gtexid = tokens[6].trim();
					int refAllele = Integer.parseInt(tokens[8]);
					int totalReads = Integer.parseInt(tokens[10]);
				//	if(sampleNames.contains(gtexid)){
						if(!reads.containsKey(gtexid)){
							reads.put(gtexid, 0);
							ref.put(gtexid, 0);
						}
						int refKey = ref.get(gtexid) + refAllele;
						int readsKey = reads.get(gtexid) + totalReads;
						ref.put(gtexid, refKey);
						reads.put(gtexid, readsKey);
				//	}
					
				} catch (Exception e){
					//do nothing
				}
			}

			br.close();
			System.out.println("reads size: "+reads.size());
			System.out.println("ref size: "+ref.size());
			BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+g.getId()+"_ase.txt")));

			
			for(String sample: reads.keySet()){
				
				ExpSample expsamp;
				double allelicRatio = 1.0*ref.get(sample) / reads.get(sample);
				int hasASE;
				if(allelicRatio>0.65){
					expsamp = new ExpSample(sample, 1);
					hasASE=1;
				}
				else if(allelicRatio<0.35){
					expsamp = new ExpSample(sample, 1);
					hasASE=1;
				}
				else{
					expsamp = new ExpSample(sample, 0);
					hasASE=0;
				}

				g.addSample(expsamp);
				
				outfile.write(sample+"\t"+allelicRatio+"\t"+hasASE+"\n");
		
			}
			g.sortSamples();
			outfile.close();
		

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return g.getExpsamples();
		
	}
	

}
