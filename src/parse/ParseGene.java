package parse;

import genome.Gene;
import genome.GenomicCoordinate;
import genome.SNP;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import sample.GenoSample;

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
		int n=0;
		try {
			while((line = reader.readLine()) != null){
				try{
					genes.add(parseGene(line));
					n++;
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
	
	/**
	public void parseExpressions(InputStream expressions) throws IOException{
		System.out.println("Reading expression data");
		BufferedReader br = new BufferedReader(new InputStreamReader(expressions));
		String line = br.readLine();
		
		String[] gtex = line.split("\\s+");
		for(int j=0; j<gtex.length;j++){
			if(!gtex[j].equals(sampleNames[j])){
				System.out.println("Samples not the same");
				return;
			}
		}
		
		
		try {
			while((line = br.readLine()) != null){
				try{
					String[] tokens = line.split("\\s+");
					
					String gene = tokens[0].trim();
					Gene g = hasASE.getGene(gene);
					
					int[] samp = new int[tokens.length-1];
					for(int i=0; i<tokens.length-1;i++){
						samp[i] = Integer.parseInt(tokens[i+1]);
					}
					esamples.put(g, samp);
					
					//System.out.println(gene+"\t"+g.getNumSamples());
					
				} catch (Exception e){
					//do nothing
				}
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	**/

}
