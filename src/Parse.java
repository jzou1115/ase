import genome.Gene;
import genome.GeneGroup;
import genome.GenomicCoordinate;
import genome.SNP;
import genome.SNPgroup;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import sample.GenoSample;


public class Parse {
	
	public static SNP parseSNP(String line){
		String[] tokens = line.split("\\s+");
//		if(Pattern.matches(SNP_REGEX, tokens[0].trim())){
			SNP s= new SNP(tokens[0].trim(), Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]));
			return s;
//		}
//		System.out.println("does not match snp regex "+tokens[0]);
//		return null;
	}
	public static SNPgroup readSNPGroup(InputStream in){
		List<SNP> snps = new ArrayList<SNP>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(in));
		String line;
		int n=0;
		try {
			while((line = reader.readLine()) != null){
				try{
					SNP s = parseSNP(line);
					if(s!=null){
						snps.add(s);
						n++;
					}
				} catch (Exception e){
					//do nothing
				}
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return new SNPgroup(snps);
	}

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
	public static GeneGroup readGeneGroup(InputStream in){
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
		return new GeneGroup(genes);
	}
	
	public static Map<Gene,List<SNP>> parseMap(InputStream in){
		Map<Gene, List<SNP>> map = new HashMap<Gene,List<SNP>>();
		
		BufferedReader reader = new BufferedReader(new InputStreamReader(in));
		String line;
		Gene g;
		List<SNP> snps= new ArrayList<SNP>();
		try{
			while((line = reader.readLine()) != null){
				if(line.charAt(0) == '>'){
					g = parseGene(line.substring(1,line.length()));
					snps = new ArrayList<SNP>();
				}
				else{
					snps.add(parseSNP(line));
				}
			}
		} catch (IOException e) {
			System.out.println("no lines");
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return map;
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
