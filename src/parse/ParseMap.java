package parse;

import genome.Gene;
import genome.SNP;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ParseMap {
	Gene g;
	//List<String> genes;
	List<SNP> snps;
	
	Map<String, SNP> snpLoc;
	Map<String, SNP> snpmap;
	
	/**
	 * Get list of SNP objects for a gene
	 * @param in Contains data in the map file
	 * @param gene Gene id
	 */
	public void parseMap(InputStream in, String gene){
		snpLoc = new HashMap<String, SNP>();
		snpmap = new HashMap<String, SNP>();
		
		g = new Gene(gene);
		snps = new ArrayList<SNP>();
		
		BufferedReader reader = new BufferedReader(new InputStreamReader(in));
		String line;
		try{
			while((line = reader.readLine()) != null){
				String[] tokens = line.split("\\s+");
				String snpid = tokens[0];
				String[] snpTokens = snpid.split("_");
				int chr = Integer.parseInt(snpTokens[0]);
				long loc = Long.parseLong(snpTokens[1]);
				
				SNP s = new SNP(snpid, chr, loc);
				snps.add(s);
				
				String key = chr+"_"+loc;
				if(!snpLoc.containsKey(key)){
					snpLoc.put(key, s);	
				}
				else{
					System.out.println("snp duplicate: "+key);
					System.exit(1);
				}
				
				if(!snpmap.containsKey(snpid)){
					snpmap.put(snpid, s);	
				}
				else{
					System.out.println("snp duplicate: "+s.getId());
					System.exit(1);
				}
				
			}
		
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	public Gene getGene(){
		return g;
	}
	
	public List<SNP> getSNPs(){
		return snps;
	}

	public Map<String, SNP> getSnpLoc() {
		return snpLoc;
	}
	
	public Map<String, SNP> getSnpMap(){
		return snpmap;
	}

}
