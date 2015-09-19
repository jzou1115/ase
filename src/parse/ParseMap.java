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
	
	Map<Gene, List<SNP>> map;
	
	Gene g;
	List<SNP> snps;
	
	Map<String, SNP> snpLoc;
	
	
	public void parseMap(InputStream in){
		map = new HashMap<Gene, List<SNP>>();
		snps = new ArrayList<SNP>();
		snpLoc = new HashMap<String, SNP>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(in));
		String line;
		Gene g=null;
		List<SNP> snps= new ArrayList<SNP>();
		try{
			while((line = reader.readLine()) != null){
				if(line.charAt(0) == '>'){
					if(g!=null){
						map.put(g, snps);
						
						g = null;
						snps.clear();
						
						g = ParseGene.parseGene(line.substring(1,line.length()));
						snps = new ArrayList<SNP>();
					}
					else{
						g = ParseGene.parseGene(line.substring(1,line.length()));
						snps = new ArrayList<SNP>();	
					}
				}
				else{
					SNP s = ParseSNP.parseSNP(line);
					snps.add(s);
					
					snpLoc.put(s.getLocation().getChromosome()+"_"+s.getLocation().getCoord(), s);
					System.out.println(s.getLocation().getChromosome()+"_"+s.getLocation().getCoord());
					
					
				}
			}
		} catch (IOException e) {
			System.out.println("no lines");
			e.printStackTrace();
		}
	}
	
	public Map<String,SNP> parseMap(InputStream in, String gene){
		map = new HashMap<Gene, List<SNP>>();
		snps = new ArrayList<SNP>();
		snpLoc = new HashMap<String, SNP>();
		
		BufferedReader reader = new BufferedReader(new InputStreamReader(in));
		String line;
		try{
			while((line = reader.readLine()) != null){
				if(line.contains(gene)){
					break;
				}
			}
			if(line.charAt(0) == '>'){
				g = ParseGene.parseGene(line.substring(1,line.length()));
				snps = new ArrayList<SNP>();
			}
			while((line = reader.readLine()) != null){
				if(line.charAt(0) != '>'){
					SNP s = ParseSNP.parseSNP(line);
					snps.add(s);
					String key = s.getLocation().getChromosome()+"_"+s.getLocation().getCoord();
					if(!snpLoc.containsKey(key)){
						snpLoc.put(key, s);	
					}
					else{
						System.out.println("snp duplicate: "+key);
					}
					
				}
				else{
					break;
				}
			}

				
		} catch (IOException e) {
			System.out.println("no lines");
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return snpLoc;
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

	public void writegene() throws IOException {
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("genesnps.txt")));

		outfile.write(g.getId()+"\n");
		
		for (SNP s:snps){
			outfile.write(s.toString()+"\n");
		}
		
		outfile.close();
	}
	
}
