package parse;

import genome.Gene;
import genome.SNP;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ParseMap {
	
	Map<Gene, List<SNP>> map;
	
	Gene g;
	List<SNP> snps;
	
	
	public void parseMap(InputStream in){
		map = new HashMap<Gene,List<SNP>>();
		
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
					snps.add(ParseSNP.parseSNP(line));
				}
			}
		} catch (IOException e) {
			System.out.println("no lines");
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void parseMap(InputStream in, String gene){
		
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
					snps.add(ParseSNP.parseSNP(line));
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
	}
	
	public Gene getGene(){
		return g;
	}
	
	public List<SNP> getSNPs(){
		return snps;
	}
	
}
