package parse;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import genome.Gene;
import genome.SNP;

public class ParseCausalVariants {
	
	public static SNP parseSNP(String line){
		String[] snpTokens = line.split("_");
		SNP s= new SNP(line, Integer.parseInt(snpTokens[0]), Integer.parseInt(snpTokens[1]));
		return s;
	}
	public static List<Gene> readVariantGroup(InputStream in){
		List<Gene> genes = new ArrayList<Gene>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(in));
		String line;
		try {
			while((line = reader.readLine()) != null){
				try{
					String[] tokens = line.split("\\s+");
					SNP s = parseSNP(tokens[1]);
					if(s!=null){
						String geneid = tokens[0];
						
						boolean addGene = true;
		
						for(Gene g:genes){
							if(g.getId().equals(geneid)){
								g.addSNP(s);
								addGene=false;
							}
						}
						
						if(addGene){
							Gene g = new Gene(tokens[0]);
							g.addSNP(s);
							genes.add(g);
						}
	
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
		return genes;
	}
}
