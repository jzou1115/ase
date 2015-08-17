package gene;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;


public class GeneGroup{

	private final HashMap<String, Gene> geneg;
	
	
	public GeneGroup(){
		geneg= new HashMap<String, Gene>();
	}
	
	public GeneGroup(Collection<Gene> genes){
		geneg= new HashMap<String, Gene>();
		for(Gene g : genes){
			geneg.put(g.id, g);
		}
	}

	public List<Gene> getGenes(){
		List<Gene> ret = new ArrayList<Gene>();
		for(String g:geneg.keySet()){
			ret.add(geneg.get(g));
		}
		return ret;
	}
	public static GeneGroup readGeneGroup(InputStream in){
		List<Gene> snps = new ArrayList<Gene>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(in));
		String line;
		try {
			while((line = reader.readLine()) != null){
				try{
					snps.add(Gene.parseGene(line));
				} catch (Exception e){
					//do nothing
				}
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return new GeneGroup(snps);
	}
	
	
	public int size(){
		return geneg.keySet().size();
	}
	
	public Gene getGene(String s){
		return geneg.get(s);
	}
	

}
