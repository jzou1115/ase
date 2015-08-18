package genome;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;


public class GeneGroup{

	private final List<Gene> geneg;
	
	
	public GeneGroup(){
		geneg= new ArrayList<Gene>();
	}
	
	public GeneGroup(List<Gene> genes){
		geneg= genes;
	}

	public List<Gene> getGenes(){
		return geneg;
	}
	public static GeneGroup readGeneGroup(InputStream in){
		List<Gene> genes = new ArrayList<Gene>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(in));
		String line;
		int n=0;
		try {
			while((line = reader.readLine()) != null){
				try{
					genes.add(Gene.parseGene(line));
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
	
	
	public int size(){
		return geneg.size();
	}
	
	public boolean contains(String s){
		if(geneg.contains(s)){
			return true;
		}
		return false;
	}
	
	public Gene getGene(String s){
		Gene nogene = new Gene("nogene");
		for(int i=0; i<geneg.size();i++){
			if(geneg.get(i).getId().equals(s)){
				return geneg.get(i);
			}
		}
		return nogene;
	}
	
	public Gene getGene(int i){
		return geneg.get(i);
	}

	public void sort(){
		Collections.sort(geneg);
	}

}
