package gene;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;


public class GeneGroup implements Iterable<Gene>{

	private final List<Gene> geneg;
	
	public GeneGroup(){
		geneg= new ArrayList<Gene>();
	}
	
	public GeneGroup(Collection<Gene> genes){
		geneg= new ArrayList<Gene>();
		for(Gene g : genes){
			geneg.add(g);
		}
	}
	
	@Override
	public Iterator<Gene> iterator() {
		// TODO Auto-generated method stub
		return geneg.iterator();
	}
	
	public void sort(){
		Collections.sort(geneg);
	}

	public static GeneGroup readGeneGroup(InputStream in){
		List<Gene> snps = new ArrayList<Gene>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(in));
		String line;
		int n=0;
		try {
			while((line = reader.readLine()) != null){
				try{
					snps.add(Gene.parseGene(line, n));
					n++;
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
	
	public List<Gene> toList(){
		return new ArrayList<Gene>(geneg);
	}
	
	public int size(){
		return geneg.size();
	}
	
	public Gene getGene(int index){
		return geneg.get(index);
	}
	
	@Override
	public String toString(){
		String s = "";
		for(int i=0; i<geneg.size(); i++){
			Gene p = geneg.get(i);
			s += p.toString() + "\n";
		}
		return s;
	}

}
