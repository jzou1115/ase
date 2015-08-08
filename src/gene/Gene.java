package gene;

import java.util.ArrayList;
import sample.*;
import snp.SNP;

public class Gene implements Comparable<Gene>{
	String id;
	int chrom;
	int location;
	String[] snpIDs;
	
	ExpSample[] esamples;
	
	Gene(String i, int c, int l){
		id=i;
		chrom=c;
		location=l;
	}
	
	@Override
	public int compareTo(Gene o) {
		return this.id.compareTo(o.id);
	}
	
	public ArrayList<String> getGenes(){
		ArrayList<String> ret = new ArrayList<String>();
		return ret;
	}

	public String toString(){
		return id;
	}
	
	/**
	public static Gene parseGene(String line, int n){
		String[] tokens = line.split("\\s+");
		
		return ;
	}
	**/
}
