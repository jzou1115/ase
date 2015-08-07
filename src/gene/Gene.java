package gene;

import java.util.ArrayList;
import sample.*;

public class Gene implements Comparable<Gene>{
	String id;
	int chrom;
	int location;
	String[] snpIDs;
	
	ExpSample[] esamples;
	
	@Override
	public int compareTo(Gene o) {
		// TODO Auto-generated method stub
		return 0;
	}
	
	public ArrayList<String> getGenes(){
		ArrayList<String> ret = new ArrayList<String>();
		return ret;
	}

}
