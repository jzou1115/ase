package gene;


import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

public class GeneGroup{

	private HashMap<String, Gene> geneMap;	
	
	public GeneGroup(){
		geneMap = new HashMap<String, Gene>();
	}
	
	public GeneGroup(Collection<Gene> genes){
		geneMap= new HashMap<String, Gene>();
		for(Gene g : genes){
			geneMap.put(g.geneId, g);
		}
	}

	public List<Gene> getGenes(){
		List<Gene> ret = new ArrayList<Gene>();
		for(String g:geneMap.keySet()){
			ret.add(geneMap.get(g));
		}
		return ret;
	}
	
	public int size(){
		return geneMap.keySet().size();
	}
	
	public Gene getGene(String s){
		return geneMap.get(s);
	}
}
