package gene;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class GeneGroup{

	private HashMap<String, Gene> geneMap;	
	
	public GeneGroup(){
		geneMap = new HashMap<String, Gene>();
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

	public void add(String geneId, Gene g) {
		geneMap.put(geneId, g);
	}
	
	public void remove(String geneId){
		geneMap.remove(geneId);
	}
}
