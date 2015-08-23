package snp;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;


public class SNPgroup{
	private final HashMap<String, SNP> snpg;
	
	public SNPgroup(Collection<SNP> snps){
		snpg= new HashMap<String, SNP>();
		for(SNP s : snps){
			snpg.put(s.snpId,s);
		}
	}

	public List<SNP> getSnps(){
		List<SNP> ret = new ArrayList<SNP>();
		for(String s:snpg.keySet()){
			ret.add(snpg.get(s));
		}
		return ret;
	}
	
	public int size(){
		return snpg.keySet().size();
	}
	
	public SNP getSNP(String s){
		return snpg.get(s);
	}
	
}
