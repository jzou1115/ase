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


public class SNPgroup{
	private final List<SNP> snpg;
	
	public SNPgroup(){
		snpg= new ArrayList<SNP>();
	}
	
	public SNPgroup(List<SNP> snps){
		snpg= snps;
	}

	public List<SNP> getSnps(){
		return snpg;
	}
	
	public int size(){
		return snpg.size();
	}
	
	public boolean contains(String s){
		for(SNP snp:snpg){
			if(snp.getId().equals(s)){
				return true;
			}
		}
		return false;
	}


	public SNP getSNP(String s){
		SNP nogene = new SNP("SNP");
		for(int i=0; i<snpg.size();i++){
			if(snpg.get(i).getId().equals(s)){
				return snpg.get(i);
			}
		}
		return nogene;
	}
	
	public SNP getSNP(int i){
		return snpg.get(i);
	}

	public void sort(){
		Collections.sort(snpg);
	}
}
