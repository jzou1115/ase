package snp;

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
	private final HashMap<String, SNP> snpg;
	
	public SNPgroup(){
		snpg= new HashMap<String,SNP>();
	}
	
	public SNPgroup(Collection<SNP> snps){
		snpg= new HashMap<String, SNP>();
		for(SNP s : snps){
			snpg.put(s.id,s);
		}
	}
	

	public static SNPgroup readSNPGroup(InputStream in){
		List<SNP> snps = new ArrayList<SNP>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(in));
		String line;
		int n=0;
		try {
			while((line = reader.readLine()) != null){
				try{
					snps.add(SNP.parseSNP(line, n));
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
		return new SNPgroup(snps);
	}
	
	public int size(){
		return snpg.keySet().size();
	}
	
	public SNP getSNP(String s){
		return snpg.get(s);
	}
	
}
