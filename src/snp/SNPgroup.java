package snp;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;


public class SNPgroup implements Iterable<SNP>{
	private final List<SNP> snpg;
	
	public SNPgroup(){
		snpg= new ArrayList<SNP>();
	}
	
	public SNPgroup(Collection<SNP> snps){
		snpg= new ArrayList<SNP>();
		for(SNP s : snps){
			snpg.add(s);
		}
	}
	
	@Override
	public Iterator<SNP> iterator() {
		// TODO Auto-generated method stub
		return snpg.iterator();
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
	
	public List<SNP> toList(){
		return new ArrayList<SNP>(snpg);
	}
	
	public int size(){
		return snpg.size();
	}
	
	public SNP getSNP(int index){
		return snpg.get(index);
	}
	
	public void sort(){
		Collections.sort(snpg);
	}
	
	@Override
	public String toString(){
		String s = "";
		for(int i=0; i<snpg.size(); i++){
			SNP p = snpg.get(i);
			s += p.toString() + "\n";
		}
		return s;
	}
}
