package genome;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import sample.GenoSample;

public class DiSNP {
	String id;
	SNP snp1;
	SNP snp2;

	List<GenoSample> gsamples;
	Map<String, GenoSample> map;
	
	public DiSNP(SNP a, SNP b){
		int snpComp = a.compareTo(b);
		if(snpComp<0){
			snp1 = a;
			snp2 = b;
		}
		else if(snpComp>0){
			snp1 = b;
			snp2 = a;
		}
		else{
			snp1 = a;
			snp2 = b;	
		}
		id = snp1.getId()+"_"+snp2.getId();
		combineGenoSamples();
	}
	
	public List<GenoSample> getGenosamples(){
		return gsamples;
	}
	public String getId(){
		return id;
	}
	public SNP getSNP1(){
		return snp1;
	}
	
	public SNP getSNP2(){
		return snp2;
	}
	
	public long getDist(){
		return snp1.getLocation().distance(snp2.getLocation());
	}
	
	public String toString(){
		return snp1.toString()+"\t"+ snp2.toString()+"\t"+ this.getDist();
	}
	
	public int compareTo(DiSNP other) {
		if(this.getSNP1().compareTo(other.getSNP1())==0 && this.getSNP2().compareTo(other.getSNP2())==0){
			return 0;
		}
		
		if(this.getSNP1().compareTo(other.getSNP1())<0){
			return -1;
		}

		else if(this.getSNP1().compareTo(other.getSNP1())==0 && this.getSNP2().compareTo(other.getSNP2())<0){
			return -1;
		}
		
		return 1;
	}
	
	public void combineGenoSamples(){
		gsamples = new ArrayList<GenoSample>();
		map = new HashMap<String, GenoSample>();
		
		List<GenoSample> g1 = snp1.getGenosamples();
		List<GenoSample> g2 = snp2.getGenosamples();
		if(g1.size() != g2.size()){
			return;
		}
		for(int i=0; i<g1.size(); i++){
			if(g1.get(i).getID().equals(g2.get(i).getID())){
				if(g1.get(i).getHetero()==1 || g2.get(i).getHetero()==1){
					GenoSample geno = new GenoSample(g1.get(i).getID(),1);
					gsamples.add(geno);
					map.put(geno.getID(),geno);
				}
				else{
					GenoSample geno = new GenoSample(g1.get(i).getID(),0);
					gsamples.add(geno);
					map.put(geno.getID(),geno);
				}	
			}
			else{
				System.out.println("Gtex ids not matching");
			}
		}
	}
}
