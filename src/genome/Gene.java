package genome;


import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import sample.*;

public class Gene implements Comparable<Gene>{
	String id;
	public GenomicRegion region;
	List<ExpSample> esamples;
	public GenomicCoordinate tss;
	Map<String, ExpSample> map;
	List<SNP> cis;
	
	public Gene(String i, GenomicCoordinate start, GenomicCoordinate end){
		id=i;
		region = new GenomicRegion(start, end);
		esamples = new ArrayList<ExpSample>();
		map = new HashMap<String, ExpSample>();
		cis = new ArrayList<SNP>();
	}
	
	public Gene(String i, GenomicCoordinate t){
		id=i;
		esamples = new ArrayList<ExpSample>();
		tss = t;
		map = new HashMap<String, ExpSample>();
		cis = new ArrayList<SNP>();
	}
	
	public Gene(String i, GenomicCoordinate start, GenomicCoordinate end, GenomicCoordinate t){
		id=i;
		region = new GenomicRegion(start, end);
		esamples = new ArrayList<ExpSample>();
		map = new HashMap<String, ExpSample>();
		tss = t;
		cis = new ArrayList<SNP>();
	}
	
	
	public Gene(String i){
		id= i;
		map = new HashMap<String, ExpSample>();
		esamples = new ArrayList<ExpSample>();
		cis = new ArrayList<SNP>();
	}
	
	public GenomicRegion copyRegion(){
		return region.copy();
	}
	
	public List<ExpSample> getExpsamples(){
		return esamples;
	}
	
	public GenomicRegion getRegion(){
		return region;
	}
	
	
	public void clearSamples(){
		esamples.clear();
		map.clear();
	}
	
	public void addSample(List<ExpSample> s){
		esamples.addAll(s);
		for(ExpSample samp:s){
			map.put(samp.getID(), samp);
		}
	}
	
	public void addSample(ExpSample e){
		esamples.add(e);
		map.put(e.getID(), e);
	}
	
	public String getId(){
		return id;
	}
	
	public int getNumSamples(){
		return esamples.size();
	}
	@Override
	public int compareTo(Gene o) {
		return this.region.compareTo(o.region);
	}
	

	public String toString(){
		if(region!=null && tss!=null){
			return id+"\t"+region.getChromosome()+"\t"+region.getStart().getCoord()+"\t"+region.getEnd().getCoord()+"\t"+tss.getCoord();
		}
		return id;
	}

	public void sortSamples() {
		Collections.sort(esamples);
	}

	public ExpSample getSample(String s){
		return map.get(s);
	}
	
	public GenomicCoordinate getTSS(){
		return tss;
	}
	
	public void addSNP(SNP s){
		cis.add(s);
	}
	
	public List<SNP> getSNPs(){
		return cis;
	}
	
	public void removeSample(String id){
		for(ExpSample e:esamples){
			if(e.getID().equals(id)){
				esamples.remove(e);
			}
		}
	}
	
	public void removeSample(int ind){
		esamples.remove(ind);
	}

	public void replaceSamples(List<ExpSample> newexpsamples) {
		esamples = newexpsamples;
		
	}
}
