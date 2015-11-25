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
	
	public Gene(String i, GenomicCoordinate start, GenomicCoordinate end){
		id=i;
		region = new GenomicRegion(start, end);
		esamples = new ArrayList<ExpSample>();
		map = new HashMap<String, ExpSample>();
	}
	
	public Gene(String i, GenomicCoordinate t){
		id=i;
		esamples = new ArrayList<ExpSample>();
		tss = t;
		map = new HashMap<String, ExpSample>();
	}
	
	public Gene(String i, GenomicCoordinate start, GenomicCoordinate end, GenomicCoordinate t){
		id=i;
		region = new GenomicRegion(start, end);
		esamples = new ArrayList<ExpSample>();
		map = new HashMap<String, ExpSample>();
		tss = t;
	}
	
	
	Gene(String i){
		id= i;
		map = new HashMap<String, ExpSample>();
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
		return id+"\t"+region.getChromosome()+"\t"+region.getStart().getCoord()+"\t"+region.getEnd().getCoord()+"\t"+tss;
	}

	public void sortSamples() {
		Collections.sort(esamples);
	}

	public ExpSample getSample(String s){
		return map.get(s);
	}
}
