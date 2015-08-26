package genome;


import java.util.List;

import sample.*;

public class Gene implements Comparable<Gene>{
	String id;
	public GenomicRegion region;
	List<ExpSample> esamples;
	
	public Gene(String i, GenomicCoordinate start, GenomicCoordinate end){
		id=i;
		region = new GenomicRegion(start, end);
	}
	
	Gene(String i){
		id= i;
	}
	public List<ExpSample> getExpsamples(){
		return esamples;
	}
	
	public GenomicRegion getRegion(){
		return region;
	}
	
	public void addSample(List<ExpSample> s){
		esamples=s;
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
		return id+"\t"+region.getChromosome()+"\t"+region.getStart().getCoord()+"\t"+region.getEnd().getCoord();
	}

	
}
