package genome;

import sample.*;

public class Gene implements Comparable<Gene>{
	String id;
	public GenomicRegion region;
	
	ExpSample[] esamples;
	
	Gene(String i, GenomicCoordinate start, GenomicCoordinate end){
		id=i;
		region = new GenomicRegion(start, end);
	}
	
	Gene(String i){
		id= i;
	}
	public ExpSample[] getExpsamples(){
		return esamples;
	}
	
	public GenomicRegion getRegion(){
		return region;
	}
	
	public void addSample(ExpSample[] s){
		esamples=s;
	}
	
	public String getId(){
		return id;
	}
	
	public int getNumSamples(){
		return esamples.length;
	}
	@Override
	public int compareTo(Gene o) {
		return this.region.compareTo(o.region);
	}
	

	public String toString(){
		return id;
	}

	public static Gene parseGene(String line){
		String[] tokens = line.split("\\s+");
		String chromosome = tokens[1];
		int chr = Integer.parseInt(chromosome);
		long s = Long.parseLong(tokens[2]);
		long e = Long.parseLong(tokens[3]);
		GenomicCoordinate start = new GenomicCoordinate(chr, s);
		GenomicCoordinate end = new GenomicCoordinate(chr, e);
		String id = tokens[0];
		return new Gene(id, start, end);
	}
}
