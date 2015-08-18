package gene;

import genome.GenomicCoordinate;
import genome.GenomicRegion;

import java.util.ArrayList;

import sample.*;

public class Gene implements Comparable<Gene>{
	String geneId;
	public GenomicRegion region;
	
	public ArrayList<ExpSample> esamples;
	
	Gene(String id, GenomicCoordinate start, GenomicCoordinate end){
		geneId=id;
		region = new GenomicRegion(start, end);
		esamples = new ArrayList<ExpSample>();
	}
	
	public ArrayList<ExpSample> getExpsamples(){
		return esamples;
	}
	
	public void addSample(String s, ExpSample e){
		esamples.add(e);
	}
	
	public String getId(){
		return geneId;
	}
	
	public int getNumSamples(){
		return esamples.size();
	}
	
	@Override
	public int compareTo(Gene o) {
		return this.region.compareTo(o.region);
	}
	

	public static Gene parseGene(String line){
		String[] tokens = line.split("\\s" , ',');
		String id = tokens[3];
		int chr = Integer.parseInt(tokens[0]);
		long s = Long.parseLong(tokens[1]);
		long e = Long.parseLong(tokens[2]);
		
		GenomicCoordinate start = new GenomicCoordinate(chr, s);
		GenomicCoordinate end = new GenomicCoordinate(chr, e);
		
		return new Gene(id, start, end);
	}
}
