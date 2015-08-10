package gene;

import genome.GenomicCoordinate;
import genome.GenomicRegion;

import java.util.ArrayList;
import java.util.HashMap;

import sample.*;
import snp.SNP;

public class Gene implements Comparable<Gene>{
	String id;
	GenomicRegion region;
	
	HashMap<String,ExpSample> esamples;
	
	Gene(String i, GenomicCoordinate start, GenomicCoordinate end){
		id=i;
		region = new GenomicRegion(start, end);
		esamples = new HashMap<String, ExpSample>();
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
		String chromosome = tokens[0];
		int chr = Integer.parseInt(chromosome.substring(3,chromosome.length()));
		long s = Long.parseLong(tokens[1]);
		long e = Long.parseLong(tokens[2]);
		GenomicCoordinate start = new GenomicCoordinate(chr, s);
		GenomicCoordinate end = new GenomicCoordinate(chr, e);
		String id = tokens[3];
		return new Gene(id, start, end);
	}
}
