package genome;
import java.util.HashMap;

import genome.GenomicCoordinate;
import sample.*;

public class SNP implements Comparable<SNP>{
	//parsed from *.map files
	String id;
	GenomicCoordinate location;
	int num;
	
	GenoSample[] gsamples;
	
	SNP(String i, int c, long l, int n){
		id=i;
		location = new GenomicCoordinate(c, l);
		num= n;
	}
	
	SNP(String i){
		id=i;
	}

	public void addSamples(GenoSample[] s){
		gsamples = s;
	}
	
	public int getNumSamples(){
		return gsamples.length;
	}
	public GenoSample[] getGenosamples(){
		return gsamples;
	}
	public String getId(){
		return id;
	}
	public GenomicCoordinate getLocation(){
		return location;
	}
	
	@Override
	public int compareTo(SNP other) {
		int chrComp = this.location.getChromosome() - other.location.getChromosome();
		if(chrComp != 0) return chrComp;
		if(this.location.getCoord() > other.location.getCoord()){
			return 1;
		}
		if(this.location.getCoord() < other.location.getCoord()){
			return -1;
		}
		return 0;
	}

	public String toString(){
		return id;
	}
	
	public static SNP parseSNP(String line, int n){
		String[] tokens = line.split("\\s+");

		SNP s= new SNP(tokens[0], Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]), n);
		
		return s;
	}
}