package genome;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Pattern;

import genome.GenomicCoordinate;
import sample.*;

public class SNP implements Comparable<SNP>{
	private static final String SNP_REGEX = "^rs.*$";
	
	//parsed from *.map files
	String id;
	GenomicCoordinate location;
	
	List<GenoSample> gsamples;
	
	public SNP(String i, int c, long l){
		id=i;
		location = new GenomicCoordinate(c, l);
		gsamples = new ArrayList<GenoSample>();
	}
	
	public SNP(String i){
		id=i;
		gsamples = new ArrayList<GenoSample>();
	}

	public void addSamples(List<GenoSample> s){
		for(GenoSample g:s){
			gsamples.add(g);
		}
	}
	
	public void addSample(GenoSample g){
		gsamples.add(g);
	}
	
	public int getNumSamples(){
		return gsamples.size();
	}
	public List<GenoSample> getGenosamples(){
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
		return id+"\t"+location.getChromosome()+"\t"+location.getCoord();
	}
	
	public String samplesToString(){
		String ret = id;
		
		for(GenoSample g:gsamples){
			ret = ret+"\t"+ g.getHetero();
		}
		
		return ret;
	}

}