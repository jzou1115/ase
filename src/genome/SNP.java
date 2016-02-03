package genome;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
//import java.util.regex.Pattern;
import java.util.Map;

import genome.GenomicCoordinate;
import sample.*;

public class SNP implements Comparable<SNP>{

	String id;
	GenomicCoordinate location;
	
	List<GenoSample> gsamples;
	
	Map<String, GenoSample> map;
	
	String chrom;
	
	public SNP(String i, int c, long l, String ch){
		id=i;
		location = new GenomicCoordinate(c, l);
		gsamples = new ArrayList<GenoSample>();
		map = new HashMap<String, GenoSample>();
		chrom = ch;
	}
	
	public SNP(String i, int c, long l){
		id=i;
		location = new GenomicCoordinate(c, l);
		gsamples = new ArrayList<GenoSample>();
		map = new HashMap<String, GenoSample>();
	}
	
	public SNP(String i){
		id=i;
		gsamples = new ArrayList<GenoSample>();
		map = new HashMap<String, GenoSample>();
	}

	public void addSamples(List<GenoSample> s){
		for(GenoSample g:s){
			gsamples.add(g);
			map.put(g.getID(), g);
		}
	}
	
	public void addSample(GenoSample g){
		gsamples.add(g);
		map.put(g.getID(), g);
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
		if(location != null && chrom!=null){
			return id+"\t"+location.getChromosome()+"\t"+location.getCoord()+"\t"+chrom;
		}
		else if(location!=null){
			return id+"\t"+location.getChromosome()+"\t"+location.getCoord();
		}
		return id;
	}
	
	public String samplesToString(){
		String ret = id;
		
		for(GenoSample g:gsamples){
			ret = ret+"\t"+ g.getHetero();
		}
		
		return ret;
	}
	
	public void sortSamples(){
		Collections.sort(gsamples);
	}

	public GenoSample getSample(String sampleID) {
		if(map.keySet().contains(sampleID)){
			return map.get(sampleID);
		}
		return null;
	}
	
	public String getChromState(){
		return chrom;
	}
	
	public void setChromState(ChromState c){
		chrom = c.getState();
	}

}