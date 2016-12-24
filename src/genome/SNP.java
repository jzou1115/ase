package genome;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import genome.GenomicCoordinate;
import sample.*;

public class SNP implements Comparable<SNP>{

	String id;
	GenomicCoordinate location;
	
	Map<String, GenoSample> map;
	
	
	public SNP(String i, int c, long l){
		id=i;
		location = new GenomicCoordinate(c, l);
		map = new HashMap<String, GenoSample>();
	}
	
	
	public SNP(String i){
		id=i;
		map = new HashMap<String, GenoSample>();
	}

	public void addSamples(List<GenoSample> s){
		for(GenoSample g:s){
			map.put(g.getID(), g);
		}
	}
	
	public void addSample(GenoSample g){
		map.put(g.getID(), g);
	}
	
	public int getNumSamples(){
		return map.size();
	}
	public List<GenoSample> getGenosamples(){
		List<GenoSample> ret = new ArrayList<GenoSample>();
		for(String key:map.keySet()){
			ret.add(map.get(key));
		}
		return ret;
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
		if(location!=null){
			return id+"\t"+location.getChromosome()+"\t"+location.getCoord();
		}
		return id;
	}

	

	public GenoSample getSample(String sampleID) {
		if(map.keySet().contains(sampleID)){
			return map.get(sampleID);
		}
		return null;
	}

	
	public void removeGenoSample(String s){
		map.remove(s);
	}

}