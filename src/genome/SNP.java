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
	int num;
	
	List<GenoSample> gsamples;
	
	SNP(String i, int c, long l, int n){
		id=i;
		location = new GenomicCoordinate(c, l);
		num= n;
		gsamples = new ArrayList<GenoSample>();
	}
	
	SNP(String i){
		id=i;
	}

	public void addSamples(List<GenoSample> s){
		for(GenoSample g:s){
			gsamples.add(g);
		}
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
		return id;
	}
	
	public static SNP parseSNP(String line, int n){
		String[] tokens = line.split("\\s+");
	//	if(Pattern.matches(SNP_REGEX, tokens[0].trim())){
			SNP s= new SNP(tokens[0].trim(), Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]), n);
			return s;
	//	}
	//	System.out.println("does not match snp regex "+tokens[0]);
	//	return null;
	}
}