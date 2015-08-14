package snp;
import java.util.ArrayList;

import com.sun.tools.javac.util.List;

import genome.GenomicCoordinate;
import sample.*;

public class SNP implements Comparable<SNP>{
	//parsed from *.map files
	String SNPid;
	GenomicCoordinate location;
	int num;
	
	ArrayList<GenoSample> gsamples;
	
	SNP(String i, int c, long l, int n){
		SNPid=i;
		location = new GenomicCoordinate(c, l);
		num= n;
		gsamples = new ArrayList<GenoSample>();
	}

	public void addSample(GenoSample g){
		gsamples.add(g);
	}
	
	public int getNumSamples(){
		return gsamples.size();
	}
	
	public ArrayList<GenoSample> getGenosamples(){
		return gsamples;
	}
	
	public String getId(){
		return SNPid;
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

	
	public static SNP parseSNP(String line, int n){
		String[] tokens = line.split("\\s+");

		SNP s= new SNP(tokens[0], Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]), n);
		
		return s;
	}
}