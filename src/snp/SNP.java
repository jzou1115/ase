package snp;
import sample.*;

public class SNP implements Comparable<SNP>{
	//parsed from *.map files
	String id;
	int chromosome;
	int location;
	int num;
	
	GenoSample[] gsamples;
	
	SNP(String i, int c, int l, int n){
		id=i;
		chromosome=c;
		location=l;
		num= n;
	}

	@Override
	public int compareTo(SNP arg0) {
		return this.id.compareTo(arg0.id);
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