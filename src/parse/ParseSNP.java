package parse;

import genome.SNP;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import sample.GenoSample;

public class ParseSNP {
	
	
	public static SNP parseSNP(String line){
		String[] tokens = line.split("\\s+");
		if(tokens.length==3){
			SNP s= new SNP(tokens[0].trim(), Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]));
			return s;
		}
		else{
			SNP s = new SNP(tokens[0]);
			return s;
		}
	}
	public static List<SNP> readSNPGroup(InputStream in){
		List<SNP> snps = new ArrayList<SNP>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(in));
		String line;
		try {
			while((line = reader.readLine()) != null){
				try{
					SNP s = parseSNP(line);
					if(s!=null){
						snps.add(s);
					}
				} catch (Exception e){
					//do nothing
				}
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return snps;
	}
	
	
	
}
