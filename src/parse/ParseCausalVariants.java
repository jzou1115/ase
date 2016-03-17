package parse;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import genome.SNP;

public class ParseCausalVariants {

	public static SNP parseSNP(String line){
		String[] tokens = line.split("\\s+");
		String[] snpTokens = tokens[1].split("_");
		return new SNP(tokens[1], Integer.parseInt(snpTokens[0]), Integer.parseInt(snpTokens[1]));
	}
	public static List<SNP> readVariantGroup(InputStream in){
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
