package parse;

import genome.SNP;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import sample.GenoSample;

public class ParseSNP {
	public static SNP parseSNP(String line){
		String[] tokens = line.split("\\s+");
//		if(Pattern.matches(SNP_REGEX, tokens[0].trim())){
			SNP s= new SNP(tokens[0].trim(), Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]));
			return s;
//		}
//		System.out.println("does not match snp regex "+tokens[0]);
//		return null;
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
	
	public static void parseGenotypes(InputStream genotypes, List<SNP> snps) throws IOException{
		System.out.println("Reading genotypes");
		BufferedReader br = new BufferedReader(new InputStreamReader(genotypes));
		String line = br.readLine();
		System.out.println(line);
		String[] sampleNames = line.split("\\s+");
		
		try {
			while((line = br.readLine()) != null){
				try{
					String[] tokens = line.split("\\s+");
					String snpid = tokens[0].trim();
					for(SNP s:snps){
						if(s.getId().equals(snpid)){
							System.out.println(s.getId());
							for(int i=1; i<tokens.length;i++){
								GenoSample genosamp = new GenoSample(sampleNames[i], (int) Math.round(Double.parseDouble(tokens[i]))%2);
								s.addSample(genosamp);
							}	
						}
					}
				} catch (Exception e){
					//do nothing
				}
			}

			br.close();
			//System.out.println("done reading genotypes");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	
}
