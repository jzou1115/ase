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
		if(tokens.length>1){
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
	
	public static void parseGenotypes(InputStream genotypes, List<SNP> snps, Map<String, SNP> snpLoc, List<String> sampleNames) throws IOException{
		System.out.println("Reading genotypes");
		BufferedReader br = new BufferedReader(new InputStreamReader(genotypes));
		String line = br.readLine();
		//System.out.println(line);
		String[] samples = line.split("\t");
		//System.out.println("snploc size: " + snpLoc.size());
		//System.out.println(samples.length);
		//System.out.println(sampleNames.size());
		try {
			while((line = br.readLine()) != null){
					String[] tokens = line.split("\\s+");
					String[] snpTokens = tokens[0].split("_");
					int chr = Integer.parseInt(snpTokens[0]);
					int pos = Integer.parseInt(snpTokens[1]);
					String snpid = chr+"_"+pos;
					SNP s = snpLoc.get(snpid);
					if(s==null){
						//System.out.println(snpid+" not found in snpLoc");
					}
					else{
						for(int i=1; i<tokens.length;i++){
							//System.out.println(i);
							if(sampleNames.contains(samples[i])){
								
								GenoSample genosamp = new GenoSample(samples[i], (int) Math.round(Double.parseDouble(tokens[i]))%2);
								s.addSample(genosamp);	
							}
						}
						s.sortSamples();
								
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
