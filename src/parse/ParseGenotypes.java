package parse;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import genome.SNP;
import sample.GenoSample;

public class ParseGenotypes {
	public static void parseGenotypes(InputStream genotypes, List<SNP> snps, Map<String, SNP> snpMap) throws IOException{
		System.out.println("Reading genotypes");
		BufferedReader br = new BufferedReader(new InputStreamReader(genotypes));
		String line = br.readLine();

		String[] samples = line.split("\t");

		try {
			while((line = br.readLine()) != null){
					String[] tokens = line.split("\\s+");
					SNP s = snpMap.get(tokens[0]);
					if(s==null){
					//	System.out.println(tokens[0] +" not found in gene to SNP map");
					}
					else{
						//skip token[0] (SNP id)
						for(int i=1; i<tokens.length;i++){	
								// 0, 2 homozygous; 1 heterozygous
								GenoSample genosamp = new GenoSample(samples[i], (int) Math.round(Double.parseDouble(tokens[i]))%2);
								s.addSample(genosamp);	
						}
					}
				
			}

			br.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
	}
	
	public static List<SNP> parseGenotypes(InputStream genotypes) throws IOException{
		System.out.println("Reading genotypes");
		BufferedReader br = new BufferedReader(new InputStreamReader(genotypes));
		String line = br.readLine();

		String[] samples = line.split("\t");

		List<SNP> snps = new ArrayList<SNP>();
		try {
			while((line = br.readLine()) != null){
				String[] tokens = line.split("\\s+");
				SNP s = new SNP(tokens[0]);
				//skip token[0] (SNP id)
				for(int i=1; i<tokens.length;i++){	
					// 0, 2 homozygous; 1 heterozygous
					GenoSample genosamp = new GenoSample(samples[i], (int) Math.round(Double.parseDouble(tokens[i]))%2);
					s.addSample(genosamp);	
				}
				
				snps.add(s);

			}

			br.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		return snps;
		
	}
	
}
