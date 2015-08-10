import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import gene.*;
import genome.*;
import run.*;
import sample.*;
import snp.*;

public class ASE {
	SNPgroup isHetero;
	GeneGroup hasASE;
	
	public void parseSnps(FileInputStream snpData){
		isHetero = SNPgroup.readSNPGroup(snpData);
	}
	
	public void parseGenes(FileInputStream geneData){
		hasASE = GeneGroup.readGeneGroup(geneData);
	}
	
	public void parseGenotypes(FileInputStream genotypes) throws IOException{
		BufferedReader br = new BufferedReader(new InputStreamReader(genotypes));
		String line = br.readLine();
		String[] sampleNames = line.split("\\s+");
		
		try {
			while((line = br.readLine()) != null){
				try{
					String[] tokens = line.split("\\s+");
					String gene = tokens[0].trim();
					Gene g = hasASE.getGene(gene);
					for(int i=1; i<tokens.length;i++){
						ExpSample e = new ExpSample(sampleNames[i], Integer.parseInt(tokens[i]));
						g.addSample(sampleNames[i], e);
					}
				} catch (Exception e){
					//do nothing
				}
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public void parseExpressions(FileInputStream expressions) throws IOException{
		BufferedReader br = new BufferedReader(new InputStreamReader(expressions));
		String line = br.readLine();
		
		String[] sampleNames = line.split("\\s+");
		
		try {
			while((line = br.readLine()) != null){
				try{
					String[] tokens = line.split("\\s+");
					String snp = tokens[0].trim();
					SNP s = isHetero.getSNP(snp);
					for(int i=1; i<tokens.length;i++){
						GenoSample g = new GenoSample(sampleNames[i], Integer.parseInt(tokens[i]));
						s.addSample(sampleNames[i], g);
					}
				} catch (Exception e){
					//do nothing
				}
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void main(String args[]) throws IOException{
		ASE a= new ASE();
		
		//String snpData = args[0];
		FileInputStream snpData = new FileInputStream(new File("/home/jennifer/ase/test/snp.map"));
		a.parseSnps(snpData);
		
		//String geneData = args[1];
		FileInputStream geneData = new FileInputStream(new File("/home/jennifer/ase/test/geneLoc.txt"));
		a.parseGenes(geneData);
		
		//String genotypeData = args[2];
		FileInputStream genotypeData = new FileInputStream(new File("/home/jennifer/ase/test/isHetero.txt"));
		a.parseGenotypes(genotypeData);
		
		//String expData = args[2];
		FileInputStream expData = new FileInputStream(new File("/home/jennifer/ase/test/hasASE.txt"));
		a.parseExpressions(expData);
	}
}
