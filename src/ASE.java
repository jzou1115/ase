import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import gene.*;
import run.*;
import sample.*;
import snp.*;

//TODO Split ASE.java into Parsing.java and ASE.java
public class ASE {

	SNPgroup isHetero;
	GeneGroup hasASE;
	Map<Gene,List<SNP>> map;
	//TODO make nSamples a variable in the future instead of global
	int nSamples = 100;
	
	public void parseSnps(FileInputStream snpData){
		List<SNP> snps = new ArrayList<SNP>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(snpData));
		String line;
		int n=0;
		try {
			while((line = reader.readLine()) != null){
				try{
					SNP snp = SNP.parseSNP(line, n);
					if (snp != null && isSNPNeeded(snp)){
						snps.add(snp);
						n++;
					}
				} catch (Exception e){
					e.printStackTrace();
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		isHetero =  new SNPgroup(snps);
	}
	
	
	//TODO factor parseGenes into this class
	public void parseGenes(FileInputStream geneData){
		hasASE = GeneGroup.readGeneGroup(geneData);
	}
	
	//TODO generate random n samples, not the 1st n samples
	public void parseGenotypes(FileInputStream genotypes) throws IOException{
		BufferedReader br = new BufferedReader(new InputStreamReader(genotypes));
		String line = br.readLine();
		String[] sampleNames = line.split("\\s+");
		
		try {
			while((line = br.readLine()) != null){
				try{
					String[] tokens = line.split("\\s+");
					String snpId = tokens[0].trim();
					SNP s = isHetero.getSNP(snpId);
					if (s != null) {
						for(int i=1; i<tokens.length && i<nSamples; i++){
							GenoSample g = new GenoSample(sampleNames[i], Math.round(Float.parseFloat(tokens[i]))%2);
							s.addSample(g);
							//System.out.println(g.toString());
						}
						//System.out.println("in parseGenotypes:  " + snpId +"\t"+s.getNumSamples());
					}
				} catch (Exception e){
					e.printStackTrace();
				}
			}
			br.close();
		}catch (IOException e) {
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
					
					String gene = tokens[0].trim();
					Gene g = hasASE.getGene(gene);
					for(int i=1; i<tokens.length;i++){
						ExpSample e = new ExpSample(sampleNames[i], Integer.parseInt(tokens[i]));
						g.addSample(sampleNames[i], e);
						//System.out.println(e.toString());
					}
					//System.out.println(gene+"\t"+g.getNumSamples());
					
				} catch (Exception e){
					e.printStackTrace();
				}
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	//TODO combine match and isSNPNeeded
	public boolean match(SNP s, Gene g){
		if(g.region.expand(100000).contains(s.getLocation())){
			return true;
		}
		else{
			return false;
		}
	}
	
	public boolean isSNPNeeded(SNP snp){
		for(Gene g : hasASE.getGenes()){
			if (match (snp, g)){
				return true;
			}
		}
		return false;
	}
	
	public void genesToSnps(){
		map = new HashMap<Gene, List<SNP>>();
	
		List<SNP> snpList = isHetero.getSnps();		
		List<Gene> geneList = hasASE.getGenes();
	
		for(SNP s:snpList){
			for(Gene g:geneList){
				if(match(s,g)){
					if(map.get(g)==null){
						List<SNP> val = new ArrayList<SNP>();
						val.add(s);
						map.put(g, val);
					}
					else{
						map.get(g).add(s);
					}
				}
			}
		}
	}
	
	public void printMapping(){
		for(Gene g:map.keySet()){
			System.out.println(g.getId());
			System.out.println(g.region.getStart().getCoord() + " - " + g.region.getEnd().getCoord());
			List<SNP> snps = map.get(g);
			for(SNP s:snps){
				System.out.println(s.getId()+"\t"+s.getLocation());
			}
		}
	}
	
	public void simulate(int errors, int reps){
		for(Gene g: hasASE.getGenes()){
			List<SNP> snps = map.get(g);
			if (snps == null){
				System.out.println(g.getId() + " has no SNPs that map to it");
			}
			else{
				int total=0;
				for(int r=0; r<reps; r++){
					Run run = new Run(g, map.get(g),errors);
					int variants = run.run();
					total = total + variants;
				}
				System.out.println(g.getId()+ " has " + snps.size() + " snps that map to it and " + 1.0*total/reps + " average variants");
			}
		}
	}
	
	public static void main(String args[]) throws IOException{
		ASE a= new ASE();
		
		
		// test data
		FileInputStream geneData = new FileInputStream(new File("./test/geneLoc.txt"));
		a.parseGenes(geneData);

		FileInputStream snpData = new FileInputStream(new File("./test/snp.map"));
		a.parseSnps(snpData);
		
		a.genesToSnps();

		FileInputStream genotypeData = new FileInputStream(new File("./test/isHetero.txt"));
		a.parseGenotypes(genotypeData);

		FileInputStream expData = new FileInputStream(new File("./test/hasASE.txt"));
		a.parseExpressions(expData);

		a.genesToSnps();
		
		
		//String geneData = args[1];
		//String snpData = args[0];
		//String genotypeData = args[2];
		//String expData = args[3];
		//a.parseExpressions(expData);

		
		/*
		FileInputStream geneData = new FileInputStream(new File("./test3/genes2.txt"));
		a.parseGenes(geneData);
		
		FileInputStream snpData = new FileInputStream(new File("./test3/ChrOne.map"));
		a.parseSnps(snpData);

		a.genesToSnps();
		
		FileInputStream genotypeData = new FileInputStream(new File("./test3/ChrOne.snps.txt"));
		a.parseGenotypes(genotypeData);
		*/
		
		a.simulate(0, 2000);
		
		//a.simulate(30, 100);
	}
	
}
