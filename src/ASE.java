import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.sun.tools.javac.util.Paths;

import gene.*;
import genome.*;
import run.*;
import sample.*;
import snp.*;

public class ASE {
	SNPgroup isHetero;
	GeneGroup hasASE;
	Map<Gene,List<SNP>> map;
	//TODO make nSamples a variable in the future instead of global
	int nSamples = 100;
	
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
					String snp = tokens[0].trim();
					SNP s = isHetero.getSNP(snp);
					//TODO: thinkaboutthis nsamples
					for(int i=1; i<tokens.length && i<nSamples; i++){
						GenoSample g = new GenoSample(sampleNames[i], Math.round(Float.parseFloat(tokens[i]))%2);
						s.addSample(g);
						//System.out.println(g.toString());
					}
					System.out.println(snp+"\t"+s.getNumSamples());
				} catch (Exception e){
					//do nothing
				}
			}
			br.close();
		}catch (IOException e) {
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
					
					String gene = tokens[0].trim();
					Gene g = hasASE.getGene(gene);
					for(int i=1; i<tokens.length;i++){
						ExpSample e = new ExpSample(sampleNames[i], Integer.parseInt(tokens[i]));
						g.addSample(sampleNames[i], e);
						//System.out.println(e.toString());
					}
					//System.out.println(gene+"\t"+g.getNumSamples());
					
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
	
	public boolean match(SNP s, Gene g){
		if(g.region.expand(100).contains(s.getLocation())){
			return true;
		}
		return false;
	}
	
	public boolean isSNPNeeded(String snpId){
		for(Gene g : hasASE.getGenes()){
			for (SNP s : map.get(g)){
				if (s.getId() == snpId){
					return true;
				}
			}
			
		}
		return false;
	}
	
	public void genesToSnps(){
		map = new HashMap<Gene, List<SNP>>();
		
		List<SNP> snpList = isHetero.getSnps();
		//Collections.sort(snpList);
		
		List<Gene> geneList = hasASE.getGenes();
		//Collections.sort(geneList);
		
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
					break;
				}
			}
		}
	}
	/**
	for(Gene g:map.keySet()){
		List<SNP> s = map.get(g);
		for(SNP x:s){
			System.out.println(x.getId()+"\t"+g.getId());
		}
	}
 	**/
	
	public void simulate(int errors, int reps){
		for(Gene g: hasASE.getGenes()){
			int total=0;
			for(int r=0; r<reps; r++){
				Run run = new Run(g, map.get(g),errors);
				int variants = run.runSim();
				total = total + variants;
			}
			System.out.println(g.getId()+"\t"+1.0*total/reps);
		}
	}

	public static void main(String args[]) throws IOException{
		ASE a= new ASE();
				
		//String geneData = args[1];
		FileInputStream geneData = new FileInputStream(new File("./test/geneLoc.txt"));
		a.parseGenes(geneData);
		
		/** Parse all data files **/
		//String snpData = args[0];
		FileInputStream snpData = new FileInputStream(new File("./test/snp.map"));
		//FileInputStream snpData = new FileInputStream(new File("./test3/ChrOne.map"));
		a.parseSnps(snpData);

		//mapping
		a.genesToSnps();
		
		//String genotypeData = args[2];
		FileInputStream genotypeData = new FileInputStream(new File("./test/isHetero.txt"));
		//FileInputStream genotypeData = new FileInputStream(new File("./test3/ChrOne.snps.txt"));
		a.parseGenotypes(genotypeData);
		
		
		//String expData = args[3];
		//FileInputStream expData = new FileInputStream(new File("./test/hasASE.txt"));
		//a.parseExpressions(expData);
		

		/** Launch simulation **/
		//int numSimulations = args[4]
		//int threshold = args[5]
		a.simulate(0, 100);
		
	}
	
}
