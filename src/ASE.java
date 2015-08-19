import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import gene.*;
import run.*;
import snp.*;

public class ASE {

	SNPgroup snpGroup;
	GeneGroup geneGroup;
	Map<Gene,List<SNP>> map;
	//TODO make nSamples a variable in the future instead of global
	int nSamples = 100;
	
	public static boolean match(SNP s, Gene g){
		if(g.region.expand(100000).contains(s.getLocation())){
			return true;
		}
		else{
			return false;
		}
	}
	
	public boolean isSNPNeeded(SNP snp){
		for(Gene g : geneGroup.getGenes()){
			if (match (snp, g)){
				return true;
			}
		}
		return false;
	}
	
	public void genesToSnps(){
		map = new HashMap<Gene, List<SNP>>();
	
		List<SNP> snpList = snpGroup.getSnps();		
		List<Gene> geneList = geneGroup.getGenes();
	
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
		for(Gene g: geneGroup.getGenes()){
			List<SNP> snps = map.get(g);
			if (snps == null){
				System.out.println(g.getId() + " has no SNPs that map to it");
			}
			else{
				
				Run firstRun = new Run(g, map.get(g), errors, 2);
				int v = firstRun.run();
				System.out.println(g.getId()+ " has " + snps.size() + " snps that map to it and " + v + " variants");
				for (String s : firstRun.variantIds){
					System.out.println(s);
				}
				
				int total=0;
				for(int r=0; r<reps; r++){
					Run run = new Run(g, map.get(g),errors, 3);
					int variants = run.run();
					total = total + variants;
					if (run.variantIds.size()!=0){
						System.out.println("Variants matched:" + run.variantIds);
					}
				}
				System.out.println(g.getId()+ " has " + snps.size() + " snps that map to it and " + 1.0*total/reps + " average variants");
			}
		}
	}
	
	public static void main(String args[]) throws IOException{
		ASE a= new ASE();
		Parse parse = new Parse(a);
		/*
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
		*/
		
		//String geneData = args[1];
		//String snpData = args[0];
		//String genotypeData = args[2];
		//String expData = args[3];
		//a.parseExpressions(expData);

		
		FileInputStream geneData = new FileInputStream(new File("./test3/genes4.txt"));
		parse.parseGenes(geneData);
		
		FileInputStream snpData = new FileInputStream(new File("./test3/ChrOne.map"));
		parse.parseSnps(snpData);

		a.genesToSnps();
		
		FileInputStream genotypeData = new FileInputStream(new File("./test3/ChrOne.snps.txt"));
		parse.parseGenotypes(genotypeData);
		
		a.simulate(0, 100);
		
		//a.simulate(30, 100);
	}
	
}
