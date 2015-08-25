import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import gene.*;
import run.*;
import snp.*;

public class ASE {

	SNPgroup snps;
	GeneGroup genes = new GeneGroup();
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
		for(Gene g : genes.getGenes()){
			if (match (snp, g)){
				return true;
			}
		}
		return false;
	}
	
	private void genesToSnps(){
		map = new HashMap<Gene, List<SNP>>();
	
		List<SNP> snpList = snps.getSnps();		
	
		for(SNP s:snpList){
			for(Gene g:genes.getGenes()){
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
	
	private void simulate(int errors){
		for (Gene g: genes.getGenes()){
			List<SNP> snps = map.get(g);
			if (snps == null){
				//genes.remove(g.geneId);
				//System.out.println(g.geneId + "," + "0");
			}
			else{
				Run r = new Run(g, snps, errors, 1);
				r.run();
				System.out.print(g.geneId + "," + snps.size() + "," + r.snp.getId() + "," + 1.0 * r.snp.numberHeterozygous / 94 + ",");
				System.out.print(r.variantIds.size() + ",");
				r.isStatisticallySignificant();
				System.out.println(r.variantIds + ",");
			}
		}
	}
	
	public static void main(String args[]) throws IOException{
			ASE a= new ASE();
			Parse parse = new Parse(a);
	
			//FileInputStream geneData = new FileInputStream(new File("./test/geneLoc.txt"));
			BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(args[0])));
			String line = null;
			do {
				try {
					for(int i = 0; i<300; i++){
						if((line = reader.readLine()) == null)
							break;
						Gene g = parse.parseGene(line);
						a.genes.add(g.geneId, g);
					}
				} catch (IOException e) {
					e.printStackTrace();
				}
	
				/*
				//test data
				FileInputStream snpData = new FileInputStream(new File("./test/snp.map"));
				FileInputStream genotypeData = new FileInputStream(new File("./test/isHetero.txt"));
				FileInputStream expData = new FileInputStream(new File("./test/hasASE.txt"));
				parse.parseExpressions(expData);
				*/
				
				parse.parseSnps("./test3/ChrOne.map");
				a.genesToSnps();
				parse.parseGenotypes("./test3/ChrOne.snps.txt");

				a.simulate(Integer.parseInt(args[1]));
				a = new ASE();
				parse = new Parse(a);
				System.gc();
			}
			while (line != null);
	}

	
}
