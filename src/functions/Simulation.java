package functions;


import genome.Gene;
import genome.SNP;

import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import parse.ParseMap;
import parse.ParseSNP;

public class Simulation {
	Gene gene;
	List<SNP> snps;

	public Simulation(InputStream map2, InputStream genotypes, String gene2,
			double threshold, int error, int perm, int n) throws IOException {
		
		setTestGene(map2, gene2);
		ParseSNP.parseGenotypes(genotypes,snps);
		startRun(threshold, error, perm, n);
	}
	
	public Simulation(){
		
	}
	
	public void setTestGene(Gene g, List<SNP> s){
		gene = g;
		snps = s;
		System.out.println(gene.toString()+"\t"+snps.size());
	}
	
	public void setTestGene(InputStream map2, String gene2){
		ParseMap parsemap = new ParseMap();
		parsemap.parseMap(map2, gene2);
		gene = parsemap.getGene();
		snps = parsemap.getSNPs();
		System.out.println(gene.toString()+"\t"+snps.size());
	}
	public void startRun(double threshold, int errors, int perms, int n) throws IOException{
		Run r = new Run(gene, snps, threshold, errors, perms, n);
		r.allSimulations();
	}
}
