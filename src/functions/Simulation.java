package functions;


import genome.Gene;
import genome.SNP;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;
import java.util.Map;

import parse.ParseMap;
import parse.ParseSNP;
import parse.ParseSamples;

public class Simulation {
	Gene gene;
	List<SNP> snps;
	Map<String,SNP> snpLoc;
	List<String> sampleNames;
	
	public void setTestGene(Gene g, List<SNP> s){
		gene = g;
		snps = s;
		//System.out.println(gene.toString()+"\t"+snps.size());
	}
	
	public void setTestGene(InputStream map2, String gene2, InputStream genotypes) throws IOException{
		ParseMap parsemap = new ParseMap();
		parsemap.parseMap(map2, gene2);
		gene = parsemap.getGene();
		snps = parsemap.getSNPs();
		snpLoc = parsemap.getSnpLoc();
		//TODO: implement -x arg
		ParseSNP.parseGenotypes(genotypes,snps, snpLoc);
		//System.out.println(gene.toString()+"\t"+snps.size());
	}
	public void startRun(double threshold, int errors, int perms, int n, int split, File outdir) throws IOException{
		Run r = new Run(gene, snps, threshold, errors, perms, n, outdir);
		r.randomSimulation(split);
	}
}
