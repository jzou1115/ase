package functions;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import genome.Gene;
import genome.SNP;
import parse.ParseGene;
import parse.ParseMap;
import parse.ParseSNP;
import parse.ParseSamples;

public class MapASE {
	Gene gene;
	List<SNP> snps;
	Map<String, SNP> snpLoc;
	
	public void setTestGene(Gene g, List<SNP> s){
		gene = g;
		snps = s;
		snpLoc = new HashMap<String, SNP>();
		//System.out.println(gene.toString()+"\t"+snps.size());
	}
	
	public void setTestGene(InputStream map, String g, InputStream genotypes) throws IOException {
		ParseMap parsemap = new ParseMap();
		parsemap.parseMap(map, g);
		gene = parsemap.getGene();
		snps = parsemap.getSNPs();
		snpLoc = parsemap.getSnpLoc();
		
	}


	public void parseGenotypes(InputStream genotypes) throws IOException{
		ParseSNP.parseGenotypes(genotypes,snps, snpLoc);
	}
	
	public void parseExpressions(InputStream expressions, File outdir) throws IOException {
		ParseGene.parseExpressions(expressions, gene, outdir);
	}

	//with perm
	public void startRun(int error, int n, int perm, File outdir) throws IOException {
		Run r = new Run(gene, snps, error, perm, n, outdir);
		r.mapASE(gene.getId());
	}
	
	
	//without perm
	public void startRun(int error, int n, int perm, File outdir, String filename) throws IOException {
		Run r = new Run(gene, snps, error, perm, n, outdir);
		r.mapASE(gene.getId(), filename);
	}

}
