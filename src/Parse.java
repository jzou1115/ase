import gene.Gene;
import genome.GenomicCoordinate;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import sample.ExpSample;
import sample.GenoSample;
import snp.SNP;
import snp.SNPgroup;

public class Parse {
	
	public ASE ase;
	public Parse (ASE a){
		ase = a;
	}
	
	public void parseSnps(FileInputStream map){
		List<SNP> snps = new ArrayList<SNP>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(map));
		String line;
		int n=0;
		try {
			while((line = reader.readLine()) != null){
				SNP snp = parseSNP(line, n);
				if (snp != null && ase.isSNPNeeded(snp)){
					snps.add(snp);
					n++;
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		ase.snps = new SNPgroup(snps);
	}
	
	private SNP parseSNP(String line, int n){
		String[] tokens = line.split("\\s+|,");
		SNP s= new SNP(tokens[0], Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]), n);
		return s;
	}
	
	public Gene parseGene(String line){
		String[] tokens = line.split("\\s|,");
		String id = tokens[0];
		int chr = Integer.parseInt(tokens[3]);
		long s = Long.parseLong(tokens[1]);
		long e = Long.parseLong(tokens[2]);
		
		GenomicCoordinate start = new GenomicCoordinate(chr, s);
		GenomicCoordinate end = new GenomicCoordinate(chr, e);
		
		return new Gene(id, start, end);
	}
	
	//TODO generate random n samples, not the 1st n samples
	/* Reads genotypes of SNPs, adds minor allele frequency to SNP */
	public void parseGenotypes(FileInputStream genotypes) throws IOException{
		BufferedReader br = new BufferedReader(new InputStreamReader(genotypes));
		String line = br.readLine();
		String[] sampleNames = line.split("\\s+");
		
		try {
			while((line = br.readLine()) != null){
				String[] tokens = line.split("\\s+");
				String snpId = tokens[0].trim();
				SNP s = ase.snps.getSNP(snpId);
				if (s != null) {
					int numberHeterozygous = 0;
					for(int i=1; i<tokens.length && i<ase.nSamples; i++){
						int isHeterozygous = Math.round(Float.parseFloat(tokens[i]))%2;
						numberHeterozygous = numberHeterozygous + isHeterozygous;
						GenoSample g = new GenoSample(sampleNames[i], isHeterozygous);
						s.addSample(g);
					}
					s.numberHeterozygous = numberHeterozygous;
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
					Gene g = ase.genes.getGene(gene);
					for(int i=1; i<tokens.length;i++){
						ExpSample e = new ExpSample(sampleNames[i], Integer.parseInt(tokens[i]));
						g.addSample(sampleNames[i], e);
					}					
				} catch (Exception e){
					e.printStackTrace();
				}
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
