package functions;

import genome.ChromState;
import genome.Gene;
import genome.GenomicCoordinate;
import genome.SNP;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import parse.ParseChromState;
import parse.ParseGene;
import parse.ParseSNP;


public class GenesToSNP {
	
	List<SNP> isHetero;
	List<Gene> hasASE;
	HashMap<Gene, List<SNP>> map;
	
	public GenesToSNP(InputStream snps, InputStream genes) throws IOException {
		parseSnps(snps);
		parseGenes(genes);
		genesToSnps();
	}

	
	public GenesToSNP(InputStream snps, InputStream genes, InputStream chrom) throws IOException {
		parseSnps(snps);
		List<ChromState> chromatin = ParseChromState.parseChromState(chrom);
		Map<SNP, ChromState> map = AssignChromState.assignStateSNP(isHetero, chromatin);
		parseGenes(genes);
		genesToSnps();
	}
	
	public void parseSnps(InputStream inputStream){
		System.out.println("Reading SNPs");
		isHetero = ParseSNP.readSNPGroup(inputStream);
		Collections.sort(isHetero);
		System.out.println("Number of SNPs: "+ isHetero.size());
	}
	
	public void parseGenes(InputStream inputStream){
		System.out.println("Reading genes");
		hasASE = ParseGene.readGeneGroup(inputStream);
		Collections.sort(hasASE);
		System.out.println("Number of genes: "+ hasASE.size());
	}

	/**
	public void genesToSnps() throws IOException{
		map = new HashMap<Gene, List<SNP>>();
	
		for(SNP s:isHetero){
			for(Gene g:hasASE){
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
	**/
	public boolean match(SNP s, Gene g){
		if(g.region.expand(100000).contains(s.getLocation())){
			return true;
		}
		return false;
	}

	public void write(File output, String filename) throws IOException {
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(output+File.separator+filename)));
		for(Gene g: map.keySet()){
			String print = ">"+g.toString();
			for(SNP s: map.get(g)){
				print = print + "\n"+s.toString();
			}
			bw.write(print+"\n");
		}
		bw.close();
	}
	
	
	public void genesToSnps() throws IOException{
		map = new HashMap<Gene, List<SNP>>();
		
		Collections.sort(hasASE);
		Collections.sort(isHetero);
		
		int snp = 0;
		for(Gene g: hasASE){
			GenomicCoordinate end = g.getRegion().getEnd();
			while(isHetero.get(snp).getLocation().compareTo(end)<=0){
				if(match(isHetero.get(snp),g)){
					if(map.get(g)==null){
						List<SNP> val = new ArrayList<SNP>();
						map.put(g, val);
					}
					List<SNP> val = map.get(g);
					val.add(isHetero.get(snp));
					map.put(g, val);
				}
				snp++;
			}
		}
		
	}


}
