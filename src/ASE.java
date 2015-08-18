import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import genome.*;
import run.*;
import sample.*;

public class ASE {
	SNPgroup isHetero;
	GeneGroup hasASE;
	Map<Gene,List<SNP>> map;
	Map<SNP, int[]> samples;
	Map<Gene, int[]> esamples;
	String[] sampleNames;
	Gene gene;
	SNPgroup snps;
	int[] ase;
	Map<SNP, List<GenoSample>> genoMap;
	
	public void parseSnps(InputStream inputStream){
		System.out.println("Reading snps");
		isHetero = SNPgroup.readSNPGroup(inputStream);
		isHetero.sort();
		System.out.println("number of snps: "+ isHetero.size());
	}
	
	public void parseGenes(InputStream inputStream){
		System.out.println("Reading genes");
		hasASE = GeneGroup.readGeneGroup(inputStream);
		hasASE.sort();
		System.out.println("number of genes: "+ hasASE.size());
	}
	

	public void parseGenotypes(InputStream genotypes) throws IOException{
		System.out.println("Reading genotypes");
		BufferedReader br = new BufferedReader(new InputStreamReader(genotypes));
		String line = br.readLine();
		
		sampleNames = line.split("\\s+");
		genoMap = new HashMap<SNP, List<GenoSample>>();
		try {
			while((line = br.readLine()) != null){
				try{
					String[] tokens = line.split("\\s+");
					String snp = tokens[0].trim();
					if(snps.contains(snp) && isHetero.contains(snp)){
						//System.out.println(snp);
						SNP s = isHetero.getSNP(snp);
						List<GenoSample> samp = new ArrayList<GenoSample>();
						for(int i=1; i<tokens.length;i++){
							GenoSample genosamp = new GenoSample(sampleNames[i], (int) Math.round(Double.parseDouble(tokens[i])%2));
							//System.out.println(genosamp.getHetero());
							s.addSamples(samp);
							samp.add(genosamp);
						}
						genoMap.put(s, samp);
					}
				} catch (Exception e){
					//do nothing
				}
			}
			for(SNP s: genoMap.keySet()){
				System.out.println("asdf: "+s.getId());
			}
			br.close();
			//System.out.println("done reading genotypes");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	
	public void parseExpressions(InputStream expressions) throws IOException{
		System.out.println("Reading expression data");
		BufferedReader br = new BufferedReader(new InputStreamReader(expressions));
		String line = br.readLine();
		
		String[] gtex = line.split("\\s+");
		for(int j=0; j<gtex.length;j++){
			if(!gtex[j].equals(sampleNames[j])){
				System.out.println("Samples not the same");
				return;
			}
		}
		
		
		try {
			while((line = br.readLine()) != null){
				try{
					String[] tokens = line.split("\\s+");
					
					String gene = tokens[0].trim();
					Gene g = hasASE.getGene(gene);
					
					int[] samp = new int[tokens.length-1];
					for(int i=0; i<tokens.length-1;i++){
						samp[i] = Integer.parseInt(tokens[i+1]);
					}
					esamples.put(g, samp);
					
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
		//System.out.println(s.getId()+"\t"+s.getLocation());
		if(g.region.expand(250000).contains(s.getLocation())){
			return true;
		}
		//System.out.println("Fail: "+g.region.expand(250000).getChromosome()+"\t"+g.region.expand(250000).getStart().getCoord()+"\t"+g.region.expand(250000).getEnd().getCoord());
		return false;
	}
	
	public void genesToSnps() throws IOException{
		map = new HashMap<Gene, List<SNP>>();
		
		List<SNP> snpList = isHetero.getSnps();
		//Collections.sort(snpList);
		
		List<Gene> geneList = hasASE.getGenes();
		//Collections.sort(geneList);
		
		
		int snp = 0;
		for(Gene g: geneList){
			GenomicCoordinate end = g.getRegion().getEnd();
			while(snpList.get(snp).getLocation().compareTo(end)<=0){
				if(match(snpList.get(snp),g)){
					if(map.get(g)==null){
						List<SNP> val = new ArrayList<SNP>();
						val.add(snpList.get(snp));
						map.put(g, val);
					}
					else{
						map.get(g).add(snpList.get(snp));
					}
				}
				snp++;
			}
		}
		
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("genesToSnps.txt")));
		for(Gene g: map.keySet()){
			String print = ">"+g.getId()+"\t"+map.get(g).size();
			for(SNP s: map.get(g)){
				print = print + "\n"+s.getId();
			}
			bw.write(print+"\n");
		}
		bw.close();
		
	}

	public void setTestGene(String s){
		gene = hasASE.getGene(s);
		System.out.println(gene.getId());
		snps = new SNPgroup(map.get(gene));
		System.out.println(snps.size());
	}
	
	public void startRun(double threshold, int errors, int perms, int n) throws IOException{
		Run r = new Run(gene, snps.getSnps(), threshold, errors, perms, n, genoMap);
		System.out.println("Gene run: "+r.getGene().getId());
		System.out.println(r.getSnps().size());
		for(SNP s: r.getGenoMap().keySet()){
			System.out.println(s.getId());
		}
		r.allSimulations();
	}

	public static void main(String args[]) throws IOException{
		
		CommandLineParams cmdArgs = new CommandLine();
		try{
			cmdArgs.parse(args);
		} catch (Exception e) {
			System.err.println(e.getMessage());
			cmdArgs.printHelp(System.err);
			System.exit(1);
		}
		if(cmdArgs.help()){
			cmdArgs.printHelp(System.err);
			System.exit(0);
		}
		
		ASE a= new ASE();
		
		a.parseSnps(cmdArgs.getSNPsInput());
	
		a.parseGenes(cmdArgs.getGenesInput());
		
		a.genesToSnps();
		
		a.setTestGene(cmdArgs.getTestGene());
		
		a.parseGenotypes(cmdArgs.getGenotypeData());
		
		//a.parseExpressions(cmdArgs.getExpressionData());
		
		a.startRun(cmdArgs.getThreshold(), cmdArgs.getErrorNum(), cmdArgs.getPermNum(), cmdArgs.getSampleNum());
		
		/** Launch simulation **/
		//int numSimulations = args[4]
		//int threshold = args[5]
		//a.simulate(0, 1000, 10);
		
	}
}
