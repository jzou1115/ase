package functions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Random;

import parse.ParseGene;
import parse.ParseMap;
import parse.ParseSNP;
import parse.ParseSamples;
import sample.ExpSample;
import sample.GenoSample;
import genome.DiSNP;
import genome.Gene;
import genome.SNP;

public class Combinations {
	Gene g;
	List<SNP> snps;
	//List<DiSNP> combs;
	List<String> gtexIds;
	File outdir;
	String filename;
	Map<String, SNP> snpLoc;
	int sampleSize;
	int errors;
	
	public Combinations(InputStream map, String gene, InputStream genotypes, InputStream expression, int n, int e, File out, String f) throws IOException {
		sampleSize = n;
		errors = e;
		outdir =out;
		filename = f;
		setTestGene(map, gene);
		ParseSNP.parseGenotypes(genotypes, snps, snpLoc);
		ParseGene.parseExpressions(expression, g, outdir);
	//	getCombinations();
	//	write("disnps.txt");
		mapASE();
		
	}

	public void setTestGene(InputStream map2, String gene2){
		ParseMap parsemap = new ParseMap();
		parsemap.parseMap(map2, gene2);
		g = parsemap.getGene();
		snps = parsemap.getSNPs();
		snpLoc = parsemap.getSnpLoc();
		//System.out.println(gene.toString()+"\t"+snps.size());
	}
	/**
	public void getCombinations(){
		combs = new ArrayList<DiSNP>();
		for(int i=0; i<snps.size(); i++){
			for(int j=i+1; j<snps.size(); j++){
				DiSNP combinedSNP = new DiSNP(snps.get(i),snps.get(j));
				combs.add(combinedSNP);
			}
		}		
	}

	public void write(String st) throws IOException{
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+st)));
		
		for(DiSNP s:combs){
			outfile.write(s.toString()+"\n");
		}
		
		outfile.close();
	}
**/
	//mapase without permutations
	public void mapASE() throws IOException{
		List<ExpSample> ase = g.getExpsamples();
		if(ase.size()<sampleSize){
			System.out.println("Not enough expsamples");
			System.exit(1);
		}
	

		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+filename)));
		BufferedWriter outfile2 = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+filename+"_significant.txt")));

		for(int i=0; i<snps.size(); i++){
			for(int j=i+1; j<snps.size(); j++){
				DiSNP s = new DiSNP(snps.get(i),snps.get(j));
				/**
				List<GenoSample> asdf = s.getGenosamples();
				String line2 = s.getId();
				for(GenoSample g:asdf){
					line2 = line2+"\t"+g.getHetero();
				}
				System.out.println(line2);
				**/
				
				//for(DiSNP s:combs){
				//subset of genosamples that have expression samples
				List<GenoSample> subset = getSubset(ase.size(), sampleSize, s, ase);
				if(subset==null){
					System.out.println("Not enough genosamples in "+s.getId());
					continue;
				}

				//m and k are used to approximate significance
				//number of ones in ase samples
				int m=0;
				//number of ones in genotype samples
				int k=0;

				int correct=0;
				int incorrect=0;
				for (int n=0; n<subset.size();n++) {
					GenoSample geno = subset.get(n);
					String sampleID = geno.getID();
					ExpSample e = g.getSample(sampleID);

					int hasASE = e.getASE();
					if(hasASE==1){
						m=m+1;
					}
					int isHetero = geno.getHetero();
					if(isHetero==1){
						k=k+1;
					}

					if(isHetero == hasASE){
						correct++;
					}
					else if(isHetero != hasASE){
						incorrect++;
					}	
					else{
						System.out.println("Missing data: "+sampleID );
					}
				}

				ComputeSig sig = new ComputeSig(sampleSize, m, k, incorrect);
				double p = sig.significance();
				String line = g.toString()+"\t"+s.toString()+"\t"+correct+"\t"+incorrect+"\t"+p;
				if(isSignificant(errors, incorrect, p)){
					line = line+"\t"+1;
					outfile2.write(line);
				}
				else{
					line = line+"\t"+0;
				}

				outfile.write(line+"\n");
			}
		}

		outfile.close();
		outfile2.close();

	}

	public boolean isSignificant(int e, int incorrect, double z){
		if(incorrect<=e){
			if(z<=.0000025){
				return true;
			}
		}
		return false;
	}
	//return subset of genosamples that also have expsamples
	public List<GenoSample> getSubset(int total, int num, DiSNP s, List<ExpSample> exp){
		List<String> sampleids = new ArrayList<String>();
		for(int j=0; j< exp.size(); j++){
			sampleids.add(exp.get(j).getID());
		}
		
		List<GenoSample> geno = new ArrayList<GenoSample>(s.getGenosamples());
		Collections.shuffle(geno, new Random(13));

		List<GenoSample> subset = new ArrayList<GenoSample>();
		int i=0; //index for geno
		int k=0; //index for retInd
	
		//iterate over genosamples and put index (i) of ones with matching expsample in retInd
		while(k<num){
			//not enough matching genosamples and expsamples
			if(i>=geno.size()){
				return null;
			}
			//matching genosample and expsample
			if(sampleids.contains(geno.get(i).getID())){
				subset.add(geno.get(i));
				k++;
			}
			i++;
		}

		//enough matching genosamples and expsamples
		return subset;
	}
	

}
