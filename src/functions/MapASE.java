package functions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import genome.Gene;
import genome.SNP;
import parse.ParseGene;
import parse.ParseMap;
import parse.ParseSNP;
import sample.ExpSample;
import sample.GenoSample;

public class MapASE {
	Gene gene;
	List<SNP> snps;
	Map<String, SNP> snpMap;
	
	public void setTestGene(Gene g, List<SNP> s){
		gene = g;
		snps = s;
		snpMap = new HashMap<String, SNP>();
		//System.out.println(gene.toString()+"\t"+snps.size());
	}
	
	public void setTestGene(InputStream map, String g, InputStream genotypes) throws IOException {
		ParseMap parsemap = new ParseMap();
		parsemap.parseMap(map, g);
		gene = parsemap.getGene();
		snps = parsemap.getSNPs();
		snpMap = parsemap.getSnpMap();
		
	}


	public void parseGenotypes(InputStream genotypes) throws IOException{
		ParseSNP.parseGenotypes(genotypes,snps, snpMap);
	}
	
	public void parseExpressions(InputStream expressions, File outdir) throws IOException {
		ParseGene.parseExpressions(expressions, gene, outdir);
	}
	

	public void startRun(int error, int n, int perm, File outdir, String filename) throws IOException {

		List<ExpSample> ase = gene.getExpsamples();
		//System.out.println("Expsamples: "+ase.size());
		List<Double> realPValues = new ArrayList<Double>();
		Map<SNP, String> lines = new HashMap<SNP, String>();
		pointwisePValue(ase, realPValues, lines);
		//System.out.println("pointwise: "+realPValues.size());
		
		double minPointwise = Collections.min(realPValues);
		System.out.println(minPointwise);
		//for adaptive permutations, notSig>25 does not meet genome wide significance threhold of 2.5*10^-6 with 10^7 permutations
		int notSig=0;
		List<Double> permPValues = new ArrayList<Double>();
		for(int p=0; p<perm; p++){
			double min = Collections.min(pointwisePValue(ase));
			permPValues.add(min);
			//System.out.println(min);
			if(min<=minPointwise){
				notSig++;
			}
			if(notSig>25){
				break;
			}
		}
		
		if(notSig>25){
			BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+filename)));
			for(int i=0; i<snps.size(); i++){
				SNP s = snps.get(i);
				String line = lines.get(s)+"\t"+1.0+"\t"+0+"\n";
				outfile.write(line);
			}
			outfile.close();
		}
		else{
			System.out.println("Perm pvalues: "+permPValues.size());
			List<Double> pValues = calcPValues(realPValues, permPValues);
			System.out.println("pValues: "+ pValues.size());
			
			BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+filename)));
			for(int i=0; i<snps.size(); i++){
				double pv = pValues.get(i);
				SNP s = snps.get(i);
				String line = lines.get(s)+"\t"+pv;

				if(isSignificant(pv)){
					line = line+"\t"+1;
					//	outfile2.write(line);
				}
				else{
					line = line+"\t"+0;
				}

				outfile.write(line+"\n");
			}
			outfile.close();

		}

		
	}



	private List<Double> calcPValues(List<Double> realPValues, List<Double> permPValues) {
		List<Double> pValues = new ArrayList<Double>();
		int perm = permPValues.size();
		for(int i=0; i<realPValues.size(); i++){
			double real = realPValues.get(i);
			int count = 0;
			for(int j=0; j<perm; j++){
				if(permPValues.get(j)<=real){
					count++;
				}
			}
			pValues.add(count*1.0/perm);
		}
		return pValues;
	}

	//real data
	public void pointwisePValue(List<ExpSample> ase, List<Double> realPValues, Map<SNP, String> lines){
		for(SNP s:snps){
			//subset of genosamples that have expression samples
			List<GenoSample> subset = getSubset(ase.size(), s, ase);
			if(subset==null){
				continue;
			}

			//m and k are used to approximate significance
			//number of ones in ase samples
			int m=0;
			//number of ones in genotype samples
			int k=0;

			int correct=0;
			int incorrect=0;
			for (int i=0; i<subset.size();i++) {
				GenoSample g = subset.get(i);
				String sampleID = g.getID();
				ExpSample e = gene.getSample(sampleID);

				int hasASE = e.getASE();
				if(hasASE==1){
					m=m+1;
				}
				int isHetero = g.getHetero();
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

			ComputeSig sig = new ComputeSig(subset.size(), m, k, incorrect);
			double p = sig.significance();
			realPValues.add(p);
			String line = gene.toString()+"\t"+s.toString()+"\t"+m+"\t"+k+"\t"+subset.size()+"\t"+incorrect+"\t"+p;
		//	System.out.println(line);
			lines.put(s, line);
		}

	}
	
	//permuted data
	public List<Double> pointwisePValue(List<ExpSample> ase){
		List<Double> ret = new ArrayList<Double>();
		for(SNP s:snps){
			//subset of genosamples that have expression samples
			List<GenoSample> subsetGeno = getSubset(ase.size(), s, ase);
			if(subsetGeno==null){
				continue;
			}
			List<ExpSample> subsetExp = new ArrayList<ExpSample>();
			for(GenoSample g: subsetGeno){
				subsetExp.add(gene.getSample(g.getID()));
			}
			Collections.shuffle(subsetGeno);
			//m and k are used to approximate significance
			//number of ones in ase samples
			int m=0;
			//number of ones in genotype samples
			int k=0;

			int correct=0;
			int incorrect=0;
			for (int i=0; i<subsetGeno.size();i++) {
				GenoSample g = subsetGeno.get(i);
				ExpSample e = subsetExp.get(i);

				int hasASE = e.getASE();
				if(hasASE==1){
					m=m+1;
				}
				int isHetero = g.getHetero();
				if(isHetero==1){
					k=k+1;
				}

				if(isHetero == hasASE){
					correct++;
				}
				else if(isHetero != hasASE){
					incorrect++;
				}	
			}
			//System.out.println(s.getId()+"\t"+subsetGeno.size()+"\t"+m+"\t"+k+"\t"+incorrect);
			ComputeSig sig = new ComputeSig(subsetGeno.size(), m, k, incorrect);
			double p = sig.significance();
			ret.add(p);
		}
		return ret;
	}
	

	//return subset of genosamples that also have expsamples
	public List<GenoSample> getSubset(int total, SNP s, List<ExpSample> exp){
		Set<String> sampleids = new HashSet<String>();
		for(ExpSample e:exp){
			sampleids.add(e.getID());
		}
		
		List<GenoSample> geno = new ArrayList<GenoSample>(s.getGenosamples());
		Set<String> sampleids2 = new HashSet<String>();
		for(int j=0; j< geno.size(); j++){
			sampleids2.add(geno.get(j).getID());
		}
		
		sampleids.retainAll(sampleids2);
		
		List<GenoSample> subset = new ArrayList<GenoSample>();
		for(GenoSample g:geno){
			if(sampleids.contains(g.getID())){
				subset.add(g);
			}
		}
		return subset;
	}
	
	public boolean isSignificant(double z){
		if(z<=.0000025){
			return true;
		}
		return false;
	}



}
