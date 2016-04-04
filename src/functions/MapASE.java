package functions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
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

import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;



import genome.Gene;
import genome.SNP;
import parse.ParseExpressions;
import parse.ParseGenotypes;
import parse.ParseMap;
import sample.ExpSample;
import sample.GenoSample;

public class MapASE {
	
	Gene gene;
	List<SNP> snps;
	//Map<String, SNP> snpMap;
	Map<SNP, double[]> pmap;
	
	int perm;
	//double adaptive;
	File outdir;
	String filename;
	
	
	public MapASE(InputStream map, InputStream genotypes,
			InputStream expressions, String g, int p, File o, String f) throws IOException{
		
		pmap = new HashMap<SNP, double[]>();
		
		perm=p; //total number of permutations
		//adaptive = .0000025*perm; //cutoff for adaptive permutations (.0000025 is Bonferonni correction for 20000 genes and alpha=.05)
		
		//output
		outdir = o;
		filename=f;
		
		ParseMap parsemap = new ParseMap();
		parsemap.parseMap(map, g);
		gene = parsemap.getGene();
		snps = parsemap.getSNPs();
		Map<String, SNP> snpMap = parsemap.getSnpMap();	
		
		ParseGenotypes.parseGenotypes(genotypes,snps, snpMap);
		ParseExpressions.parseExpressions(expressions, gene, outdir);
	}
	
	public MapASE(Gene g, List<SNP> s, int p, File o, String f){
		gene = g;
		snps = s;
		perm = p;
		outdir = o;
		filename = f;
		
		pmap = new HashMap<SNP, double[]>();
	}


	public void mapase() throws IOException {

		List<ExpSample> ase = gene.getExpsamples();

		//Calculate pointwise p-values using real genotype and ase data
		//also fills pmap
		List<Double> pointwisePValues = new ArrayList<Double>(); //list of p-values for SNPs in snps
		Map<SNP, String> lines = new HashMap<SNP, String>(); //lines for output
		
		pointwisePValue(ase, pointwisePValues, lines);
		
		writePointwise(lines);
		
		//When there are more than $adaptive permutations that are lower than $minPointwise, no SNPs are significant.  Stop permutations when $notSig>$adaptive 
		double minPointwise = Collections.min(pointwisePValues);
		int notSig=0;
		List<Double> permPValues = new ArrayList<Double>(); //list of minimum p-value for each permutation

		boolean adaptive=false;
		for(int p=1; p<=perm; p++){
			double min = permutationPValue(ase);
			permPValues.add(min);
			
			if(min>=minPointwise){
				notSig++;
			}
			//New adaptive permutations
			BinomialTest bt = new BinomialTest();
			//System.out.println("bt param: "+p+"\t"+notSig);
			if(bt.binomialTest(p, notSig, .0000025, AlternativeHypothesis.GREATER_THAN, .0000025)){
				System.out.println("Permutations: "+ p);
				adaptive=true;
				break;
			}
		}
		
		//nothing is significant
		if(adaptive){
			BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+filename)));
			for(int i=0; i<snps.size(); i++){
				SNP s = snps.get(i);
				String line = lines.get(s)+"\t"+1.0+"\t"+0+"\n";
				outfile.write(line);
			}
			outfile.close();
		}
		else{
			System.out.println("Permutations: "+ perm);
			
			List<Double> pValues = calcPValues(pointwisePValues, permPValues); //p-values from permutation test
			
			BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+filename)));
			for(int i=0; i<snps.size(); i++){
				double pv = pValues.get(i);
				SNP s = snps.get(i);
				String line = lines.get(s)+"\t"+pv;

				if(isSignificant(pv)){
					line = line+"\t"+1;
				}
				else{
					line = line+"\t"+0;
				}

				outfile.write(line+"\n");
			}
			outfile.close();

		}

		
	}


	private void writePointwise(Map<SNP, String> lines) throws IOException {
		System.out.println("Writing pointwise p-values to file");
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+gene.getId()+"_pointwise.txt")));
		for(SNP s: lines.keySet()){
			outfile.write(lines.get(s));
		}
		outfile.close();
	}

	//Calculates pointwise p-values and populates $realPValues and $lines
	public void pointwisePValue(List<ExpSample> ase, List<Double> realPValues, Map<SNP, String> lines){
		for(SNP s:snps){
			//subset of genosamples that have expression samples
			List<GenoSample> subset = getSubset(ase.size(), s, ase);
		
			if(subset==null){
				continue;
			}

			//number of ones in ase samples
			int m=0;
			//number of ones in genotype samples
			int k=0;

			int incorrect=0;
			for (int i=0; i<subset.size();i++) {
				GenoSample g = subset.get(i);
				String sampleID = g.getID();
				ExpSample e = gene.getSample(sampleID);
				
				if(g==null || e==null){
					System.out.println("Subset function does not work.  Either missing GenoSample or missing ExpSample for individual "+sampleID);
					System.exit(1);
				}

				int hasASE = e.getASE();
				if(hasASE==1){
					m=m+1; //increment count of number of individuals with ASE
				}
				int isHetero = g.getHetero();
				if(isHetero==1){
					k=k+1; //increment count of number of individuals that are heterozygous
				}

				if(isHetero != hasASE){
					incorrect++; //increment counter for number of mismatches
				}	

			}

			//calculate p-value based on hypergeometric distribution
			ComputeSig sig = new ComputeSig(subset.size(), m, k, incorrect);
			double p = sig.significance();
			realPValues.add(p);
			
			double[] poss = sig.possiblePValues();
			pmap.put(s, poss);
			/**
			int j = (incorrect - Math.abs(m-k))/2;
			if(poss[j]!=p){
				System.out.println("Bug in cache hypergeometric");
				System.exit(1);
			}
			**/
			String line = gene.toString()+"\t"+s.toString()+"\t"+m+"\t"+k+"\t"+subset.size()+"\t"+incorrect+"\t"+p;
			lines.put(s, line);
		}

	}
	
	//returns minimum pointwise p-value from permuted data
	public double permutationPValue(List<ExpSample> ase){
		List<Double> pvals = new ArrayList<Double>();
		for(SNP s:snps){
			List<GenoSample> subsetGeno = s.getGenosamples();
			/**
			//subset of genosamples that have expression samples
			List<GenoSample> subsetGeno = getSubset(ase.size(), s, ase);
			if(subsetGeno==null){
				continue;
			}
			**/
			List<ExpSample> subsetExp = new ArrayList<ExpSample>();
			for(GenoSample g: subsetGeno){
				ExpSample e = gene.getSample(g.getID());
				if(e==null){
					System.out.println("ERROR");
					System.exit(1);
				}
				
				subsetExp.add(e);
			}
			
			//permute genotypes
			Collections.shuffle(subsetGeno);
			
			//number of ones in ase samples
			int m=0;
			//number of ones in genotype samples
			int k=0;

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

				else if(isHetero != hasASE){
					incorrect++;
				}	
			}
			int j = (incorrect - Math.abs(m-k))/2;
			double p = pmap.get(s)[j];
		
			/**
			//System.out.println(s.getId()+"\t"+subsetGeno.size()+"\t"+m+"\t"+k+"\t"+incorrect);
			ComputeSig sig = new ComputeSig(subsetGeno.size(), m, k, incorrect);
			double p = sig.significance();
		**/	
			pvals.add(p);
		}
		return Collections.min(pvals);
	}
	

	//return subset of genosamples that also have expsamples
	public List<GenoSample> getSubset(int total, SNP s, List<ExpSample> exp){
		//sampleids of expsamples
		Set<String> sampleids = new HashSet<String>();
		for(ExpSample e:exp){
			sampleids.add(e.getID());
		}

		//find subset of genosamples that have matching expsamples
		List<GenoSample> subset = new ArrayList<GenoSample>();
		List<String> remove = new ArrayList<String>();
		for(GenoSample g: s.getGenosamples()){
			if(sampleids.contains(g.getID())){
				subset.add(g);
			}
			else{
				remove.add(g.getID());
			}
		}
		
		//remove genosamples that do not have expsample
		//System.out.println("SNPs before: "+s.getGenosamples().size());
		for(String id:remove){
			s.removeGenoSample(id);
		}
		
		//System.out.println("SNPs after: "+s.getGenosamples().size());
		return subset;
	}
	
	


	private List<Double> calcPValues(List<Double> pointwisePValues, List<Double> permPValues) {
		if(permPValues.size()!=perm){
			System.out.println("Not all permutations completed");
			System.exit(1);
		}
		
		List<Double> pValues = new ArrayList<Double>();
		
		//iterate over pointwise p-values for SNPs
		for(int i=0; i<pointwisePValues.size(); i++){
			double p = pointwisePValues.get(i); //pointwise p-value for ith SNP in $snps
			int count = 0; //number of permuted p-values that are equal to or lower than $p
			for(int j=0; j<perm; j++){
				if(permPValues.get(j)<=p){
					count++;
				}
			}
			pValues.add(count*1.0/perm); //add permutation p-value for SNP i to $pValues
		}
		return pValues;
	}
	
	//test whether pvalue passes genome-wide significance threshold
	public boolean isSignificant(double z){
		if(z<=.0000025){
			return true;
		}
		return false;
	}



}
