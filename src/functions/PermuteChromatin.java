package functions;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import genome.ChromState;
import genome.GenomicCoordinate;
import genome.GenomicRegion;
import genome.SNP;

import org.apache.commons.math3.util.CombinatoricsUtils;// MathUtils.binomialCoefficient;

public class PermuteChromatin {
	
	List<ChromState> chromatin;
	List<SNP> snps;
	List<SNP> variants;
	
	GenomicRegion chromosome;
	int chromNum;
	long min;
	long max;

	int size;
	File outdir;
	public PermuteChromatin(List<ChromState> chrom, List<SNP> snp, List<SNP> variant, File out){
		chromatin = chrom;
		snps = snp;
		variants = variant;
		
		chromNum = snps.get(0).getLocation().getChromosome();
		min = chromatin.get(0).getRegion().getStart().getCoord();
		max = chromatin.get(chromatin.size()-1).getRegion().getEnd().getCoord();
		
		chromosome = new GenomicRegion(new GenomicCoordinate(chromNum,min), new GenomicCoordinate(chromNum,max));
	
		outdir=out;
	}
	
	public void filterChromatin(){
		List<ChromState> chrom = new ArrayList<ChromState>();
		for(ChromState c:chromatin){
			if(c.getRegion().getChromosome()==chromNum){
				chrom.add(c);
			}
		}
		//System.out.println("filter chromatin: "+chromatin.size()+"\t"+chrom.size());
		chromatin = chrom;
		size = chrom.size();
	}
	
	public int getSize(){
		return size;
	}
	public void permute(int n) throws IOException{
	//	System.out.println("Starting permutation with "+chromatin.size());
	//	BufferedWriter file = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+"permutation"+n)));
		
		List<ChromState> chromatin2 = new ArrayList<ChromState>();
		for(int i=0; i<chromatin.size();i++){
			int ind;
			if((i+n)<chromatin.size()){
				ind = i+n;
			}
			else{
				//System.out.println(i+"+"+n+"="+(i+n));
				int offset = (i+n)-chromatin.size();
				//System.out.println(offset);
				ind = offset;
			}
			
			ChromState c = chromatin.get(i);
			ChromState c2 = chromatin.get(ind);
			chromatin2.add(new ChromState(c2.getState(), c.copyRegion()));
			
		}
		
		if(chromatin2.size()!=chromatin.size()){
			System.out.print("errorerrorerror");
		}

		chromatin = chromatin2;
		Collections.sort(chromatin);
		/**
		for(ChromState c:chromatin){
			file.write(c.toString()+"\n");
		}
		file.close();
**/
	}
	
	public int mapSnps(String state){
		Map<SNP,ChromState> map = AssignChromState.assignStateSNP(snps, chromatin);
		int num=0;
		for(SNP s:snps){
			if(map.containsKey(s)){
				//System.out.print(s.getChromState()+"_"+state);
				if(s.getChromState().equals(state) || s.getChromState().equals("5_TxWk") || s.getChromState().equals("6_EnhG") || s.getChromState().equals("7_Enh")){
					//System.out.print("YES\n");
					num++;
				}
				else{
					//System.out.print("\n");
				}
			}
		}
		
		return num;
	}
	
	public int mapVar(String state){
		Map<SNP,ChromState> map = AssignChromState.assignStateSNP(variants, chromatin);
		int num=0;
		for(SNP s:variants){
			if(map.containsKey(s)){
				//System.out.print(s.getChromState()+"_"+state);
				if(s.getChromState().equals(state)){
				//	System.out.print("YES\n");
					num++;
				}
			}
		}
		return num;
		
	}

	public double calculateT(int snpsSize, int varSize, int s, int v) {
		//System.out.println(snpsSize+"\t"+varSize+"\t"+s+"\t"+v);
		double a = CombinatoricsUtils.binomialCoefficientDouble(snpsSize,s);
		double b = CombinatoricsUtils.binomialCoefficientDouble(varSize,v);
		double c = CombinatoricsUtils.binomialCoefficientDouble(snpsSize+varSize, s+v);
		
		return 1.0*a*b/c;
	}
	
	public double calculateT(int total, int s){
		double a = CombinatoricsUtils.binomialCoefficientDouble(total,s);
		double p = 1.0*s/total;
		double q  = 1-p;
		int not = total-s;
		
		return a*Math.pow(p, s)*Math.pow(q,not);
	}

	public int getVarSize() {
		// TODO Auto-generated method stub
		return variants.size();
	}

	public int getSnpsSize() {
		// TODO Auto-generated method stub
		return snps.size();
	}
	
}
