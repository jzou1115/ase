package sample;

import java.util.List;

public class ExpMatrix {
	int[] expressions;
	String[] sampleids;
	
	public ExpMatrix(List<ExpSample> expr, List<String> samples){
		sampleids = new String[samples.size()];
		for(int j=0; j<samples.size(); j++){
			sampleids[j] = samples.get(j);
		}
		
		expressions = new int[samples.size()];
		int i=0; //index for sample in expressions
		for(ExpSample e:expr){
			if(samples.contains(e.getID())){
				expressions[i] = e.getASE();
				i++;
			}
		}
		
	}
	
	public ExpMatrix(List<GenoSample> geno){
		sampleids = new String[geno.size()];
		expressions = new int[geno.size()];
		for(int i=0; i< geno.size(); i++){
			GenoSample g = geno.get(i);
			sampleids[i]=g.getID();
			expressions[i] = g.getHetero();
		}
		
	}
	
	public ExpMatrix(int[] expr, String[] sampleids2) {
		expressions = expr;
		sampleids = sampleids2;
	}

	public int[] getExpressions(){
		return expressions;
	}
	
	public String[] getSampleids(){
		return sampleids;
	}
}
