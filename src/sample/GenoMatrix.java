package sample;

import java.util.List;

import genome.SNP;

public class GenoMatrix {
	int[][] genotypes;
	String[] sampleids;
	String[] snpids;
	
	public GenoMatrix(List<SNP> snps, List<String> samples){
		sampleids = new String[samples.size()];
		for(int j=0; j<samples.size(); j++){
			sampleids[j] = samples.get(j);
		}
		
		snpids = new String[snps.size()];
		genotypes = new int[snps.size()][samples.size()];
		for(int i=0; i<snps.size(); i++){
			SNP s = snps.get(i);
			String sid = s.getId();
			snpids[i] = sid;
			
			List<GenoSample> genos = s.getGenosamples();
			int j=0; //index for sample
			for(GenoSample g:genos){
				String id = g.getID();
				if(samples.contains(id)){
					int isHetero = g.getHetero();
					genotypes[i][j] = isHetero;
					j++;
				}
			}
		}
	}
	
	public int[][] getGenotypes(){
		return genotypes;
	}
	
	public String[] getSampleids(){
		return sampleids;
	}
	
	public String[] getSnpids(){
		return snpids;
	}
}
