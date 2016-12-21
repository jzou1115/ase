package sample;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
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
	
	public GenoMatrix(int[][] geno, String[] sampleids2, String[] snpids2) {
		genotypes = geno;
		sampleids = sampleids2;
		snpids = snpids2;
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
	
	public void write(File out) throws IOException{
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(out)));
		String header= "";
		for(String sample: sampleids){
			header = header+ sample+"\t";
		}
		outfile.write(header.trim()+"\n");
		
		for(int i=0; i< snpids.length; i++){
			String snp = snpids[i];
			String line = snp+"\t";
			for(int j=0; j< genotypes[i].length; j++){
				line = line+genotypes[i][j]+"\t";
			}
			outfile.write(line.trim()+"\n");
		}
		
		outfile.close();
	}
}
