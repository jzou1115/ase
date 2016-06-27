package functions;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;

import parse.ParseMatrix;

public class ComputeFDR {
	public ComputeFDR(InputStream genotypes, InputStream expressions, int geneNum, int perm, File outfile, String filename) throws IOException{
		int[][] geno = ParseMatrix.parseMatrix(genotypes);
		int[][] expr = ParseMatrix.parseMatrix(expressions);
		
		String geneid = "gene"+geneNum;
		if(filename==null){
			filename= geneid+"_mapase.txt";
		}
		
		int[] snp1 = geno[0];
		String[] sampleids = new String[snp1.length];
		for(int i=0; i<snp1.length; i++){
			sampleids[i] = "sample"+i;
		}
		
		String[] snpids = new String[geno.length];
		for(int i=0; i<geno.length; i++){
			snpids[i] = "snp"+i;
		}
		
		MapASE ase = new MapASE(geno, expr[geneNum-1], sampleids, snpids, geneid, perm, outfile, filename);
		ase.mapase();
	}
}
