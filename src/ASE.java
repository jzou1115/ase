import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import gene.*;
import genome.*;
import run.*;
import sample.*;
import snp.*;

public class ASE {
	SNPgroup isHetero;
	GeneGroup hasASE;
	
	public void parseSnps(FileInputStream snpData){
		isHetero = SNPgroup.readSNPGroup(snpData);
	}
	
	public void parseGenes(FileInputStream geneData){
		hasASE = GeneGroup.readGeneGroup(geneData);
	}
	
	public void parseGenotypes(FileInputStream genotypes) throws IOException{
		BufferedReader br = new BufferedReader(new InputStreamReader(genotypes));
		String line = br.readLine();
		
		String[] sampleNames = line.split("\\s+");
		
		
	}
	
	public void parseExpressions(FileInputStream expressions){
		
	}
	
	public static void main(String args[]) throws FileNotFoundException{
		ASE a= new ASE();
		
		//String snpData = args[0];
		FileInputStream snpData = new FileInputStream(new File("/home/david/Documents/Jennifer/workspace/ase/test/snp.map"));
		a.parseSnps(snpData);
		
		//String geneData = args[1];
		FileInputStream geneData = new FileInputStream(new File("/home/david/Documents/Jennifer/workspace/ase/test/geneLoc.txt"));
		a.parseGenes(geneData);
		
		//String genotypeData = args[2];
		FileInputStream genotypeData = new FileInputStream(new File("/home/david/Documents/Jennifer/workspace/ase/test/isHetero.txt"));
		
		//String expData = args[2];
		FileInputStream expData = new FileInputStream(new File("/home/david/Documents/Jennifer/workspace/ase/test/hasASE.txt"));
	}
}
