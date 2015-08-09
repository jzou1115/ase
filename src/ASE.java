import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import gene.*;
import genome.*;
import run.*;
import sample.*;
import snp.*;

public class ASE {
	SNPgroup snps;
	GeneGroup genes;
	
	public void parseSnps(FileInputStream snpData){
		snps = SNPgroup.readSNPGroup(snpData);
		snps.sort();
		System.out.println(snps.toString());
	}
	
	public void parseGenes(FileInputStream geneData){
		genes = GeneGroup.readGeneGroup(geneData);
		genes.sort();
		System.out.println(genes.toString());
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
