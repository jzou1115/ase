package functions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;



import parse.ParseMatrix;


public class MapASE {
	int[][] genotypes;
	int[][] expressions;
	
	String[] snpids;
	int numSNPs;
	int numSamples;
	
	

	double[] pointwise;
	
	File outdir;
	String filename;
	
	
	public MapASE(InputStream genotypesInput,
			InputStream expressionsInput, File o, String f) throws IOException{

		outdir = o;
		filename=f;
		
		//parse genotype and ASE data
		ParseMatrix parseGeno = new ParseMatrix(genotypesInput);
		genotypes = parseGeno.getMat();
		snpids = parseGeno.getRowIDs();
		
		ParseMatrix parseASE = new ParseMatrix(expressionsInput);
		expressions = parseASE.getMat();
		
		int numGenoSamples = parseGeno.getNumCols();
		int numASESamples = parseASE.getNumCols();

		if(numGenoSamples==numASESamples){
			numSamples = numGenoSamples;
		}
		else{
			System.out.println("Unequal number of samples in genotype and ASE input");
			System.exit(1);
		}
		
		int numGenoSNPs = parseGeno.getNumRows();
		int numExpSNPs = parseASE.getNumRows();
		if(numGenoSNPs == numExpSNPs){
			numSNPs = numGenoSNPs;
		}
		else{
			System.out.println("Unequal number of SNPs in genotype and ASE input");
			System.exit(1);
		}
		

	}

	
	public void mapase() throws IOException {
	
		pointwise = new double[numSNPs]; //array to store pointwise p-values
		String[] line = pointwisePValue(); //populate pointwise array and return array of output for each snp

		
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+filename)));
		for(int i=0; i<numSNPs; i++){
			outfile.write(line[i]+"\n");
		}
		outfile.close();

	}



	//Calculates pointwise p-values and populates $realPValues and $lines
	public String[] pointwisePValue(){
		String[] line=new String[numSNPs];
		for(int i=0; i<numSNPs; i++){
			//number of ones in ase samples
			int m=0;
			//number of ones in genotype samples
			int k=0;
			
			//number of ASE and heterozygous
			int a = 0;
			//number of not ASE and heterozygous
			int b = 0;
			
			for (int j=0; j<numSamples;j++) {
				int hasASE = expressions[j][j];
				int isHetero = genotypes[i][j];
				
				if(hasASE==1){
					m=m+1; //increment count of number of individuals with ASE
				}
				if(isHetero==1){
					k=k+1; //increment count of number of individuals that are heterozygous
				}


				if(hasASE==1 && isHetero==1){
					a=a+1;
				}
				
				if(hasASE==0 && isHetero==1){
					b=b+1;
				}
				
			}

			//proportion of heterozygous individuals w/i ASE subset
			double p1;
			if(m==0){
				p1=0;
			}
			else{
				p1 = a*1.0/m;	
			}
			int n1 = m*2;
			
			//proportion of heterozygous individuals w/i balanced subset
			double p2;
			if(numSamples-m==0){
				p2=0;
			}
			else{
				p2 = b*1.0/(numSamples - m);	
			}
			int n2 = 2*(numSamples-m);
			
			//calculate statistic
			double p3 = (p1+p2)/2;
			if(p3==0){
				p3=.0000001;
			}
			if(p3==1){
				p3=.9999999;
			}
			double s = (p1 - p2) / Math.sqrt((p3*(1-p3))*(n1+n2)/(n1*n2));
			
			//SNPID and ASE summary statistic
			line[i] = snpids[i]+"\t"+s;
		}
		
		return line;

	}
	
	
	


}
