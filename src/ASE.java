import java.util.HashMap;

import gene.*;
import genome.*;
import run.*;
import sample.*;
import snp.*;

public class ASE {
	
	int[][] hasASE;
	int[][] isHetero;
	int threshold;
	HashMap<String, String> map;
	HashMap<Integer, String> snpName;
	HashMap<String, Integer> geneNum;
	int numIncorrect[];
	int numSamples;
	
	public void setThreshold(int i){
		threshold=i;
	}
	
	public void parseASE(String s){
		
	}
	
	public void parseGeno(String s){
		
	}
	
	public void parseMap(String s){
		
	}
	
	public int findGene(String snp){
		String ens= map.get(snp);
		return geneNum.get(ens);	
	}
	
	public void findVariants(){
		for(int i=0; i<isHetero.length; i++){
			String snp= snpName.get(i);
			int gene= findGene(snp);
			
			int correct=0;
			int incorrect=0;
			for(int j=0; j<isHetero[0].length; j++){
				if(isHetero[i][j] == hasASE[gene][j]){
					correct++;
				}
				else{
					incorrect++;
				}
			}
			if(correct+incorrect!=numSamples){
				throw new RuntimeException();
			}
			
			numIncorrect[i]= incorrect;
		}
	}
	
	public void printVariants(){
		for(int i=0; i<numIncorrect.length;i++){
			if(i<=threshold){
				System.out.println(snpName.get(i));
			}
		}
	}
	
	
	
	public static void main(String args[]){
		ASE a= new ASE();
		
		String aseFile= args[0];
		a.parseASE(aseFile);
		String genoFile= args[1];
		a.parseGeno(genoFile);
		
		a.setThreshold(Integer.parseInt(args[2]));
		
		String map= args[3];
		a.parseMap(map);
		
		a.findVariants();
		
		a.printVariants();
	}
}
