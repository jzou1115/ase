package functions;

import java.util.Collections;
import java.util.List;

import genome.ChromState;
import genome.GenomicRegion;
import genome.SNP;
import genome.GenomicCoordinate;

public class AssignChromState {
	
	public static void assignStateSNP(List<SNP> snps, List<ChromState> chrom){
		int ind=0;
		
		Collections.sort(snps);
		Collections.sort(chrom);
		
		for(ChromState c:chrom){
			GenomicRegion region = c.getRegion();
			GenomicCoordinate end = region.getEnd();
			
			SNP snp = snps.get(ind);
			
			if(region.contains(snp.getLocation())){
				while(region.contains(snp.getLocation())){
					snp.setChromState(c);
					ind++;
					snp = snps.get(ind);
				}	
			}
		}
		
		System.out.println("Chromatin state assigned for snp index: "+ind);
	}
}