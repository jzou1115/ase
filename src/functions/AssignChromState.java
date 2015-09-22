package functions;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import genome.ChromState;
import genome.GenomicRegion;
import genome.SNP;
import genome.GenomicCoordinate;

public class AssignChromState {
	
	public static Map<SNP, ChromState> assignStateSNP(List<SNP> snps, List<ChromState> chrom){
		Map<SNP, ChromState> map = new HashMap<SNP, ChromState>();
		int ind=0;
		
		Collections.sort(snps);
		Collections.sort(chrom);
		
		for(ChromState c:chrom){
			//System.out.println(c.getState()+"\t"+chrom.size());
			GenomicRegion region = c.getRegion();
			GenomicCoordinate end = region.getEnd();
			
			SNP snp = snps.get(ind);
			//System.out.println(snp.getId()+"\t"+snps.size());
			if(region.contains(snp.getLocation())){
				while(region.contains(snp.getLocation())){
					map.put(snp, c);
					snp.setChromState(c);
					//System.out.println(snp.getId());
					if(ind>snps.size()-1){
						ind++;
						snp = snps.get(ind);	
					}
				}	
			}
		}
		
		System.out.println("Chromatin state assigned for snp index: "+ind);
		
		return map;
	}
}
