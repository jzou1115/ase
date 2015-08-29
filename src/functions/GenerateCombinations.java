package functions;

import java.util.ArrayList;
import java.util.List;

import sample.GenoSample;

import genome.SNP;

public class GenerateCombinations {

	public List<SNP> getCombinations(List<SNP> snps){
		List<SNP> ret = new ArrayList<SNP>();
		for(int i=0; i<snps.size(); i++){
			for(int j=i+1; j<snps.size(); j++){
				SNP combinedSNP = combineGenotypes(snps.get(i),snps.get(j));
				ret.add(combinedSNP);
			}
		}		
		return ret;
	}

	private SNP combineGenotypes(SNP snp, SNP snp2) {
		String id = snp.getId()+"_"+snp2.getId();
		SNP ret= new SNP(id);

		List<GenoSample> g1 = snp.getGenosamples();
		List<GenoSample> g2 = snp2.getGenosamples();
		if(g1.size() != g2.size()){
			return null;
		}
		for(int i=0; i<g1.size(); i++){
			if(g1.get(i).getID().equals(g2.get(i).getID()))
			if(g1.get(i).getHetero()==1 && g2.get(i).getHetero()==1){
				ret.addSample(new GenoSample(g1.get(i).getID(),1));
			}
			else{
				ret.addSample(new GenoSample(g1.get(i).getID(),0));
			}
		}
		
		return ret;
	}

}
