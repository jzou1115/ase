package functions;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;

import parse.ParseMap;
import parse.ParseSNP;

import sample.GenoSample;

import genome.Gene;
import genome.SNP;

public class Combinations {
	Gene gene;
	List<SNP> snps;
	List<SNP> combs;
	List<String> gtexIds;
	
	public Combinations(InputStream map, String gene2,
			InputStream genotypes) throws IOException {
		setTestGene(map, gene2);
		ParseSNP.parseGenotypes(genotypes,snps);
		setSampIDs();
		getCombinations();
	}

	private void setSampIDs() {
		gtexIds = new ArrayList<String>();
		for(GenoSample g: snps.get(0).getGenosamples()){
			gtexIds.add(g.getID());
		}
	}

	public void setTestGene(InputStream map2, String gene2){
		ParseMap parsemap = new ParseMap();
		parsemap.parseMap(map2, gene2);
		gene = parsemap.getGene();
		snps = parsemap.getSNPs();
		//System.out.println(gene.toString()+"\t"+snps.size());
	}
	public List<SNP> getCombinations(){
		combs = new ArrayList<SNP>();
		for(int i=0; i<snps.size(); i++){
			for(int j=i+1; j<snps.size(); j++){
				SNP combinedSNP = combineGenotypes(snps.get(i),snps.get(j));
				combs.add(combinedSNP);
			}
		}		
		return combs;
	}
	
	public void write() throws IOException{
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("combined.txt")));
		
		String labels="";
		for(String id: gtexIds){
			labels = labels+"\t"+id;
		}
		
		outfile.write(labels.trim()+"\n");
		
		for(SNP s:combs){
			outfile.write(s.samplesToString()+"\n");
		}
		
		outfile.close();
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
			if(g1.get(i).getID().equals(g2.get(i).getID())){
				if(g1.get(i).getHetero()==1 && g2.get(i).getHetero()==1){
					ret.addSample(new GenoSample(g1.get(i).getID(),1));
				}
				else{
					ret.addSample(new GenoSample(g1.get(i).getID(),0));
				}	
			}
			else{
				System.out.println("Gtex ids not matching");
			}
		}
		
		return ret;
	}
	

	public void simulate(double threshold, int error, int perm, int n) throws IOException {
		Simulation sim = new Simulation();
		sim.setTestGene(gene, combs);
		sim.startRun(threshold, error, perm, n);
	}

}
