import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.util.List;
import java.util.Map;

//import parse.ParseSNP;
import functions.*;
import genome.ChromState;
import genome.SNP;
import parse.ParseChromState;
import parse.ParseMap;
import parse.ParseSNP;


public class ASE {
	
	private void assignChromatin(InputStream states, InputStream variants, InputStream genesmap, String gene, File outdir, String filename) throws IOException{
		ParseMap parsemap = new ParseMap();
		parsemap.parseMap(genesmap, gene);
		List<SNP> s = parsemap.getSNPs();
		
		List<SNP> var = ParseSNP.readVariantGroup(variants);
		
		List<ChromState> chrom = ParseChromState.parseChromState(states);
		Map<SNP, ChromState> map = AssignChromState.assignStateSNP(s, chrom);
		
		if(filename==null){
			filename = "assignChromatinOutput.txt";
		}
		BufferedWriter outfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+filename)));

		for(SNP snp: s){
			outfile.write(snp.toString()+"\n");
		}
		outfile.close();
		
		BufferedWriter outfile2 = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outdir+File.separator+filename+"_variants.txt")));

		for(SNP snp:var){
			String line = snp.toString()+"\t"+map.get(snp).getState()+"\n";
			outfile2.write(line);
		}
		outfile2.close();
	}

	
	private void createMap(InputStream snps, InputStream genes, File outdir, String filename) throws IOException {
		GenesToSNP map = new GenesToSNP(snps, genes);
		if(filename!=null){
			map.write(outdir, filename);	
		} else{
			filename = "genestosnps.txt";
			map.write(outdir, filename);	
		}
	}
	

	private void startSimulation(InputStream map2, InputStream genotypes,
			String gene2, double threshold, int error, int perm, int n, File outdir, InputStream samples, int splits) throws IOException {
		Simulation sim = new Simulation();
		sim.setTestGene(map2, gene2, genotypes, samples);
		sim.startRun(threshold, error, perm, n, splits, outdir);
		
	}
	
	

	private void mapASE(InputStream map, InputStream genotypes,
			InputStream expressions, String gene, double threshold, int error, int n, File outdir, InputStream samples) throws IOException {
		MapASE ase = new MapASE();
		ase.setTestGene(map, gene, genotypes, samples);
		ase.parseGenotypes(genotypes);
		ase.parseExpressions(expressions, outdir, samples);
		ase.startRun(threshold, error, n, outdir);
		
	}
	
	private void generateCombinations(InputStream map, String gene,
			InputStream genotypes, double threshold, int error, int perm, int n, File outdir, String filename, InputStream samples) throws IOException {
		if(filename!=null){
			Combinations combinations = new Combinations(map, gene, genotypes, outdir, samples);
			combinations.write(filename);
			combinations.simulate(threshold, error, perm, n);
		}
		else{
			filename = "combinations.txt";
			Combinations combinations = new Combinations(map, gene, genotypes, outdir, samples);
			combinations.write(filename);
			combinations.simulate(threshold, error, perm, n);
		}
	}
	
	
	public static void main(String args[]) throws IOException{
		
		CommandLineParams cmdArgs = new CommandLine();
		try{
			cmdArgs.parse(args);
		} catch (Exception e) {
			System.err.println(e.getMessage());
			cmdArgs.printHelp(System.err);
			System.exit(1);
		}
		if(cmdArgs.help()){
			cmdArgs.printHelp(System.err);
			System.exit(0);
		}
		
		ASE a= new ASE();
		String fcn = cmdArgs.getFunction();
		
		if(fcn.equals("genestosnps")){
			InputStream snps = cmdArgs.getSNPsInput();
			InputStream genes = cmdArgs.getGenesInput();
			File outdir = cmdArgs.getOutputDir();
			String filename = cmdArgs.getFilename();
			if(snps!=null && genes!=null){
				a.createMap(snps, genes, outdir, filename);	
			}	
			else{
				cmdArgs.printHelp(System.err);
				System.exit(0);
			}
		}
		
		else if(fcn.equals("simulation")){
			InputStream map = cmdArgs.getMap();
			InputStream genotypes = cmdArgs.getGenotypeData();
			String gene = cmdArgs.getTestGene();
			double threshold = cmdArgs.getThreshold();
			int error = cmdArgs.getErrorNum();
			int perm = cmdArgs.getPermNum();
			int n = cmdArgs.getSampleNum();
			File outdir = cmdArgs.getOutputDir();
			InputStream samples = cmdArgs.getSamples();
			int split = cmdArgs.getSplits();
			
			if(map!=null && genotypes!=null && gene!=null){
				a.startSimulation(map, genotypes, gene, threshold, error, perm, n, outdir, samples, split);

			} else{
				cmdArgs.printHelp(System.err);
				System.exit(0);
			}


		}
		
		else if(fcn.equals("mapase")){
			InputStream map = cmdArgs.getMap();
			InputStream genotypes = cmdArgs.getGenotypeData();
			InputStream expressions = cmdArgs.getExpressionData();
			String gene = cmdArgs.getTestGene();
			double threshold = cmdArgs.getThreshold();
			int error = cmdArgs.getErrorNum();
			int n = cmdArgs.getSampleNum();
			File outdir = cmdArgs.getOutputDir();
			InputStream samples = cmdArgs.getSamples();
			
			if(map!=null && genotypes!=null && gene!=null){
				a.mapASE(map, genotypes, expressions, gene, threshold, error, n, outdir, samples);

			} else{
				cmdArgs.printHelp(System.err);
				System.exit(0);
			}
		}
		
		else if(fcn.equals("chromatin")){
			InputStream chrom = cmdArgs.getChrom();
			InputStream variants = cmdArgs.getVariants();
			InputStream genesmap = cmdArgs.getMap();
			String gene = cmdArgs.getTestGene();
			
			File outdir = cmdArgs.getOutputDir();
			String filename = cmdArgs.getFilename();
			
			if(genesmap!=null && chrom!=null){
				a.assignChromatin(chrom, variants,genesmap, gene, outdir, filename);
			}
			
		}
		/**
		else if(fcn.equals("combinations")){
			InputStream map = cmdArgs.getMap();
			String gene = cmdArgs.getTestGene();
			InputStream genotypes = cmdArgs.getGenotypeData();
			File outdir = cmdArgs.getOutputDir();
			String outfile = cmdArgs.getFilename();
			//InputStream expressions = cmdArgs.getExpressionData();
			 InputStream samples = cmdArgs.getSamples();
			
			double threshold = cmdArgs.getThreshold();
			int error = cmdArgs.getErrorNum();
			int perm = cmdArgs.getPermNum();
			int n = cmdArgs.getSampleNum();
			
			if(map!=null && gene!=null && genotypes!=null){
				a.generateCombinations(map, gene, genotypes, threshold, error, perm, n, outdir, outfile, samples);
			} else{
				cmdArgs.printHelp(System.err);
				System.exit(0);
			}
			
		}
		**/
		else{
			System.out.println("Function not recognized");
			cmdArgs.printHelp(System.err);
			System.exit(0);
		}

		
	}





}
