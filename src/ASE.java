import java.io.File;
import java.io.IOException;
import java.io.InputStream;

import parse.ParseSNP;
import functions.*;


public class ASE {
	
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
			String gene2, double threshold, int error, int perm, int n, File outdir) throws IOException {
		Simulation sim = new Simulation();
		sim.setTestGene(map2, gene2, genotypes);
		sim.startRun(threshold, error, perm, n, outdir);
		
	}
	
	

	private void startRun(InputStream map, InputStream genotypes,
			InputStream expressions, String gene, double threshold, int perm, int n) {
		// start run
		
	}
	
	private void generateCombinations(InputStream map, String gene,
			InputStream genotypes, double threshold, int error, int perm, int n, File outdir, String filename) throws IOException {
		if(filename!=null){
			Combinations combinations = new Combinations(map, gene, genotypes, outdir);
			combinations.write(filename);
			combinations.simulate(threshold, error, perm, n);
		}
		else{
			filename = "combinations.txt";
			Combinations combinations = new Combinations(map, gene, genotypes, outdir);
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
			
			if(map!=null && genotypes!=null && gene!=null){
				a.startSimulation(map, genotypes, gene, threshold, error, perm, n, outdir);

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
			
			if(map!=null && genotypes!=null && gene!=null){
				a.startRun(map, genotypes, expressions, gene, threshold, error, n);

			} else{
				cmdArgs.printHelp(System.err);
				System.exit(0);
			}
		}
		
		else if(fcn.equals("combinations")){
			InputStream map = cmdArgs.getMap();
			String gene = cmdArgs.getTestGene();
			InputStream genotypes = cmdArgs.getGenotypeData();
			File outdir = cmdArgs.getOutputDir();
			String outfile = cmdArgs.getFilename();
			//InputStream expressions = cmdArgs.getExpressionData();
			
			double threshold = cmdArgs.getThreshold();
			int error = cmdArgs.getErrorNum();
			int perm = cmdArgs.getPermNum();
			int n = cmdArgs.getSampleNum();
			
			if(map!=null && gene!=null && genotypes!=null){
				a.generateCombinations(map, gene, genotypes, threshold, error, perm, n, outdir, outfile);
			} else{
				cmdArgs.printHelp(System.err);
				System.exit(0);
			}
			
		}
		else{
			System.out.println("Function not recognized");
			cmdArgs.printHelp(System.err);
			System.exit(0);
		}

		
	}





}
