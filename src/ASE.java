import java.io.IOException;
import java.io.InputStream;

import functions.*;


public class ASE {
	
	private void createMap(InputStream snps, InputStream genes) throws IOException {
		GenesToSNP map = new GenesToSNP(snps, genes);
		map.write();
	}
	

	private void startSimulation(InputStream map2, InputStream genotypes,
			String gene2, double threshold, int error, int perm, int n) throws IOException {
		Simulation sim = new Simulation(map2, genotypes, gene2, threshold, error, perm, n);
		//sim.setTestGene();
		//sim.parseGenotypes();
		
	}
	

	private void startRun(InputStream map, InputStream genotypes,
			InputStream expressions, String gene, double threshold, int perm, int n) {
		// start run
		
	}
	
	private void generateCombinations(InputStream map, String gene,
			InputStream genotypes, InputStream expressions) {
		// TODO Auto-generated method stub
		
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
			if(snps!=null && genes!=null){
				a.createMap(snps, genes);	
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
			
			if(map!=null && genotypes!=null && gene!=null){
				a.startSimulation(map, genotypes, gene, threshold, error, perm, n);

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
			int perm = cmdArgs.getPermNum();
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
			InputStream expressions = cmdArgs.getExpressionData();
			
			if(map!=null && gene!=null && genotypes!=null && expressions!=null){
				a.generateCombinations(map, gene, genotypes, expressions);
			} else{
				cmdArgs.printHelp(System.err);
				System.exit(0);
			}
		}
		else{
			//a.parseExpressions(cmdArgs.getExpressionData());
			;
		}

		
	}



}
