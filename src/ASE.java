import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import functions.*;


public class ASE {
	
	

	private void startSimulation(InputStream map, InputStream genotypes,
			String gene, int n, File outdir, String filename) throws IOException {
		if(filename==null){
			filename=gene+"_simulation.txt";
		}
		Simulation sim = new Simulation();
		sim.setTestGene(map, gene, genotypes, outdir);
		sim.mapASE(filename, n);
	}
	
	void mapase(InputStream map, InputStream genotypes,
			InputStream expressions, String gene, int perm, File outdir, String filename) throws IOException {
		if(filename==null){
			MapASE ase = new MapASE( map,  genotypes, expressions,  gene,  perm,  outdir,  gene+"_mapase.txt");
			ase.mapase();
		}
		else{
			MapASE ase = new MapASE( map,  genotypes, expressions,  gene,  perm,  outdir,  filename);
			ase.mapase();
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
		
	
		if(fcn.equals("simulation")){
			InputStream map = cmdArgs.getMap();
			InputStream genotypes = cmdArgs.getGenotypeData();
			String gene = cmdArgs.getGene();
			int sampleSize = cmdArgs.getSampleSize();
			File outdir = cmdArgs.getOutputDir();
			String filename = cmdArgs.getFilename();
			
			if(map!=null && genotypes!=null && gene!=null){
				a.startSimulation(map, genotypes, gene, sampleSize, outdir, filename);

			} else{
				cmdArgs.printHelp(System.err);
				System.exit(0);
			}


		}

		else if(fcn.equals("mapase")){
			InputStream map = cmdArgs.getMap();
			InputStream genotypes = cmdArgs.getGenotypeData();
			InputStream expressions = cmdArgs.getExpressionData();
			String gene = cmdArgs.getGene();
			int perm = cmdArgs.getPermNum();
			File outdir = cmdArgs.getOutputDir();
			String filename = cmdArgs.getFilename();
			
			if(map!=null && genotypes!=null && expressions!=null && gene!=null){
				a.mapase(map, genotypes, expressions, gene, perm, outdir, filename);

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
