import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

//import parse.ParseSNP;
import functions.*;
import genome.ChromState;
import genome.Gene;
import genome.SNP;
import parse.ParseCausalVariants;
import parse.ParseChromState;
import parse.ParseMap;
import parse.ParseSNP;


public class ASE {
	
	void mapase(InputStream map, InputStream genotypes,
			InputStream expressions, String gene, int perm, File outdir, String filename) throws IOException {
		if(filename==null){
			MapASE ase = new MapASE( map,  genotypes, expressions,  gene,  perm,  outdir,  gene+"_mapase");
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
		

		if(fcn.equals("mapase")){
			InputStream map = cmdArgs.getMap();
			InputStream genotypes = cmdArgs.getGenotypeData();
			InputStream expressions = cmdArgs.getExpressionData();
			String gene = cmdArgs.getTestGene();
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
