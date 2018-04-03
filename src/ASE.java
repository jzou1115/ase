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

		if(cmdArgs.isMatrixFormat()){
			InputStream genotypes = cmdArgs.getGenotypeData();
			InputStream expressions = cmdArgs.getExpressionData();

			File outdir = cmdArgs.getOutputDir();
			String filename = cmdArgs.getFilename();
			
			a.mapase(genotypes, expressions, outdir, filename);
		}


		else{
			System.out.println("Function not recognized");
			cmdArgs.printHelp(System.err);
			System.exit(0);
		}

		
	}

	//for matrix format
	private void mapase(InputStream genotypes, InputStream expressions, File outdir, String filename) {
		try {
			MapASE ase = new MapASE(genotypes, expressions, outdir, filename);
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}







}
