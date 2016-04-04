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
	
	private void assignChromatin(File states, InputStream genesmap, InputStream variants, int p, String gene, File outdir, String filename) throws IOException{

		List<Gene> genes = ParseCausalVariants.readVariantGroup(variants);
		List<String> geneids = new ArrayList<String>();
		List<SNP> var = new ArrayList<SNP>();
		for(Gene g:genes){
			var.addAll(g.getSNPs());
			geneids.add(g.getId());
		}
	
		
		ParseMap map = new ParseMap();
		map.parseMap(genesmap, geneids);
		List<SNP> snps = map.getSNPs();
		
		/**
		Map<String,SNP> snpmap = parsemap.getSnpMap();
		List<SNP> var2 = new ArrayList<SNP>();
		for(SNP v: var){
			if(snpmap.containsKey(v.getId())){
				var2.add(snpmap.get(v.getId()));
			}
		}
		**/
		PermuteChromatin perm = new PermuteChromatin(states, snps, var, genes, outdir, filename);
		
		perm.testEnrichmentGene(p);
		
	}

	
	private void createMap(InputStream snps, InputStream genes, InputStream chrom, File outdir, String filename) throws IOException {
		GenesToSNP map;
		if(chrom==null){
			map = new GenesToSNP(snps, genes);
		}
		else{
			map = new GenesToSNP(snps, genes, chrom);
		}
		
		if(filename!=null){
			map.write(outdir, filename);	
		} else{
			filename = "genestosnps.txt";
			map.write(outdir, filename);	
		}
	}
	

	private void startSimulation(InputStream map2, InputStream genotypes,
			String gene2, int error, int perm, int n, File outdir, String filename) throws IOException {
		Simulation sim = new Simulation();
		sim.setTestGene(map2, gene2, genotypes);
		sim.startRun(error, perm,n, outdir);
		
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
	
	private void generateCombinations(InputStream map, String gene,
			InputStream genotypes, InputStream expression, int p, File out, String f) throws IOException {
		if(f!=null){
			Combinations combinations = new Combinations(map, gene, genotypes, expression, p, out, f);
		}
		else{
			f = gene+"_combASE.txt";
			Combinations combinations = new Combinations(map, gene, genotypes, expression, p, out, f);
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
			InputStream chrom = cmdArgs.getChromFile();
			File outdir = cmdArgs.getOutputDir();
			String filename = cmdArgs.getFilename();
			if(snps!=null && genes!=null){
				a.createMap(snps, genes, chrom, outdir, filename);	
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
			int error = cmdArgs.getErrorNum();
			int perm = cmdArgs.getPermNum();
			int n = cmdArgs.getSampleNum();
			File outdir = cmdArgs.getOutputDir();
			//InputStream samples = cmdArgs.getSamples();
			String filename = cmdArgs.getFilename();
			
			if(map!=null && genotypes!=null && gene!=null){
				a.startSimulation(map, genotypes, gene, error, perm, n, outdir, filename);

			} else{
				cmdArgs.printHelp(System.err);
				System.exit(0);
			}


		}
		/**
		else if(fcn.equals("mapase")){
			InputStream map = cmdArgs.getMap();
			InputStream genotypes = cmdArgs.getGenotypeData();
			InputStream expressions = cmdArgs.getExpressionData();
			String gene = cmdArgs.getTestGene();
			int error = cmdArgs.getErrorNum();
			int n = cmdArgs.getSampleNum();
			int perm = cmdArgs.getPermNum();
			File outdir = cmdArgs.getOutputDir();
			String filename = cmdArgs.getFilename();
			
			if(map!=null && genotypes!=null && gene!=null){
				a.mapASE(map, genotypes, expressions, gene, error, n, perm, outdir, filename);

			} else{
				cmdArgs.printHelp(System.err);
				System.exit(0);
			}
		}
		**/
		else if(fcn.equals("mapase")){
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
		
		else if(fcn.equals("chromatin")){
			File chrom = cmdArgs.getChrom();
			InputStream genesmap = cmdArgs.getMap();
			InputStream variants = cmdArgs.getVariants();
			String gene = cmdArgs.getTestGene();
			int p = cmdArgs.getPermNum();
			
			File outdir = cmdArgs.getOutputDir();
			String filename = cmdArgs.getFilename();
			
			if(genesmap!=null && chrom!=null){
				a.assignChromatin(chrom,genesmap, variants, p, gene, outdir, filename);
			}
			
		}
		else if(fcn.equals("combinations")){
			InputStream map = cmdArgs.getMap();
			InputStream genotypes = cmdArgs.getGenotypeData();
			InputStream expression = cmdArgs.getExpressionData();
			String gene = cmdArgs.getTestGene();
			int p = cmdArgs.getPermNum();
			File out = cmdArgs.getOutputDir();
			String f = cmdArgs.getFilename();
			
			if(map!=null && gene!=null && genotypes!=null && expression!=null){
				a.generateCombinations(map, gene, genotypes, expression, p, out, f);
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
