import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.io.PrintStream;

public class CommandLine implements CommandLineParams{
	
	private static final String DEFAULT_OUTPUT_DIR = "ase_output";
	
	private static final String SNP_TAG = "-s";
	private static final String GENE_TAG = "-g";
	private static final String GENOTYPES_TAG = "-a";
	private static final String EXPRESSIONS_TAG = "-b";
	private static final String PERM_TAG = "-p";
	private static final String ERROR_TAG = "-e";
	private static final String HELP_TAG = "-h";
	private static final String OUTPUT_TAG = "-o";
	private static final String FILE_TAG = "-f";
	private static final String TEST_TAG = "-t";
	private static final String SAMPLE_TAG = "-n";
	private static final String MAP_TAG = "-m";
	private static final String CHROM_TAG = "-c";
	private static final String VAR_TAG = "-v";
	private static final String TSS_TAG ="-z";
	private static final String MAP_FCN = "genestosnps";
	private static final String SIM_FCN = "simulation";
	private static final String ASE_FCN = "mapase";
	private static final String COMB_FCN = "combinations";
	private static final String CHROM_FCN = "chromatin";
	
	
	
	
	private InputStream snps;
	private InputStream genes;
	private InputStream genotypes;
	private InputStream expressions;
	private InputStream map;
	private InputStream samples;
	private InputStream chrom;
	private InputStream variant;
	private InputStream tss;
	private boolean help = false;
	private int perm;
	private int errors;
	private String test;
	private int sampleNum;
	private String function;
	private File output = new File(DEFAULT_OUTPUT_DIR);
	private String outfile;
	
	@Override
	public void parse(String ... args) throws Exception{
		if(args.length==0){
			printHelp(System.err);
			System.exit(0);
		}
		for( int i = 0 ; i < args.length ; ++i ){
			String cur = args[i];
			switch(cur){
			case MAP_FCN: 
				function = MAP_FCN;
				break;
			case SIM_FCN: 
				function = SIM_FCN;
				break;
			case ASE_FCN:
				function = ASE_FCN;
				break;
			case COMB_FCN:
				function = COMB_FCN;
				break;
			case CHROM_FCN:
				function = CHROM_FCN;
				break;
			case SNP_TAG: 
				assertNextArg(SNP_TAG, i, args);
				snps = parseDataFile(args[++i]);
				break;
			case GENE_TAG: 
				assertNextArg(GENE_TAG, i, args);
				genes = parseDataFile(args[++i]);
				break;
			case MAP_TAG: 
				assertNextArg(MAP_TAG, i, args);
				map = parseDataFile(args[++i]);
				break;
			case PERM_TAG: 
				assertNextArg(PERM_TAG, i, args);
				perm = parseIntArg(args[++i]);
				break;
			case ERROR_TAG: 
				assertNextArg(ERROR_TAG, i, args);
				errors = parseIntArg(args[++i]);
				break;
			case SAMPLE_TAG: 
				assertNextArg(SAMPLE_TAG, i, args);
				sampleNum = parseIntArg(args[++i]);
				break;
			case TEST_TAG: 
				assertNextArg(TEST_TAG, i, args);
				test = parseStringArg(args[++i]);
				break;
			case HELP_TAG:
				help = true;
				return;
			case OUTPUT_TAG:
				assertNextArg(OUTPUT_TAG, i, args);
				output = new File(args[++i]);
				break;
			case FILE_TAG:
				assertNextArg(FILE_TAG, i, args);
				outfile = parseStringArg(args[++i]);
				break;
			case GENOTYPES_TAG:
				assertNextArg(GENOTYPES_TAG, i, args);
				genotypes = parseDataFile(args[++i]);
				break;
			case EXPRESSIONS_TAG:
				assertNextArg(EXPRESSIONS_TAG, i, args);
				expressions = parseDataFile(args[++i]);
				break;
			case CHROM_TAG:
				assertNextArg(CHROM_TAG, i, args);
				chrom = parseDataFile(args[++i]);
				break;
			case VAR_TAG:
				assertNextArg(VAR_TAG,i,args);
				variant = parseDataFile(args[++i]);
				break;
			case TSS_TAG:
				assertNextArg(TSS_TAG,i,args);
				tss = parseDataFile(args[++i]);
			default:
				throw new Exception("Unrecognized flag: "+cur);
			}
		}
	}
	
	private InputStream parseDataFile(String string) throws FileNotFoundException {
		return new BufferedInputStream( new FileInputStream(new File(string)));
	}


	private String parseStringArg(String string) {
		return string.trim();
	}

	private int parseIntArg(String string) {
		return Integer.parseInt(string);
	}

	private void assertNextArg(String flag, int index, String ... args) throws Exception{
		if( index + 1 >= args.length){
			throw new Exception("Flag `"+flag+"' requires an argument.");
		}
	}
	
	@Override
	public void printHelp(PrintStream out){
		out.println("Usage: ase <function> [<args>] \n");
		out.println("Functions:");
		out.println("<genestosnps>\tThis function takes a list of genes and a list of snps.  It creates a map from each gene to a list of snps.");
		out.println("<simulation>\tThis function starts a simulation for one gene.");
		out.println("<mapase>\tThis function maps variants to ASE.");
		out.println("\nOptions:");
		out.println("-s\tSNP Locations");
		out.println("-g\tGene Locations");
		out.println("-m\tMap from gene to SNPs");
		out.println("-a\tGenotype file");
		out.println("-b\tExpression file");
		out.println("-p\tNumber of permutations");
		out.println("-e\tMaximum number of errors allowed");
		out.println("-v\tVariants");
		out.println("-n\tNumber of samples used");
		out.println("-o\tOutput directory");
		out.println("-f\tFilename");
		out.println("-h\thelp statement");

	}
	

	@Override
	public InputStream getSNPsInput() {
		return snps;
	}
	
	@Override
	public InputStream getGenesInput() {
		return genes;
	}
	
	@Override
	public File getOutputDir() {
		if(!output.exists()){
			output.mkdirs();
		}
		return output;
	}
	
	@Override
	public boolean help() {
		return help;
	}

	@Override
	public int getPermNum() {
		return perm;
	}

	@Override
	public int getErrorNum() {
		return errors;
	}

	@Override
	public String getTestGene() {
		return test;
	}

	@Override
	public int getSampleNum() {
		return sampleNum;
	}


	@Override
	public InputStream getGenotypeData() {
		return genotypes;
	}

	@Override
	public InputStream getExpressionData() {
		return expressions;
	}

	@Override
	public InputStream getMap() {
		return map;
	}

	@Override
	public String getFunction() {
		return function;
	}

	@Override
	public String getFilename() {
		return outfile;
	}

	@Override
	public InputStream getSamples() {
		return samples;
	}


	@Override
	public InputStream getChrom() {
		return chrom;
	}

	@Override
	public InputStream getVariants() {
		return variant;
	}

	@Override
	public InputStream getTss() {
		return tss;
	}



	
}
