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
	private static final String TEST_TAG = "-t";
	private static final String SAMPLE_TAG = "-n";
	private static final String THRESHOLD_TAG = "-z";
	
	private File output = new File(DEFAULT_OUTPUT_DIR);
	private InputStream snps = System.in;
	private InputStream genes = System.in;
	//private InputStream map = System.in;
	private boolean help = false;
	private int perm;
	private int errors;
	private String test;
	private int sampleNum;
	private double threshold;
	private InputStream genotypes;
	private InputStream expressions;
	
	@Override
	public void parse(String ... args) throws Exception{
		for( int i = 0 ; i < args.length ; ++i ){
			String cur = args[i];
			switch(cur){
			case SNP_TAG: 
				assertNextArg(SNP_TAG, i, args);
				snps = parseSNPArg(args[++i]);
				break;
			case GENE_TAG: 
				assertNextArg(GENE_TAG, i, args);
				genes = parseGeneArg(args[++i]);
				break;
			case PERM_TAG: 
				assertNextArg(PERM_TAG, i, args);
				perm = parsePermArg(args[++i]);
				break;
			case ERROR_TAG: 
				assertNextArg(ERROR_TAG, i, args);
				errors = parseERRORArg(args[++i]);
				break;
			case SAMPLE_TAG: 
				assertNextArg(SAMPLE_TAG, i, args);
				sampleNum = parseSAMPLEArg(args[++i]);
				break;
			case TEST_TAG: 
				assertNextArg(TEST_TAG, i, args);
				test = parseTESTArg(args[++i]);
				break;
			case THRESHOLD_TAG: 
				assertNextArg(THRESHOLD_TAG, i, args);
				threshold = parseTHRESHOLDArg(args[++i]);
				break;
			case HELP_TAG:
				help = true;
				return;
			case OUTPUT_TAG:
				assertNextArg(OUTPUT_TAG, i, args);
				output = new File(args[++i]);
				break;
			case GENOTYPES_TAG:
				assertNextArg(GENOTYPES_TAG, i, args);
				genotypes = parseGenotypesArg(args[++i]);
				break;
			case EXPRESSIONS_TAG:
				assertNextArg(EXPRESSIONS_TAG, i, args);
				expressions = parseExpressionArg(args[++i]);
				break;
			default:
				throw new Exception("Unrecognized flag: "+cur);
			}
		}
	}
	
	private InputStream parseExpressionArg(String string) throws FileNotFoundException {
		return new BufferedInputStream( new FileInputStream(new File(string)));
	}

	private InputStream parseGenotypesArg(String string) throws FileNotFoundException {
		return new BufferedInputStream( new FileInputStream(new File(string)));
	}

	private double parseTHRESHOLDArg(String string) {
		return Double.parseDouble(string);
	}

	private int parseSAMPLEArg(String string) {
		return Integer.parseInt(string);
	}

	private String parseTESTArg(String string) {
		return string.trim();
	}

	private int parseERRORArg(String string) {
		return Integer.parseInt(string);
	}

	private int parsePermArg(String string) {
		return Integer.parseInt(string);
	}

	private void assertNextArg(String flag, int index, String ... args) throws Exception{
		if( index + 1 >= args.length){
			throw new Exception("Flag `"+flag+"' requires an argument.");
		}
	}
	
	@Override
	public void printHelp(PrintStream out){
		out.println("Usage:\nASEsimulator [-s] [-g] [-p] [-e] [-t] [-n] [-z] [-o] [-h] \n");
		out.println("Options:");
		out.println("-s\tSNP file");
		out.println("-g\tGene file");
		out.println("-a\tGenotype data");
		out.println("-b\tExpression data");
		out.println("-p\tNumber of permutations");
		out.println("-e\tMaximum number of errors allowed");
		out.println("-t\tTest gene ID");
		out.println("-n\tNumber of samples");
		out.println("-z\tSignificance threshold");
		out.println("-o\tOutput file");
		out.println("-h\thelp statement");
	}
	
	private InputStream parseSNPArg(String s) throws FileNotFoundException {
		return new BufferedInputStream( new FileInputStream(new File(s)));
	}
	
	private InputStream parseGeneArg(String s) throws FileNotFoundException {
		return new BufferedInputStream( new FileInputStream(new File(s)));
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
	public double getThreshold() {
		return threshold;
	}

	@Override
	public InputStream getGenotypeData() {
		return genotypes;
	}

	@Override
	public InputStream getExpressionData() {
		return expressions;
	}


	
}
