import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.io.PrintStream;

public class CommandLine implements CommandLineParams{
	
	private static final String DEFAULT_OUTPUT_DIR = "ase_output";
	
	private static final String GENOTYPES_TAG = "-a";
	private static final String ASE_TAG = "-b";
	private static final String HELP_TAG = "-h";
	private static final String OUTPUT_TAG = "-o";
	private static final String FILE_TAG = "-f";
	private static final String INPUT_TAG = "-i";
	private static final String GTEX_FORMAT = "gtex";
	private static final String MATRIX_FORMAT = "matrix";
	
	private InputStream genotypes;
	private InputStream expressions;
	private File output = new File(DEFAULT_OUTPUT_DIR);
	private String outfile;
	private Boolean help;
	private String inputFormat;
	
	@Override
	public void parse(String ... args) throws Exception{
		if(args.length==0){
			printHelp(System.err);
			System.exit(0);
		}
		for( int i = 0 ; i < args.length ; ++i ){
			String cur = args[i];
			switch(cur){
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
			case ASE_TAG:
				assertNextArg(ASE_TAG, i, args);
				expressions = parseDataFile(args[++i]);
				break;
			case INPUT_TAG:
				assertNextArg(INPUT_TAG, i, args);
				inputFormat = parseStringArg(args[++i]);
				break;
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
		out.println("Usage: ase.jar [<args>] \n");
		out.println("-a\tGenotype file");
		out.println("-b\tExpression file");
		out.println("-o\tOutput directory");
		out.println("-f\tFilename");
		out.println("-h\thelp statement");
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
	public InputStream getGenotypeData() {
		return genotypes;
	}

	@Override
	public InputStream getExpressionData() {
		return expressions;
	}

	@Override
	public String getFilename() {
		return outfile;
	}

	@Override
	public String getASEDataFormat(){
		if(inputFormat.equals(GTEX_FORMAT)){
			return GTEX_FORMAT;
		}
		else if(inputFormat.equals(MATRIX_FORMAT)){
			return MATRIX_FORMAT;
		}

		System.out.println(inputFormat);
		System.exit(1);
		
		return inputFormat;
	}



	
}
