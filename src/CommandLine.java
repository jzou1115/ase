import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.io.PrintStream;

public class CommandLine implements CommandLineParams{
	
	private static final String DEFAULT_OUTPUT_DIR = "ase_output";
//	private static final String DEFAULT_PACKAGE = "";
	
	private static final String SNP_TAG = "-s";
	private static final String GENE_TAG = "-g";
//	private static final String PACKAGE_TAG = "-p";
	private static final String HELP_TAG = "-h";
	private static final String OUTPUT_TAG = "-o";

	//private String root = DEFAULT_PACKAGE;
	private File output = new File(DEFAULT_OUTPUT_DIR);
	//private InputStream grammar = System.in;
	private InputStream snps = System.in;
	private InputStream genes = System.in;
	//private InputStream map = System.in;
	private boolean help = false;
	
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
	//		case PACKAGE_TAG:
		//		assertNextArg(PACKAGE_TAG, i, args);
			//	root = args[++i];
				//break;
			case HELP_TAG:
				help = true;
				return;
			case OUTPUT_TAG:
				assertNextArg(OUTPUT_TAG, i, args);
				output = new File(args[++i]);
				break;
			default:
				throw new Exception("Unrecognized flag: "+cur);
			}
		}
	}
	
	private void assertNextArg(String flag, int index, String ... args) throws Exception{
		if( index + 1 >= args.length){
			throw new Exception("Flag `"+flag+"' requires an argument.");
		}
	}
	
	@Override
	public void printHelp(PrintStream out){
		out.println("Usage:\nASEsimulator [-s] [-g] [-p] [-h] [-o] \n");
		out.println("Options:");
		out.println("-s\tSNP file");
		out.println("-g\tGene file");
		out.println("-p\tPackage");
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
		return snps;
	}
	
	@Override
	public File getOutputDir() {
		if(!output.exists()){
			output.mkdirs();
		}
		return output;
	}

	//@Override
	//public String getRootPackage() {
	//	return root;
	//}

	@Override
	public boolean help() {
		return help;
	}


	
}
