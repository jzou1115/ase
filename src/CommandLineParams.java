
import java.io.File;
import java.io.InputStream;
import java.io.PrintStream;

public interface CommandLineParams {
	
	public void parse(String ... args) throws Exception;
	
	public void printHelp(PrintStream out);

	public File getOutputDir();
	public boolean help();
	public InputStream getGenotypeData();
	public InputStream getExpressionData();
	public String getFilename();
	public String getASEDataFormat();
	
}
