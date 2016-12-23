
import java.io.File;
import java.io.InputStream;
import java.io.PrintStream;

public interface CommandLineParams {
	
	public void parse(String ... args) throws Exception;
	
	public void printHelp(PrintStream out);


	public File getOutputDir();
	public boolean help();
	public String getGene();
	public InputStream getGenotypeData();
	public InputStream getExpressionData();
	public InputStream getMap();
	public String getFunction();
	public String getFilename();
	public int getPermNum();

	public int getSampleSize();
}
