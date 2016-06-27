package parse;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;


public class ParseMatrix {
	public static int[][] parseMatrix(InputStream genotypes) throws IOException{
		BufferedReader br = new BufferedReader(new InputStreamReader(genotypes));
		List<String> lines = new ArrayList<String>();
	
		String line;

		try {
			while((line = br.readLine()) != null){
				lines.add(line);
			}

			br.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		String[] linetokens = lines.get(0).split("\\s+");
		System.out.println("Reading matrix: "+lines.size()+"x"+linetokens.length);
		int[][] ret = new int[lines.size()][linetokens.length];
		for(int i=0; i<lines.size(); i++){
			String l = lines.get(i);
			String[] tokens = l.split("\\s+");
			for(int j=0; j<tokens.length;j++){
				ret[i][j] = Integer.parseInt(tokens[j]);
			}
		}
		
		return ret;
	}
}
