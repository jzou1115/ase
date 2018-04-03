package parse;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;


public class ParseMatrix {
	private String[] rowids; 
	private int[][] mat;
	private int numRows;
	private int numCols;
	
	public ParseMatrix(InputStream genotypes) {
		
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
		numRows = lines.size();
		numCols = linetokens.length - 1;
		
		int[][] ret = new int[numRows][numCols];
		rowids = new String[numRows];
		for(int i=0; i<numRows; i++){
			String l = lines.get(i);
			
			String[] tokens = l.split("\\s+");
			//add rowid to list
			rowids[i] = tokens[0];
			//
			for(int j=1; j<numCols ;j++){
				ret[i][j] = Integer.parseInt(tokens[j]);
			}
		}
		
		mat = ret;
		
		
	}
	
	public String[] getRowIDs(){
		return rowids;
	}
	
	public int[][] getMat(){
		return mat;
	}
	
	public int getNumRows(){
		return numRows;
	}
	
	public int getNumCols(){
		return numCols;
	}
}
