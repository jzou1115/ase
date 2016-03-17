package parse;

import genome.Gene;
import genome.GenomicCoordinate;
import genome.SNP;
import sample.ExpSample;
import sample.GenoSample;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class ParseGene {

	public static Gene parseGene(String line){
		String[] tokens = line.trim().split("\\s+");
		if(tokens.length==6){
			String id = tokens[0];
			int chr = Integer.parseInt(tokens[1]);
			long s = Long.parseLong(tokens[2]);
			long e = Long.parseLong(tokens[3]);
			long t = Long.parseLong(tokens[5]);
			GenomicCoordinate start = new GenomicCoordinate(chr, s);
			GenomicCoordinate end = new GenomicCoordinate(chr, e);
			GenomicCoordinate tss = new GenomicCoordinate(chr, t);
			return new Gene(id, start, end, tss);
		}
		if(tokens.length==5){
			String id = tokens[0];
			int chr = Integer.parseInt(tokens[1]);
			long s = Long.parseLong(tokens[2]);
			long e = Long.parseLong(tokens[3]);
			long t = Long.parseLong(tokens[4]);
			GenomicCoordinate start = new GenomicCoordinate(chr, s);
			GenomicCoordinate end = new GenomicCoordinate(chr, e);
			GenomicCoordinate tss = new GenomicCoordinate(chr, t);
			return new Gene(id, start, end, tss);
		}
		String id = tokens[0];
		int chr = Integer.parseInt(tokens[1]);
		long s = Long.parseLong(tokens[2]);
		long e = Long.parseLong(tokens[3]);
		GenomicCoordinate start = new GenomicCoordinate(chr, s);
		GenomicCoordinate end = new GenomicCoordinate(chr, e);
		return new Gene(id, start, end);
	}
	public static List<Gene> readGeneGroup(InputStream in){
		List<Gene> genes = new ArrayList<Gene>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(in));
		String line;
		try {
			while((line = reader.readLine()) != null){
				try{
					genes.add(parseGene(line));
				} catch (Exception e){
					//do nothing
					System.out.println("cannot add Gene");
				}
			}
			reader.close();
		} catch (IOException e) {
			System.out.println("no lines");
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return genes;
	}
	

	

}
