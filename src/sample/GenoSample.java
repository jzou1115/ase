package sample;

public class GenoSample extends Sample{
	int isHetero;
	
	public GenoSample(String s, int b){
		id = s;
		isHetero = b;
	}
	
	public String toString(){
		return id+"\t"+isHetero;
	}
	
	public int getHetero(){
		return isHetero;
	}
	
	public String getID(){
		return id;
	}
}
