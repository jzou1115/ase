package sample;

public class ExpSample extends Sample{
	public int hasASE;
	
	public ExpSample(String s, int b){
		id = s;
		hasASE = b;
	}
	
	public String toString(){
		return id+"\t"+hasASE;
	}
	
	public int getASE(){
		return hasASE;
	}
}
