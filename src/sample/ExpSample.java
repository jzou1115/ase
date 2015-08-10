package sample;

public class ExpSample extends Sample{
	int hasASE;
	
	public ExpSample(String s, int b){
		id = s;
		hasASE = b;
	}
	
	public String toString(){
		return id+"\t"+hasASE;
	}
}
