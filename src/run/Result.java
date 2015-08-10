package run;

public class Result {
	int numCorrect;
	int numIncorrect;
	String snpID;
	
	Result(int c, int i, String snp){
		numCorrect=c;
		numIncorrect=i;
		snpID=snp;
	}
	
	public String toString(){
		return snpID+"\t"+numCorrect+"\t"+numIncorrect;
	}
}
