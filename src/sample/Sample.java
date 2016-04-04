package sample;

public class Sample implements Comparable<Sample>{
	String id;

	@Override
	public int compareTo(Sample o) {
		return this.getID().compareTo(o.getID());
	}
	
	public String getID(){
		return id;
	}
	
}
