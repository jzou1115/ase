package sample;

public class Sample implements Comparable<Sample>{
	String id;

	@Override
	public int compareTo(Sample o) {
		return id.compareTo(o.id);
	}
	
}
