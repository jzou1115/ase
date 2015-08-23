package sample;

public class Sample implements Comparable<Object>{
	public String id;

	public int compareTo(Sample s) {
		return id.compareTo(s.id);
	}

	@Override
	public int compareTo(Object o) {
		// TODO Auto-generated method stub
		return 0;
	}
	
}
