package genome;

public class ChromState implements Comparable<ChromState>{
	String state;
	GenomicRegion region;
	
	public ChromState(String s, GenomicRegion r){
		state = s;
		region = r;
	}
	
	public GenomicRegion getRegion(){
		return region;
	}
	
	public String getState(){
		return state;
	}
	
	@Override
	public int compareTo(ChromState c) {
		return this.region.compareTo(c.region);
	}

	@Override
	public boolean equals(Object o){
		if(o == null) return false;
		if(this == o) return true;
		if(o instanceof ChromState){
			ChromState other = (ChromState) o;
			return state.equals(other.state) && region.equals(other.region);
		}
		return false;
	}
	
	@Override
	public int hashCode(){
		return (int) (7*state.hashCode()+13*region.hashCode());
	}
}
