package sample;

public class Sample implements Comparable<Sample>{
	String id;

	public String getID(){
		return id;
	}
	
	@Override
	public int compareTo(Sample o) {
		return this.getID().compareTo(o.getID());
	}
	
    @Override
    public boolean equals(Object obj) {
    	if(obj instanceof Sample){
    		Sample other = (Sample) obj;
    		if(this.getID().equals(other.getID())){
    			return true;
    		}
    	}
        return false;
    }

    @Override
    public int hashCode() {
        return id.hashCode();
    }
	
}
