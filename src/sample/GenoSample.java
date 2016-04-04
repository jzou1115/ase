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
