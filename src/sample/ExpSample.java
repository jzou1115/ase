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
	
	public int getASE(){
		return hasASE;
	}
	
	
    @Override
    public boolean equals(Object obj) {
    	if(obj instanceof GenoSample){
    		GenoSample other = (GenoSample) obj;
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
