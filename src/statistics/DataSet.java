package statistics;

public class DataSet {

	int[] dataSet;
	
	public double getMean() {
		int sum = 0;
		for (int i : dataSet){
			sum = sum + i;
		}
		return 1.0 * sum/dataSet.length;
	}
	
	public double getStandardDeviation(){
		double meanSquared = Math.pow(this.getMean(), 2);
		double expSquared = 0;
		for (int i: dataSet){
			expSquared = expSquared + Math.pow(i,2);
		}
		expSquared = expSquared / dataSet.length;
		return Math.sqrt(expSquared - meanSquared);
	}

}
