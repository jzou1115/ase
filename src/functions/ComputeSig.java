package functions;
import org.apache.commons.math3.util.CombinatoricsUtils;

public class ComputeSig {
	//Number of samples = length of w and v
	int n;
	//Number of 1's in vector w
	int k;
	//Number of 1's in vector v
	int m;
	//Maximum number of errors allowed
	int e;
	
	double[] sig;
	
	ComputeSig(int[] w, int[] v, int errors){
		n = w.length;
		k = ones(w);
		m = ones(v);
		e = errors;
	}
	
	public int ones(int[] a){
		int ret = 0;
		for(int i=0; i<a.length; i++){
			ret = ret + a[i];
		}
		return ret;
	}
	
	public double[] significance(){
		sig = new double[e];
		for(int err=0; err<=e; e++){
			long w1 = CombinatoricsUtils.binomialCoefficient(m,k-((err-m+k)/2));
			long w2 = CombinatoricsUtils.binomialCoefficient(n-m, (err-m+k)/2);
			long possible = CombinatoricsUtils.binomialCoefficient(n, k);
			long p = w1*w2/possible;
			sig[err] = p;
		}
		return sig;
	}
	
}
