package functions;
import java.util.List;

import org.apache.commons.math3.util.CombinatoricsUtils;

import sample.ExpSample;

public class ComputeSig {
	//Number of samples = length of w and v
	int n;
	//Number of 1's in vector w
	int k;
	//Number of 1's in ase
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
	
	ComputeSig(int w,List<ExpSample> ase, int errors){
		n = ase.size();
		k = w;
		m = ones(ase);
		e = errors;
	}
	
	public int ones(int[] a){
		int ret = 0;
		for(int i=0; i<a.length; i++){
			ret = ret + a[i];
		}
		return ret;
	}
	
	public int ones(List<ExpSample> ase){
		int ret =0;
		for(int i=0; i< ase.size();i++){
			ret=ret+ase.get(i).getASE();
		}
		return ret;
	}
	public double[] significance(){
		//significance for each level of errors
		sig = new double[e+1];
		
		//probability that a permutation will have hamming distance <= err
		double cumulative=0;
		for(int err=0; err<e+1; err++){
			int alpha = err+ Math.abs(m+k);
			//alpha must be even and positive
			if((alpha & 1) == 0  && alpha>=0){
				System.out.println("e="+err);
				System.out.println("n="+n);
				System.out.println("m="+m);
				System.out.println("k="+k);
				System.out.println("alpha="+alpha);
				
				long w1 = CombinatoricsUtils.binomialCoefficient(m,k-(alpha/2));
				long w2 = CombinatoricsUtils.binomialCoefficient(n-m, alpha/2);
				long possible = CombinatoricsUtils.binomialCoefficient(n, k);
				//probability that a permutation of w will have hamming distance of err
				double p = (double) w1*w2*1.0/possible;
				System.out.println("p="+p);
				cumulative = cumulative +p;
			}
			else{
				System.out.println("e="+err);
				System.out.println("n="+n);
				System.out.println("m="+m);
				System.out.println("k="+k);
				System.out.println("alpha="+alpha);
			}
			
			sig[err] = cumulative;
		}
		return sig;
	}
	
}
