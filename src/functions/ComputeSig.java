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
	
	//this is currently the right constructor
	ComputeSig(int size, int ase, int geno, int errors){
		n=size;
		m=ase;
		k=geno;
		e=errors;
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
	public double significance(){
		double sig=0.0;
		int min = Math.min(m, k);
		int max = Math.max(m, k);
		int upper1= (e-max+min)/2;
		int upper2 = Math.min(min, upper1);
		int upper3 = Math.min(n- max, upper2);
		
		
		for(int i=0; i<=upper3;i++){
			//System.out.println("min="+min);
			//System.out.println("max="+max);
			//System.out.println("i="+i);
			long w1 = CombinatoricsUtils.binomialCoefficient(max, min-i);
			long w2 = CombinatoricsUtils.binomialCoefficient(n-max, i);
			long total = CombinatoricsUtils.binomialCoefficient(n, min);
			double p = w1*w2*1.0/total;
			//System.out.println(p);
			sig= sig+p;
		}

		return sig;
	}
	
}
