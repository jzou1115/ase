package functions;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.List;

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
			double logw1 = logchoose(max, min-i);
			double logw2 = logchoose(n-max, i);
			double logtotal = logchoose(n, min);
			
			double logp = logw1+logw2-logtotal;
			double p = Math.exp(logp);
			sig= sig+p;
		}
		
		//roundoff error
		if(sig>1){
			sig = 1.0;
		}

		return sig;
	}
	
	public double logchoose(int a, int b){
		double lc = 0.0;
	
		if(a==b){
			return 0.0;
		}
		if(b==1 || b==(a-1)){
			return Math.log(a);
		}
		
		if(b > (a-b)){
			for(int i=b+1; i<=a; i++){
				lc = lc + Math.log(i);
			}
			for(int j=1; j<=(a-b); j++){
				lc = lc - Math.log(j);
			}
		}
		else{
			for(int i=(a-b+1); i<=a; i++){
				lc = lc + Math.log(i);
			}
			for(int j=1; j<=b; j++){
				lc = lc - Math.log(j);
			}
		}
		
		return lc;
	}
	
	/**
	public BigDecimal choose(int a, int b){
		BigDecimal bc;
		if(a==b){
			return BigDecimal.valueOf(1);
		}
		if(b==1 || b==a-1){
			return BigDecimal.valueOf(a);
		}
		BigDecimal numerator = BigDecimal.valueOf(1);
		BigDecimal denominator = BigDecimal.valueOf(1);
		if(b > (a-b)){
			for(int i=b+1; i<=a;i++){
				numerator = numerator.multiply(BigDecimal.valueOf(i));
			}
			for(int j=1; j<=(a-b); j++){
				denominator = denominator.multiply(BigDecimal.valueOf(j));
			}
			bc = numerator.divide(denominator);
		}
		
		else{
			for(int i=(a-b+1); i<=a;i++){
				numerator = numerator.multiply(BigDecimal.valueOf(i));
			}
			for(int j=1; j<=b; j++){
				denominator = denominator.multiply(BigDecimal.valueOf(j));
			}
			bc = numerator.divide(denominator);
		}
		return bc;
	}
	**/
	
}
