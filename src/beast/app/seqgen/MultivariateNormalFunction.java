package beast.app.seqgen;

import java.util.Arrays;

import org.apache.commons.math3.distribution.AbstractMultivariateRealDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;

public class MultivariateNormalFunction {
	
	public static void example(){
	
		double mu [] = {0.0,0.0};
		System.out.println(Arrays.toString(mu));
		double sigma [][]= {{1.0,0.0},{0.0,1.0}};
		System.out.println(Arrays.deepToString(sigma));
		//double vals [] = new double [1000];
		
		double lTimeElapsed = 0.3;
		double [][]l_sigma = new double [sigma.length][sigma[0].length];
		for (int i = 0; i < sigma.length; i++){
			for (int j = 0; j < sigma[i].length; j++){
				l_sigma[i][j] = lTimeElapsed * sigma[i][j];
			}
		}
		
		System.out.println(Arrays.deepToString(sigma));
		System.out.println(Arrays.deepToString(l_sigma));
		
		for(int i = 0; i < 1; i++){
			MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(mu, sigma);
			double vals[] = mnd.sample();
			System.out.println(Arrays.toString(vals));
		}
		

	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		example();
		
	}


}
