package beast.app.seqgen;

import sim.util.distribution.VonMises;
import ec.util.MersenneTwisterFast;

public class VonMisesTest{
	double k = 1.0;
	MersenneTwisterFast mtf = new MersenneTwisterFast();
	VonMises vm = new VonMises(k, mtf);
	int numberOfDraws = 1000;
	double [] vmDraws = new double[numberOfDraws];
	double sum;
	double mean;
	double e; 
	double total_e;
	double std;
	
	public void example(){
		for(int i = 0; i < numberOfDraws; i++){
			vmDraws[i]= vm.nextDouble();
			System.out.println(vmDraws[i]);
		}
		sum = 0.0;
		for(int i = 0; i < numberOfDraws; i ++){
			sum = sum + vmDraws[i];
		}
		mean = sum/numberOfDraws;
		System.out.println("Mean :" + mean);
		
		
		for(int i = 0; i < numberOfDraws; i++){
			e = Math.pow(vmDraws[i] - mean, 2);
			total_e = total_e + e;
		}
		total_e = total_e/numberOfDraws;
		std = Math.sqrt(total_e);
		System.out.println("Standard deviation: " + std);
		
	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		(new VonMisesTest()).example();
	}
}
