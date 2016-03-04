package beast.app.seqgen;

import ec.util.MersenneTwisterFast;
import sim.util.distribution.VonMises;

public class TrigTests {
	double angle;
	double lat;
	double lon;
	double k = 1.0;
	MersenneTwisterFast mtf = new MersenneTwisterFast();
	VonMises vm = new VonMises(k, mtf);
	

	
	
	public void example(){
		angle = vm.nextDouble();
		lat = Math.sin(angle);
		lon = Math.cos(angle);
		System.out.println("Angle: " + angle);
		System.out.println("Latitude: " + lat);
		System.out.println("Longitude: " + lon);
		
	}	

	public static void main(String[]args){
		(new TrigTests()).example();
	}
	
}	