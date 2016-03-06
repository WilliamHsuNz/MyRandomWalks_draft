package beast.app.seqgen;

import beast.util.Randomizer;

public class SRWDirectionTest {

	double angle;
	boolean clockwiseOrAnticlockwise;

	public void example(){

		for(int i = 0; i <1000; i++){
			angle = Randomizer.nextDouble()*Math.PI;
			
			angle = (Randomizer.nextBoolean() == true)? angle: angle *-1;
			//if(Randomizer.nextDouble() >= 0.5)
				//clockwiseOrAnticlockwise = true;
			//else
				//clockwiseOrAnticlockwise = false;
			//if(clockwiseOrAnticlockwise == true)
				//angle = angle;
			//else
				//angle = -1 *angle;
			System.out.println(angle);
		}
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		(new SRWDirectionTest()).example();

	}

}
