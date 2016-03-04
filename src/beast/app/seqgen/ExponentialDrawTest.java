package beast.app.seqgen;
import beast.util.Randomizer;


public class ExponentialDrawTest {
	public void example(){
		double lTimeElapsed = 20;
		double m_lambda = 1.0;
		int num_direction_change = 0;
		int count = 0;
		double tSum = 0;
		double exp_Draw;
		while(tSum < lTimeElapsed){
			exp_Draw = Randomizer.nextExponential(m_lambda);
			tSum += exp_Draw; 
			System.out.println("Exponential Draws: "+ exp_Draw);
			count += 1;
		}
		System.out.println(Math.PI);
	}
	
	

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		(new ExponentialDrawTest()).example();
	}

	
}
