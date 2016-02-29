package beast.app.seqgen;

import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.BivariateRealFunction;

public class BivariateRealFunctionTest implements BivariateRealFunction {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		//double value;
		BivariateRealFunction bv = new BivariateRealFunctionTest();
		double cur;
		try{
		cur = bv.value(0.1, 0.2);
		System.out.println(cur);
		}catch(Exception e){
			System.out.print(e);
		}
	}

	@Override
	public double value(double x, double y) throws FunctionEvaluationException {
		// TODO Auto-generated method stub
		double value = x * y;
		return value;
	}

}
