
public class TestICA 
{

	/**
	 * @param args
	 */
	public static void main(String[] args) 
	{
		int problemDimension = 2; 
		double[] minBounds = {-10,10};
		double[] maxBounds = {-10,10};
		
		Object[] argsICA ={problemDimension,minBounds,maxBounds};
		ImperialistCompetitiveAlgorithm ica = new ImperialistCompetitiveAlgorithm(argsICA);

		ica.runICA();
		
	}

}
