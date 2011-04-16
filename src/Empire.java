
/**
 * Class creating the empire type, with its imperialists, colonies 
 * and their respective positions and costs
 * @author Robin Roche
 */
public class Empire
{

	// Empire variables
	int problemDimension;
	double[] imperialistPosition = new double[problemDimension];
	double imperialistCost;
	double[][] coloniesPosition;
	double[] coloniesCost;
	double totalCost;
	
	
	/**
	 * Constructor
	 * @param problemDimension
	 */
	public Empire(int problemDimension)
	{
		this.problemDimension = problemDimension;
	}

	
	/**
	 * To initialize the colonies
	 * @param numOfColonies
	 */
	public void init(int numOfColonies)
	{
		coloniesPosition = new double[numOfColonies][problemDimension];
		coloniesCost = new double[numOfColonies];
	}

	
	// Getters and setters
	public double[] getImperialistPosition() {
		return imperialistPosition;
	}
	public void setImperialistPosition(double[] imperialistPosition) {
		this.imperialistPosition = imperialistPosition;
	}
	public double getImperialistCost() {
		return imperialistCost;
	}
	public void setImperialistCost(double imperialistCost) {
		this.imperialistCost = imperialistCost;
	}
	public double[][] getColoniesPosition() {
		return coloniesPosition;
	}
	public void setColoniesPosition(double[][] coloniesPosition) {
		this.coloniesPosition = coloniesPosition;
	}
	public double[] getColoniesCost() {
		return coloniesCost;
	}
	public void setColoniesCost(double[] coloniesCost) {
		this.coloniesCost = coloniesCost;
	}
	public double getTotalCost() {
		return totalCost;
	}
	public void setTotalCost(double totalCost) {
		this.totalCost = totalCost;
	}
	public void setColoniesPosition(int bestColonyInd, double[] oldImperialistPosition) 
	{
		this.coloniesPosition[bestColonyInd] = oldImperialistPosition;
	}
	public void setColoniesCost(int bestColonyInd, double oldImperialistCost) 
	{
		this.coloniesCost[bestColonyInd] = oldImperialistCost;
	}
	public void setColonyPosition(int i, double[] ds) 
	{
		this.coloniesPosition[i] = ds;			
	}

}

