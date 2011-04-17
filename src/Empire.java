/*	
 * Copyright 2011, Robin Roche
 * This file is part of jica.

    jica is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    jica is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with jica.  If not, see <http://www.gnu.org/licenses/>.
*/


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
	public double[] getImperialistPosition() 
	{
		return imperialistPosition;
	}
	public void setImperialistPosition(double[] imperialistPosition) 
	{
		this.imperialistPosition = imperialistPosition;
	}
	public double getImperialistCost() 
	{
		return imperialistCost;
	}
	public void setImperialistCost(double imperialistCost) 
	{
		this.imperialistCost = imperialistCost;
	}
	public double[][] getColoniesPosition() 
	{
		return coloniesPosition;
	}
	public void setColoniesPosition(double[][] coloniesPosition) 
	{
		this.coloniesPosition = coloniesPosition;
	}
	public double[] getColoniesCost() 
	{
		return coloniesCost;
	}
	public void setColoniesCost(double[] coloniesCost) 
	{
		this.coloniesCost = coloniesCost;
	}
	public double getTotalCost() 
	{
		return totalCost;
	}
	public void setTotalCost(double totalCost) 
	{
		this.totalCost = totalCost;
	}
	public void setColonyPosition(int colonyIndex, double[] position) 
	{
		this.coloniesPosition[colonyIndex] = position;
	}
	public void setColonyCost(int colonyIndex, double cost) 
	{
		this.coloniesCost[colonyIndex] = cost;
	}


}

