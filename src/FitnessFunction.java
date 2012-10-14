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
 * This class contains the fitness function
 * @author Robin Roche
 */
public class FitnessFunction 
{

	public int dimension = 20;
	private double minVal = -2.048;
	private double maxVal = 2.048;
	public double[] minBounds;	// The minimum bounds for each dimension
	public double[] maxBounds;	// The maximum bounds for each dimension
	public int nbEvals = 0;
	
	
	// Constructor
	public FitnessFunction()
	{
		// Initialize boundaries
		minBounds = new double[dimension];
		maxBounds = new double[dimension];
		for(int i=0;i<dimension;i++)
		{
			minBounds[i] = minVal;
			maxBounds[i] = maxVal;
		}
	}
	
	
	/**
	 * Returns the fitness of one country
	 * @param individual the solution to evaluate
	 * @return the fitness
	 */
	public double getFitnessValue(double[] individual) 
	{
		double fitness = 0;		
		
		// Sphere function	
		for(int i=0; i<individual.length; i++)
		{
			fitness = fitness + Math.pow(individual[i],2);
		}
		
//		// Rastrigin function	
//		for(int i=0; i<individual.length; i++)
//		{
//			fitness = fitness + (Math.pow(individual[i],2)-10*Math.cos(2*Math.PI*individual[i]));
//		}
//		fitness = 10*dimension + fitness;
		
//		// Rosenbrock function
//		for(int i=0; i<individual.length-1; i++)
//		{
//			fitness = fitness + 100*Math.pow((Math.pow(individual[i],2)-individual[i+1]),2) + Math.pow((individual[i]-1),2);
//		}
		
//		// Ackley function
//		double a = 20; 
//		double b = 0.2; 
//		double c = 2*Math.PI;
//		double s1 = 0; 
//		double s2 = 0;
//		for(int i=0; i<individual.length; i++)
//		{
//			s1 = s1 + Math.pow(individual[i],2);
//			s2 = s2 + Math.cos(c*individual[i]);
//		}
//		fitness = -a * Math.exp( -b * Math.sqrt(1/individual.length*s1)) - Math.exp(1/individual.length*s2) + a + Math.exp(1);

		nbEvals++;
		return fitness;
	}
	
}
