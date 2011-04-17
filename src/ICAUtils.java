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


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;


/**
 * Contains several methos on arrays used by the ICA algorithm
 * @author Robin Roche
 */
public class ICAUtils 
{

	/**
	 * Returns the index of the max value contained in the vector
	 * @param vector values
	 * @return the index of the max value
	 */
	public int getMaxIndex(double[] vector)
	{
		double max = Double.MIN_VALUE;
		int i;
		int bestIndex = 0;
		for(i=0; i<vector.length; i++) 
		{
			if(vector[i] > max)
			{
				max = vector[i];
				bestIndex = i;
			}
		}
		return bestIndex;
	}



	/**
	 * Returns the index of the min value contained in the vector
	 * @param vector values
	 * @return the index of the min value
	 */
	public int getMinIndex(double[] vector)
	{
		double min = Double.MAX_VALUE;
		int i;
		int bestIndex = 0;
		for(i=0; i<vector.length; i++) 
		{
			if(vector[i] < min)
			{
				min = vector[i];
				bestIndex = i;
			}
		}
		return bestIndex;
	}



	/**
	 * Returns the mean value of a vector
	 * @param vector the vector
	 * @return the mean value
	 */
	public double getMean(double[] vector) 
	{
		double sum = 0;
		for (int i=0; i<vector.length; i++) 
		{
			sum += vector[i];
		}
		return sum / vector.length;
	}



	/**
	 * Returns the norm of a vector
	 * @param vector the vector
	 * @return the norm
	 */
	public double getNorm(double[] vector) 
	{
		double sum = 0;		
		for(int i=0; i<vector.length; i++)
		{
			sum = sum + Math.pow(vector[i],2);
		}
		return Math.sqrt(sum);
	}



	/**
	 * Returns the sum of the elements on a vector
	 * @param vector the vector
	 * @return the sum
	 */
	public double getSum(double[] vector)
	{
		double sum = 0;
		for (double i : vector) 
		{
			sum += i;
		}
		return sum;
	}



	/**
	 * Returns the sum of the elements on a vector
	 * @param vector the vector
	 * @return the sum
	 */
	public int getSum(int[] vector)
	{
		int sum = 0;
		for (int i : vector) 
		{
			sum += i;
		}
		return sum;
	}



	/**
	 * Returns an array with the maximum values of two two-dimensional arrays
	 * @param array1 the first array
	 * @param array2 the second array
	 * @return an array with the max values
	 */
	public double[][] max(double[][] array1, double[][] array2) 
	{
		double maxArray[][] = new double[array1.length][array1[0].length];

		for(int i=0; i<array1.length; i++)
		{
			for(int j=0; j<array1[0].length; j++)
			{
				maxArray[i][j] = Math.max(array1[i][j], array2[i][j]);
			}
		}
		return maxArray;
	}



	/**
	 * Returns an array with the minimum values of two two-dimensional arrays
	 * @param array1 the first array
	 * @param array2 the second array
	 * @return an array with the min values
	 */
	public double[][] min(double[][] array1, double[][] array2) 
	{
		double minArray[][] = new double[array1.length][array1[0].length];

		for(int i=0; i<array1.length; i++)
		{
			for(int j=0; j<array1[0].length; j++)
			{
				minArray[i][j] = Math.min(array1[i][j], array2[i][j]);
			}
		}
		return minArray;
	}



	/**
	 * Returns a vector with the n integers (each appearing once) in a random order
	 * This is equivalent to the Matlab function randperm()
	 * @param n the number of values
	 * @param r the random generator
	 * @return the vector of integers
	 */
	public int[] randperm(int n, Random r) 
	{
		ArrayList<Integer> nVector = new ArrayList<Integer>();
		for(int i=0; i<n; i++)
		{
			nVector.add(i);
		}

		int[] outputVector = new int[n];
		int outputIndex = 0;

		while(nVector.size()>0)
		{
			int position = r.nextInt(nVector.size());
			outputVector[outputIndex] = nVector.get(position);
			outputIndex++;
			nVector.remove(position);
		}

		return outputVector;
	}



	/**
	 * Returns an array with a copy of a pattern, done n times
	 * This is equivalent to the Matlab function repmat()
	 * @param pattern the pattern to to repeat
	 * @param n the number of times
	 * @return the array with the copies
	 */
	public double[][] repmat(double[] pattern, int n) 
	{
		double[][] outputArray = new double[n][pattern.length];
		for(int i=0; i<n; i++)
		{
			outputArray[i] = pattern;
		}
		return outputArray;
	}



	/**
	 * Returns an extract of an array between two indexes (i.e. a range of it)
	 * @param array the array
	 * @param startIndex the index where to start
	 * @param endIndex the index where to stop
	 * @return the array extract
	 */
	public double[][] extractArrayRange(double[][] array, int startIndex, int endIndex)
	{
		double[][] arrayExtract = new double[endIndex-startIndex][array[0].length];
		int newIndex = 0;
		for(int i=startIndex; i<endIndex ; i++)
		{
			arrayExtract[newIndex] = array[i];
			newIndex++;
		}
		return arrayExtract;
	}



	/**
	 * Returns an extract of an array, with only selected indexes extracted
	 * @param array the array to be extracted from
	 * @param selectedIndexes the indexes to extract data from
	 * @return the array with the extracted data
	 */
	public double[][] extractGivenArrayParts(double[][] array, int[] selectedIndexes) 
	{
		double[][] arrayExtract = new double[selectedIndexes.length][array[0].length];
		int index;
		for(int i=0; i<selectedIndexes.length; i++)
		{
			index = selectedIndexes[i];
			arrayExtract[i] = array[index];
		}
		return arrayExtract;
	}



	/**
	 * Returns an extract of an array, with only selected indexes extracted
	 * @param array the array to be extracted from
	 * @param selectedIndexes the indexes to extract data from
	 * @return the array with the extracted data
	 */
	public double[] extractGivenArrayParts(double[] array, int[] selectedIndexes) 
	{
		double[] arrayExtract = new double[selectedIndexes.length];
		int index;
		for(int i=0; i<selectedIndexes.length; i++)
		{
			index = selectedIndexes[i];
			arrayExtract[i] = array[index];
		}
		return arrayExtract;
	}



	/**
	 * Prints the values of an array
	 * @param arrayName
	 * @param array
	 */
	public void printArray(String arrayName, double[] array)
	{
		System.out.println(arrayName + ": " + Arrays.toString(array));
	}



	/**
	 * Prints the values of an array
	 * @param arrayName
	 * @param array
	 */
	public void printArray(String arrayName, double[][] array)
	{
		for(int i=0; i<array.length; i++)
		{
			System.out.println(arrayName + "[" + i + "]: " + Arrays.toString(array[i]));
		}
	}



	/**
	 * Prints the values of an array
	 * @param arrayName
	 * @param array
	 */
	public void printArray(String arrayName, int[] array) 
	{
		System.out.println(arrayName + ": " + Arrays.toString(array));
	}



	/**
	 * Prints the characteristics of an empire
	 * @param empire
	 */
	public void printEmpire(Empire empire, int empireIndex) 
	{
		System.out.println("Empire " + empireIndex);
		System.out.println("Number of colonies: " + empire.getColoniesCost().length);
		printArray("imperialistPosition", empire.getImperialistPosition());
		System.out.println("imperialistCost: " + empire.getImperialistCost());
		printArray("coloniesPosition", empire.getColoniesPosition());
		printArray("coloniesCost", empire.getColoniesCost());
		System.out.println("totalCost: " + empire.getTotalCost());
	}


	
	/**
	 * Returns the total number of colonies in the empires list
	 * @param empiresList
	 * @return the number of colonies
	 */
	public int getTotalColoniesCount(Empire[] empiresList) 
	{
		int currentNumOfColonies = 0;
		for(int i=0; i<empiresList.length; i++)
		{
			currentNumOfColonies = currentNumOfColonies + empiresList[i].getColoniesPosition().length;
		}
		return currentNumOfColonies;
	}

}
