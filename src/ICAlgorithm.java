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


import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;


/**
 * This class contains the ICA algorithm methods
 * @author Robin Roche
 */
public class ICAlgorithm
{

	// ICA parameters
	int numOfCountries = 80;               		// Number of initial countries
	int numOfInitialImperialists = 8;      		// Number of initial imperialists
	int numOfAllColonies = numOfCountries - numOfInitialImperialists;
	int numOfDecades = 1000;					// Number of decades (generations)
	double revolutionRate = 0.1;               	// Revolution is the process in which the socio-political characteristics of a country change suddenly
	double assimilationCoefficient = 2;        	// In the original paper assimilation coefficient is shown by "beta"
	double assimilationAngleCoefficient = .5; 	// In the original paper assimilation angle coefficient is shown by "gama"
	double zeta = 0.02;							// Total Cost of Empire = Cost of Imperialist + Zeta * mean(Cost of All Colonies)
	double dampRatio = 0.99;					// The damp ratio
	boolean stopIfJustOneEmpire = false;        // Use "true" to stop the algorithm when just one empire is remaining. Use "false" to continue the algorithm
	double unitingThreshold = 0.02;           	// The percent of search space size, which enables the uniting process of two empires


	// Variables
	protected long seed = System.currentTimeMillis();       
	Random r = new Random(seed);
	int problemDimension;						// Problem dimension
	double[] minBounds;							// Minimum bounds
	double[] maxBounds;							// Maximum bounds
	Empire[] empiresList = new Empire[numOfInitialImperialists];	// List of Empires
	double[][] initialCountries;				// The initial countries with their positions
	double[] initialCosts;						// The costs of the initial countries
	double[][] bestDecadePosition; 				// The best found position for each decade
	double[] minimumCost = new double[numOfDecades];	// The minimum cost of each decade
	double[] meanCost = new double[numOfDecades];		// The mean cost of each decade
	double[] searchSpaceSize;					// The search space size (between the min and max bounds)
	ICAUtils utils = new ICAUtils();			// A class with useful methods for array operations

	

	/**
	 * Constructor of the class
	 * @param args an object array with the parameters
	 * The first is the problem dimension
	 * The second a vector with the min bounds
	 * The third a vector with the max bounds
	 */
	public ICAlgorithm(Object[] args) 
	{
		problemDimension = (Integer) args[0]; 
		minBounds = (double[]) args[1];
		maxBounds = (double[]) args[2];
	}

	

	/**
	 * Runs the ICA algorithm
	 * @return the best found solution
	 */
	protected double[] runICA()
	{
		
		// Initialize variables
		bestDecadePosition = new double[numOfDecades][problemDimension];
		searchSpaceSize = new double[problemDimension];
		
		// Compute the problem search space, between the min and max bounds
		for(int i=0; i<problemDimension; i++)
		{
			searchSpaceSize[i] = maxBounds[i] - minBounds[i];
		}
		//utils.printArray("searchSpaceSize", searchSpaceSize);

		// Create an initial population of individuals (countries)
		initialCountries = generateNewCountries(numOfCountries, problemDimension, minBounds, maxBounds, r);
		//utils.printArray("initialCountries", initialCountries);
		
		// Compute the cost of each country: the lesser the cost, the more powerful the country is
		initialCosts = getCountriesCosts(initialCountries);
		//utils.printArray("initialCosts", initialCosts);

		// Sort the costs and the corresponding countries in assending order. The best countries will be in higher places.
		sortArray(initialCosts, initialCountries);
		//utils.printArray("initialCountries", initialCountries);
		//utils.printArray("initialCosts", initialCosts);
		
		// Create the initial empires
		createInitialEmpires();

		
		// While no stopping condition is met
		for(int decade=0; decade<numOfDecades; decade++)
		{
			
			//System.out.println("-------------");
			//System.out.println("Decade: " + decade);
			//System.out.println("Number of empires: " + empiresList.length);
			//System.out.println("Total number of colonies: " + utils.getTotalColoniesCount(empiresList));
			
			// Update the revolution rate
			revolutionRate = dampRatio * revolutionRate;
			
			// For each empire
			for (int i=0; i<empiresList.length; i++)
			{
				// Assimilation: movement of colonies toward their imperialist
				assimilateColonies(empiresList[i]);
				
				// Revolution: change the position of some colonies, to avoid getting stuck into local minima
				revolveColonies(empiresList[i]);
				empiresList[i].setColoniesCost( getCountriesCosts( empiresList[i].getColoniesPosition()) );

				// Empire possession: a strong colony can become the imperialist
				possesEmpire(empiresList[i]);
				
				// Update the empire's total cost
				empiresList[i].setTotalCost( empiresList[i].getImperialistCost() + zeta * utils.getMean(empiresList[i].getColoniesCost()) );
			}

			// Uniting Similiar Empires
			uniteSimilarEmpires();
			
			// Imperialistic Competition
			imperialisticCompetition();
			
			// If the user wants it, the algorithm can stop when only one empire remains
			if (empiresList.length == 1 && stopIfJustOneEmpire)
			{
				break;
			}

			// Extract and save results of the round
			double[] imperialistCosts = new double[empiresList.length];
			for(int i=0; i<empiresList.length ; i++)
			{
				imperialistCosts[i] = empiresList[i].getImperialistCost();
			}
			minimumCost[decade] = imperialistCosts[utils.getMinIndex(imperialistCosts)];
			meanCost[decade] = utils.getMean(imperialistCosts);
			bestDecadePosition[decade] = empiresList[utils.getMinIndex(imperialistCosts)].getImperialistPosition();

		}


		// Return results of the optimization
		double bestCost = minimumCost[utils.getMinIndex(minimumCost)];
		int bestIndex = utils.getMinIndex(minimumCost);
		double[] bestSolution = bestDecadePosition[bestIndex];
		System.out.println("Best solution: " + Arrays.toString(bestSolution));
		System.out.println("Best fitness: " + bestCost);
		
		return bestSolution;
	}



	/** Generates new countries
	 * @param numberOfCountries the number of countries to generate
	 * @param dimension the dimension of each country
	 * @param minVector the minimum values for each dimension 
	 * @param maxVector the maximum values for each dimension
	 * @param rand a random number generator
	 * @return an array with the new countries positions
	 */
	private double[][] generateNewCountries(int numberOfCountries, int dimension, double[] minVector, double[] maxVector, Random rand) 
	{
		double[][] countriesArray = new double[numberOfCountries][dimension];  
		for(int i=0; i<numberOfCountries; i++)
		{
			for(int j=0; j<dimension; j++)
			{
				countriesArray[i][j] = (maxVector[j] - minVector[j]) * rand.nextDouble() + minVector[j];
			}
		}
		return countriesArray;
	}



	/**
	 * Returns a vector with the cost of all countries computed according to their position
	 * @param numberOfCountries the number of countries
	 * @param countriesArray
	 * @return a vector with the costs of all countries
	 */
	private double[] getCountriesCosts(double[][] countriesArray) 
	{
		double[] costsVector = new double[countriesArray.length];
		for(int i=0; i<countriesArray.length; i++)
		{
			costsVector[i] = getCountryCost(countriesArray[i]);
		}
		return costsVector;
	}

	

	/**
	 * Returns the cost of one country
	 * @param country the country
	 * @return the cost
	 */
	private double getCountryCost(double[] country) 
	{
		double cost = 0;		
		for(int i=0; i<country.length; i++)
		{
			cost = cost + Math.pow(country[i],2);
		}
		return cost;
	}



	// TODO: directly operate on the arrays, not copies?
	/**
	 * Sorts an array according to its values and sorts another array in the same order
	 * The lowest value is put first
	 * @param arrayToSort the array to sort according to its values
	 * @param matchingArray the array that should be sorted according to the first one
	 */
	private void sortArray(final double[] arrayToSort, double[][] matchingArray) 
	{
		// Create an index array
		Integer[] sortOrder = new Integer[arrayToSort.length];
		for(int i=0; i<sortOrder.length; i++)
		{
			sortOrder[i] = i;
		}

		// Sort the array using the index array, thanks to a custom comparator
		Arrays.sort(sortOrder,new Comparator<Integer>() 
			{   
				public int compare(Integer a, Integer b)
				{
					double delta = arrayToSort[b]-arrayToSort[a];
					if(delta < 0) return 1;
					if(delta > 0) return -1;
					return 0;
				}
			}
		);
		
		// Create copies of the arrays
		double[] arrayToSortCopy = arrayToSort.clone();
        double[][] matchingArrayCopy = matchingArray.clone();
            
        // Output the values using the sorted indexes
		for(int i=0;i<sortOrder.length;i++)
		{
			initialCosts[i] = arrayToSortCopy[sortOrder[i]];
			initialCountries[i] = matchingArrayCopy[sortOrder[i]];
        }
	}


	
	/**
	 * Generates the initial empires
	 */
	private void createInitialEmpires()
	{

		// Extract the best countries to create empires
		double[][] allImperialistsPosition = utils.extractArrayRange(initialCountries, 0, numOfInitialImperialists);

		// Extract their costs
		double[] allImperialistsCost = new double[numOfInitialImperialists];
		System.arraycopy(initialCosts, 0, allImperialistsCost, 0, numOfInitialImperialists);

		// Extract the rest to create colonies
		double[][] allColoniesPosition = utils.extractArrayRange(initialCountries, numOfInitialImperialists, initialCountries.length);

		// Extract their costs
		double[] allColoniesCost = new double[initialCosts.length-numOfInitialImperialists]; 
		System.arraycopy(initialCosts, numOfInitialImperialists, allColoniesCost, 0, initialCosts.length-numOfInitialImperialists);	

		//utils.printArray("allImperialistsPosition", allImperialistsPosition);
		//utils.printArray("allImperialistsCost", allImperialistsCost);
		//utils.printArray("allColoniesPosition", allColoniesPosition);
		//utils.printArray("allColoniesCost", allColoniesCost);
		
		// Compute the power of imperialists
		double[] allImperialistsPower = new double[numOfInitialImperialists];
		if(allImperialistsCost[utils.getMaxIndex(allImperialistsCost)]>0)
		{
			for(int i=0; i<allImperialistsCost.length; i++)
			{
				allImperialistsPower[i] = 1.3 * allImperialistsCost[utils.getMaxIndex(allImperialistsCost)] - allImperialistsCost[i];
			}
		}
		else
		{
			for(int i=0; i<allImperialistsCost.length; i++)
			{
				allImperialistsPower[i] = 0.7 * allImperialistsCost[utils.getMaxIndex(allImperialistsCost)] - allImperialistsCost[i];
			}
		}
		//utils.printArray("allImperialistsPower", allImperialistsPower);
		
		// Set the number of colonies for the imperialists 
		int[] allImperialistNumOfColonies = new int[numOfInitialImperialists];
		for(int i=0; i<allImperialistsPower.length; i++)
		{
			allImperialistNumOfColonies[i] = (int) Math.round(allImperialistsPower[i]/utils.getSum(allImperialistsPower) * numOfAllColonies);
		}
		allImperialistNumOfColonies[allImperialistNumOfColonies.length-1] = 
			Math.max(
					numOfAllColonies - 
					utils.getSum(Arrays.copyOfRange(allImperialistNumOfColonies, 0, allImperialistNumOfColonies.length-1)),
					0
					);
		//utils.printArray("allImperialistNumOfColonies", allImperialistNumOfColonies);

		// Initialize the empires
		for(int i=0; i<numOfInitialImperialists; i++)
		{
			empiresList[i] = new Empire(problemDimension);
		}
		
		// Create a random permutation of integers
		int[] randomIndex = utils.randperm(numOfAllColonies, r);
		
		// Create the empires and attribute them their colonies
		for(int i=0; i<numOfInitialImperialists; i++)
		{
			int[] R = Arrays.copyOfRange(randomIndex, 0, allImperialistNumOfColonies[i]);
			empiresList[i].init(R.length);
			randomIndex = Arrays.copyOfRange(randomIndex, allImperialistNumOfColonies[i], randomIndex.length);
			
			empiresList[i].setImperialistPosition(allImperialistsPosition[i]);
			empiresList[i].setImperialistCost(allImperialistsCost[i]);
			empiresList[i].setColoniesPosition(utils.extractGivenArrayParts(allColoniesPosition, R));
			empiresList[i].setColoniesCost(utils.extractGivenArrayParts(allColoniesCost, R));
			empiresList[i].setTotalCost(empiresList[i].getImperialistCost() + zeta * utils.getMean(empiresList[i].getColoniesCost()));
			
			//utils.printEmpire(empiresList[i], i);
		}

		// If an empire has no colony, give it one
		for(int i=0; i<empiresList.length; i++)
		{
			if(empiresList[i].getColoniesPosition().length == 0)
			{
				empiresList[i].setColoniesPosition(generateNewCountries(1, problemDimension, minBounds, maxBounds, r));
				empiresList[i].setColoniesCost(getCountriesCosts(empiresList[i].getColoniesPosition()));
			}
		}

	}



	/**
	 * Assimilates the colonies of an empire: move their positions
	 * @param theEmpire
	 */
	private void assimilateColonies(Empire theEmpire)
	{

		// Get the number of colonies of the empire
		int numOfColonies = theEmpire.getColoniesPosition().length;

		// Create an array containing the distance between the imperialist and the colonies
		double[][] repmatArray = utils.repmat(theEmpire.getImperialistPosition(),numOfColonies);
		double[][] array = new double[numOfColonies][problemDimension];
		for(int i=0; i<numOfColonies; i++)
		{
			for(int j=0; j<problemDimension; j++)
			{
				array[i][j] = repmatArray[i][j] - theEmpire.getColoniesPosition()[i][j];
			}
		}
		//utils.printArray("array",array);

		// Create a new array to store the updated colonies positions
		double[][] coloniesPosition = new double[numOfColonies][problemDimension];

		// Fill the array
		for(int i=0; i<array.length; i++)
		{
			for(int j=0; j<array[0].length; j++)
			{
				coloniesPosition[i][j] = theEmpire.getColoniesPosition()[i][j] + 2 * assimilationCoefficient * r.nextDouble() * array[i][j];
			}
		}
		
		// Update the positions
		theEmpire.setColoniesPosition(coloniesPosition);

		// Bound the values with the min and max bounds
		double[][] minVarMatrix = utils.repmat(minBounds,numOfColonies);
		double[][] maxVarMatrix = utils.repmat(maxBounds,numOfColonies);
		theEmpire.setColoniesPosition(utils.max(theEmpire.getColoniesPosition(),minVarMatrix));
		theEmpire.setColoniesPosition(utils.min(theEmpire.getColoniesPosition(),maxVarMatrix));
		//utils.printArray("theEmpire.getColoniesPosition()", theEmpire.getColoniesPosition());
	}



	/**
	 * Make colonies of an empire revolve.
	 * This is equivalent to a perturbation which avoid getting stuck 
	 * into local minima.
	 * @param theEmpire to revolve
	 */
	private void revolveColonies(Empire theEmpire)
	{
		// Get the number of colonies to revolve
		int numOfRevolvingColonies = (int) Math.round((revolutionRate * theEmpire.getColoniesCost().length));

		// Create a new array with new random positions for the revolving colonies
		double[][] revolvedPosition = generateNewCountries(numOfRevolvingColonies, problemDimension, minBounds, maxBounds, r);

		// Generate a vector with integer values in a random order
		int[] R = utils.randperm(theEmpire.getColoniesCost().length, r);
		R = Arrays.copyOfRange(R, 0, numOfRevolvingColonies);

		// Update the positions of the revolved colonies of the empire
		for(int i=0; i<R.length; i++)
		{
			theEmpire.setColonyPosition(R[i], revolvedPosition[i]);
		}
	}



	/**
	 * Can make a colony become the imperialist 
	 * if it is more powerful than the imperialist.
	 * @param theEmpire
	 */
	private void possesEmpire(Empire theEmpire)
	{
		// Get the costs of the colonies
		double[] coloniesCost = theEmpire.getColoniesCost();

		// Get the cost of the best colony (the lowest cost)
		int bestColonyInd = utils.getMinIndex(coloniesCost);
		double minColoniesCost = coloniesCost[bestColonyInd]; 

		// If this cost is lower than the one of the imperialist
		if(minColoniesCost < theEmpire.getImperialistCost())
		{
			// Backup the position and cost of the former imperialist
			double[] oldImperialistPosition = theEmpire.getImperialistPosition();
			double oldImperialistCost = theEmpire.getImperialistCost();

			// Update the position and cost of the imperialist with the ones of the colony
			theEmpire.setImperialistPosition(theEmpire.getColoniesPosition()[bestColonyInd]);
			theEmpire.setImperialistCost(theEmpire.getColoniesCost()[bestColonyInd]);

			// Update the position and cost of the former colony with the ones of the former imperialist
			theEmpire.setColonyPosition(bestColonyInd, oldImperialistPosition);
			theEmpire.setColonyCost(bestColonyInd, oldImperialistCost);
		}
	}



	/**
	 * Unites imperialists that are close to each other
	 */
	private void uniteSimilarEmpires()
	{
		// Get the threshold distance between two empires
		double thresholdDistance = unitingThreshold * utils.getNorm(searchSpaceSize);
		
		// Get the number of empires
		int numOfEmpires = empiresList.length;

		// Compare each empire with the other ones
		for(int i=0; i<(numOfEmpires-1); i++)
		{
			for(int j=i+1; j<numOfEmpires; j++)
			{
				// Compute the distance between the two empires i and j
				double[] distanceVector = new double[empiresList[i].getImperialistPosition().length];
				for(int k=0; k<empiresList[i].getImperialistPosition().length; k++)
				{
					distanceVector[k] = empiresList[i].getImperialistPosition()[k] - empiresList[j].getImperialistPosition()[k];
				}
				double distance = utils.getNorm(distanceVector);

				// If the empires are too close
				if(distance<=thresholdDistance)
				{
					// Get the best and worst empires of the two
					int betterEmpireInd;
					int worseEmpireInd;
					if(empiresList[i].getImperialistCost() < empiresList[j].getImperialistCost())
					{
						betterEmpireInd=i;
						worseEmpireInd=j;
					}
					else
					{
						betterEmpireInd=j;
						worseEmpireInd=i;
					}

					// Update the positions
					double[][] newColoniesPosition = getColonyPositionsOfUnitedEmpire(betterEmpireInd, worseEmpireInd);
					empiresList[betterEmpireInd].setColoniesPosition(newColoniesPosition);

					// Update the costs
					double[] newColoniesCost = getColonyCostsOfUnitedEmpire(betterEmpireInd, worseEmpireInd);
					empiresList[betterEmpireInd].setColoniesCost(newColoniesCost);

					// Update the total cost of the united empire                                     
					empiresList[betterEmpireInd].setTotalCost(
							empiresList[betterEmpireInd].getImperialistCost() +
							zeta * utils.getMean(empiresList[betterEmpireInd].getColoniesCost())
					);

					// Update the empires list
					deleteAnEmpire(worseEmpireInd);
					//System.out.println("New number of empires: " + empiresList.length);

					return;
				}

			}
		}  
	}



	/**
	 * Returns the costs of the colonies of the united empire (after two empires merge) 
	 * @param betterEmpireInd the best empire
	 * @param worseEmpireInd the worst empire
	 * @return the corresponding colony costs
	 */
	private double[] getColonyCostsOfUnitedEmpire(int betterEmpireInd, int worseEmpireInd) 
	{
		// Get the new number of colonies of the united empire
		int newColoniesCount = 
			empiresList[betterEmpireInd].getColoniesCost().length + 
			1 + 
			empiresList[worseEmpireInd].getColoniesCost().length;

		// Create a new vector to store the costs of the colonies
		double[] newColoniesCost = new double[newColoniesCount];

		// At first, store the costs of the colonies of the best empire in the vector
		int i;
		for(i=0; i<empiresList[betterEmpireInd].getColoniesCost().length; i++)
		{
			newColoniesCost[i] = empiresList[betterEmpireInd].getColoniesCost()[i];
		}

		// Then add the cost of the former worst imperialist
		newColoniesCost[i] = empiresList[worseEmpireInd].getImperialistCost();
		
		// Finally, add the costs of the colonies of the worst empire
		int i2;	
		for(i2=i+1; i2<newColoniesCount; i2++)
		{
			newColoniesCost[i2] = empiresList[worseEmpireInd].getColoniesCost()[i2-empiresList[betterEmpireInd].getColoniesCost().length-1];
		}
		
		// Return the vector with the updated costs
		return newColoniesCost;
	}



	/**
	 * Returns the positions of the colonies of the united empire (after two empires merge) 
	 * @param betterEmpireInd the best empire
	 * @param worseEmpireInd the worst empire
	 * @return the corresponding colony positions
	 */
	private double[][] getColonyPositionsOfUnitedEmpire(int betterEmpireInd, int worseEmpireInd) 
	{
		// Get the new number of colonies of the united empire
		int newSize = 
			empiresList[betterEmpireInd].getColoniesPosition().length + 
			1 + 
			empiresList[worseEmpireInd].getColoniesPosition().length;
		
		// Create a new array to store the positions of the colonies
		double[][] newColoniesPosition = new double[newSize][problemDimension];

		// At first, store the positions of the colonies of the best empire in the array
		int i;
		for(i=0; i<empiresList[betterEmpireInd].getColoniesPosition().length; i++)
		{
			newColoniesPosition[i] = empiresList[betterEmpireInd].getColoniesPosition()[i];
		}
		
		// Then add the position of the former worst imperialist
		newColoniesPosition[i] = empiresList[worseEmpireInd].getImperialistPosition();
		
		// Finally, add the costs of the colonies of the worst empire
		int i2;	
		for(i2=i+1; i2<newSize; i2++)
		{
			newColoniesPosition[i2] = empiresList[worseEmpireInd].getColoniesPosition()[i2-empiresList[betterEmpireInd].getColoniesPosition().length-1];
		}

		// Return the array with the updated positions
		return newColoniesPosition;
	
	}



	/**
	 * Runs the competition between empires
	 */
	private void imperialisticCompetition()
	{
		
		// Generate a random number, and return if this number is too high
		double rand = r.nextDouble();
		if(rand > .11)
		{
			return;
		}
		
		// Idem if their is only one empire
		if(empiresList.length<=1)
		{
			return;
		}

		// Get the total cost of each empire
		double[] totalCosts = new double[empiresList.length]; 
		for(int i=0; i<empiresList.length; i++)
		{
			totalCosts[i] = empiresList[i].getTotalCost();
		}
		//utils.printArray("totalCosts", totalCosts);
		
		// Get the weakest empire (the one with the highest cost) and its cost
		int weakestEmpireInd = utils.getMaxIndex(totalCosts);
		double maxTotalCost = totalCosts[weakestEmpireInd]; 
		
		// Get the power of each empire
		double[] totalPowers = new double[empiresList.length];
		for(int i=0; i<empiresList.length; i++)
		{
			totalPowers[i] = maxTotalCost - totalCosts[i];
		}
		//utils.printArray("totalPowers", totalPowers);

		// Get the possession probability of each empire
		double[] possessionProbability = new double[empiresList.length];
		for(int i=0; i<empiresList.length; i++)
		{
			possessionProbability[i] = totalPowers[i] / utils.getSum(totalPowers);
		}
		//utils.printArray("possessionProbability", possessionProbability);
		
		// Select an empire according to their probabilities
		int selectedEmpireInd = selectAnEmpire(possessionProbability);
		//System.out.println("selectedind: " + selectedEmpireInd);
		
		// Generate a random integer
		int numOfColoniesOfWeakestEmpire = empiresList[weakestEmpireInd].getColoniesCost().length;
		int indexOfSelectedColony = r.nextInt(numOfColoniesOfWeakestEmpire);

		// Update the positions of the colonies of the selected empire 
		// by adding the position of the randomly selected colony of the weakest empire
		empiresList[selectedEmpireInd].setColoniesPosition(	
				concatenatePositions(
						empiresList[selectedEmpireInd].getColoniesPosition(), 
						empiresList[weakestEmpireInd].getColoniesPosition()[indexOfSelectedColony]
				)
		);

		// Idem for the costs
		empiresList[selectedEmpireInd].setColoniesCost( 
				concatenateCosts(
						empiresList[selectedEmpireInd].getColoniesCost(), 
						empiresList[weakestEmpireInd].getColoniesCost()[indexOfSelectedColony]
				) 
		);

		// Update the positions of the colonies of the weakest empire 
		// by removing the position of the randomly selected colony of the empire
		empiresList[weakestEmpireInd].setColoniesPosition( 
				removeColonyPosition(empiresList[weakestEmpireInd].getColoniesPosition(), indexOfSelectedColony) 
		);

		// Idem for the costs
		empiresList[weakestEmpireInd].setColoniesCost( 
				removeColonyCost(empiresList[weakestEmpireInd].getColoniesCost(), indexOfSelectedColony) 
		);

		// Get the number of colonies of the weakest empire
		numOfColoniesOfWeakestEmpire = empiresList[weakestEmpireInd].getColoniesCost().length;
		
		// If it has not more than 1 colony, then make it disapear/collapse
		// It is then absorbed by the selected empire
		if(numOfColoniesOfWeakestEmpire<=1)
		{
			// Update the positions of the colonies by adding the collapsed imperialist
			empiresList[selectedEmpireInd].setColoniesPosition( 
					concatenatePositions(
							empiresList[selectedEmpireInd].getColoniesPosition(), 
							empiresList[weakestEmpireInd].getImperialistPosition()
					) 
			);
			
			// Idem for the costs
			empiresList[selectedEmpireInd].setColoniesCost( 
					concatenateCosts(
							empiresList[selectedEmpireInd].getColoniesCost(), 
							empiresList[weakestEmpireInd].getImperialistCost()
					) 
			);
			
			// Erase the collapsed empire from the empires list
			deleteAnEmpire(weakestEmpireInd);
			//System.out.println("An empire has collapsed");
			//System.out.println("New number of empires: " + empiresList.length);
		}

	}


	
	/**
	 * Removes a position from the colony positions array
	 * @param colonyPositions
	 * @param indexToRemove
	 * @return the updated positions
	 */
	private double[][] removeColonyPosition(double[][] colonyPositions, int indexToRemove)
	{
		// Create a new array to store the updated positions
		double[][] newColonyPositions = new double[colonyPositions.length-1][colonyPositions[0].length];

		// Copy in it the positions before the colony to remove
		for(int i=0; i<indexToRemove; i++)
		{
			newColonyPositions[i] = colonyPositions[i];
		}

		// Then copy the rest of the positions, without including to colony to remove
		for(int j=indexToRemove; j<newColonyPositions.length; j++)
		{
			newColonyPositions[j] = colonyPositions[j+1];
		}

		// Return the updated positions
		return newColonyPositions;
	}

	

	/**
	 * Removes a position from the colony costs vector
	 * @param colonyCosts
	 * @param indexToRemove
	 * @return the updated costs vector
	 */
	private double[] removeColonyCost(double[] colonyCosts, int indexToRemove)
	{
		// Create a new vector to store the updated costs
		double[] newColonyCosts = new double[colonyCosts.length-1];

		// Copy in it the costs before the colony to remove
		for(int i=0; i<indexToRemove; i++)
		{
			newColonyCosts[i] = colonyCosts[i];
		}

		// Then copy the rest of the costs, without including to colony to remove
		for(int j=indexToRemove; j<newColonyCosts.length; j++)
		{
			newColonyCosts[j] = colonyCosts[j+1];
		}

		// Return the updated costs
		return newColonyCosts;
	}


	
	/**
	 * Concatenates the positions of an empire with an additionnal one
	 * @param positions1 the positions array of the empire
	 * @param position2 the position to add
	 * @return the updated positions
	 */
	private double[][] concatenatePositions(double[][] positions1, double[] position2)
	{
		// Create a new array to store the updated positions 
		double[][] newPositions = new double[positions1.length+1][positions1[0].length];

		// Add the positions of the existing empire in the array
		int i;
		for(i=0; i<positions1.length; i++)
		{
			newPositions[i] = positions1[i];
		}

		// Then add the position to add at the end
		newPositions[i] = position2;

		// Return the updated positions
		return newPositions;
	}

	

	/**
	 * Concatenates the costs of an empire with an additionnal one
	 * @param costs1 the costs vector of the empire
	 * @param cost2 the cost to add
	 * @return the updated costs
	 */
	private double[] concatenateCosts(double[] costs1, double cost2)
	{
		// Create a new vector to store the updated costs
		double[] newCosts = new double[costs1.length+1];

		// Add the costs of the existing empire in the array
		int i;
		for(i=0; i<costs1.length; i++)
		{
			newCosts[i] = costs1[i];
		}

		// Then add the cost to add at the end
		newCosts[i] = cost2;

		// Return the updated costs
		return newCosts;
	}


	
	/**
	 * Deletes an empire from the empires list
	 * @param indexToDelete
	 */
	private void deleteAnEmpire(int indexToDelete)
	{
		// Split the empires list into two sub lists, before and after the empire to remove
		Empire[] empiresList1 = Arrays.copyOfRange(empiresList, 0, indexToDelete);
		Empire[] empiresList2 = Arrays.copyOfRange(empiresList, indexToDelete+1, empiresList.length);

		// Create a new list with the updated size
		empiresList = new Empire[empiresList1.length+empiresList2.length];
		
		// Copy the empires of the sub lists into the new one
		for(int n=0; n<(empiresList1.length + empiresList2.length); n++)
		{
			if(n<empiresList1.length)
			{
				empiresList[n] = empiresList1[n];
			}

			if(n>= empiresList1.length)
			{
				empiresList[n] = empiresList2[n-empiresList1.length];
			}
		}
	}


	
	/**
	 * Selects an empire according to their probabilities 
	 * @param probability the probability vector
	 * @return the selected empire index
	 */
	private int selectAnEmpire(double[] probability)
	{
		// Create a vector of random numbers
		double[] randVector = new double[probability.length];
		for(int i=0; i<probability.length; i++)
		{
			randVector[i] = r.nextDouble();
		}
		
		// Substract to each element of this vector the corresponding 
		// value of the probability vector
		double[] dVector = new double[probability.length];
		for(int i=0; i<probability.length; i++)
		{
			dVector[i] = probability[i] - randVector[i];
		}
		
		// Return the index of the maximum value of the vector
		return utils.getMaxIndex(dVector);
	}



	/**
	 * Returns a string with informtion about the algorithm
	 * @return the string
	 */
	public String getDetails()
	{
		return "Imperialist Competitive Algorithm (ICA): " + 
		"as described in: Atashpaz-Gargari, E. and Lucas, C., Imperialist Competitive Algorithm: An Algorithm for Optimization Inspired by Imperialistic Competition, IEEE Congress on Evolutionary Computation, 2007, pp. 4661-4667." +
		" Adapted from: http://www.mathworks.com/matlabcentral/fileexchange/22046-imperialist-competitive-algorithm-ica";
	}



	/**
	 * Returns a string with the name of the algorithm
	 * @return the string
	 */
	public String getName() 
	{
		return "Imperialist Competitive Algorithm";
	}

}
