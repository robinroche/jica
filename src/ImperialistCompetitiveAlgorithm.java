

import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;


/**
 * This class contains the ICA algorithm methods
 * @author Robin Roche
 */
public class ImperialistCompetitiveAlgorithm
{

	// ICA parameters
	int numOfCountries = 30;               		// Number of initial countries
	int numOfInitialImperialists = 8;      		// Number of initial imperialists
	int numOfAllColonies = numOfCountries - numOfInitialImperialists;
	int numOfDecades = 100;					// Number of decades (generations)
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
	public ImperialistCompetitiveAlgorithm(Object[] args) 
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
		
		// Create initial empires
		createInitialEmpires();

		
		// While no stopping condition is met
		for(int decade=0; decade<numOfDecades; decade++)
		{
			System.out.println("Decade: " + decade);
			
			revolutionRate = dampRatio * revolutionRate;

			System.out.println("Number of empires: " + empiresList.length);
			
			int currentNumOfCountries = 0;
			for(int i=0; i<empiresList.length; i++)
			{
				currentNumOfCountries = currentNumOfCountries + empiresList[i].getColoniesPosition().length;
				System.out.println("Empire[" + i + "]   Colonies: " + empiresList[i].getColoniesPosition().length);
			}
			System.out.println("Total number of colonies: " + currentNumOfCountries);
			System.out.println("Total number of countries: " + (currentNumOfCountries + empiresList.length));
			
			
			for (int i=0; i<empiresList.length; i++)
			{
				
				// Assimilation;  Movement of Colonies Toward Imperialists (Assimilation Policy)
				assimilateColonies(empiresList[i]);
				//System.out.println("assimilate colonies done");
				
				// Revolution;  A Sudden Change in the Socio-Political Characteristics
				revolveColonies(empiresList[i]);
				//System.out.println("revolve colonies done");
				
				// New Cost Evaluation
				empiresList[i].setColoniesCost( getCountriesCosts( empiresList[i].getColoniesPosition()) );
				//System.out.println("empire " + i + " colonies cost empire = " + mean(empiresList[i].getColoniesCost()));
				
				// Empire Possession (a strong colony becomes the imperialist)
				possesEmpire(empiresList[i]);
				//System.out.println("possess empires done");
				
				// Computation of Total Cost for Empires
				empiresList[i].setTotalCost( empiresList[i].getImperialistCost() + zeta * utils.getMean(empiresList[i].getColoniesCost()) );
				//System.out.println("empire " + i + " totalcost of empire = " + empiresList[i].getTotalCost());
			}

			// Uniting Similiar Empires
			uniteSimilarEmpires();
			//System.out.println("unite empires done");
			
			// Imperialistic Competition
			imperialisticCompetition();
			//System.out.println("competition done");
			
			
			if (empiresList.length == 1 && stopIfJustOneEmpire)
			{
				break;
			}

			// Extract results of the round
			double[] imperialistCosts = new double[empiresList.length];
			for(int i=0; i<empiresList.length ; i++)
			{
				imperialistCosts[i] = empiresList[i].getImperialistCost();
			}
			minimumCost[decade] = imperialistCosts[utils.getMinIndex(imperialistCosts)];
			//System.out.println("Minimum decade cost: " + minimumCost[decade]);
			meanCost[decade] = utils.getMean(imperialistCosts);
			
			
			bestDecadePosition[decade] = empiresList[utils.getMinIndex(imperialistCosts)].getImperialistPosition();

		}


		double bestCost = minimumCost[utils.getMinIndex(minimumCost)];
		int bestIndex = utils.getMinIndex(minimumCost);
		//double[] bestSolution = empiresList[bestIndex].getImperialistPosition();
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
		utils.printArray("allImperialistNumOfColonies", allImperialistNumOfColonies);

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
			
			utils.printEmpire(empiresList[i], i);
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
				array[i][j] = repmatArray[i][j] - theEmpire.getColoniesPosition()[i][j];//TODO: problem here
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
		
		utils.printArray("theEmpire.getColoniesPosition()", theEmpire.getColoniesPosition());

		// Bound the values with the min and max bounds
		double[][] minVarMatrix = utils.repmat(minBounds,numOfColonies);
		double[][] maxVarMatrix = utils.repmat(maxBounds,numOfColonies);
		theEmpire.setColoniesPosition(utils.max(theEmpire.getColoniesPosition(),minVarMatrix));
		theEmpire.setColoniesPosition(utils.min(theEmpire.getColoniesPosition(),maxVarMatrix));
		utils.printArray("theEmpire.getColoniesPosition()", theEmpire.getColoniesPosition());
	}







	private void revolveColonies(Empire theEmpire)
	{
		int numOfRevolvingColonies = (int) Math.round((revolutionRate * theEmpire.getColoniesCost().length));

		double[][] revolvedPosition = generateNewCountries(numOfRevolvingColonies, problemDimension, minBounds, maxBounds, r);

		//System.out.println("about to run randperm : size = " + theEmpire.getColoniesCost().length + " values = " + Arrays.toString(theEmpire.getColoniesCost()));
		int[] R = utils.randperm(theEmpire.getColoniesCost().length, r);
		R = Arrays.copyOfRange(R, 0, numOfRevolvingColonies);
		//System.out.println("R built : " + Arrays.toString(R));

		for(int i=0; i<R.length; i++)
		{
			theEmpire.setColonyPosition(R[i], revolvedPosition[i]);
		}
	}



	private void possesEmpire(Empire theEmpire)
	{
		double[] coloniesCost = theEmpire.getColoniesCost();

		int bestColonyInd = utils.getMinIndex(coloniesCost);
		double minColoniesCost = coloniesCost[bestColonyInd]; 

		if(minColoniesCost < theEmpire.getImperialistCost())
		{
			double[] oldImperialistPosition = theEmpire.getImperialistPosition();
			double oldImperialistCost = theEmpire.getImperialistCost();

			theEmpire.setImperialistPosition(theEmpire.getColoniesPosition()[bestColonyInd]);
			theEmpire.setImperialistCost(theEmpire.getColoniesCost()[bestColonyInd]);

			theEmpire.setColoniesPosition(bestColonyInd, oldImperialistPosition);
			theEmpire.setColoniesCost(bestColonyInd, oldImperialistCost);
		}
	}



	private void uniteSimilarEmpires()
	{
		double theresholdDistance = unitingThreshold * utils.getNorm(searchSpaceSize);
		int numOfEmpires = empiresList.length;

		for(int i=0; i<(numOfEmpires-1); i++)
		{
			for(int j=i+1; j<numOfEmpires; j++)
			{

				// Compute the distance between two empires
				double[] distanceVector = new double[empiresList[i].getImperialistPosition().length];
				for(int k=0; k<empiresList[i].getImperialistPosition().length; k++)
				{
					distanceVector[k] = empiresList[i].getImperialistPosition()[k] - empiresList[j].getImperialistPosition()[k];
				}
				double distance = utils.getNorm(distanceVector);

				// If the empires are too close
				if(distance<=theresholdDistance)
				{
					System.out.println("distance = " + distance + " thresesholdDistance = " + theresholdDistance + "fusion of empires " + i + "+" + j);
					
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


					// Update TotalCost for new United Empire                                     
					empiresList[betterEmpireInd].setTotalCost(
							empiresList[betterEmpireInd].getImperialistCost() +
							zeta * utils.getMean(empiresList[betterEmpireInd].getColoniesCost())
					);

					// Update the empires list
					deleteAnEmpire(worseEmpireInd);
					System.out.println("New number of empires: " + empiresList.length);

					return;
				}

			}
		}  
	}



	private double[] getColonyCostsOfUnitedEmpire(int betterEmpireInd, int worseEmpireInd) 
	{
		int newSize2 = 
			empiresList[betterEmpireInd].getColoniesCost().length + 
			1 + 
			empiresList[worseEmpireInd].getColoniesCost().length;

		double[] newColoniesCost = new double[newSize2];

		int m2;

		for(m2=0; m2<empiresList[betterEmpireInd].getColoniesCost().length; m2++)
		{
			newColoniesCost[m2] = empiresList[betterEmpireInd].getColoniesCost()[m2];
		}

		newColoniesCost[m2] = empiresList[worseEmpireInd].getImperialistCost();
		
		int index2;	
		for(index2=m2+1; index2<newSize2; index2++)
		{
			newColoniesCost[index2] = empiresList[worseEmpireInd].getColoniesCost()[index2-empiresList[betterEmpireInd].getColoniesCost().length-1];
		}
		return newColoniesCost;
	}


	private double[][] getColonyPositionsOfUnitedEmpire(int betterEmpireInd, int worseEmpireInd) 
	{
		
		System.out.println("Uniting empire colonies positions:");
		System.out.println("best colonies: " + empiresList[betterEmpireInd].getColoniesPosition().length + " worst colonies = " + empiresList[worseEmpireInd].getColoniesPosition().length);
		
		// The new size = the best empire's colony count + the weakest empire's colony's count + its imperialist
		int newSize = 
			empiresList[betterEmpireInd].getColoniesPosition().length + 
			1 + 
			empiresList[worseEmpireInd].getColoniesPosition().length;

		double[][] newColoniesPosition = new double[newSize][problemDimension];

		int m;
		for(m=0; m<empiresList[betterEmpireInd].getColoniesPosition().length; m++)
		{
			newColoniesPosition[m] = empiresList[betterEmpireInd].getColoniesPosition()[m];
		}
		
		newColoniesPosition[m] = empiresList[worseEmpireInd].getImperialistPosition();
		
		int index;	
		for(index=m+1; index<newSize; index++)
		{
			newColoniesPosition[index] = empiresList[worseEmpireInd].getColoniesPosition()[index-empiresList[betterEmpireInd].getColoniesPosition().length-1];
		}
		
		System.out.println("united empire colonies: " + newColoniesPosition.length);
		
		return newColoniesPosition;
	
	}





	private void imperialisticCompetition()
	{
		double rand = r.nextDouble();
		if(rand > .11)
		{
			return;
		}
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
		System.out.println("totalcosts = " + Arrays.toString(totalCosts));//TODO
		
		// Get the weakest empire and its cost
		int weakestEmpireInd = utils.getMaxIndex(totalCosts);
		double maxTotalCost = totalCosts[weakestEmpireInd]; 

		System.out.println("weakestind = " + weakestEmpireInd + " cost = " + maxTotalCost);//TODO
		
		// Get the power of each empire
		double[] totalPowers = new double[empiresList.length];
		for(int i=0; i<empiresList.length; i++)
		{
			totalPowers[i] = maxTotalCost - totalCosts[i];
		}

		// Get the possession probability of each empire
		double[] possessionProbability = new double[empiresList.length];
		for(int i=0; i<empiresList.length; i++)
		{
			possessionProbability[i] = totalPowers[i] / utils.getSum(totalPowers);
		}
		System.out.println("totalpowers = " + Arrays.toString(totalPowers));//TODO
		System.out.println("proba = " + Arrays.toString(possessionProbability));//TODO
		
		// Select an empire according to their probabilities
		int selectedEmpireInd = selectAnEmpire(possessionProbability);
		System.out.println("selectedind = " + selectedEmpireInd);//TODO
		
		// Generate a random integer
		int nn = empiresList[weakestEmpireInd].getColoniesCost().length;
		int jj = r.nextInt(nn);

		System.out.println("jj = " + jj);
		System.out.println("colonies length = " + empiresList[weakestEmpireInd].getColoniesPosition().length);
		
		empiresList[selectedEmpireInd].setColoniesPosition(	concetenatePositions(empiresList[selectedEmpireInd].getColoniesPosition(), empiresList[weakestEmpireInd].getColoniesPosition()[jj]) );

		empiresList[selectedEmpireInd].setColoniesCost( concatenateCosts(empiresList[selectedEmpireInd].getColoniesCost(), empiresList[weakestEmpireInd].getColoniesCost()[jj]) );

		empiresList[weakestEmpireInd].setColoniesPosition( removeColonyPosition(empiresList[weakestEmpireInd].getColoniesPosition(), jj) );

		empiresList[weakestEmpireInd].setColoniesCost( removeColonyCost(empiresList[weakestEmpireInd].getColoniesCost(), jj) );

		// Collapse of the the weakest colony-less Empire
		nn = empiresList[weakestEmpireInd].getColoniesCost().length;
		System.out.println("weakest empire has: " + nn + " colonies");
		if(nn<=1)
		{
			empiresList[selectedEmpireInd].setColoniesPosition( concetenatePositions(empiresList[selectedEmpireInd].getColoniesPosition(), empiresList[weakestEmpireInd].getImperialistPosition()) );

			empiresList[selectedEmpireInd].setColoniesCost( concatenateCosts(empiresList[selectedEmpireInd].getColoniesCost(), empiresList[weakestEmpireInd].getImperialistCost()) );

			// Update the empires list
			deleteAnEmpire(weakestEmpireInd);
			System.out.println("New number of empires: " + empiresList.length);
		}

	}


	private double[][] removeColonyPosition(double[][] colonyPositions, int indexToRemove)
	{
		double[][] newColonyPositions = new double[colonyPositions.length-1][colonyPositions[0].length];

		for(int i=0; i<indexToRemove; i++)
		{
			newColonyPositions[i] = colonyPositions[i];
		}

		for(int j=indexToRemove; j<newColonyPositions.length; j++)
		{
			newColonyPositions[j] = colonyPositions[j+1];
		}

		return newColonyPositions;
	}


	private double[] removeColonyCost(double[] colonyCosts, int indexToRemove)
	{
		double[] newColonyCosts = new double[colonyCosts.length-1];

		for(int i=0; i<indexToRemove; i++)
		{
			newColonyCosts[i] = colonyCosts[i];
		}

		for(int j=indexToRemove; j<newColonyCosts.length; j++)
		{
			newColonyCosts[j] = colonyCosts[j+1];
		}

		return newColonyCosts;
	}


	private double[][] concetenatePositions(double[][] positions1, double[] position2)
	{
		System.out.println("conatenating colony positions: " + positions1.length);
		
		int newSize = positions1.length+1;
		double[][] newPositions = new double[newSize][positions1[0].length];

		int i;

		for(i=0; i<positions1.length; i++)
		{
			newPositions[i] = positions1[i];
		}

		newPositions[i] = position2;

		return newPositions;
	}


	private double[] concatenateCosts(double[] costs1, double cost2)
	{
		int newSize = costs1.length+1;
		double[] newCosts = new double[newSize];

		int i;

		for(i=0; i<costs1.length; i++)
		{
			newCosts[i] = costs1[i];
		}

		newCosts[i] = cost2;

		return newCosts;
	}


	private void deleteAnEmpire(int indexToDelete)
	{
		System.out.println("Deleting an empire");
		
		Empire[] empiresList1 = Arrays.copyOfRange(empiresList, 0, indexToDelete);
		Empire[] empiresList2 = Arrays.copyOfRange(empiresList, indexToDelete+1, empiresList.length);

		empiresList = new Empire[empiresList1.length+empiresList2.length];
		
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


	private int selectAnEmpire(double[] probability)
	{
		int index;
		
		//System.out.println("proba = " + Arrays.toString(probability));//TODO
		
		double[] randVector = new double[probability.length];
		for(int i=0; i<probability.length; i++)
		{
			randVector[i] = r.nextDouble();
		}
		System.out.println("randvector = " + Arrays.toString(randVector));//TODO
		
		double[] dVector = new double[probability.length];
		for(int i=0; i<probability.length; i++)
		{
			dVector[i] = probability[i] - randVector[i];
		}
		System.out.println("dvector = " + Arrays.toString(dVector));//TODO
		
		index =  utils.getMaxIndex(dVector);
		return index;
	}


	// This functions creates random numbers between [1 MaxInt] (1 itself and MaxInt itself)



	public String getDetails()
	{
		return "Imperialist Competitive Algorithm (ICA): " + 
		"as described in: Atashpaz-Gargari, E. and Lucas, C., Imperialist Competitive Algorithm: An Algorithm for Optimization Inspired by Imperialistic Competition, IEEE Congress on Evolutionary Computation, 2007, pp. 4661-4667." +
		" Adapted from: http://www.mathworks.com/matlabcentral/fileexchange/22046-imperialist-competitive-algorithm-ica";
	}



	public String getName() 
	{
		return "Imperialist Competitive Algorithm";
	}

}
