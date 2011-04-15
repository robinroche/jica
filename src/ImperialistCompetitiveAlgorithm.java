

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;


/**
 * @author Robin
 *
 */
public class ImperialistCompetitiveAlgorithm
{

	// ICA parameters
	int numOfCountries = 20;               		// Number of initial countries
	int numOfInitialImperialists = 8;      		// Number of initial imperialists
	int numOfAllColonies = numOfCountries - numOfInitialImperialists;
	int numOfDecades = 10;					// Number of decades (generations)
	double revolutionRate = 0.1;               	// Revolution is the process in which the socio-political characteristics of a country change suddenly
	double assimilationCoefficient = 2;        	// In the original paper assimilation coefficient is shown by "beta"
	double assimilationAngleCoefficient = .785; // In the original paper assimilation angle coefficient is shown by "gama"
	double zeta = 0.02;							// Total Cost of Empire = Cost of Imperialist + Zeta * mean(Cost of All Colonies)
	double dampRatio = 0.99;
	boolean stopIfJustOneEmpire = false;        // Use "true" to stop the algorithm when just one empire is remaining. Use "false" to continue the algorithm
	double unitingThreshold = 0.02;           	// The percent of search space size, which enables the uniting process of two empires


	// Variables
	protected long seed = System.currentTimeMillis();       
	Random r = new Random(seed);
	int problemDimension;						// Problem dimension
	double[] minBounds;							// Minimum bounds
	double[] maxBounds;							// Maximum bounds
	Empire[] empiresList = new Empire[numOfInitialImperialists];						// List of Empires
	double[][] initialCountries;				// List of countries
	double[] initialCosts;						// Costs of the countries
	double[] minimumCost = new double[numOfDecades];
	double[] meanCost = new double[numOfDecades];
	double[] searchSpaceSize;



	public ImperialistCompetitiveAlgorithm(Object[] args) 
	{
		problemDimension = (Integer) args[0]; 
		minBounds = (double[]) args[1];
		maxBounds = (double[]) args[2];
	}


	protected double[] runICA()
	{

		System.out.println("Algorithm launched");
		
		searchSpaceSize = new double[problemDimension];
		for(int i=0; i<problemDimension; i++)
		{
			searchSpaceSize[i] = maxBounds[i] - minBounds[i];
		}

	
		// Create an initial population of individuals (countries)
		generateNewCountries();
		System.out.println("New countries generated");
		
		// Compute the cost of each country: the lesser the cost, the more powerful the country is
		getCountriesCosts();
		System.out.println("Countries costs computed");

		// Sort the costs and the corresponding countries in assending order. The best countries will be in higher places.
		sortArray(initialCosts, initialCountries);
		System.out.println("Arrays sorted");

		// Create initial empires
		createInitialEmpires();
		System.out.println("Empires created");
		
		// While no stopping condition is met
		for(int decade=0; decade<numOfDecades; decade++)
		{
			System.out.println("Decade: " + decade);
			
			revolutionRate = dampRatio * revolutionRate;

			for (int i=0; i<empiresList.length; i++)
			{
				// Assimilation;  Movement of Colonies Toward Imperialists (Assimilation Policy)
				assimilateColonies(empiresList[i]);

				// Revolution;  A Sudden Change in the Socio-Political Characteristics
				revolveColonies(empiresList[i]);

				// New Cost Evaluation
				empiresList[i].setColoniesCost( getCountriesCosts( empiresList[i].getColoniesPosition()) );
				System.out.println("empire " + i + " colonies cost empire = " + mean(empiresList[i].getColoniesCost()));
				
				// Empire Possession
				possesEmpire(empiresList[i]);

				// Computation of Total Cost for Empires
				empiresList[i].setTotalCost( empiresList[i].getImperialistCost() + zeta * mean(empiresList[i].getColoniesCost()) );
				System.out.println("empire " + i + " totalcost of empire = " + empiresList[i].getTotalCost());
			}

			// Uniting Similiar Empires
			uniteSimilarEmpires();

			// Imperialistic Competition
			imperialisticCompetition();

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
			minimumCost[decade] = imperialistCosts[min(imperialistCosts)];
			System.out.println("Minimum decade cost: " + minimumCost[decade]);
			meanCost[decade] = mean(imperialistCosts);

		}


		double bestCost = minimumCost[min(minimumCost)];
		int bestIndex = min(minimumCost);
		double[] bestSolution = empiresList[bestIndex].getImperialistPosition();

		System.out.println("Best solution: ");
		for(int i=0; i<bestSolution.length; i++)
		{
			 System.out.println(bestSolution[i]);
		}
		System.out.println("Best fitness: " + bestCost);
		
		return bestSolution;

	}



	private double[] getCountriesCosts(double[][] coloniesPosition) 
	{
		double[] costs = new double[coloniesPosition.length];
		
		for(int i=0; i<coloniesPosition.length; i++)
		{
			costs[i] = getCountryCost(coloniesPosition[i]);
		}
		return costs;
	}


	private void sortArray(final double[] arrayToSort, double[][] matchingArray) 
	{
		Integer[] sortOrder = new Integer[arrayToSort.length];

		// Create index array.
		for(int i=0; i<sortOrder.length; i++)
		{
			sortOrder[i] = i;
		}

		Arrays.sort(sortOrder,new Comparator<Integer>() 
				{   
			public int compare(Integer a, Integer b)
			{
				double delta = arrayToSort[b]-arrayToSort[a];
			    if(delta > 0) return 1;
			    if(delta < 0) return -1;
			    return 0;
			}});
		
		double[] arrayToSortCopy = arrayToSort.clone();
        double[][] matchingArrayCopy = matchingArray.clone();
            
		for(int i=0;i<sortOrder.length;i++)
		{
			initialCosts[i] = arrayToSortCopy[sortOrder[i]];
			initialCountries[i] = matchingArrayCopy[sortOrder[i]];
        }
	}


	private void getCountriesCosts() 
	{
		initialCosts = new double[numOfCountries];

		for(int i=0; i<numOfCountries; i++)
		{
			initialCosts[i] = getCountryCost(initialCountries[i]);
		}
	}



	private double getCountryCost(double[] country) 
	{
		// TODO Auto-generated method stub
		double fitness = 0;
		
		for(int i=0; i<country.length; i++)
		{
			fitness = fitness + country[i]*country[i];
		}
		
		return fitness;
	}



	private void generateNewCountries() 
	{
		initialCountries = new double[numOfCountries][problemDimension];  

		for(int i=0; i<numOfCountries; i++)
		{
			for(int j=0; j<problemDimension; j++)
			{
				initialCountries[i][j] = (maxBounds[j] - minBounds[j]) * r.nextDouble() + minBounds[j];//TODO
				//initialCountries[i][j] = (maxBounds[j] - minBounds[j]) * Math.random() + minBounds[j];
			}
		}
		
		for(int i=0; i<numOfCountries; i++)
		{
			System.out.println("initialcountries i = " + i + " = " + Arrays.toString(initialCountries[i]));
		}
	}



	private double[][] extractPartOfArray(double[][] array, int startIndex, int endIndex)
	{
		int secondDimSize = array[0].length;
		double[][] arrayExtract = new double[endIndex-startIndex][secondDimSize];

		int newIndex = 0;
		for(int i=startIndex; i<endIndex ; i++)
		{
			for(int j=0; j<secondDimSize; j++)
			{
				arrayExtract[newIndex][j] = array[i][j];
			}
			newIndex++;
		}

		return arrayExtract;
	}


	private void createInitialEmpires()
	{

		// Extract the best countries to create empires
		double[][] allImperialistsPosition = extractPartOfArray(initialCountries, 0, numOfInitialImperialists);

		// Extract their costs
		double[] allImperialistsCost = new double[numOfInitialImperialists];
		System.arraycopy(initialCosts, 0, allImperialistsCost, 0, numOfInitialImperialists);

		// Extract the rest to create colonies
		double[][] allColoniesPosition = extractPartOfArray(initialCountries, numOfInitialImperialists, initialCountries.length);

		// Extract their costs
		double[] allColoniesCost = new double[initialCosts.length-numOfInitialImperialists]; 
		System.arraycopy(initialCosts, numOfInitialImperialists, allColoniesCost, 0, initialCosts.length-numOfInitialImperialists);	

		System.out.println("Extractions done");
		
		// Compute the power of imperialists
		double[] allImperialistsPower = new double[numOfInitialImperialists];

		if(allImperialistsCost[max(allImperialistsCost)]>0)
		{
			for(int i=0; i<allImperialistsCost.length; i++)
			{
				allImperialistsPower[i] = 1.3 * allImperialistsCost[max(allImperialistsCost)] - allImperialistsCost[i];
			}
		}
		else
		{
			for(int i=0; i<allImperialistsCost.length; i++)
			{
				allImperialistsPower[i] = 0.7 * allImperialistsCost[max(allImperialistsCost)] - allImperialistsCost[i];
			}
		}

		System.out.println("Powers computed");
		
		// Set the number of colonies for the imperialists 
		int[] allImperialistNumOfColonies = new int[numOfInitialImperialists];
		for(int i=0; i<allImperialistsPower.length; i++)
		{
			allImperialistNumOfColonies[i] = (int) Math.round(allImperialistsPower[i]/sum(allImperialistsPower) * numOfAllColonies);
		}
		allImperialistNumOfColonies[allImperialistNumOfColonies.length-1] = Math.max(
			numOfAllColonies - 
			sum(Arrays.copyOfRange(allImperialistNumOfColonies, 0, allImperialistNumOfColonies.length-1)),
			0);//TODO: à revoir

		System.out.println("Colonies created");
		
		// Create a random permutation of integers
		int[] randomIndex = randperm(numOfAllColonies);

		// Initialize empires
		for(int i=0; i<numOfInitialImperialists; i++)
		{
			empiresList[i] = new Empire();
		}
		
		// Set the position of the last imperialist
		double[] zeros = new double[problemDimension];
		Arrays.fill(zeros, 0);

		empiresList[numOfInitialImperialists-1].setImperialistPosition(zeros);

		System.out.println("Last imperialist set");
		
		// Create imperialist empires
		for(int i=0; i<numOfInitialImperialists; i++)
		{

			empiresList[i].setImperialistPosition(allImperialistsPosition[i]);
			empiresList[i].setImperialistCost(allImperialistsCost[i]);

			int[] R = Arrays.copyOfRange(randomIndex, 0, allImperialistNumOfColonies[i]);
			// randomIndex(allImperialistNumOfColonies(i)+1:end);

			empiresList[i].init(R.length);

			empiresList[i].setColoniesPosition(extractGivenArrayParts(allColoniesPosition, R));
			empiresList[i].setColoniesCost(extractGivenArrayParts(allColoniesCost, R));
			empiresList[i].setTotalCost(empiresList[i].getImperialistCost() + zeta * mean(empiresList[i].getColoniesCost()));
		}

		System.out.println("Imperialists created");
		
		for(int i=0; i<empiresList.length; i++)
		{
			if(empiresList[i].getColoniesPosition().length == 0)
			{
				empiresList[i].setColoniesPosition(generateNewCountry(1));
				// TODO not always necessary?
				//empiresList[i].setColoniesCost( getCountriesCosts(empiresList[i].getColoniesPosition()));
			}
		}

	}


	private double[][] generateNewCountry(int number) 
	{
		double[][] countries = new double[number][problemDimension];  

		for(int i=0; i<number; i++)
		{
			for(int j=0; j<problemDimension; j++)
			{
				countries[i][j] = (maxBounds[j] - minBounds[j]) * r.nextDouble() + minBounds[j];
			}
		}
		return countries;
	}


	public static double mean(double[] p) 
	{
		double sum = 0;
		for (int i=0; i<p.length; i++) 
		{
			sum += p[i];
		}
		return sum / p.length;
	}


	private double[][] extractGivenArrayParts(double[][] array, int[] selectedIndexes) 
	{
		double[][] arrayExtract = new double[selectedIndexes.length][array[0].length];

		for(int i=0; i<selectedIndexes.length; i++)
		{
			int index = selectedIndexes[i];
			for(int j=0; j<array[0].length; j++)
			{
				arrayExtract[i][j] = array[index][j];
			}
		}

		return arrayExtract;
	}


	private double[] extractGivenArrayParts(double[] array, int[] selectedIndexes) 
	{
		double[] arrayExtract = new double[selectedIndexes.length];

		for(int i=0; i<selectedIndexes.length; i++)
		{
			int index = selectedIndexes[i];
			arrayExtract[i] = array[index];
		}

		return arrayExtract;
	}


	private int[] randperm(int size) 
	{
		ArrayList<Integer> nVector = new ArrayList<Integer>();
		for(int i=0; i<size; i++)
		{
			nVector.add(i);
		}
		
		int[] outputVector = new int[size];
		
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



	private double sum(double[] array)
	{
		double sum = 0;
		for (double i : array) 
		{
			sum += i;
		}
		return sum;
	}


	private int sum(int[] array)
	{
		int sum = 0;
		for (int i : array) 
		{
			sum += i;
		}
		return sum;
	}


	// Empire class
	private class Empire
	{

		public Empire()
		{
			super();
		}
		
		// Empire characteristics
		double[] imperialistPosition = new double[problemDimension];
		double imperialistCost;
		double[][] coloniesPosition;
		double[] coloniesCost;
		double totalCost;

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


	private int max(double[] values)
	{
		double max = Double.MIN_VALUE;
		int i;
		int bestIndex = 0;
		for(i=0; i<values.length; i++) 
		{
			if(values[i] > max)
			{
				max = values[i];
				bestIndex = i;
			}
		}
		return bestIndex;
	}



	private int min(double[] values)
	{
		double min = Double.MAX_VALUE;
		int i;
		int bestIndex = 0;
		for(i=0; i<values.length; i++) 
		{
			if(values[i] < min)
			{
				min = values[i];
				bestIndex = i;
			}
		}
		return bestIndex;
	}



	private void assimilateColonies(Empire theEmpire)
	{

		int numOfColonies = theEmpire.getColoniesPosition().length;

		double[][] vector = repmat(theEmpire.getImperialistPosition(),numOfColonies);
		for(int i=0; i<vector.length; i++)
		{
			for(int j=0; j<vector[0].length; j++)
			{
				vector[i][j] = vector[i][j] - theEmpire.getColoniesPosition()[i][j];
			}
		}

		double[][] coloniesPosition = new double[vector.length][vector[0].length];

		for(int i=0; i<vector.length; i++)
		{
			for(int j=0; j<vector[0].length; j++)
			{
				coloniesPosition[i][j] = theEmpire.getColoniesPosition()[i][j] + 2 * assimilationCoefficient * r.nextInt(vector.length) * vector[i][j];
			}
		}

		theEmpire.setColoniesPosition(coloniesPosition);

		double[][] minVarMatrix = repmat(minBounds,numOfColonies);
		double[][] maxVarMatrix = repmat(maxBounds,numOfColonies);

		theEmpire.setColoniesPosition(max(theEmpire.getColoniesPosition(),minVarMatrix));
		theEmpire.setColoniesPosition(min(theEmpire.getColoniesPosition(),maxVarMatrix));
	}



	private double[][] max(double[][] array1, double[][] array2) 
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


	private double[][] min(double[][] array1, double[][] array2) 
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


	private double[][] repmat(double[] motif, int number) 
	{
		double[][] patchwork = new double[number][motif.length];

		for(int i=0; i<number; i++)
		{
			for(int j=0; j<motif.length; j++)
			{
				patchwork[i][j] = motif[j];
			}
		}
		return patchwork;
	}


	private void revolveColonies(Empire theEmpire)
	{
		int numOfRevolvingColonies = (int) Math.round((revolutionRate * theEmpire.getColoniesCost().length));

		double[][] revolvedPosition = generateNewCountry(numOfRevolvingColonies);

		int[] R = randperm(theEmpire.getColoniesCost().length);
		R = Arrays.copyOfRange(R, 0, numOfRevolvingColonies);

		for(int i=0; i<R.length; i++)
		{
			theEmpire.setColonyPosition(R[i], revolvedPosition[i]);
		}
	}



	private void possesEmpire(Empire theEmpire)
	{
		double[] coloniesCost = theEmpire.getColoniesCost();

		int bestColonyInd = min(coloniesCost);
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
		double theresholdDistance = unitingThreshold * norm(searchSpaceSize);
		int numOfEmpires = empiresList.length;

		for(int i=0; i<(numOfEmpires-1); i++)
		{
			for(int j=i; j<numOfEmpires; j++)
			{

				// Compute the distance between two empires
				double[] distanceVector = new double[empiresList[i].getImperialistPosition().length];
				for(int k=0; k<empiresList[i].getImperialistPosition().length; k++)
				{
					distanceVector[k] = empiresList[i].getImperialistPosition()[k] - empiresList[j].getImperialistPosition()[k];
				}
				double distance = norm(distanceVector);

				// If the empires are too close
				if(distance<=theresholdDistance)
				{
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
					int newSize = 
						empiresList[betterEmpireInd].getColoniesPosition().length + 
						1 + 
						empiresList[worseEmpireInd].getColoniesPosition().length;

					double[][] newColoniesPosition = new double[newSize][empiresList[betterEmpireInd].getColoniesPosition()[0].length];

					int m;

					for(m=0; m<empiresList[betterEmpireInd].getColoniesPosition().length; m++)
					{
						newColoniesPosition[m] = empiresList[betterEmpireInd].getColoniesPosition()[m];
					}
					
					newColoniesPosition[m] = empiresList[worseEmpireInd].getImperialistPosition();
					int index=m+1;	
					for(m=index; m<newSize; m++)
					{
						newColoniesPosition[m] = empiresList[worseEmpireInd].getColoniesPosition()[m-empiresList[betterEmpireInd].getColoniesPosition().length-1];
					}

					empiresList[betterEmpireInd].setColoniesPosition(newColoniesPosition);

					// Update the costs
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
					int index2=m2+1;	
					for(m2=index2; m2<newSize2; m2++)
					{
						newColoniesCost[m2] = empiresList[worseEmpireInd].getColoniesCost()[m2-empiresList[betterEmpireInd].getColoniesCost().length-1];
					}

					empiresList[betterEmpireInd].setColoniesCost(newColoniesCost);


					// Update TotalCost for new United Empire                                     
					empiresList[betterEmpireInd].setTotalCost(
							empiresList[betterEmpireInd].getImperialistCost() +
							zeta * mean(empiresList[betterEmpireInd].getColoniesCost())
					);

					// Update the empires list
					deleteAnEmpire(worseEmpireInd);

					return;
				}

			}
		}  
	}



	private double norm(double[] vector) 
	{
		double sum = 0;		
		for(int i=0; i<vector.length; i++)
		{
			sum = sum + Math.pow(vector[i],2);
		}
		return Math.sqrt(sum);
	}


	private void imperialisticCompetition()
	{
		if(r.nextDouble() > .11)
		{
			return;//TODO
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
		int weakestEmpireInd = max(totalCosts);
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
			possessionProbability[i] = totalPowers[i] / sum(totalPowers);
		}
		System.out.println("totalpowers = " + Arrays.toString(totalPowers));//TODO
		System.out.println("proba = " + Arrays.toString(possessionProbability));//TODO
		
		// Select an empire according to their probabilities
		int selectedEmpireInd = selectAnEmpire(possessionProbability);
System.out.println("selectedind = " + selectedEmpireInd);
		
		// Generate a random integer
		int nn = empiresList[weakestEmpireInd].getColoniesCost().length;
		int jj = myRandInt(nn);

		empiresList[selectedEmpireInd].setColoniesPosition(	concetenatePositions(empiresList[selectedEmpireInd].getColoniesPosition(), empiresList[weakestEmpireInd].getColoniesPosition()[jj]) );

		empiresList[selectedEmpireInd].setColoniesCost( concatenateCosts(empiresList[selectedEmpireInd].getColoniesCost(), empiresList[weakestEmpireInd].getColoniesCost()[jj]) );

		empiresList[weakestEmpireInd].setColoniesPosition( removeColonyPosition(empiresList[weakestEmpireInd].getColoniesPosition(), jj) );

		empiresList[weakestEmpireInd].setColoniesCost( removeColonyCost(empiresList[weakestEmpireInd].getColoniesCost(), jj) );

		// Collapse of the the weakest colony-less Empire
		nn = empiresList[weakestEmpireInd].getColoniesCost().length;
		if(nn<=1)
		{
			empiresList[selectedEmpireInd].setColoniesPosition( concetenatePositions(empiresList[selectedEmpireInd].getColoniesPosition(), empiresList[weakestEmpireInd].getImperialistPosition()) );

			empiresList[selectedEmpireInd].setColoniesCost( concatenateCosts(empiresList[selectedEmpireInd].getColoniesCost(), empiresList[weakestEmpireInd].getImperialistCost()) );

			// Update the empires list
			deleteAnEmpire(weakestEmpireInd);
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
		Empire[] empiresList1 = Arrays.copyOfRange(empiresList, 0, indexToDelete);
		Empire[] empiresList2 = Arrays.copyOfRange(empiresList, indexToDelete+1, empiresList.length);

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
		
		index =  max(dVector);
		return index;
	}


	// This functions creates random numbers between [1 MaxInt] (1 itself and MaxInt itself)
	private int myRandInt(int maxInt)
	{
		int y = (int) Math.ceil(r.nextDouble()*maxInt);
		return y;
	}



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
