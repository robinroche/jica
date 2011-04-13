

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.Random;


/**
 * @author Robin
 *
 */
public class ImperialistCompetitiveAlgorithm
{

	// Variables
	protected long seed = System.currentTimeMillis();       
	Random r = new Random(seed);
	int problemDimension;						// Problem dimension
	double[] minBounds;							// Minimum bounds
	double[] maxBounds;							// Maximum bounds
	 
	
	// ICA parameters
	int numOfCountries = 80;               		// Number of initial countries
	int numOfInitialImperialists = 8;      		// Number of initial imperialists
	int numOfAllColonies = numOfCountries - numOfInitialImperialists;
	int numOfDecades = 1500;					// Number of decades (generations)
	double revolutionRate = 0.1;               	// Revolution is the process in which the socio-political characteristics of a country change suddenly
	double assimilationCoefficient = 2;        	// In the original paper assimilation coefficient is shown by "beta"
	double assimilationAngleCoefficient = .785; // In the original paper assimilation angle coefficient is shown by "gama"
	double zeta = 0.02;							// Total Cost of Empire = Cost of Imperialist + Zeta * mean(Cost of All Colonies)
	double dampRatio = 0.99;
	boolean stopIfJustOneEmpire = false;        // Use "true" to stop the algorithm when just one empire is remaining. Use "false" to continue the algorithm
	double unitingThreshold = 0.02;           	// The percent of search space size, which enables the uniting process of two empires


	
	public ImperialistCompetitiveAlgorithm(Object[] args) 
	{
		problemDimension = (Integer) args[0]; 
		minBounds = (double[]) args[1];
		maxBounds = (double[]) args[2];
	}
	
	
	
	protected void runICA()
	{

		// Create an initial population of individuals (countries)
		ArrayList<double[]> initialCountries = generateNewCountries(problemDimension, minBounds, maxBounds, r);
		
		// Compute the cost of each country: the lesser the cost, the more powerful the country is
		double[] initialCosts = getCountriesCosts(initialCountries);
		
		// Sort the costs in assending order. The best countries will be in higher places.
		java.util.Arrays.sort(initialCosts);
		
		// Then sort the population with respect to their cost.
		//TODO: [initialCost,sortInd] = sort(initialCost);
		//TODO: initialCountries = initialCountries(sortInd,:);

		// Create initial empires
		ArrayList<Empire> empiresList = createInitialEmpires(initialCountries, initialCosts);


		// While no stopping condition is met
		for(int decade=0; decade<numOfDecades; decade++)
		{
			revolutionRate = dampRatio * revolutionRate;

			for (int i=0; i<empiresList.size(); i++)
			{
				// Assimilation;  Movement of Colonies Toward Imperialists (Assimilation Policy)
				empiresList(ii) = assimilateColonies(empires(ii));

				// Revolution;  A Sudden Change in the Socio-Political Characteristics
				empiresList(ii) = revolveColonies(empires(ii));

				// New Cost Evaluation
				empiresList(i).ColoniesCost = feval(ProblemParams.CostFuncName,empiresList(i).ColoniesPosition,ProblemParams.CostFuncExtraParams);

				// Empire Possession
				empiresList(ii) = possesEmpire(empiresList(ii));

				// Computation of Total Cost for Empires
				empiresList(i).TotalCost = empiresList(i).ImperialistCost + zeta * mean(empiresList(i).ColoniesCost);
			}

			// Uniting Similiar Empires
			empiresList = uniteSimilarEmpires(empiresList);

			// Imperialistic Competition
			empiresList = imperialisticCompetition(empiresList);

			if (empiresList.size() == 1 && stopIfJustOneEmpire)
			{
				break;
			}

		} // End of Algorithm

		
		
		BestCost = MinimumCost(end);
		BestIndex = find(ImerialistCosts == min(ImerialistCosts)); BestIndex = BestIndex(1);
		BestSolution = Empires(BestIndex).ImperialistPosition;

		y=BestSolution;

	}



	private double[] getCountriesCosts(ArrayList<double[]> initialCountries) 
	{
		double[] costs = new double[initialCountries.size()];
		
		for(int i=0; i<initialCountries.size(); i++)
		{
			costs[i] = getCountryCost(initialCountries.get(i));
		}
		return costs;
	}



	private double getCountryCost(double[] ds) 
	{
		// TODO Auto-generated method stub
		return 0;
	}



	private ArrayList<double[]> generateNewCountries(int dimension, double[] minBounds2, double[] maxBounds2, Random r2) 
	{
		ArrayList<double[]> initialCountries = new ArrayList<double[]>();  
		
		for(int i=0; i<initialCountries.size(); i++)
		{
			double[] newRandomCountry = new double[dimension];
			
			for(int j=0; j<dimension; j++)
			{
				newRandomCountry[j] = (maxBounds[j] - minBounds[j]) * r2.nextDouble() + minBounds[j];
			}
			initialCountries.add(newRandomCountry);
		}
		return initialCountries;
	}



	private ArrayList<double[]> createInitialEmpires(ArrayList<double[]> initialCountries, double[] initialCosts)
	{

		// Extract the best countries and their costs to create empires
		ArrayList<double[]> allImperialistsPosition = (ArrayList<double[]>) initialCountries.subList(0, numOfInitialImperialists);
		double[] allImperialistsCost = new double[numOfInitialImperialists];
		System.arraycopy(initialCosts, 0, allImperialistsCost, 0, numOfInitialImperialists);

		// Extract the rest to create colonies
		ArrayList<double[]> allColoniesPosition = (ArrayList<double[]>) initialCountries.subList(numOfInitialImperialists,initialCountries.size());
		double[] allColoniesCost = new double[initialCosts.length-numOfInitialImperialists]; 
		System.arraycopy(initialCosts, numOfInitialImperialists, allColoniesCost, 0, initialCosts.length-numOfInitialImperialists);	

		// Compute the power of imperialists
		double[] allImperialistsPower;
		if(max(allImperialistsCost)>0)
		{
			for(int i=0; i<allImperialistsCost.length; i++)
			{
				allImperialistsPower[i] = 1.3 * max(allImperialistsCost) - allImperialistsCost[i];
			}
		}
		else
		{
			for(int i=0; i<allImperialistsCost.length; i++)
			{
				allImperialistsPower[i] = 0.7 * max(allImperialistsCost) - allImperialistsCost[i];
			}
		}

		// 
		double[] allImperialistNumOfColonies = round(allImperialistsPower/sum(allImperialistsPower) * numOfAllColonies);
		double[] allImperialistNumOfColonies(allImperialistNumOfColonies.length-1) = numOfAllColonies - sum(allImperialistNumOfColonies(1:end-1));
		randomIndex = randperm(numOfAllColonies);

		
		empires(numOfInitialImperialists).ImperialistPosition = 0;

		// Create imperialist empires
		for(ii=1; ii<numOfInitialImperialists; ii++)
		{
			Empires(ii).ImperialistPosition = AllImperialistsPosition(ii,:);
			Empires(ii).ImperialistCost = AllImperialistsCost(ii,:);
			R = RandomIndex(1:AllImperialistNumOfColonies(ii)); 
			RandomIndex(AllImperialistNumOfColonies(ii)+1:end);
			Empires(ii).ColoniesPosition = AllColoniesPosition(R,:);
			Empires(ii).ColoniesCost = AllColoniesCost(R,:);
			Empires(ii).TotalCost = Empires(ii).ImperialistCost + zeta * mean(Empires(ii).ColoniesCost);
		}

		for(ii=1; ii<=numel(Empires); ii++)
		{
			if(numel(Empires(ii).ColoniesPosition) == 0)
			{
				Empires(ii).ColoniesPosition = GenerateNewCountry(1,ProblemParams);
				Empires.ColoniesCost = feval(ProblemParams.FunctionName,Empires.ColoniesPosition);
			}
		}

		return empires;
	}

	
	
	private class Empire
	{
		
		//TODO: tableau/liste d'empires
		int imperialistPosition;
		double imperialistCost;
		int[] coloniesPosition;
		double[] coloniesCost;
		double totalCost;
		
		public int getImperialistPosition() {
			return imperialistPosition;
		}
		public void setImperialistPosition(int imperialistPosition) {
			this.imperialistPosition = imperialistPosition;
		}
		public double getImperialistCost() {
			return imperialistCost;
		}
		public void setImperialistCost(double imperialistCost) {
			this.imperialistCost = imperialistCost;
		}
		public int[] getColoniesPosition() {
			return coloniesPosition;
		}
		public void setColoniesPosition(int[] coloniesPosition) {
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

	}
	

	 private double max(double[] values) {
         double max = Double.MIN_VALUE;
         for(double value : values) 
         {
                 if(value > max)
                         max = value;
         }
         return max;
	 }
	
	

	private Empire assimilateColonies(Empire theEmpire)
	{

		numOfColonies = size(TheEmpire.ColoniesPosition,1);

		Vector = repmat(TheEmpire.ImperialistPosition,numOfColonies,1)-theEmpire.ColoniesPosition;

		TheEmpire.ColoniesPosition = TheEmpire.ColoniesPosition + 2 * assimilationCoefficient * rand(size(Vector)) .* Vector;

		MinVarMatrix = repmat(varMin,numOfColonies,1);
		MaxVarMatrix = repmat(varMax,numOfColonies,1);

		TheEmpire.ColoniesPosition=max(TheEmpire.ColoniesPosition,MinVarMatrix);
		TheEmpire.ColoniesPosition=min(TheEmpire.ColoniesPosition,MaxVarMatrix);

		return TheEmpire;
	}



	private TheEmpire revolveColonies(TheEmpire)
	{
		numOfRevolvingColonies = round(revolutionRate * numel(TheEmpire.ColoniesCost));
		revolvedPosition = GenerateNewCountry(numOfRevolvingColonies);
		R = randperm(numel(TheEmpire.ColoniesCost));
		R = R(1:NumOfRevolvingColonies);
		TheEmpire.ColoniesPosition(R,:) = RevolvedPosition;
		return TheEmpire;
	}



	private TheEmpire possesEmpire(TheEmpire)
	{
		ColoniesCost = TheEmpire.ColoniesCost;

		[MinColoniesCost BestColonyInd]=min(ColoniesCost);

		if(MinColoniesCost < TheEmpire.ImperialistCost)
		{
			OldImperialistPosition = TheEmpire.ImperialistPosition;
			OldImperialistCost = TheEmpire.ImperialistCost;

			TheEmpire.ImperialistPosition = TheEmpire.ColoniesPosition(BestColonyInd,:);
			TheEmpire.ImperialistCost = TheEmpire.ColoniesCost(BestColonyInd);

			TheEmpire.ColoniesPosition(BestColonyInd,:) = OldImperialistPosition;
			TheEmpire.ColoniesCost(BestColonyInd) = OldImperialistCost;
		}
		return TheEmpire;
	}



	private Empires uniteSimilarEmpires(Empires)
	{

		TheresholdDistance = AlgorithmParams.UnitingThreshold * norm(ProblemParams.SearchSpaceSize);
		NumOfEmpires = numel(Empires);

		for(ii=1; ii<=(NumOfEmpires-1); ii++)
		{
			for(jj=ii+1; jj<=NumOfEmpires; jj++)
			{
				DistanceVector = Empires(ii).ImperialistPosition - Empires(jj).ImperialistPosition;
				Distance = norm(DistanceVector);
				if(Distance<=TheresholdDistance)
				{
					if(Empires(ii).ImperialistCost < Empires(jj).ImperialistCost)
					{
						BetterEmpireInd=ii;
						WorseEmpireInd=jj;
					}
					else
					{
						BetterEmpireInd=jj;
						WorseEmpireInd=ii;
					}

					Empires(BetterEmpireInd).ColoniesPosition = [Empires(BetterEmpireInd).ColoniesPosition
					                                             Empires(WorseEmpireInd).ImperialistPosition
					                                             Empires(WorseEmpireInd).ColoniesPosition];

					Empires(BetterEmpireInd).ColoniesCost = [Empires(BetterEmpireInd).ColoniesCost
					                                         Empires(WorseEmpireInd).ImperialistCost
					                                         Empires(WorseEmpireInd).ColoniesCost];

					// Update TotalCost for new United Empire                                     
					Empires(BetterEmpireInd).TotalCost = Empires(BetterEmpireInd).ImperialistCost + AlgorithmParams.Zeta * mean(Empires(BetterEmpireInd).ColoniesCost);

					Empires = Empires([1:WorseEmpireInd-1 WorseEmpireInd+1:end]);
					return;
				}

			}
		}
		return Empires;   
	}



	private Empires imperialisticCompetition(Empires)
	{
		if(rand > .11)
		{
			return;
		}
		if(numel(Empires)<=1)
		{
			return;
		}

		TotalCosts = [Empires.TotalCost];
		[MaxTotalCost WeakestEmpireInd] = max(TotalCosts);
		TotalPowers = MaxTotalCost - TotalCosts;
		PossessionProbability = TotalPowers / sum(TotalPowers);

		SelectedEmpireInd = SelectAnEmpire(PossessionProbability);

		nn = numel(Empires(WeakestEmpireInd).ColoniesCost);
		jj = myrandint(nn,1,1);

		Empires(SelectedEmpireInd).ColoniesPosition = [Empires(SelectedEmpireInd).ColoniesPosition
		                                               Empires(WeakestEmpireInd).ColoniesPosition(jj,:)];

		Empires(SelectedEmpireInd).ColoniesCost = [Empires(SelectedEmpireInd).ColoniesCost
		                                           Empires(WeakestEmpireInd).ColoniesCost(jj)];

		Empires(WeakestEmpireInd).ColoniesPosition = Empires(WeakestEmpireInd).ColoniesPosition([1:jj-1 jj+1:end],:);
		Empires(WeakestEmpireInd).ColoniesCost = Empires(WeakestEmpireInd).ColoniesCost([1:jj-1 jj+1:end],:);

		// Collapse of the the weakest colony-less Empire
		nn = numel(Empires(WeakestEmpireInd).ColoniesCost);
		if(nn<=1)
		{
			Empires(SelectedEmpireInd).ColoniesPosition = [Empires(SelectedEmpireInd).ColoniesPosition
			                                               Empires(WeakestEmpireInd).ImperialistPosition];

			Empires(SelectedEmpireInd).ColoniesCost = [Empires(SelectedEmpireInd).ColoniesCost
			                                           Empires(WeakestEmpireInd).ImperialistCost];

			Empires=Empires([1:WeakestEmpireInd-1 WeakestEmpireInd+1:end]);
		}

		return Empires;

	}



	private index selectAnEmpire(Probability)
	{
		R = rand(size(Probability));
		D = Probability - R;
		[MaxD Index] = max(D);
		return index;
	}

	

	private y myrandint(MaxInt,m,n)
	{
		// This functions creates random numbers between [1 MaxInt] (1 itself and MaxInt itself)
		if(nargin == 1)
		{
			y = ceil(rand*MaxInt);
		}
		else
		{
			if(nargin == 3)
			{
				y = ceil(rand(m,n)*MaxInt);
			}
			else
			{
				//warning('Incorrect Number of Inputs');
			}
		}

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
