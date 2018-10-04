import java.util.Random;

public class Utils {
	
	/*
	 * Sums the values of fitnesses of a pop.
	 */
	public static double sumFitness(Individual[] pop)
	{
		double sumFitness = 0;
		
		for(int i = 0; i < pop.length; i++)
		{
			sumFitness += pop[i].fitness;
		}
		
		return sumFitness;
	}
	
	/*
	 * Sums the values of ranks of a pop.
	 */
	public static int sumRanks(Individual[] pop)
	{
		int sumRanks = 0;
		
		for(int i = 0; i < pop.length; i++)
		{
			sumRanks += pop[i].rank;
		}
		
		return sumRanks;
	}
	
	/*
	 * A random number between 0 and sumFitness.
	 */
	public static double randomNumber(double sumFitness)
	{	
		return sumFitness * player72.rnd_.nextDouble();	
	}
	
	/*
	 * Produces a random value between upper and lower if both possitive.
	 */
	public static double double_in_range(double upper, double lower)
	{
    	double randomValue = (upper - lower) * player72.rnd_.nextDouble() + lower;
    	return randomValue;
	}
	
	/*
	 * Returns -1 or 1.
	 */
	public static double getSign()
	{
		double chanse = player72.rnd_.nextDouble();
		double sign = 1;
		
		if(chanse > 0.5)
		{
			sign = 1;
		}
		else
		{
			sign = -1;
		}
			
		return sign;
	}
	
	/*
	 * Returns the rank of an individual given the participants and the number of different ranks.
	 */
	public static int getRank(int index, int number_participants)
	{
		double individual_rank = (number_participants - index) / player72.ranks;
		int rank_as_int = (int) Math.ceil(individual_rank);
		
		// 1 is the minimum rank so everybody got chances in the roulette. 
		int rank = Math.max(1, rank_as_int);
		return rank;
	}
}
