
public class FitnessSelection implements IndividualSelection {
	
	/*
	 * It returns 1 individual out of interval number of individuals proportionally to their fitness.
	 */
	public Individual partial_rouletteWheel(Individual[] pop, int interval, int participants)
	{
		Individual[] pop_interval = new Individual[interval];

		for(int i = 0; i < interval; i++)
		{
			int choosen = player72.rnd_.nextInt(participants);
			
			pop_interval[i] = pop[choosen];
		}
		
		Individual choosen_indiv = rouletteWheel(pop_interval, interval);
	
		return choosen_indiv;
	}
	
	/*
	 * Returns 1 individual from the list fitnesses with chanse proportionally to their fitness.
	 */
	public Individual rouletteWheel(Individual[] pop, int participants)
	{
		double sumFitness = Utils.sumFitness(pop, participants);
		
		double rand = Utils.randomNumber(sumFitness);
		
		for(int i = 0; i < participants; i++)
		{
			rand -= pop[i].fitness;
			
			if(rand < 0)
			{
				return pop[i];
			}
		}
		
		System.out.println("error----------------------------------------------------------------------------------------------------------------------------------");
		
		return pop[player72.rnd_.nextInt(pop.length)];
	}
	
	/*
	 * used for island models
	 */
//	public int island_wheel(FitnessIndex[] fi, int start, int interval)
//	{
//		double[] fitnesses_interval = new double[interval];
//		// Gather interval consecutive fitnesses.
//		for(int i = 0; i < interval; i++)
//		{
//			fitnesses_interval[i] = fi[i+start].fitness;
//		}
//		
//		int choosen_indiv = rouletteWheel(fitnesses_interval);
//		
//		int parent = fi[choosen_indiv + start].index;
//		return parent;
//	}
}
