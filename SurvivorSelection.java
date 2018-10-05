public class SurvivorSelection 
{
	
	private int pop_size = player72.pop_size;
	/*
	 * Chooses the individuals to go to the next generation probabilistically, proportionally to their fitness. 
	 */
	public void stochastic(Individual[] pop, Individual[] next_gen, int elites, double randomly_chosen_individual, IndividualSelection individual_selection)
	{	
		//Individual[] temp = new Individual[pop_size];
		int chosen_individual = 0;
	
		for(int i = 0; i < pop_size; i++)
		{
			int participants = next_gen.length-i;
			boolean randomly_chosen_indiv = (player72.rnd_.nextDouble() > randomly_chosen_individual);
			
			// Elitism.
			if(i < elites) 
			{
				chosen_individual = i;
			}
			if(randomly_chosen_indiv)
			{
				chosen_individual = player72.rnd_.nextInt(participants);
			}
			else 
			{
				chosen_individual = individual_selection.rouletteWheel(next_gen, participants).index;
			}
			
			pop[i] = new Individual(next_gen[chosen_individual]);
			pop[i].index = i;
			
			// Move the picked individual to the end so that he cannot be picked again.
			next_gen[i] = next_gen[participants - 1];
		}
		//copy_populations(pop, temp);
	}
	
	public void copy_populations(Individual[] pop, Individual[] temp)
	{
		for(int i = 0; i < pop_size; i++)
		{
			pop[i] = new Individual(temp[i]);
		}
	}
}
