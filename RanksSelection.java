
public class RanksSelection implements IndividualSelection{

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
	
	public Individual rouletteWheel(Individual[] pop, int participants)
	{
		double sumRanks = Utils.sumRanks(pop, participants);
		
		double rand = Utils.randomNumber(sumRanks);
		
		for(int i = 0; i < pop.length; i++)
		{
			rand -= pop[i].rank;
			
			if(rand < 0)
			{
				return pop[i];
			}
		}
		
		System.out.println("error----------------------------------------------------------------------------------------------------------------------------------");
		
		return pop[player72.rnd_.nextInt(pop.length)];
	}
	
}
