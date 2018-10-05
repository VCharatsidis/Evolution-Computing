
public class Printer {
	
	public Printer()
	{
		
	}
	public void printRanks(Individual[] pop, int individuals_to_display) 
	{
		for(int i = 0; i < individuals_to_display; i++)
		{
			System.out.println("indiv "+i+" ranks "+pop[i].rank);
		}
	}

	public void printFitnesses(Individual[] pop, int individuals_to_display, int rank_populations)
	{
		for(int i = 0; i < individuals_to_display; i++)
		{
			System.out.println("indiv "+i+" fitness "+pop[i].fitness);
		}
		
		System.out.println();
		
		System.out.println("pop_size "+ pop.length+" num_indiv_per_rank "+rank_populations);
	}

	public void printPop(Individual[] pop, int dimensions)
	{
		for(int i = 0; i < pop.length; i++)
		{
		    System.out.println();
		    System.out.println("indiv "+i+" dims: ");
		    
			for(int dim = 0; dim < dimensions; dim++) 
			{
				System.out.print(pop[i].genome[dim] + " , ");
			}
			
			System.out.println();
		}
	}
	
	public void printInfo( double sumFitness, int evals_left, double fitness, int pop_size)
	{
		System.out.println();
		System.out.println(" avg fitness : "+(sumFitness/pop_size) +" best indiv : "+fitness+" evals left "+evals_left);
		System.out.println();
	}
	
	public void printBestIndivValues(Individual[] pop) {
		System.out.println();
		System.out.println("fitness of best indiv"+ pop[0].fitness);
		
		for(int dim = 0; dim < 10; dim++)
		{
			System.out.print("genome dim "+dim+" "+pop[0].genome[dim]);
		}
		
		for(int dim = 0; dim < 10; dim++)
		{
			System.out.print("mutation_step dim "+dim+" "+pop[0].mutation_steps[dim]);
		}
		
		System.out.println();
	}
	
}
