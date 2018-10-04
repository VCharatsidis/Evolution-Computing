import java.util.ArrayList;
import java.util.Collections;
import java.util.PriorityQueue;
import java.util.Properties;
import java.util.Random;

import org.vu.contest.ContestEvaluation;
import org.vu.contest.ContestSubmission;

public class player72 implements ContestSubmission
{
	public static double uper_bound = 5;
	public static double lower_bound = -5;
	
	public static int rank_populations = 2200;
	public static int ranks = 1;
	public static int pop_size = ranks * rank_populations;
	public static int offsprings = (int)(pop_size * 5);
	public static int dimensions = 10;
	public static int evals_left;

	private boolean shock = false;
	public static Random rnd_;
	
	ContestEvaluation evaluation_;
    private int evaluations_limit_;
	
	public player72()
	{
		dimensions = 10;
		rnd_ = new Random();
	}
	
	public void setSeed(long seed)
	{
		// Set seed of algortihms random process
		rnd_.setSeed(seed);
	}

	public void setEvaluation(ContestEvaluation evaluation)
	{
		// Set evaluation problem used in the run
		evaluation_ = evaluation;
		
		// Get evaluation properties
		Properties props = evaluation.getProperties();
        // Get evaluation limit
        evaluations_limit_ = Integer.parseInt(props.getProperty("Evaluations"));
        evals_left = evaluations_limit_;
		// Property keys depend on specific evaluation
		// E.g. double param = Double.parseDouble(props.getProperty("property_name"));
        boolean isMultimodal = Boolean.parseBoolean(props.getProperty("Multimodal"));
        boolean hasStructure = Boolean.parseBoolean(props.getProperty("Regular"));
        boolean isSeparable = Boolean.parseBoolean(props.getProperty("Separable"));

		// Do sth with property values, e.g. specify relevant settings of your algorithm
        if(isMultimodal){
            // Do sth
        }else{
            // Do sth else
        }
    }
    
	public void run()
	{
		// Run your algorithm here
		System.out.println("eval limit "+evaluations_limit_);
		
        // init population
        Individual[] pop = new Individual[pop_size];
        initPop(pop);
        
        // calculate fitness
       
        double sumFitness = 0;
        double childrenSumFitness = 0;
       
        //int offsprings = 800;
        System.out.println("offsprings "+offsprings);
        
        PriorityQueue<Individual> sorted_pop = new PriorityQueue<Individual>(new IndividualComparator());
        evaluate_pop(pop, sorted_pop, sumFitness);
        
        pop = pq_to_array(sorted_pop);
    	
        fill_ranks(pop);
  
        IndividualSelection fitnessSelection = new FitnessSelection();
        IndividualSelection rankSelection = new RanksSelection();
        ArithmeticCrossOver arithmeticCrossover = new ArithmeticCrossOver();
        SurvivorSelection survivorSelection = new SurvivorSelection();
       
        while(evals_left > 0)
        {
        	
        	double last_avg_children_ftiness = childrenSumFitness/offsprings;
        	//print(pop);
        	Individual[] next_gen = new Individual[offsprings];
        
        	//print(pop);
        	if(evals_left == 1)
    		{
        		printFitnesses(pop, 1);
        		printRanks(pop, 1);
        		printInfo(sumFitness, evals_left, pop[0].fitness);
    		}
        	
        	sumFitness = Utils.sumFitness(pop, pop.length);
    		
    		Individual[] copy_pop = new Individual[pop_size];
//        	for(int j  = 0; j < pop_size; j++)
//        	{
//        		copy_pop[j] = new Individual(pop[j]);
//        	}
        	
        	for(int i = 0; i < offsprings; i++)
        	{
        		Individual child = produce_child(pop, sumFitness, fitnessSelection, rankSelection, arithmeticCrossover,
						copy_pop, i);	
        		//System.out.println("parent a "+parent_a + " parent b "+parent_b);
            	next_gen[i] = child;
        	}
        	
        	//MUTATE
        	int individuals_to_mutate = offsprings/10;
        	int dims_to_mutate = 1;
        	if(evals_left < evaluations_limit_/30)
        	{
        		individuals_to_mutate = offsprings/4;
        	}
        	MutateChildren(next_gen, individuals_to_mutate, dims_to_mutate);
        	
        	sorted_pop.clear();
        	
        	for(int indiv = 0; indiv < pop_size; indiv++)
        	{
                sorted_pop.add(pop[indiv]);
            }
            
            childrenSumFitness = 0;
            
            evaluate_pop(next_gen, sorted_pop, childrenSumFitness);
    		
        	next_gen = pq_to_array(sorted_pop);
        	fill_ranks(next_gen);
        	
        	survivorSelection(pop, fitnessSelection, survivorSelection, next_gen);
        	
        	// Defines epsilon given how small is the average fitness
        	double epsilon = define_epsilon(last_avg_children_ftiness);
        	double avg_children_fitness = childrenSumFitness/offsprings;
        	
    		diversity_check(pop, sumFitness, sorted_pop, last_avg_children_ftiness, epsilon, avg_children_fitness);
        
        }  
	}

	private Individual produce_child(Individual[] pop, double sumFitness, IndividualSelection fitnessSelection,
			IndividualSelection rankSelection, ArithmeticCrossOver arithmeticCrossover, Individual[] copy_pop, int i) 
	{
		Individual child = new Individual();
		
		Individual parent_a = new Individual();
		Individual parent_b = new Individual();
		
		//boolean partial = (runs <= total_runs || runs % 20 ==0);
		boolean partial = true;
         
		if(partial) {
			//Partial wheel
			int interval = 2;
		
			//if(true)
			// Parent selection with no replacement
			if(false)
			{
				int participants = pop.length - i;
				parent_a = rankSelection.partial_rouletteWheel(copy_pop, interval, participants);
				int last_individual = participants - 1;
				
				copy_pop[parent_a.index] = copy_pop[last_individual];
				
				parent_b = rankSelection.partial_rouletteWheel(copy_pop, interval, participants -1);	
				last_individual = participants - 2;
			
				copy_pop[parent_b.index] = copy_pop[last_individual];
				
			}
			else
			{
				int random_a = rnd_.nextInt(pop_size);
				while(pop[random_a].fitness > 0.5 && (evals_left > (evaluations_limit_/2)))
				{
					random_a = rnd_.nextInt(pop_size);
				}
				
				int random_b = rnd_.nextInt(pop_size);
				while(pop[random_b].fitness > 0.5 && (evals_left > (evaluations_limit_/2)))
				{
					random_b = rnd_.nextInt(pop_size);
				}
				
				while(random_a == random_b)
				{
					random_b = rnd_.nextInt(pop_size);
				}
				
				
				parent_a = pop[random_a];
				parent_b = pop[random_b];
					
			}
			
			if(i ==0) {
				printInfo(sumFitness, evals_left, pop[0].fitness);
				System.out.print(" partial ");
			}	
		}
		
		boolean roulette = !partial;
		   	
		if(roulette)
		{
			parent_a = fitnessSelection.rouletteWheel(pop, pop.length);
			parent_b = fitnessSelection.rouletteWheel(pop, pop.length);
			
			if(i ==0)
			{
				printInfo(sumFitness, evals_left, pop[0].fitness);
				System.out.print(" roulette ");
			}
		}
		
		child = arithmeticCrossover.cross_over(parent_a, parent_b);
		return child;
	}

	private void diversity_check(Individual[] pop, double sumFitness, PriorityQueue<Individual> sorted_pop,
			double last_avg_children_ftiness, double epsilon, double avg_children_fitness)
	{
		if((Math.abs(last_avg_children_ftiness - avg_children_fitness) < epsilon) && !shock)
		{
			if(last_avg_children_ftiness > 0.1)
			{
				initPop(pop);
				evaluate_pop(pop, sorted_pop, sumFitness);
				System.out.println("Pop Init -------------------------------------------------------------------------------------------------");
			}
		}
	}

	//TODO not working properly
	private void survivorSelection(Individual[] pop, IndividualSelection fitnessSelection,
			SurvivorSelection survivorSelection, Individual[] next_gen)
	{
//		if(evals_left < evaluations_limit_/100)
//		{
//			int elites = 0;
//			double chanse_to_randomly_choose_individual = 0;
//			survivorSelection.stochastic(pop, next_gen, elites, chanse_to_randomly_choose_individual, fitnessSelection);
//		}
//		else
		{
			survivorSelection.copy_populations(pop, next_gen);
		}
	}

	private void printBestIndivValues(Individual[] pop) {
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

	private double define_epsilon(double last_avg_children_ftiness)
	{
		double epsilon = 0;
		
		if(last_avg_children_ftiness > 1)
		{
			epsilon = 0.0001;		
		}
		else if((last_avg_children_ftiness * 10) > 1)
		{
			epsilon = 0.00001;
		}
		else if((last_avg_children_ftiness * 100) > 1)
		{
			epsilon = 0.000001;
		}
		else if(last_avg_children_ftiness * 1000 > 1)
		{
			epsilon = 0.0000001;
		}
		else if(last_avg_children_ftiness * 10000 > 1)
		{
			epsilon = 0.00000001;
		}else if(last_avg_children_ftiness * 100000 > 1)
		{
			epsilon = 0.000000001;
		}
		
		return epsilon;
	}
	
	/*
	 * This method fills the individuals of the population with their corresponding rank.
	 */
	private void fill_ranks(Individual[] pop) 
	{
		for(int i = 0; i < pop.length; i++)
		{
			pop[i].rank = Utils.getRank(i, pop.length);
		}
	}

	private void evaluate_pop(Individual[] pop, PriorityQueue<Individual> sorted_pop, double sumFitness)
	{
		for(int i = 0; i < pop.length; i++)
        {
            Double fitness = (double) evaluation_.evaluate(pop[i].genome);
            
            pop[i].fitness = fitness;
            pop[i].index = i;
            
            sorted_pop.add(pop[i]);
            
            sumFitness += fitness;
            evals_left--;
        }
	}

	private void printRanks(Individual[] pop, int individuals_to_display) 
	{
		for(int i = 0; i < individuals_to_display; i++)
		{
			System.out.println("indiv "+i+" ranks "+pop[i].rank);
		}
	}

	private void printFitnesses(Individual[] pop, int individuals_to_display)
	{
		for(int i = 0; i < individuals_to_display; i++)
		{
			System.out.println("indiv "+i+" fitness "+pop[i].fitness);
		}
		
		System.out.println();
		
		System.out.println("pop_size "+pop_size+" num_indiv_per_rank "+rank_populations);
	}

	private void MutateChildren(Individual[] next_gen, int number_individuals_to_mutate, int number_dims_to_mutate) 
	{
		ArrayList<Integer> list_of_dimensions = new ArrayList<Integer>();
		
		for(int dim = 0; dim < dimensions; dim++)
		{
			list_of_dimensions.add(dim);
		}
		
		Collections.shuffle(list_of_dimensions);
		
		// Randomly pick number_individuals_to_mutate individuals and store their index in the to_mutate array.
		int[] to_mutate = new int[number_individuals_to_mutate];
		int max = offsprings;
		int min = 0;
		for(int i = 0; i < number_individuals_to_mutate; i++)
		{
			to_mutate[i] = rnd_.nextInt(max-min) +min;
		}
		
		for(int i = 0; i < number_individuals_to_mutate; i++)
		{
			int indiv = to_mutate[i];
			
			for(int dim = 0; dim < number_dims_to_mutate; dim++)
			{
				int dim_to_mutate = list_of_dimensions.get(dim);
				mutate_individual(next_gen[indiv], dim_to_mutate);
			}
		}
	}
	
	/*
	 * Fills the individual array with random values from -5 to 5 for genome and 0.05 to 0.25 for mutation_steps.
	 */
	private void initPop(Individual[] pop)
	{
		double upper_mutation_step_bound = 0.25;
		double lower_mutation_step_bound = 0.05;
		
		for(int individual = 0; individual < pop_size; individual++)
        {
			pop[individual] = new Individual(uper_bound, lower_bound, upper_mutation_step_bound, lower_mutation_step_bound);
        }
	}
	
	private void printInfo( double sumFitness, int evals_left, double fitness)
	{
		System.out.println();
		System.out.println(" avg fitness : "+(sumFitness/pop_size) +" best indiv : "+fitness+" evals left "+evals_left);
		System.out.println();
	}

	
	private void printPop(Individual[] pop)
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
	
	public Individual[] pq_to_array(PriorityQueue<Individual> pq)
	{
		int size = pq.size();
		Individual[] fi = new Individual[size];
		int counter = 0;
		
		while(counter < size)
		{
			
			fi[counter] = pq.poll();
//			if(counter < 100) {
//				//System.out.println("index "+ fi[counter].index +" fitness "+fi[counter].fitness);
//			}
				
			counter++;
		}
		
		return fi;
	}
	
	public void mutate_individual(Individual indiv, int dim_to_mutate)
	{
		double chanse = rnd_.nextDouble();
		double sign = Utils.getSign();
		double change = sign * indiv.mutation_steps[dim_to_mutate];
		
		if(chanse > 0.5)
		{
			if(change < 0)
			{
				indiv.genome[dim_to_mutate] = Math.max(-5, indiv.genome[dim_to_mutate] + change);
			}
			else
			{
				indiv.genome[dim_to_mutate] = Math.min(5, indiv.genome[dim_to_mutate] + change);
			}
		}
	}
		
	public static void main(String args[])
	{
		System.out.println("start!");
	}
	
	
}