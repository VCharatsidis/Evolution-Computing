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
	public static int dimensions = 10;
	
	public static int rank_populations;
	public static int ranks;
	public static int pop_size;
	public static int offsprings;
	
	public static int evals_left;
	public static CrossOver crossOver;

	private boolean shock = false;
	public static Random rnd_;
	public Printer printer;
	
	ContestEvaluation evaluation_;
    private int evaluations_limit_;
	
	public player72()
	{
		printer = new Printer();
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
        	rank_populations = 2200;
        	ranks = 1;
        	pop_size = ranks * rank_populations;
        	offsprings = (int)(pop_size * 5);
        	crossOver = new ArithmeticCrossOver();
        }else{
            // Do sth else
        }
    }
    
	public void run()
	{
		
		Katsuura();
	}
	
	private void Katsuura()
	{
		System.out.println("Katsuura");
		System.out.println("eval limit "+evaluations_limit_);
		System.out.println("pop_size "+pop_size);
		System.out.println("offsprings "+offsprings);
		//Printer printer = new Printer();
		
        // init population
        Individual[] pop = new Individual[pop_size];
        initPop(pop);
        PriorityQueue<Individual> sorted_pop = new PriorityQueue<Individual>(new IndividualComparator());
        
        double sumFitness = 0;
        double childrenSumFitness = 0;
       
        evaluate_pop(pop, sorted_pop, sumFitness);
        pop = pq_to_array(sorted_pop, pop_size);
        fill_ranks_on_sorted_pop(pop);
  
        FitnessSelection fitnessSelection = new FitnessSelection();
        RanksSelection rankSelection = new RanksSelection();
        SurvivorSelection survivorSelection = new SurvivorSelection();
       
        while(evals_left > 0)
        {
        	double last_avg_children_ftiness = childrenSumFitness/offsprings;
        	//print(pop);
        	Individual[] next_gen = new Individual[offsprings];
        
        	//print(pop);
        	if(evals_left == 1)
    		{
        		printer.printFitnesses(pop, 1, rank_populations);
        		printer.printRanks(pop, 1);
        		printer.printInfo(sumFitness, evals_left, pop[0].fitness, pop_size);
    		}
        	
        	sumFitness = Utils.sumFitness(pop, pop.length);
        	
        	//REPRODUCTION
        	reproduction(pop, next_gen, sumFitness, fitnessSelection, rankSelection, crossOver);	

        	//MUTATE
        	int individuals_to_mutate = offsprings/20;
        	int dims_to_mutate = 1;
        	MutateChildren(next_gen, individuals_to_mutate, dims_to_mutate);
        	
        	sorted_pop.clear();
        	
        	for(int indiv = 0; indiv < pop_size; indiv++)
        	{
                sorted_pop.add(pop[indiv]);
            }
            
        	//CHILDREN EVALUATION
        	childrenSumFitness = 0;
            evaluate_pop(next_gen, sorted_pop, childrenSumFitness);
    		
            //SURVIVOR SELECTION
            pop = pq_to_array(sorted_pop, pop.length);
            fill_ranks_on_sorted_pop(pop);
        	
            //DIVERSITY CHECK
    		diversity_check(pop, sumFitness, sorted_pop, last_avg_children_ftiness, childrenSumFitness);
        }  
	}
	
	private void reproduction(Individual[] pop, Individual[] next_gen, double sumFitness, FitnessSelection fitnessSelection, RanksSelection rankSelection, CrossOver crossOver)
	{
		randomReproduction(pop, next_gen, crossOver, offsprings, sumFitness);
	}
	
	private void randomReproduction(Individual[] pop, Individual[] next_gen, CrossOver crossOver, int offsprings, double sumFitness)
	{
		System.out.println(" random reproduction ");
		printer.printInfo(sumFitness, evals_left, pop[0].fitness, pop_size);
		
		for(int i = 0; i < offsprings; i++)
    	{
			Individual child = new Individual();
			
			Individual parent_a = new Individual();
			Individual parent_b = new Individual();
			
			int random_a = rnd_.nextInt(pop.length);
			
			int random_b = rnd_.nextInt(pop.length);
			
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
			
			//randomParentSelection(pop, parent_a, parent_b);	
			
			child = crossOver.cross_over(parent_a, parent_b);
			next_gen[i] = child;
    	}
		
	}
	
	private void rank_partial_roulette_with_No_replacement_reproduction(
			Individual[] pop, RanksSelection rankSelection, Individual[] next_gen, int offsprings)
	{
		System.out.print(" partial no replacement ");
		Individual[] copy_pop = new Individual[pop_size];
		
    	for(int j  = 0; j < pop_size; j++)
    	{
    		copy_pop[j] = new Individual(pop[j]);
    	}
    	
		for(int i = 0; i < offsprings; i++)
    	{
			Individual child = new Individual();
			
			Individual parent_a = new Individual();
			Individual parent_b = new Individual();
			
			rankParentSelection_NoReplacement(copy_pop, i, parent_a, parent_b, rankSelection);
			
			child = crossOver.cross_over(parent_a, parent_b);
			next_gen[i] = child;
    	}
	}

	private void rank_rouletteWheel_reproduction(Individual[] pop, double sumFitness, FitnessSelection fitnessSelection,
			RanksSelection rankSelection, CrossOver arithmeticCrossover, Individual[] next_gen, int offsprings) 
	{
	
		printer.printInfo(sumFitness, evals_left, pop[0].fitness, pop_size);
		System.out.print(" roulette ");
		

		for(int i = 0; i < offsprings; i++)
    	{
			Individual child = new Individual();
			
			Individual parent_a = new Individual();
			Individual parent_b = new Individual();
			
			rankRouletteWheel(parent_a, parent_b, rankSelection,pop);
			
			child = arithmeticCrossover.cross_over(parent_a, parent_b);
			next_gen[i] = child;
    	}
	}
	
	public void rankParentSelection_NoReplacement(Individual[] copy_pop, int i, Individual parent_a, Individual parent_b, RanksSelection rankSelection )
	{
		int interval = 2;
		int participants = copy_pop.length - i;
		parent_a = rankSelection.partial_rouletteWheel(copy_pop, interval, participants);
		int last_individual = participants - 1;
		
		copy_pop[parent_a.index] = copy_pop[last_individual];
		
		parent_b = rankSelection.partial_rouletteWheel(copy_pop, interval, participants -1);	
		last_individual = participants - 2;
	
		copy_pop[parent_b.index] = copy_pop[last_individual];
	}
	
	public void randomParentSelection(Individual[] pop, Individual parent_a, Individual parent_b)
	{
		
		int random_a = rnd_.nextInt(pop.length);
		
		int random_b = rnd_.nextInt(pop.length);
		
//		while(pop[random_b].fitness > 0.5 && (evals_left > (evaluations_limit_/2)))
//		{
//			random_b = rnd_.nextInt(pop_size);
//		}
//		
//		while(random_a == random_b)
//		{
//			random_b = rnd_.nextInt(pop_size);
//		}
		
		parent_a = pop[random_a];
		parent_b = pop[random_b];
	}
	
	public void rankRouletteWheel(Individual parent_a, Individual parent_b, RanksSelection rankSelection, Individual[] pop)
	{
		parent_a = rankSelection.rouletteWheel(pop, pop.length);
		parent_b = rankSelection.rouletteWheel(pop, pop.length);
			
	}
	

	private void diversity_check(Individual[] pop, double sumFitness, PriorityQueue<Individual> sorted_pop,
			double last_avg_children_ftiness, double childrenSumFitness)
	{
		// Defines epsilon given how small is the average fitness
    	double epsilon = define_epsilon(last_avg_children_ftiness);
    	double avg_children_fitness = childrenSumFitness/offsprings;
    	
		sorted_pop.clear();
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
			SurvivorSelection survivorSelection, Individual[] next_gen, PriorityQueue<Individual> sorted_pop)
	{
//		next_gen = pq_to_array(sorted_pop);
//		fill_ranks_on_sorted_pop(next_gen);
//		if(evals_left < evaluations_limit_/100)
//		{
//			int elites = 0;
//			double chanse_to_randomly_choose_individual = 0;
//			survivorSelection.stochastic(pop, next_gen, elites, chanse_to_randomly_choose_individual, fitnessSelection);
//		}
//		else
//		{
			pop = pq_to_array(sorted_pop, pop.length);
			//survivorSelection.copy_populations(pop, next_gen);
//		}
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
	 * Pop must be sorted.
	 */
	private void fill_ranks_on_sorted_pop(Individual[] pop) 
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

	public Individual[] pq_to_array(PriorityQueue<Individual> pq, int size)
	{
		Individual[] fi = new Individual[size];
		int counter = 0;
		
		while(counter < size)
		{	
			fi[counter] = pq.poll();				
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