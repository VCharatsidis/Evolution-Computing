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
        	offsprings = (int)(pop_size * 4);
        	
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
		int individuals_to_mutate = (int)((double)offsprings/4);
		
		System.out.println("Katsuura");
		System.out.println("eval limit "+evaluations_limit_);
		System.out.println("pop_size "+pop_size);
		System.out.println("offsprings "+offsprings);
		System.out.println("arithmeticCrossOver");
		System.out.println("mutation rate "+ ((double)individuals_to_mutate / (double)offsprings));
		
		//FITNESS SHARING SIGMA
		boolean fitness_sharing = true;
		double sigma = 10;
		
		if(fitness_sharing)
		{
			System.out.println("fs sigma " + sigma);
		}
		
        // init population
        Individual[] pop = new Individual[pop_size];
        initPop(pop);
        PriorityQueue<Individual> sorted_pop = new PriorityQueue<Individual>(new IndividualComparator());
        
        double sumFitness = 0;
        double childrenSumFitness = 0;
        
        evaluate_pop(pop, sorted_pop, sumFitness);
        pq_to_array(sorted_pop, pop);
        
        fill_ranks_on_sorted_pop(pop);
  
        FitnessSelection fitnessSelection = new FitnessSelection();
        RanksSelection rankSelection = new RanksSelection();
        SurvivorSelection survivorSelection = new SurvivorSelection();
        CrossOver arithmeticCrossOver = new ArithmeticCrossOver();
        CrossOver uniCrossOver = new UniformCrossOver();
        
       
        while(evals_left > 0)
        {
        	double last_avg_children_ftiness = childrenSumFitness/offsprings;
        	//print(pop);
        	Individual[] next_gen = new Individual[offsprings];
            
        	sumFitness = Utils.sumFitness(pop, pop.length);
        	
        	double[] sh = new double[pop_size];
        	
        	//FITNESS SHARING
        	if(fitness_sharing)
        	{
            	sh_values(sh, pop, next_gen, sigma);
            	
//            	for(int i  = 0; i <10000;i++) {
//            		System.out.println(sh[0]);
//            	}
            	
            	for(int indiv = 0; indiv < pop_size; indiv++)
            	{
            		pop[indiv].fitness = pop[indiv].fitness / (1 + sh[indiv]);
                }
        	}
        	
        	//REPRODUCTION
        	double chanse = rnd_.nextDouble();
        	if(chanse >0.9)
        	{
        		System.out.println("Uniform CrossOver");
        		reproduction(pop, next_gen, sumFitness, fitnessSelection, rankSelection, uniCrossOver);	
        	}
        	else
        	{
        		System.out.println("Arithemtic CrossOver");
        		reproduction(pop, next_gen, sumFitness, fitnessSelection, rankSelection, arithmeticCrossOver);
        	}
        		
        	//MUTATE
        	int dims_to_mutate = 1;
        	MutateChildren(next_gen, individuals_to_mutate, dims_to_mutate);
        	
        	sorted_pop.clear();
        	fill_sorted_pop(fitness_sharing, pop, sorted_pop, sh);
    		
        	//CHILDREN EVALUATION
        	childrenSumFitness = 0;
            evaluate_pop(next_gen, sorted_pop, childrenSumFitness);
            
            //SURVIVOR SELECTION
            pq_to_array(sorted_pop, pop);
            //fill_ranks_on_sorted_pop(pop);
        	
            //DIVERSITY CHECK
    		diversity_check(pop, sumFitness, sorted_pop, last_avg_children_ftiness, childrenSumFitness);
        }  
	}

	private void fill_sorted_pop(boolean fitness_sharing, Individual[] pop, PriorityQueue<Individual> sorted_pop,
			double[] sh) {
		for(int indiv = 0; indiv < pop_size; indiv++)
		{
			// restore fitness for survivor selection.
			if(fitness_sharing)
			{
				pop[indiv].fitness = pop[indiv].fitness * (1 + sh[indiv]);
			}
			
			Individual individual = new Individual(pop[indiv]);
			sorted_pop.add(individual);
		}
	}
	
	/*
	 * Fills the distances array which has length pop_size + offsprings.
	 */
	private void sh_values(double[] sh, Individual[] pop, Individual[] next_gen, double sigma)
	{
		for(int i = 0; i < pop.length; i++)
		{
			for(int j = i+1; j  < pop.length; j ++)
			{
				calculate_sh(pop[i], pop[j], i, j, sh, sigma);
			}
			
//			for(int j = 0; j < next_gen.length; j ++)
//			{
//				calculate_sh(pop[i], next_gen[j], i, j + pop.length, sh, sigma);
//			}
		}
		
//		for(int i = 0; i < next_gen.length; i ++)
//		{
//			for(int j = i + 1; j < next_gen.length; j ++)
//			{
//				calculate_sh(next_gen[i], next_gen[j], i + pop.length, j + pop.length, sh, sigma);
//			}
//		}
	}
	
	private void calculate_sh(Individual a, Individual b, int index_a, int index_b, double[] sh, double sigma)
	{
		double total_distance = 0;
		
		for(int dim = 0; dim < dimensions; dim++)
		{
			if(a.genome[dim] > 0 && b.genome[dim] > 0)
			{
				total_distance += Math.abs(a.genome[dim] - b.genome[dim]);
			}
			else if(a.genome[dim] < 0 && b.genome[dim] < 0)
			{
				total_distance += Math.abs(Math.abs(a.genome[dim]) - Math.abs(b.genome[dim]));
			}
			else
			{
				total_distance += Math.abs(a.genome[dim] + b.genome[dim]);
			}	
		}
		
		double sh_value = 1 - total_distance/sigma;
		sh_value = Math.max(sh_value, 0);
		//System.out.println(" td "+total_distance+" index a "+index_a+" index b "+index_b + " sh_value "+sh_value) ;
		
		sh[index_a] += sh_value;
		sh[index_b] += sh_value;
	}
	
	private void reproduction(Individual[] pop, Individual[] next_gen, double sumFitness, FitnessSelection fitnessSelection, RanksSelection rankSelection, CrossOver crossOver)
	{
		//randomReproduction(pop, next_gen, crossOver, offsprings, sumFitness);
		fitness_rouletteWheel_reproduction(pop, next_gen, crossOver, offsprings, sumFitness, fitnessSelection);
	}
	
	private void fitness_rouletteWheel_reproduction(Individual[] pop, Individual[] next_gen, CrossOver crossOver, int offsprings, double sumFitness, FitnessSelection fitnessSelection )
	{
		System.out.println("fitness roulette Wheel ");
		printer.printInfo(sumFitness, evals_left, pop[0].fitness, pop_size);
		
		for(int i = 0; i < offsprings; i++)
    	{
			Individual[] indivs = rouletteWheel_fitness_parentSelection(pop, fitnessSelection);
			
			next_gen[i] = crossOver.cross_over(indivs[0], indivs[1]);
			//System.out.println(next_gen[i].genome[0] +" "+ next_gen[i].genome[9]);
    	}
	}
	
	private Individual[] rouletteWheel_fitness_parentSelection(Individual[] pop, FitnessSelection fitnessSelection)
	{
		Individual parent_a = fitnessSelection.partial_rouletteWheel(pop, 2, pop.length);
		Individual parent_b = fitnessSelection.partial_rouletteWheel(pop, 2, pop.length);
		
		while( parent_a.index == parent_b.index)
		{
			parent_b = fitnessSelection.partial_rouletteWheel(pop, 2, pop.length);
		}
		
		Individual[] idnivs = {parent_a, parent_b};
		return idnivs;
	}
	
	private void randomReproduction(Individual[] pop, Individual[] next_gen, CrossOver crossOver, int offsprings, double sumFitness)
	{
		System.out.println("random reproduction ");
		printer.printInfo(sumFitness, evals_left, pop[0].fitness, pop_size);
		
		for(int i = 0; i < offsprings; i++)
    	{
			int[] indexes = randomParentSelection(pop.length);
			
			next_gen[i] = crossOver.cross_over(pop[indexes[0]], pop[indexes[1]]);
			//System.out.println(next_gen[i].genome[0] +" "+ next_gen[i].genome[9]);
    	}
	}
	
	private void rank_partial_roulette_with_No_replacement_reproduction(
			Individual[] pop, RanksSelection rankSelection, Individual[] next_gen, int offsprings, CrossOver crossOver)
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
	
	public int[] randomParentSelection(int pop_length)
	{
		
		int index_parent_a = rnd_.nextInt(pop_length);
		
		int index_parent_b = rnd_.nextInt(pop_length);
		
		int[] indexes = {index_parent_a, index_parent_b};
		return indexes;
		
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
			 pq_to_array(sorted_pop, pop);
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
            
            pop[i].fitness = fitness ;
            pop[i].index = i;
            
            Individual indiv = new Individual(pop[i]);
            
            sorted_pop.add(indiv);
            
            sumFitness += fitness;
            evals_left--;
            
            if(evals_left ==1 )
    		{
        		printer.printFitnesses(pop, pop.length, rank_populations);
        		//printer.printRanks(pop, 1);
        		printer.printInfo(sumFitness, evals_left, pop[0].fitness, pop_size);
        		printer.printBestIndivValues(pop, 100);
    		}
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

	public void pq_to_array(PriorityQueue<Individual> pq, Individual[] pop)
	{
		int counter = 0;
		
		while(counter < pop.length)
		{	
			pop[counter] = pq.poll();				
			counter++;
		}	
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