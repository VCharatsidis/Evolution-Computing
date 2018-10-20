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
	public boolean katsuura;
	public boolean schaffers;
	public boolean bentCigar;
	
	public static int evals_left;
	public int infinityProtection = 1;

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
        
        if(isMultimodal && hasStructure){
            // Do sth
        	schaffers = true;
        	rank_populations = 100;
        	ranks = 1;
        	pop_size = ranks * rank_populations;
        	offsprings = (int)((double)pop_size * 2);
        	
        }else if(isMultimodal){
            // Do sth else
        	katsuura = true;
        	
        	rank_populations = 100;
        	ranks = 1;
        	pop_size = ranks * rank_populations;
        	offsprings = (int)((double)pop_size * 4);
        }
        if(!isMultimodal)
        {
        	katsuura =false;
        	schaffers =false;
        	bentCigar = true;
        }
    }
    
	public void run()
	{
//		if(schaffers)
//		{
//			BentCigar();
//			
//		}
//		else if(katsuura)
//		{
//			Katsuura();
//		}else {
//			BentCigar();
//		}
		long start = System.currentTimeMillis();    
		
		
		double[] center = {0,0,0,0,0 ,0,0,0,0,0};
		double[] best_indiv = {0,0,0,0,0 ,0,0,0,0,0};
		double prec = 0.01;
		int epoch = 0;
		double best_score = 0;
		int runs = 150;
		double gamma = 0.80;
		double max_prec = 5*0.0000000000001;
		
		if(schaffers)
		{
			runs = 200;
			gamma = 0.95;
			max_prec = 5*0.0000000000001;
		}
		else if(katsuura)
		{

			runs = 580;
			gamma = 0.98;
			max_prec = 5*0.00000001;
			
		}
		else if(bentCigar)
		{
			runs = 90;
			gamma = 0.9;
			max_prec = 5*0.00000000000001;
		}
		
		
		System.out.println("max_prec "+max_prec+" runs "+ runs+" gamma "+gamma);
		while(true)
		{
			
			double current_score = 0;
//			runs = runs - 100;
//			runs = Math.max(200, runs);
			
			current_score = public_wellfare(runs, prec, center, 1000, true, best_score);
		
			if(current_score <-500)
			{
				break;
			}

			prec = Math.max(max_prec, prec);	
//			epoch++;
			prec *= gamma;
//			//beta *= gamma;
//			
//			if(katsuura)
//			{  
//				if(epoch>20)
//				{
//					runs = 500;
//					
//					
//				}
//				
//			}
//			
//			if(schaffers)
//			{
//				prec = Math.max(5*0.00001, prec);
//			}
			
			
//			if(epoch >10 )
//			{
//				runs = 3000;
//			}
			
			if(current_score > best_score)
			{
				best_score = current_score;
				
				for(int dim = 0; dim < dimensions; dim++)
				{
					best_indiv[dim] = center[dim];
				}
			}
			else if(bentCigar)
			{
				for(int dim = 0; dim < dimensions; dim++)
				{
					center[dim] = best_indiv[dim];
				}
			}
			
			
			

//		    System.out.println("precision "+prec+" gamma "+gamma +" runs "+runs + " best_score "+best_score);
//			System.out.println(best_score);
//			System.out.println("------------------------------------------------------------------");
		    long elapsedTime = System.currentTimeMillis() - start;
		    if(elapsedTime > 12000)
		    {
		    	break;
		    }
		    epoch++;
		}
		
		//Katsuura();
	}
	
	
	
//	private void Katsuura()
//	{
//		int individuals_to_mutate = (int)((double)offsprings/15);
//	
//        ArithmeticCrossOver arithmeticCrossOver = new ArithmeticCrossOver();
//        CrossOver uniCrossOver = new UniformCrossOver();
//        
//		System.out.println("Katsuura");
//		System.out.println("eval limit "+evaluations_limit_);
//		System.out.println("pop_size "+pop_size);
//		System.out.println("offsprings "+offsprings);
//		System.out.println("arithmeticCrossOver distance "+ arithmeticCrossOver.distance + " genome increment "+arithmeticCrossOver.genome_increment);
//		System.out.println("mutation rate "+ ((double)individuals_to_mutate / (double)offsprings));
//		
//		
//		//FITNESS SHARING SIGMA
//		boolean fitness_sharing = false;
//		double sigma = 100;
//		
//		if(sigma == 100)
//		{
//			infinityProtection = 0;
//		}
//		
//		if(fitness_sharing)
//		{
//			System.out.println("fs sigma " + sigma);
//		}
//		
//        // init population
//        Individual[] pop = new Individual[pop_size];
//        initPop(pop);
//        PriorityQueue<Individual> sorted_pop = new PriorityQueue<Individual>(new IndividualComparator());
//        
//        double sumFitness = 0;
//        double childrenSumFitness = 0;
//        double last_pop_fitness = 0;
//        double previous_top_idniv_fitness = 0;
//        
//        double last_best_fitness = 0;
//        int inti_counter = 0;
//        int no_progress = 0;
//        double score = 0;
//        
//        evaluate_pop(pop, sorted_pop, sumFitness);
//        pq_to_array(sorted_pop, pop);
//        
//        fill_ranks_on_sorted_pop(pop);
//  
//        while(evals_left > 0)
//        {
//        	
//       
//        	Individual[] next_gen = new Individual[offsprings];
//            
//        	sumFitness = Utils.sumFitness(pop, pop.length);
//        	
//        	double[] sh = new double[pop_size];
//        	
//        	//FITNESS SHARING
//        	if(fitness_sharing)
//        	{
//            	adjust_fitnesses(sigma, pop, next_gen, sh);
//        	}
//        	
//        	//REPRODUCTION
//        	double chanse = rnd_.nextDouble();
//        	if(chanse > 0.9)
//        	{
//        		//System.out.println("Uniform CrossOver");
//        		randomReproduction(pop, next_gen, uniCrossOver, offsprings, sumFitness);
//        		//(pop, next_gen, sumFitness, fitnessSelection, rankSelection, uniCrossOver);	
//        	}
//        	else
//        	{
//        		//System.out.println("Arithemtic CrossOver");
//        		randomReproduction(pop, next_gen, arithmeticCrossOver, offsprings, sumFitness);
//        		//reproduction(pop, next_gen, sumFitness, fitnessSelection, rankSelection, arithmeticCrossOver);
//        	}
//        		
//        	//MUTATE
//        	int dims_to_mutate = 1;
//        	MutateChildren(next_gen, individuals_to_mutate, dims_to_mutate);
//        	
//        	sorted_pop.clear();
//        	fill_sorted_pop(fitness_sharing, pop, sorted_pop, sh);
//    		
//        	//CHILDREN EVALUATION
//        	childrenSumFitness = 0;
//            evaluate_pop(next_gen, sorted_pop, childrenSumFitness);
//           
//            if(score < pop[0].fitness)
//            {
//            	score = pop[0].fitness;
//            }
//           // System.out.println("score "+score);
//            //System.out.println("pop[0].fitness "+pop[0].fitness);
//            sumFitness = Utils.sumFitness(pop, pop.length);
//           // System.out.println("SumFitness "+sumFitness+" last_pop_fitness "+last_pop_fitness);
//           // System.out.println("difference "+ Math.abs(sumFitness - last_pop_fitness));
//           // printer.printInfo(sumFitness, evals_left, pop[0].fitness, pop.length);
//            
//           
//            double epsilon = define_epsilon(sumFitness);
////            boolean noDifference = ((Math.abs(sumFitness - last_pop_fitness) < epsilon));
////            boolean pop_degenerated = (sumFitness > (pop[0].fitness * (pop_size/2 + pop_size/3)));
//            
////        	if(noDifference || pop_degenerated)
////            {
//             double difference = 80;
//             if(sumFitness < 1000)
//             {
//            	 difference = 8;
//             }
//             if(sumFitness >10000)
//             {
//            	 difference = 800;
//             }
//           
//      
//             
//             double best_fitness = pop[0].fitness;
//             
////             System.out.println("total distance "+total_distance);
////             System.out.println("condition 1 "+(Math.abs(sumFitness - last_pop_fitness) < difference));
////             System.out.println("conditions 2 "+ (sumFitness > (best_fitness * (pop_size/10))));
////             System.out.println("condition 3 "+ (Math.abs(best_fitness - last_best_fitness) < 0.001));
////             System.out.println("init counter "+inti_counter);
////             System.out.println();
//             
//    		 if((Math.abs(sumFitness - last_pop_fitness) < difference) && sumFitness > (best_fitness * (pop_size/10)) &&
//    				 Math.abs(best_fitness - last_best_fitness) < 0.001)
//             {
//            	//System.out.println("childrenSumFitness " + (sumFitness/offsprings));
//             	//System.out.println("epsilon "+ (1 - (sumFitness/offsprings)/pop[0].fitness));
//    			inti_counter++;
//            	shock = true;
//            	Individual center_indiv = new Individual(pop[0]);
//            	double radious = 0.05;
//            	
//            	epsilon = define_epsilon(center_indiv.fitness);
//            	
//            	if(Math.abs(center_indiv.fitness - previous_top_idniv_fitness) < epsilon)
//            	{
//            		//center_indiv = new Individual(pop[10]);
//            		
//            		radious = 0.05;
//            		
//            		center_indiv = new Individual(pop[no_progress]);
//            		no_progress++;
//            	}
//            	else
//            	{
//            		no_progress = 0;
//            	}
//            	
//            	previous_top_idniv_fitness = center_indiv.fitness;
//            	//printer.printFitnesses(pop, pop.length, rank_populations);
//            	initPop_around_individual(pop, center_indiv, radious);
//            	
//            	sorted_pop.clear();
//            	evaluate_pop(pop, sorted_pop, sumFitness);
//            	
//    			///System.out.println("Pop Init -------------------------------------------------------------------------------------");
//            }
//    		last_best_fitness = best_fitness;
//    		
//           // System.out.println("kill first "+kill_first);
//           // System.out.println();
//            last_pop_fitness = sumFitness;
//            //SURVIVOR SELECTION
//            pq_to_array(sorted_pop, pop);
//            //fill_ranks_on_sorted_pop(pop);
//           
//            //DIVERSITY CHECK
//    		//diversity_check(pop, sumFitness, sorted_pop, last_avg_children_ftiness, childrenSumFitness);
//        }  
//	}
//	
//	private void BentCigar()
//	{
//		int individuals_to_mutate = (int)((double)offsprings/15);
//		
//		FitnessSelection fitnessSelection = new FitnessSelection();
//
//        ArithmeticCrossOver arithmeticCrossOver = new ArithmeticCrossOver();
//        CrossOver uniCrossOver = new UniformCrossOver();
//        
//		System.out.println("BentCigar");
//		System.out.println("eval limit "+evaluations_limit_);
//		System.out.println("pop_size "+pop_size);
//		System.out.println("offsprings "+offsprings);
//		System.out.println("arithmeticCrossOver distance "+ arithmeticCrossOver.distance + " genome increment "+arithmeticCrossOver.genome_increment);
//		System.out.println("mutation rate "+ ((double)individuals_to_mutate / (double)offsprings));
//		
//		
//		//FITNESS SHARING SIGMA
//		boolean fitness_sharing = false;
//		double sigma = 100;
//		
//		if(sigma == 100)
//		{
//			infinityProtection = 0;
//		}
//		
//		if(fitness_sharing)
//		{
//			System.out.println("fs sigma " + sigma);
//		}
//		
//        // init population
//        Individual[] pop = new Individual[pop_size];
//        initPop(pop);
//        PriorityQueue<Individual> sorted_pop = new PriorityQueue<Individual>(new IndividualComparator());
//        
//        double sumFitness = 0;
//        double childrenSumFitness = 0;
//        double last_pop_fitness = 0;
//        double previous_top_idniv_fitness = 0;
//        boolean kill_first = false;
//        double last_best_fitness = 0;
//        int inti_counter = 0;
//        int no_progress = 0;
//        double score = 0;
//        
//        evaluate_pop(pop, sorted_pop, sumFitness);
//        pq_to_array(sorted_pop, pop);
//        
//        fill_ranks_on_sorted_pop(pop);
//  
//        while(evals_left > 0)
//        {
//        	
//        	
//        	//print(pop);
//        	Individual[] next_gen = new Individual[offsprings];
//            
//        	sumFitness = Utils.sumFitness(pop, pop.length);
//        	
//        	double[] sh = new double[pop_size];
//        	
//        	//FITNESS SHARING
//        	if(fitness_sharing)
//        	{
//            	adjust_fitnesses(sigma, pop, next_gen, sh);
//        	}
//        	
//        	//REPRODUCTION
//        	double chanse = rnd_.nextDouble();
//        	if(chanse > 0.25)
//        	{
//        		System.out.println("Uniform CrossOver");
//        		fitness_rouletteWheel_reproduction(pop, next_gen, uniCrossOver, offsprings, sumFitness, fitnessSelection);
//        		
//        	}
//        	else
//        	{
//        		System.out.println("Arithemtic CrossOver");
//        		fitness_rouletteWheel_reproduction(pop, next_gen, arithmeticCrossOver, offsprings, sumFitness, fitnessSelection);
//        		
//        	}
//        		
//        	//MUTATE
//        	int dims_to_mutate = 1;
//        	MutateChildren(next_gen, individuals_to_mutate, dims_to_mutate);
//        	
//        	sorted_pop.clear();
//        	fill_sorted_pop(fitness_sharing, pop, sorted_pop, sh);
//    		
//        	//CHILDREN EVALUATION
//        	childrenSumFitness = 0;
//            evaluate_pop(next_gen, sorted_pop, childrenSumFitness);
//           
//            if(score < pop[0].fitness)
//            {
//            	score = pop[0].fitness;
//            }
//            System.out.println("score "+score);
//            System.out.println("pop[0].fitness "+pop[0].fitness);
//            sumFitness = Utils.sumFitness(pop, pop.length);
//            System.out.println("SumFitness "+sumFitness+" last_pop_fitness "+last_pop_fitness);
//            System.out.println("difference "+ Math.abs(sumFitness - last_pop_fitness));
//            printer.printInfo(sumFitness, evals_left, pop[0].fitness, pop.length);
//            
//           
//            double epsilon = define_epsilon(sumFitness);
//            boolean noDifference = ((Math.abs(sumFitness - last_pop_fitness) < epsilon));
//            boolean pop_degenerated = (sumFitness > (pop[0].fitness * (pop_size/2 + pop_size/3)));
//            
////        	if(noDifference || pop_degenerated)
////            {
//             double difference = 120;
//             if(sumFitness < 1000)
//             {
//            	 difference = 8;
//             }
//             if(sumFitness >10000)
//             {
//            	 difference = 800;
//             }
//           
//             double total_distance = 0;
//             for(int i = 0; i < sh.length; i++)
//             {
//            	 total_distance += sh[i];
//             }
//             
//             total_distance /= (pop_size + offsprings);
//             double best_fitness = pop[0].fitness;
//             
//             System.out.println("total distance "+total_distance);
//             System.out.println("condition 1 "+(Math.abs(sumFitness - last_pop_fitness) < difference));
//             System.out.println("conditions 2 "+ (sumFitness > (best_fitness * (pop_size/2))));
//             System.out.println("condition 3 "+ (Math.abs(best_fitness - last_best_fitness) < 0.001));
//             System.out.println("init counter "+inti_counter);
//             System.out.println();
//             
//    		 if((Math.abs(sumFitness - last_pop_fitness) < difference) && sumFitness > (best_fitness * (pop_size/2)) &&
//    				 Math.abs(best_fitness - last_best_fitness) < (define_epsilon(last_best_fitness)/100))
//             {
//            	//System.out.println("childrenSumFitness " + (sumFitness/offsprings));
//             	//System.out.println("epsilon "+ (1 - (sumFitness/offsprings)/pop[0].fitness));
//    			inti_counter++;
//            	shock = true;
//            	Individual center_indiv = new Individual(pop[0]);
//            	double radious = 0.25;
//            	kill_first = false;
//            	epsilon = define_epsilon(center_indiv.fitness);
//            	if(Math.abs(center_indiv.fitness - previous_top_idniv_fitness) < epsilon)
//            	{
//            		radious = 0.25;
//            		kill_first = true;
//            		center_indiv = new Individual(pop[no_progress]);
//            		no_progress++;
//            	}
//            	else
//            	{
//            		no_progress = 0;
//            	}
//            	
//            	previous_top_idniv_fitness = center_indiv.fitness;
//            	printer.printFitnesses(pop, pop.length, rank_populations);
//            	initPop_around_individual(pop, center_indiv, radious);
//            	
//            	sorted_pop.clear();
//            	evaluate_pop(pop, sorted_pop, sumFitness);
//            	
//    			System.out.println("Pop Init -------------------------------------------------------------------------------------");
//            }
//    		last_best_fitness = best_fitness;
//    		
//        	System.out.println("noDifference "+noDifference);
//            System.out.println("pop_degenerated "+pop_degenerated);
//            System.out.println("kill first "+kill_first);
//            System.out.println();
//            last_pop_fitness = sumFitness;
//            //SURVIVOR SELECTION
//            pq_to_array(sorted_pop, pop);
//            //fill_ranks_on_sorted_pop(pop);
//           
//            //DIVERSITY CHECK
//    		//diversity_check(pop, sumFitness, sorted_pop, last_avg_children_ftiness, childrenSumFitness);
//        }  
//	}
//
//	private void adjust_fitnesses(double sigma, Individual[] pop, Individual[] next_gen, double[] sh) {
//		sh_values(sh, pop, next_gen, sigma);
//		
////            	for(int i  = 0; i <10000;i++) {
////            		System.out.println(sh[0]);
////            	}
//		
//		for(int indiv = 0; indiv < pop_size; indiv++)
//		{
//			pop[indiv].fitness = pop[indiv].fitness / (infinityProtection + sh[indiv]);
//		}
//	}
//
//	private void fill_sorted_pop(boolean fitness_sharing, Individual[] pop, PriorityQueue<Individual> sorted_pop,
//			double[] sh) {
//		for(int indiv = 0; indiv < pop_size; indiv++)
//		{
//			// restore fitness for survivor selection.
//			if(fitness_sharing)
//			{
//				pop[indiv].fitness = pop[indiv].fitness * (infinityProtection + sh[indiv]);
//			}
//			
//			Individual individual = new Individual(pop[indiv]);
//			sorted_pop.add(individual);
//		}
//	}
//	
//	/*
//	 * Fills the distances array which has length pop_size + offsprings.
//	 */
//	private void sh_values(double[] sh, Individual[] pop, Individual[] next_gen, double sigma)
//	{
//		for(int i = 0; i < pop.length; i++)
//		{
//			for(int j = i+1; j  < pop.length; j ++)
//			{
//				calculate_sh(pop[i], pop[j], i, j, sh, sigma);
//			}
//			
////			for(int j = 0; j < next_gen.length; j ++)
////			{
////				calculate_sh(pop[i], next_gen[j], i, j + pop.length, sh, sigma);
////			}
//		}
//		
////		for(int i = 0; i < next_gen.length; i ++)
////		{
////			for(int j = i + 1; j < next_gen.length; j ++)
////			{
////				calculate_sh(next_gen[i], next_gen[j], i + pop.length, j + pop.length, sh, sigma);
////			}
////		}
//	}
//	
//	private void calculate_sh(Individual a, Individual b, int index_a, int index_b, double[] sh, double sigma)
//	{
//		double total_distance = 0;
//		
//		for(int dim = 0; dim < dimensions; dim++)
//		{
//			if(a.genome[dim] > 0 && b.genome[dim] > 0)
//			{
//				total_distance += Math.abs(a.genome[dim] - b.genome[dim]);
//			}
//			else if(a.genome[dim] < 0 && b.genome[dim] < 0)
//			{
//				total_distance += Math.abs(Math.abs(a.genome[dim]) - Math.abs(b.genome[dim]));
//			}
//			else
//			{
//				total_distance += Math.abs(a.genome[dim] + b.genome[dim]);
//			}	
//		}
//		
//		double sh_value = 1 - total_distance/sigma;
//		sh_value = Math.max(sh_value, 0);
//		//System.out.println(" td "+total_distance+" index a "+index_a+" index b "+index_b + " sh_value "+sh_value) ;
//		
//		sh[index_a] += sh_value;
//		sh[index_b] += sh_value;
//	}
//	
//	private void fitness_rouletteWheel_reproduction(Individual[] pop, Individual[] next_gen, CrossOver crossOver, int offsprings, double sumFitness, FitnessSelection fitnessSelection )
//	{
//		System.out.println("fitness roulette Wheel ");
//		
//		for(int i = 0; i < offsprings; i++)
//    	{
//			int interval = 2;
////			if(evals_left < evaluations_limit_/20)
////			{
////				interval = pop.length;
////			}
//			Individual[] indivs = rouletteWheel_fitness_parentSelection(pop, fitnessSelection, interval);
//			
//			next_gen[i] = crossOver.cross_over(indivs[0], indivs[1]);
//			//System.out.println(next_gen[i].genome[0] +" "+ next_gen[i].genome[9]);
//    	}
//	}
//	
//	private Individual[] rouletteWheel_fitness_parentSelection(Individual[] pop, FitnessSelection fitnessSelection, int interval)
//	{
//		Individual parent_a = fitnessSelection.partial_rouletteWheel(pop, interval, pop.length);
//		Individual parent_b = fitnessSelection.partial_rouletteWheel(pop, interval, pop.length);
//		
//		while( parent_a.index == parent_b.index)
//		{
//			parent_b = fitnessSelection.partial_rouletteWheel(pop, 2, pop.length);
//		}
//		
//		Individual[] idnivs = {parent_a, parent_b};
//		return idnivs;
//	}
//	
//	private void randomReproduction(Individual[] pop, Individual[] next_gen, CrossOver crossOver, int offsprings, double sumFitness)
//	{
//		//System.out.println("random reproduction ");
//		//printer.printInfo(sumFitness, evals_left, pop[0].fitness, pop_size);
//		
//		for(int i = 0; i < offsprings; i++)
//    	{
//			int[] indexes = randomParentSelection(pop.length);
//			
//			next_gen[i] = crossOver.cross_over(pop[indexes[0]], pop[indexes[1]]);
//			//System.out.println(next_gen[i].genome[0] +" "+ next_gen[i].genome[9]);
//    	}
//	}
//	
//	private void rank_partial_roulette_with_No_replacement_reproduction(
//			Individual[] pop, RanksSelection rankSelection, Individual[] next_gen, int offsprings, CrossOver crossOver)
//	{
//		System.out.print(" partial no replacement ");
//		Individual[] copy_pop = new Individual[pop_size];
//		
//    	for(int j  = 0; j < pop_size; j++)
//    	{
//    		copy_pop[j] = new Individual(pop[j]);
//    	}
//    	
//		for(int i = 0; i < offsprings; i++)
//    	{
//			Individual child = new Individual();
//			
//			Individual parent_a = new Individual();
//			Individual parent_b = new Individual();
//			
//			rankParentSelection_NoReplacement(copy_pop, i, parent_a, parent_b, rankSelection);
//			
//			child = crossOver.cross_over(parent_a, parent_b);
//			next_gen[i] = child;
//    	}
//	}
//
//	private void rank_rouletteWheel_reproduction(Individual[] pop, double sumFitness, FitnessSelection fitnessSelection,
//			RanksSelection rankSelection, CrossOver arithmeticCrossover, Individual[] next_gen, int offsprings) 
//	{
//	
//		printer.printInfo(sumFitness, evals_left, pop[0].fitness, pop_size);
//		System.out.print(" roulette ");
//		
//
//		for(int i = 0; i < offsprings; i++)
//    	{
//			Individual child = new Individual();
//			
//			Individual parent_a = new Individual();
//			Individual parent_b = new Individual();
//			
//			rankRouletteWheel(parent_a, parent_b, rankSelection,pop);
//			
//			child = arithmeticCrossover.cross_over(parent_a, parent_b);
//			next_gen[i] = child;
//    	}
//	}
//	
//	public void rankParentSelection_NoReplacement(Individual[] copy_pop, int i, Individual parent_a, Individual parent_b, RanksSelection rankSelection )
//	{
//		int interval = 2;
//		int participants = copy_pop.length - i;
//		parent_a = rankSelection.partial_rouletteWheel(copy_pop, interval, participants);
//		int last_individual = participants - 1;
//		
//		copy_pop[parent_a.index] = copy_pop[last_individual];
//		
//		parent_b = rankSelection.partial_rouletteWheel(copy_pop, interval, participants -1);	
//		last_individual = participants - 2;
//	
//		copy_pop[parent_b.index] = copy_pop[last_individual];
//	}
//	
//	public int[] randomParentSelection(int pop_length)
//	{
//		
//		int index_parent_a = rnd_.nextInt(pop_length);
//		
//		int index_parent_b = rnd_.nextInt(pop_length);
//		
//		int[] indexes = {index_parent_a, index_parent_b};
//		return indexes;
//		
//	}
//	
//	public void rankRouletteWheel(Individual parent_a, Individual parent_b, RanksSelection rankSelection, Individual[] pop)
//	{
//		parent_a = rankSelection.rouletteWheel(pop, pop.length);
//		parent_b = rankSelection.rouletteWheel(pop, pop.length);
//			
//	}
//	
//
//	private void diversity_check(Individual[] pop, double sumFitness, PriorityQueue<Individual> sorted_pop,
//			double last_avg_children_ftiness, double childrenSumFitness)
//	{
//		// Defines epsilon given how small is the average fitness
//    	double epsilon = define_epsilon(last_avg_children_ftiness);
//    	double avg_children_fitness = childrenSumFitness/offsprings;
//    	System.out.println("last avg "+last_avg_children_ftiness+" current avg "+avg_children_fitness+" epsilon "+epsilon);
//    	Individual center_indiv = new Individual(sorted_pop.poll());
//		sorted_pop.clear();
//		if((Math.abs(pop[0].fitness - avg_children_fitness) < epsilon) && !shock)
//		{
//			if(last_avg_children_ftiness > 0.1)
//			{
//				initPop_around_individual(pop, center_indiv, 0.75);
//				evaluate_pop(pop, sorted_pop, sumFitness);
//				System.out.println("Pop Init -------------------------------------------------------------------------------------------------");
//			}
//		}
//	}
//
//	//TODO not working properly
//	private void survivorSelection(Individual[] pop, IndividualSelection fitnessSelection,
//			SurvivorSelection survivorSelection, Individual[] next_gen, PriorityQueue<Individual> sorted_pop)
//	{
////		next_gen = pq_to_array(sorted_pop);
////		fill_ranks_on_sorted_pop(next_gen);
////		if(evals_left < evaluations_limit_/100)
////		{
////			int elites = 0;
////			double chanse_to_randomly_choose_individual = 0;
////			survivorSelection.stochastic(pop, next_gen, elites, chanse_to_randomly_choose_individual, fitnessSelection);
////		}
////		else
////		{
//			 pq_to_array(sorted_pop, pop);
//			//survivorSelection.copy_populations(pop, next_gen);
////		}
//	}
//
//	private double define_epsilon(double last_avg_children_ftiness)
//	{
//		double epsilon = 0;
//		
//		if(last_avg_children_ftiness > 1)
//		{
//			epsilon = 0.01;		
//		}
//		else if((last_avg_children_ftiness * 10) > 1)
//		{
//			epsilon = 0.001;
//		}
//		else if((last_avg_children_ftiness * 100) > 1)
//		{
//			epsilon = 0.0001;
//		}
//		else if(last_avg_children_ftiness * 1000 > 1)
//		{
//			epsilon = 0.00001;
//		}
//		else if(last_avg_children_ftiness * 10000 > 1)
//		{
//			epsilon = 0.000001;
//		}else if(last_avg_children_ftiness * 100000 > 1)
//		{
//			epsilon = 0.0000001;
//		}
//		
//		return epsilon;
//	}
//	
//	/*
//	 * This method fills the individuals of the population with their corresponding rank.
//	 * Pop must be sorted.
//	 */
//	private void fill_ranks_on_sorted_pop(Individual[] pop) 
//	{
//		for(int i = 0; i < pop.length; i++)
//		{
//			pop[i].rank = Utils.getRank(i, pop.length);
//		}
//	}
//
//	private void evaluate_pop(Individual[] pop, PriorityQueue<Individual> sorted_pop, double sumFitness)
//	{
//		for(int i = 0; i < pop.length; i++)
//        {
//            Double fitness = (double) evaluation_.evaluate(pop[i].genome);
//            
//            pop[i].fitness = fitness;
//            pop[i].index = i;
//            
//            Individual indiv = new Individual(pop[i]);
//            
//            sorted_pop.add(indiv);
//            
//            sumFitness += fitness;
//            evals_left--;
//            
//            if(evals_left ==1 )
//    		{
//        		printer.printFitnesses(pop, pop.length, rank_populations);
//        		//printer.printRanks(pop, 1);
//        		printer.printInfo(sumFitness, evals_left, pop[0].fitness, pop_size);
//        		//printer.printBestIndivValues(pop, 100);
//    		}
//        }
//	}
//
//	
//	private void MutateChildren(Individual[] next_gen, int number_individuals_to_mutate, int number_dims_to_mutate) 
//	{
//		ArrayList<Integer> list_of_dimensions = new ArrayList<Integer>();
//		
//		for(int dim = 0; dim < dimensions; dim++)
//		{
//			list_of_dimensions.add(dim);
//		}
//		
//		Collections.shuffle(list_of_dimensions);
//		
//		// Randomly pick number_individuals_to_mutate individuals and store their index in the to_mutate array.
//		int[] to_mutate = new int[number_individuals_to_mutate];
//		int max = offsprings;
//		int min = 0;
//		for(int i = 0; i < number_individuals_to_mutate; i++)
//		{
//			to_mutate[i] = rnd_.nextInt(max-min) +min;
//		}
//		
//		for(int i = 0; i < number_individuals_to_mutate; i++)
//		{
//			int indiv = to_mutate[i];
//			
//			for(int dim = 0; dim < number_dims_to_mutate; dim++)
//			{
//				int dim_to_mutate = list_of_dimensions.get(dim);
//				mutate_individual(next_gen[indiv], dim_to_mutate);
//			}
//		}
//	}
//	
//	private void initPop_around_individual(Individual[] pop, Individual center_indiv, double radious)
//	{
//		//System.out.println("Reform pop around best Individual");
//		
//		pop[0] = new Individual(center_indiv);
//		int ignore_first = 1;
//		
//		for(int individual = ignore_first; individual < pop_size; individual++)
//        {	
//			
//			for(int dim = 0; dim < dimensions; dim++)
//			{
//				double center = center_indiv.genome[dim] ;
//				pop[individual].genome[dim] =  Utils.double_in_range(Math.min(center+radious, 5), Math.max(-5, center-radious));
//				//pop[individual].mutation_steps[dim] = Utils.double_in_range(upper_mutation_step, lower_mutation_step);
//			}
//			pop[individual].fitness = 0;
//			pop[individual].index = individual;
//			pop[individual].rank = 0;
//			
//        }
//	}
//	
//	/*
//	 * Fills the individual array with random values from -5 to 5 for genome and 0.05 to 0.25 for mutation_steps.
//	 */
//	private void initPop(Individual[] pop)
//	{
//		double upper_mutation_step_bound = 0.25;
//		double lower_mutation_step_bound = 0.05;
//		
//		for(int individual = 0; individual < pop_size; individual++)
//        {
//			pop[individual] = new Individual(uper_bound, lower_bound, upper_mutation_step_bound, lower_mutation_step_bound);
//        }
//	}
//
//	public void pq_to_array(PriorityQueue<Individual> pq, Individual[] pop)
//	{
//		int counter = 0;
//		
//		while(counter < pop.length)
//		{	
//			pop[counter] = pq.poll();				
//			counter++;
//		}	
//	}
//	
//	public void mutate_individual(Individual indiv, int dim_to_mutate)
//	{
//		double chanse = rnd_.nextDouble();
//		double sign = Utils.getSign();
//		double change = sign * indiv.mutation_steps[dim_to_mutate];
//		
//		if(chanse > 0.5)
//		{
//			if(change < 0)
//			{
//				indiv.genome[dim_to_mutate] = Math.max(-5, indiv.genome[dim_to_mutate] + change);
//			}
//			else
//			{
//				indiv.genome[dim_to_mutate] = Math.min(5, indiv.genome[dim_to_mutate] + change);
//			}
//		}
//	}
		
	public static void main(String args[])
	{
		System.out.println("start!");
	}
	
	public double public_wellfare(int jobs, double precision, double[] center, int participants, boolean directed, double center_fitness)
	{
		double[][] community = new double[10][participants];
		//double[] start = new double[10];
		boolean[] neighbours = new boolean[dimensions];
		
//		if(directed)
//		{
//			for(int i = 0; i < dimensions; i++)
//			{
//				double[] neighbour = copy_center(center);
//				neighbour[i] = neighbour[i]-precision;
//
//				Double fitness = (double) evaluation_.evaluate(neighbour);
//				
//				neighbours[i] = (fitness > center_fitness);
//			}
			gather_participants(precision, center, community);
			//directed_gather_participants(precision, start, community);
//		}
//		
//		else
//		{
//			gather_participants(precision, center, community);
//		}
		
		double[][] fitnesses = new double[10][participants];
		init_fitnesses(0.00000000000000000000000000000000000000000000000000000001, fitnesses, community);
		
		int evals = 0;
		double best_score = 0;
		while(evals < jobs)
		{
//			System.out.println();
//			System.out.println("evals "+evals);
//			System.out.println();
			double[] job = new double[dimensions];
			int[] currents = new int[dimensions];
			
			for(int dim = 0; dim < dimensions; dim++)
			{
				int participant = 0;
				//participant = rnd_.nextInt(community[dim].length);
				if(evals>(9*(jobs/10)))
				{
					participant = rouletteWheel(fitnesses[dim]);
					while(fitnesses[dim][participant] < 0)
					{
						participant = rnd_.nextInt(community[dim].length);
					}
				}
				else if(evals>(4*(jobs/5)))
				{
					participant = partial_rouletteWheel(fitnesses[dim], 2);
					while(fitnesses[dim][participant] < 0)
					{
						participant = rnd_.nextInt(community[dim].length);
					}
				}
				else
				{
					participant = rnd_.nextInt(community[dim].length);
					
					while(fitnesses[dim][participant] < 0)
					{
						participant = rnd_.nextInt(community[dim].length);
					}
				}
//				
				currents[dim] = participant;
				job[dim] = community[dim][participant];
				
			}
			
			if(evals_left ==0)
			{
				return -1000;
			}
			
			
			Double fitness = (double) evaluation_.evaluate(job);
			evals_left--;
			
			if(fitness > best_score)
			{
				center_fitness = fitness;
				best_score = fitness;
				for(int i = 0; i< dimensions; i++)
				{
					center[i] = job[i];
				}
			}
			// Reword workers!
			for(int dim = 0; dim < dimensions; dim++)
			{
				
				int participant = currents[dim];
				
				fitnesses[dim][participant] += fitness;
				
				//one step neighbours
//				if(participant > 0)
//				{
//					if(neighbours[dim])
//					{
//						fitnesses[dim][participant-1] += (fitness + fitness/8);
//					}
//					else
//					{
//						fitnesses[dim][participant-1] += (fitness - fitness/8);
//					}
//				}
//				else if(participant<(participants-1))
//				{
//					if(neighbours[dim])
//					{
//						fitnesses[dim][participant+1] += (fitness - fitness/8);
//					}
//					else
//					{
//						fitnesses[dim][participant+1] += (fitness + fitness/8);
//					}
//				}
//				
//				//two step neighbours
//				if(participant>1)
//				{
//					if(neighbours[dim])
//					{
//						fitnesses[dim][participant-2] += (fitness + fitness/4);
//					}
//					else
//					{
//						fitnesses[dim][participant-2] += (fitness - fitness/4);
//					}
//				}
//			    if(participant<(participants-2))
//				{
//			    	if(neighbours[dim])
//					{
//						fitnesses[dim][participant+2] += (fitness - fitness/4);
//					}
//					else
//					{
//						fitnesses[dim][participant+2] += (fitness + fitness/4);
//					}
//				}
////				
//				//three step neighbours
//				if(participant>2)
//				{
//					if(neighbours[dim])
//					{
//						fitnesses[dim][participant-3] += (fitness + fitness/2);
//					}
//					else
//					{
//						fitnesses[dim][participant-3] += (fitness - fitness/2);
//					}
//					
//				}
//				if(participant<(participants-3))
//				{
//					if(neighbours[dim])
//					{
//						fitnesses[dim][participant+3] += (fitness - fitness/2);
//					}
//					else
//					{
//						fitnesses[dim][participant+3] += (fitness + fitness/2);
//					}
//				}
				if(participant > 0)
				{
					fitnesses[dim][participant-1] += (fitness/2);
					
				}
				else if(participant<(participants-1))
				{
					fitnesses[dim][participant+1] += ( fitness/2);
					
				}
				
				//two step neighbours
				if(participant>1)
				{
					fitnesses[dim][participant-2] += ( fitness/4);
					
				}
			    if(participant<(participants-2))
				{
			    	fitnesses[dim][participant+2] += ( fitness/4);
			    	
				}
//				
				//three step neighbours
				if(participant>2)
				{
					fitnesses[dim][participant-3] += ( fitness/8);
					
					
				}
				if(participant<(participants-3))
				{
					fitnesses[dim][participant+3] += (fitness/8);
					
				}
				
				// four step neighbours
				if(participant>3)
				{
					fitnesses[dim][participant-4] += ( fitness/16);
					
					
				}
				if(participant<(participants-4))
				{
					fitnesses[dim][participant+4] += (fitness/16);
					
				}

				
				//System.out.println("dim "+dim+" value "+job[dim] +" participant "+participant+" fitness "+fitness);
			}
			
//			if(evals % 100000 == 0)
//			{
//				System.out.println("best score "+bestscore);
//			}
			evals++;
		}
		
		
//		for(int dim = 0; dim < dimensions; dim++)
//		{
//			double max_fitness = 0;
//			int max_index = 0;
//			for(int i = 0; i < fitnesses[dim].length; i++)
//			{
//				if(fitnesses[dim][i] > max_fitness)
//				{
//					max_fitness = fitnesses[dim][i];
//					max_index = i;
//				}
//			}
//			center[dim] = community[dim][max_index];
//		}
		
		
		return center_fitness;
	}
	
	public void directed_gather_participants(double precision, double[] start, double[][] community) 
	{   
		
		for(int dimension = 0; dimension < dimensions; dimension++)
		{
			double dim_start = start[dimension];
			
			for(int participant = 0; participant < community[0].length; participant ++)
			{
				if(dim_start >5 || dim_start <-5)
				{
					community[dimension][participant] = 1000;
				}
				else
				{
					community[dimension][participant] = dim_start;
				}
				
				dim_start += precision;
			}
		}
	}
	
	public void gather_participants(double precision, double[] start, double[][] community) 
	{   
		
		for(int dimension = 0; dimension < dimensions; dimension++)
		{
			double dim_start = start[dimension] - ((community[0].length/2) * precision);
			
			for(int participant = 0; participant < community[0].length; participant ++)
			{
				if(dim_start >5 || dim_start <-5)
				{
					community[dimension][participant] = 1000;
				}
				else
				{
					community[dimension][participant] = dim_start;
				}
				
				dim_start += precision;
			}
		}
	}
	
	public void init_fitnesses(double competence_pressure, double[][] fitnesses, double[][] community)
	{
		for(int dimension = 0; dimension < dimensions; dimension++)
		{
			for(int participant = 0; participant < fitnesses[0].length; participant ++)
			{
				if(community[dimension][participant] > 500)
				{
					fitnesses[dimension][participant] = -5;
				}
				else
				{
					fitnesses[dimension][participant] = competence_pressure;
				}
			}
		}
	}
	
	public int rouletteWheel(double[] fitnesses)
	{
		double sumFitness = sumFitness(fitnesses);
		
		double rand = Utils.randomNumber(sumFitness);
		
		for(int i = 0; i < fitnesses.length; i++)
		{
			rand -= fitnesses[i];
			
			if(rand < 0)
			{
				return i;
			}
		}
		
		return rnd_.nextInt(fitnesses.length);
	}
	
	public double sumFitness(double[] fitnesses)
	{
		double sumFitness = 0;
		
		for(int i = 0; i < fitnesses.length; i++)
		{
			sumFitness += fitnesses[i];
		}
		
		return sumFitness;
	}
	
	public int partial_rouletteWheel(double[] fitnesses, int interval)
	{
		double[] pop_interval = new double[interval];
		int[] indexes = new int[interval];

		for(int i = 0; i < interval; i++)
		{
			int choosen = player72.rnd_.nextInt(fitnesses.length);
			
			pop_interval[i] = fitnesses[choosen];
			indexes[i] = choosen;
		}
		
		int choosen_indiv = rouletteWheel(pop_interval);
	
		return indexes[choosen_indiv];
	}
	
	public double[] copy_center(double[] center)
	{
		double[] neighbour = new double[10];
		for(int dim = 0; dim < dimensions; dim++)
		{
			neighbour[dim] = center[dim];
		}
		
		return neighbour;
	}
	
	
}