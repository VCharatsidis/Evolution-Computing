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
		
	}
	
		
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
//				if(participant>3)
//				{
//					fitnesses[dim][participant-4] += ( fitness/16);
//					
//					
//				}
//				if(participant<(participants-4))
//				{
//					fitnesses[dim][participant+4] += (fitness/16);
//					
//				}

				
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