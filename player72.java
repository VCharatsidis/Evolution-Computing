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
        }
    }
    
	public void run()
	{

		long start = System.currentTimeMillis();    
		
		
		double[] center = {0,0,0,0,0 ,0,0,0,0,0};
		double prec = 0.05;
		int epoch = 0;
		double best_score = 0;
		
		int runs = 400;
		double gamma = 0.95;
		
//		if(schaffers)
//		{
//			runs = 1200;
//			gamma = 0.75;
//		}
//		else if(katsuura)
//		{
//			runs = 450;
//			gamma = 0.98;
//			//beta = 0.05;
//			
//		}
//		
		while(true)
		{
			
			double current_score = 0;

			current_score = public_wellfare(runs, prec, center);

			epoch++;
			prec *= gamma;
			
			prec = Math.max(5*0.000001, prec);
			

			if(current_score > best_score)
			{
				best_score = current_score;
			}


//		    System.out.println("precision "+prec+" gamma "+gamma +" runs "+runs + " best_score "+best_score);
//			System.out.println(best_score);
//			System.out.println("------------------------------------------------------------------");
		    long elapsedTime = System.currentTimeMillis() - start;
		    if(elapsedTime > 12000)
		    {
		    	break;
		    }
		}
		
		//Katsuura();
	}
	

		
	public static void main(String args[])
	{
		System.out.println("start!");
	}
	
	public double public_wellfare(int jobs, double precision, double[] center)
	{
		double[][] community = new double[10][200];
		gather_participants(precision, center, community);
		
		double[][] fitnesses = new double[10][200];
		//init_fitnesses(0.00000000000000000000000000000000000000000000000000000001, fitnesses, community);
		
		int evals = 0;
		double bestscore = 0;
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
					while(community[dim][participant] > 500)
					{
						participant = rouletteWheel(fitnesses[dim]);
					}
				}
				else if(evals>(2*(jobs/3)))
				{
					participant = partial_rouletteWheel(fitnesses[dim], 2);
					while(community[dim][participant] > 500)
					{
						participant = partial_rouletteWheel(fitnesses[dim], 2);
					}
					
				}
				else
				{
					participant = rnd_.nextInt(community[dim].length);
					while(community[dim][participant] > 500)
					{
						participant = rnd_.nextInt(community[dim].length);
					}
					
				}
//				
				currents[dim] = participant;
				job[dim] = community[dim][participant];
				
			}
			
			Double fitness = (double) evaluation_.evaluate(job);
			
			if(fitness > bestscore)
			{
				bestscore = fitness;
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
				//fitness /= 2;
				if(participant > 0)
				{
					fitnesses[dim][participant-1] += fitness/4;
					//fitness /= 3;
				}
				if(participant<199)
				{
					fitnesses[dim][participant+1]+=fitness/4;
					//fitness /= 3;
				}
				
				if(participant>1)
				{
					fitnesses[dim][participant-2] += fitness/8;
				}
				
				if(participant<198)
				{
					fitnesses[dim][participant+2]+= fitness/8;
				}
				
				//System.out.println("dim "+dim+" value "+job[dim] +" participant "+participant+" fitness "+fitness);
			}
			
			
			evals++;
		}
		
		
		return bestscore;
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
					fitnesses[dimension][participant] = 0;
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
}	
	