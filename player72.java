import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import java.util.Random;

import javax.naming.PartialResultException;

import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.Properties;

public class player72 implements ContestSubmission
{
	Random random = new Random();
	int pop_size = 100;
	Random rnd_;
	ContestEvaluation evaluation_;
    private int evaluations_limit_;
	
	public player72()
	{
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
        
		// Run your algorithm here
		System.out.println("eval limit "+evaluations_limit_);
		
        int evals = 0;
        // init population
        double[][] pop = new double[pop_size][10];
        
        for(int individual = 0; individual < pop_size; individual++)
        {
        	pop[individual] = new double[10];
        	
        	for(int dim = 0; dim < 10; dim++)
        	{
        		int rangeMax = 10;
            	double randomValue = (rangeMax * random.nextDouble()) - 5;
            	
        		pop[individual][dim] = randomValue;
        	}
        	
        }
        // init population
        // calculate fitness
        PriorityQueue<FitnessIndex> fitnesses_pq = new PriorityQueue<FitnessIndex>(new FitnessComparator());
        double sumFitness = 0;
        int runs = 0 ;
        int total_runs = evaluations_limit_/pop_size;
        while(runs < total_runs)
        {
        	fitnesses_pq.clear();
        	sumFitness = 0;
        	evals = 0;
        	double[] fitnesses = new double[pop_size];
        	double[][] next_gen = new double[pop_size][10];
        	
        	//if(runs ==100) {
//        	for(int j = 0; j<1;j++)
//			{
//				for(int i = 0; i < 10; i++)
//    			{
//    				System.out.print(pop[evals][i] +" ");
//    			}
//    			System.out.println("indiv  "+j);
//			}
        	//}
        	while(evals < pop_size)
            {

                Double fitness = (double) evaluation_.evaluate(pop[evals]);
               
                FitnessIndex fitnessesIndexes = new FitnessIndex(evals, fitness);
                fitnesses_pq.add(fitnessesIndexes);
                
                fitnesses[evals] = fitness;
                sumFitness += fitness;
                evals++;
                // Select survivors
            }
        	System.out.println("run "+runs);
        	
        	FitnessIndex[] fi = pq_to_array(fitnesses_pq);
        	
//        	System.out.println("fitness "+fi[0].fitness);
//        	for(int dim = 0; dim < 10; dim++)
//        	{
//        		System.out.print("dim "+dim+" "+pop[fi[0].index][dim]);
//        	}
        	
        	for(int i = 0; i < pop_size; i++)
        	{
        		int parent_a = 0;
        		int parent_b = 0;
        		double[] child = new double[10];
        		
        		
//        		if(runs < 5)
//        		{
//        			//indexes when shorted according to fitness.
//        		    int shorted_index_a = random.nextInt(20);
//        		    int shorted_index_b = random.nextInt(20);
//        			
//        			parent_a = fi[shorted_index_a].index;
//        			parent_b = fi[shorted_index_b].index;
//        			
//        			child = new double[10];
//            		uniformCrossOver(pop[parent_a], pop[parent_b], child);
//            		for(int j = 0; j < 10; j++)
//            		{
//            			next_gen[i][j] = child[j];
//            		}	
//        			
//        		}else {
//        			 parent_a = rouletteWheel(sumFitness, fitnesses);
//            		 parent_b = rouletteWheel(sumFitness, fitnesses);
//            		 child = new double[10];
//             		 uniformCrossOver(pop[parent_a], pop[parent_b], child);
//        		}
//        		
        		if(runs > (total_runs * 1))
				{
        			parent_a = rouletteWheel(sumFitness, fitnesses);
        			parent_b = rouletteWheel(sumFitness, fitnesses);
				}
        		else {
        			//Partial wheel
            		int interval = 2;
            		
            		int start_a = random.nextInt(pop_size-interval+1);
            		parent_a = partial_rouletteWheel(fitnesses, start_a, interval);
            		
            		int start_b = random.nextInt(pop_size-interval+1);
            		parent_b = partial_rouletteWheel(fitnesses, start_b, interval);
            		
				}
        		
        		uniformCrossOver(pop[parent_a], pop[parent_b], child);
       
        		//System.out.println("parent a "+parent_a + " parent b "+parent_b);
        		
        		
        		for(int j = 0; j < 10; j++)
        		{
        			next_gen[i][j] = child[j];
        		}	
        	}
        	
        	//MUTATE
        	int individuals_to_mutate = 10;
        	int[] to_mutate = new int[individuals_to_mutate];
        	
        	for(int i = 0; i <individuals_to_mutate; i++)
        	{
        		to_mutate[i] = random.nextInt(100);
        	}
        	
        	for(int i = 0; i < individuals_to_mutate; i++)
        	{
        		int indiv = to_mutate[i];

        		mutate(next_gen[indiv]);

        	}
        	
        	pop = transfer_gen(next_gen);

        	runs++;
        	if(runs % 5 ==0)
        	{
        		System.out.println("run "+runs +" avg fitness : "+sumFitness/pop_size);
        	}
        	
        }  
	}
	
	public double[][] transfer_gen(double[][] next_gen)
	{
		double[][] pop = new double[pop_size][10];
		
		for(int i = 0; i < pop_size; i++)
    	{
    		for(int j = 0; j < 10; j++)
    		{
    			pop[i][j] = next_gen[i][j];
    		}	
    	}
		return next_gen;
	}
	
	public FitnessIndex[] pq_to_array(PriorityQueue<FitnessIndex> pq)
	{
		//Iterator<FitnessIndex> fitness_iterator = pq.iterator();
		FitnessIndex[] fi = new FitnessIndex[pop_size];
		int counter = 0;
		
		while(!pq.isEmpty())
		{
			
			fi[counter] = pq.poll();
			if(counter < 100) {
				//System.out.println("index "+ fi[counter].index +" fitness "+fi[counter].fitness);
			}
				
			counter++;
		}
		
		return fi;
	}
	
	public void mutate(double[] indiv)
	{
		int dim = random.nextInt(10);
		double chanse = random.nextDouble();
		double change = 0.05;
		
		if(chanse >0.5 && ((indiv[dim]+change) <5))
		{
			indiv[dim] = indiv[dim] + change;
		}
		else if(((indiv[dim]-change) > -5)){
			indiv[dim] = indiv[dim] - change;
		}
	}
	
	public void uniformCrossOver(double[] parent_a, double[] parent_b, double[] child)
	{
			
			
			for(int i = 0; i < 10; i++)
			{
				double chanse = random.nextDouble();
				if(chanse > 0.5 )
				{
					child[i] = parent_a[i];
				}
				else {
					child[i] = parent_b[i];
				}
			}
				
	}
	
	public double[] crossOver(double[] parent_a, double[] parent_b)
	{
		double child[] = new double[10];
		int crossOverPoint = 4;
		
		for(int i = 0; i < crossOverPoint; i++)
		{
			child[i] = parent_a[i];
		}
		for(int i = child.length-crossOverPoint; i < child.length; i++)
		{
			child[i] = parent_b[i];
		}
		
		return child;
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
	/*
	 * It will return an individual from start to start+interval proportionally to his fitness.
	 * The fitest individual in that interval has higher chanses to be picked.
	 */
	public int partial_rouletteWheel(double[] fitnesses, int start, int interval)
	{
		double[] fitnesses_interval = new double[interval];
		//gather interval consecutive fitnesses
		for(int i = 0; i < interval; i++)
		{
			fitnesses_interval[i] = fitnesses[start+i];
		}
		
		double sumFitness_interval = sumFitness(fitnesses_interval);
		int choosen_indiv = rouletteWheel(sumFitness_interval, fitnesses_interval);
	
		int parent = start+choosen_indiv;
		return parent;
	}
	
	public int rouletteWheel(double sumFitness,  double[] fitnesses)
	{
		
		double rand = randomNumber(sumFitness);
		
		for(int i = 0; i < fitnesses.length; i++)
		{
			rand -= fitnesses[i];
			
			if(rand < 0)
			{
				return i;
			}
		}
		return 0;
	}
	
//	public int rouletteWheel(int sumFitness,  Iterator<FitnessIndex> population)
//	{
//		
//		int chosen_indiv = 0;
//		int rand = random.nextInt(50);
//		
//		int counter =0;
//		while(population.hasNext())
//		{
//			if(counter == rand)
//			{
//				
//			}
//			FitnessIndex fitnessIndex = population.next();
//			double fitness = fitnessIndex.fitness;
//			rand -= fitness;
//			if(rand < 0)
//			{
//				return fitnessIndex.index;
//			}
//		}
		
		//return chosen_indiv;
	//}
	
	public double randomNumber(double sumFitness)
	{	
		return sumFitness * random.nextDouble();	
	}
	
	public static void main(String args[])
	{
		System.out.println("start!");
	}
	
	
}


