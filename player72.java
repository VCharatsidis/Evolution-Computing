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
	int pop_size = 1000;
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
		
       
        // init population
        double[][] pop = new double[pop_size][10];
        
        for(int individual = 0; individual < pop_size; individual++)
        {
        	pop[individual] = new double[10];
        	
        	for(int dim = 0; dim < 10; dim++)
        	{
        		double rangeMax = 10;
        		double rangeMin = -5;
//        		if(dim == 0)
//        		{
//        			rangeMax = 1.5;
//        			rangeMin = -2;
//        		}
//        		else if(dim == 1)
//        		{
//        			rangeMax = 1.5;
//        			rangeMin = 3;
//        		}
//        		else if(dim == 2)
//        		{
//        			rangeMax = 2;
//        			rangeMin = -2;
//        		}
//        		else if(dim == 3)
//        		{
//        			rangeMax = 2;
//        			rangeMin = -4;
//        		}
//        		else if(dim == 4)
//        		{
//        			rangeMax = 3;
//        			rangeMin = 0;
//        		}
//        		else if(dim == 5)
//        		{
//        			rangeMax = 3;
//        			rangeMin = -3;
//        		}
//        		else if(dim == 6)
//        		{
//        			rangeMax = 3;
//        			rangeMin = 0;
//        		}
//        		else if(dim == 7)
//        		{
//        			rangeMax = 3;
//        			rangeMin = -2;
//        		}
//        		else if(dim == 8)
//        		{
//        			rangeMax = 2;
//        			rangeMin = 1;
//        		}
//        		else if(dim == 9)
//        		{
//        			rangeMax = 1.5;
//        			rangeMin = -1;
//        		}
            	double randomValue = (rangeMax * random.nextDouble()) + rangeMin;
            	
        		pop[individual][dim] = randomValue;
        	}
        	
        }
        // init population
        // calculate fitness
      
        PriorityQueue<FitnessIndex> fitnesses_pq = new PriorityQueue<FitnessIndex>(new FitnessComparator());
        double sumFitness = 0;
        int runs = 0 ;
        int offsprings = pop_size/2;
        int total_runs = (evaluations_limit_ - pop_size) / (offsprings);
        double[] fitnesses = new double[pop_size];
        
        int evals = 0;
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
        
        while(runs < total_runs)
        {

        	double[][] next_gen = new double[offsprings][10];
        	
        	for(int i = 0; i < offsprings; i++)
        	{
        		double[] child = new double[10];
        		
        		int parent_a = 0;
        		int parent_b = 0;
            	int island_length = pop_size/4;
            	
            	int ten_percent = total_runs/10;
            	boolean partial_border = (runs <= ten_percent || runs % 20 ==0);
            	boolean partial_rouletteWheel = (partial_border);
            	
            	if(true) {
            		if( i == 0) {
            			System.out.println();
                		System.out.println("run "+runs +" avg fitness : "+sumFitness/pop_size);
            		}
            		
        			//Partial wheel
            		int interval = 2;
            		
            		parent_a = partial_rouletteWheel(fitnesses, interval);
            		parent_b = partial_rouletteWheel(fitnesses, interval);	
            		if(i ==0) {
            			System.out.print(" partial ");
            		}
            		
				}
            	
            	boolean island_border = (runs < ten_percent * 9);
            	boolean island_wheel = (!partial_rouletteWheel && island_border);
            	
            	if(false) {
            		if( i == 0) {
            			System.out.println();
                		System.out.println("run "+runs +" avg fitness : "+sumFitness/pop_size);
            		}
            		int start = ((int)(i/island_length))*island_length;
            		System.out.println("start "+start);
            		//if(runs % 5 ==0) {
//            			System.out.println();
//            			System.out.println(" 0 ");
//            			for(int ij =0; ij<10; ij++) {
//            				System.out.print(pop[0][ij]+" ");
//            				
//        
//            			}
//            			System.out.println();
//            			System.out.println(" isnaldn_length ");
//            			for(int ij =0; ij<10; ij++) {
//            			
//            				System.out.print(pop[island_length][ij]+ " ");
//
//            			}
//            			System.out.println();
//            			System.out.println(" 2 isnaldn_length ");
//            			for(int ij =0; ij<10; ij++) {
//            		
//            			
//            				System.out.print(pop[island_length*2][ij]+" ");
//            				
//            				
//            			}
//            			System.out.println();
//            			System.out.println(" 3 isnaldn_length ");
//            			for(int ij =0; ij<10; ij++) {
//            		
//            				System.out.print(pop[island_length*3][ij]+" ");
//            			}
//            			System.out.println();
            		if(runs %2 ==0) 
            		{
            			double[] temp = pop[0];
            			double fitness_temp = fitnesses[0];
            			
            			pop[0] = pop[island_length];
            			fitnesses[0] = fitnesses[island_length];
            			
            			pop[island_length] = temp;
            			fitnesses[island_length] = fitness_temp;
            			
            			temp = pop[island_length * 2];
            			fitness_temp = fitnesses[island_length * 2];
            			
            			pop[island_length * 2] = pop[island_length * 3];
            			fitnesses[island_length * 2] = fitnesses[island_length * 3];
            			
            			pop[island_length * 3] = temp;
            			fitnesses[island_length * 3] = fitness_temp;
            		}else 
            		{
            			double[] temp = pop[0];
            			double fitness_temp = fitnesses[0];
            			
            			pop[0] = pop[island_length * 2];
            			fitnesses[0] = fitnesses[island_length * 2];
            			
            			pop[island_length * 2] = temp;
            			fitnesses[island_length * 2] = fitness_temp;
            			
            			temp = pop[island_length];
            			fitness_temp = fitnesses[island_length];
            			
            			pop[island_length] = pop[island_length * 3];
            			fitnesses[island_length] = fitnesses[island_length * 3];
            			
            			pop[island_length * 3] = temp;
            			fitnesses[island_length * 3] = fitness_temp;
            		}
            			
//            			System.out.println();
//            			System.out.println(" 0 ");
//            			for(int ij =0; ij<10; ij++) {
//            				System.out.print(pop[0][ij]+" ");
//            				
//        
//            			}
//            			System.out.println();
//            			System.out.println(" isnaldn_length ");
//            			for(int ij =0; ij<10; ij++) {
//            			
//            				System.out.print(pop[island_length][ij]+ " ");
//
//            			}
//            			System.out.println();
//            			System.out.println(" 2 isnaldn_length ");
//            			for(int ij =0; ij<10; ij++) {
//            		
//            			
//            				System.out.print(pop[island_length*2][ij]+" ");
//            				
//            				
//            			}
//            			System.out.println();
//            			System.out.println(" 3 isnaldn_length ");
//            			for(int ij =0; ij<10; ij++) {
//            		
//            				System.out.print(pop[island_length*3][ij]+" ");
//            			}
//            			System.out.println();
            			
            		//}
            		parent_a = island_wheel(fitnesses, start, island_length);
            		parent_b = island_wheel(fitnesses, start, island_length);
            		if(i ==0) {
            			System.out.print(" island ");
            		}
            		
            	}
            	
        		if(!island_border)
				{
        			if( i == 0) {
            			System.out.println();
                		System.out.println("run "+runs +" avg fitness : "+ sumFitness/pop_size);
            		}
        			parent_a = rouletteWheel(sumFitness, fitnesses);
        			parent_b = rouletteWheel(sumFitness, fitnesses);
        			if(i ==0) {
        				System.out.print(" roulette ");
            		}
        			
				}
        		
//        		if(runs>50) {
//        			avg_CrossOver(pop[parent_a], pop[parent_b], child);
//        		}
//        		else {
//        			uniformCrossOver(pop[parent_a], pop[parent_b], child);
//        		}
        		//avg_CrossOver(pop[parent_a], pop[parent_b], child);
        		if(runs % 3 == 0)
        		{
        			avg_CrossOver(pop[parent_a], pop[parent_b], child);
        		}else {
        			uniformCrossOver(pop[parent_a], pop[parent_b], child);
        		}
        		
        		//System.out.println("parent a "+parent_a + " parent b "+parent_b);
        		
        		for(int j = 0; j < 10; j++)
        		{
        			next_gen[i][j] = child[j];
        		}	
        	}
        	
        	//MUTATE
        	int individuals_to_mutate = offsprings/5;
        	int[] to_mutate = new int[individuals_to_mutate];
        	
        	for(int i = 0; i <individuals_to_mutate; i++)
        	{
        		//to_mutate[i] = random.nextInt(pop_size);
        		int max = offsprings;
        		int min = 0;
        		
        		to_mutate[i] = random.nextInt(max-min) +min;
        	}
        	
        	for(int i = 0; i < individuals_to_mutate; i++)
        	{
        		int indiv = to_mutate[i];
        	    mutateAllDims(next_gen[indiv]);
        		//mutate(pop[indiv]);
        		//mutate(next_gen[indiv]);
        		//mutate(next_gen[indiv]);
        	}
        	
          	fitnesses_pq.clear();
        	evals = 0;
            while(evals < pop_size)
            {
                FitnessIndex fitnessesIndexes = new FitnessIndex(evals,  fitnesses[evals] );
                fitnesses_pq.add(fitnessesIndexes);
                evals++;
                // Select survivors
            }
            
            int children = pop_size;
    		while(children < (pop_size + offsprings))
            {
                Double fitness = (double) evaluation_.evaluate(next_gen[children - pop_size]);
               
                FitnessIndex fitnessesIndexes = new FitnessIndex(children, fitness);
                fitnesses_pq.add(fitnessesIndexes);
                
                children++;
                // Select survivors
            }
        
        	System.out.println("run "+runs);
        	FitnessIndex[] fi = new FitnessIndex[pop_size];
        	
        	fi = pq_to_array(fitnesses_pq);
        	
        	double[][] sorted = new double[pop_size][10];
        	transfer_gen(pop, next_gen, fi, fitnesses, sorted);
        	
        	for(int ind = 0; ind < pop_size; ind++)
        	{
        		for(int j = 0; j < 10; j++) {
        			pop[ind][j] = sorted[ind][j];
        		}
        	}
        	
        	runs++;
        	if(runs % 50 ==0)
        	{
        		System.out.println();
        		sumFitness = 0;
        		for( int fit = 0; fit < pop_size; fit++)
        		{
        			sumFitness += fitnesses[fit];
        		}
        		
        		System.out.println("run "+ runs +" avg fitness : "+ (sumFitness/pop_size));
            	System.out.println("fitness "+fi[0].fitness);
            	
            	for(int dim = 0; dim < 10; dim++)
            	{
            		int index = fi[0].index;
            		if(index >= pop_size)
            		{
            			System.out.print("dim "+dim+" "+pop[index - pop_size][dim]);
            		}
            		else
            		{
            			System.out.print("dim "+dim+" "+pop[index][dim]);
            		}
            		
            	}
            	System.out.println();
        	}
        	
        }  
	}
	
	public void transfer_gen(double[][] pop, double[][] next_gen,
			FitnessIndex[] fi, double[] fitnesses, double[][] sorted)
	{
		for(int i = 0; i < pop_size; i++)
    	{
			fitnesses[i] = fi[i].fitness;
			int index = fi[i].getIndex();
			
    		for(int j = 0; j < 10; j++)
    		{
    			if(index >= pop_size)
    			{
    				sorted[i][j] = next_gen[index - pop_size][j];
    			}
    			else
    			{
    				sorted[i][j] = pop[index][j];
    			}	
    		}	
    	}
	
		//return sorted;
	}
	
	public FitnessIndex[] pq_to_array(PriorityQueue<FitnessIndex> pq)
	{
		//Iterator<FitnessIndex> fitness_iterator = pq.iterator();
		FitnessIndex[] fi = new FitnessIndex[pop_size];
		int counter = 0;
		
		while(counter<pop_size)
		{
			
			fi[counter] = pq.poll();
			if(counter < 100) {
				//System.out.println("index "+ fi[counter].index +" fitness "+fi[counter].fitness);
			}
				
			counter++;
		}
		
		return fi;
	}
	
	public void mutateAllDims(double[] indiv)
	{
		for(int dim = 0; dim < 10; dim++)
		{
			double chanse = random.nextDouble();
			double change = 0.15;
			
			if(chanse > 0.5)
			{
				if(indiv[dim]+change <5)
				{
					indiv[dim] = indiv[dim] + change;
				}
				else {
					indiv[dim] = indiv[dim] - change;
				}
			}
			else if(((indiv[dim]-change) > -5))
			{
				indiv[dim] = indiv[dim] - change;
			}else 
			{
				indiv[dim] = indiv[dim] + change;
			}
		}
		
	}
	
	public void mutate(double[] indiv)
	{
		int dim = random.nextInt(10);
		double chanse = random.nextDouble();
		double change = 0.1;
		
		if(chanse >0.5 && ((indiv[dim]+change) <5))
		{
			indiv[dim] = indiv[dim] + change;
		}
		else if(((indiv[dim]-change) > -5)){
			indiv[dim] = indiv[dim] - change;
		}
	}
	
	public void avg_CrossOver(double[] parent_a, double[] parent_b, double[] child)
	{
		for(int i = 0; i < 10; i++)
		{
			child[i] = (parent_a[i] + parent_b[i])/2;
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
	
	public int island_wheel(double[] fi, int start, int interval)
	{
		
		double[] fitnesses_interval = new double[interval];
		//gather interval consecutive fitnesses
		for(int i = 0; i < interval; i++)
		{
			fitnesses_interval[i] = fi[i+start];
		}
		
		double sumFitness_interval = sumFitness(fitnesses_interval);
		//int choosen_indiv = rouletteWheel(sumFitness_interval, fitnesses_interval);
		int choosen_indiv = partial_rouletteWheel(fitnesses_interval, 2);
		int parent = choosen_indiv+start;
		return parent;
	}
	
	
	public int island_wheel(FitnessIndex[] fi, int start, int interval)
	{
		double[] fitnesses_interval = new double[interval];
		//gather interval consecutive fitnesses
		for(int i = 0; i < interval; i++)
		{
			fitnesses_interval[i] = fi[i+start].fitness;
		}
		
		double sumFitness_interval = sumFitness(fitnesses_interval);
		int choosen_indiv = rouletteWheel(sumFitness_interval, fitnesses_interval);
		
		int parent = fi[choosen_indiv+start].index;
		return parent;
	}
	/*
	 * It will return an individual from start to start+interval proportionally to his fitness.
	 * The fitest individual in that interval has higher chanses to be picked.
	 */
	public int partial_rouletteWheel(double[] fitnesses, int interval)
	{
		double[] fitnesses_interval = new double[interval];
		int[] indexes_interval = new int[interval];
		//gather interval consecutive fitnesses
		for(int i = 0; i < interval; i++)
		{
			int choosen = random.nextInt(fitnesses.length);
			indexes_interval[i] = choosen;
			fitnesses_interval[i] = fitnesses[choosen];
		}
		
		double sumFitness_interval = sumFitness(fitnesses_interval);
		int choosen_indiv = rouletteWheel(sumFitness_interval, fitnesses_interval);
	
		int parent = indexes_interval[choosen_indiv];
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
	
	public double randomNumber(double sumFitness)
	{	
		return sumFitness * random.nextDouble();	
	}
	
	public static void main(String args[])
	{
		System.out.println("start!");
	}
	
	
}


