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
	int offsprings = (int)(pop_size * 1);
	int dimensions = 10;
	int attributes = 3;
	boolean shock =  false;
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
        double[][][] pop = new double[pop_size][3][dimensions];
        initPop(pop);
                // calculate fitness
      
        PriorityQueue<FitnessIndex> fitnesses_pq = new PriorityQueue<FitnessIndex>(new FitnessComparator());
        double sumFitness = 0;
        double childrenSumFitness = 0;
        int runs = 0 ;
       
        //int offsprings = 800;
        System.out.println("offsprings "+offsprings);
        int total_runs = (evaluations_limit_ - pop_size) / (offsprings);
        double[] fitnesses = new double[pop_size];
        double[] ranks = new double[pop_size];
        
        int evals = 0;
        while(evals < pop_size)
        {
            Double fitness = (double) evaluation_.evaluate(pop[evals][0]);
           
            FitnessIndex fitnessesIndexes = new FitnessIndex(evals, fitness);
            fitnesses_pq.add(fitnessesIndexes);
            
            fitnesses[evals] = fitness;
            ranks[evals] = fitness;
            sumFitness += fitness;
            evals++;
            // Select survivors
        }
        int afterShock = 0;
        
        
        while(runs < total_runs)
        {
        	double last_avg_children_ftiness = childrenSumFitness/offsprings;
        	//print(pop);
        	double[][][] next_gen = new double[offsprings][3][dimensions];
        	int ten_percent = total_runs/10;
        	
        	//print(pop);
        	if(runs == total_runs-1)
    		{
        		for(int i = 0; i < fitnesses.length; i++)
            	{
            		System.out.println("indiv "+i+" fitness "+fitnesses[i]);
            	}
        		
        		System.out.println();
        		
        		for(int i = 0; i < ranks.length; i++)
            	{
            		System.out.println("indiv "+i+" fitness "+ranks[i]);
            	}
    		}
        	sumFitness = sumFitness(fitnesses);
    		
    		if(runs % 10 == 0)
    		{
    			System.out.println();
        		System.out.println("run "+runs +" avg fitness : "+(sumFitness/pop_size) +" best indiv : "+fitnesses[0]);
        		System.out.println();
    		}
        	
    		
    		
        	for(int i = 0; i < offsprings; i++)
        	{
        		double[][] child = new double[3][dimensions];
        		
        		int parent_a = 0;
        		int parent_b = 0;
            	
            	boolean partial = (runs <= (ten_percent * 9) || runs % 20 ==0);
            	
            	double[] copy_fitness = new double[pop_size];
            	int[] indexes = new int[pop_size];
            	
            	for(int j  = 0; j < pop_size; j++)
            	{
            		copy_fitness[j] = fitnesses[j];
            		
            		
            		indexes[j] = j;
            		
            	}
            	
            	
            	if(partial) {
        			//Partial wheel
            		int interval = 2;
            	
            		if(offsprings < (pop_size/2))
    				{
            			int participants = ranks.length - i;
            			parent_a = partial_rouletteWheel(ranks, interval, participants);
            			int last = participants - 1;
            			ranks[parent_a] = ranks[last];
                		indexes[parent_a] = indexes[last];
                		
                		parent_b = partial_rouletteWheel(ranks, interval, participants -1);	
                		last = participants - 2;
                		ranks[parent_b] = ranks[last];
                		indexes[parent_b] = indexes[last];
                		
                		arithmeticCrossOver(pop[indexes[parent_a]], pop[indexes[parent_b]], child);
                		
//                		if(runs % 2 == 0)
//                		{
//                			arithmeticCrossOver(pop[indexes[parent_a]], pop[indexes[parent_b]], child);
//                		}
//                		else 
//                		{
//                			uniformCrossOver(pop[indexes[parent_a]], pop[indexes[parent_b]], child);
//                		}
    				}
            		else
            		{
            			parent_a = partial_rouletteWheel(ranks, interval, ranks.length);
                		parent_b = partial_rouletteWheel(ranks, interval, ranks.length);	
            			arithmeticCrossOver(pop[parent_a], pop[parent_b], child);
            		}
            		
            		if(i ==0) {
            			System.out.print(" partial ");
            		}
            		
				}
            	
            	boolean roulette = !partial;
            	   	
            	if(roulette)
				{
        			parent_a = rouletteWheel(ranks);
        			parent_b = rouletteWheel(ranks);
        			
        			if(i ==0) {
        				System.out.print(" roulette ");
            		}
        			arithmeticCrossOver(pop[parent_a], pop[parent_b], child);
				}
        	
        		
        		//System.out.println("parent a "+parent_a + " parent b "+parent_b);
        		
        		for(int attr = 0; attr < attributes; attr++)
        		{
        			for(int j = 0; j < dimensions; j++)
            		{
            			next_gen[i][attr][j] = child[attr][j];
            		}	
        		}
        		
        	}
        	
        	//MUTATE
        	
        	int individuals_to_mutate = offsprings/10;
        	int[] to_mutate = new int[individuals_to_mutate];
        	
        	for(int i = 0; i <individuals_to_mutate; i++)
        	{
        		//to_mutate[i] = random.nextInt(100);
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
               
            }
            
            int children = 0;
            childrenSumFitness = 0;
    		while(children < offsprings)
            {
                Double fitness = (double) evaluation_.evaluate(next_gen[children][0]);
               //System.out.println("fitnesses[evals] "+fitness);
                FitnessIndex fitnessesIndexes = new FitnessIndex(children + pop_size, fitness);
                fitnesses_pq.add(fitnessesIndexes);
                childrenSumFitness += fitness;
                children++;
               
            }
   
        	FitnessIndex[] fi = new FitnessIndex[pop_size+offsprings];
        	
        	fi = pq_to_array(fitnesses_pq, offsprings);
        	
        	if(runs == total_runs-1) {
        		for(int i = 0; i < fi.length; i++)
            	{
            		System.out.println("indiv "+i+" fi "+fi[i].fitness);
            	}
        		System.out.println();
        		
        		for(int attr = 1; attr < attributes; attr++)
        		{
        			System.out.println();
        			for(int dim = 0; dim < 10; dim++)
        			{
        				System.out.println("pop[0][+"+attr+"]"+"["+ dim+"] : "+pop[0][attr][dim]);
        			}
        			System.out.println();
        		}
        	}
        	// Select survivors	
        	//stochastic_transfer_gen(pop, next_gen, fi, fitnesses,ranks);
        	//transfer_gen(pop, next_gen, fi, fitnesses, ranks);
        	if(runs % 10 == 0)
        	{
        		stochastic_transfer_gen(pop, next_gen, fi, fitnesses, ranks);
        	}
        	else
        	{
        		transfer_gen(pop, next_gen, fi, fitnesses, ranks);
        	}
        	
        	double epsilon = 0;
        	if(last_avg_children_ftiness > 1)
        	{
        		epsilon = 0.001;
        				
        	}
        	else if((last_avg_children_ftiness * 10) > 1)
        	{
        		epsilon = 0.0001;
        	}
        	else if((last_avg_children_ftiness * 100) > 1)
        	{
        		epsilon = 0.0001;
        	}else if(last_avg_children_ftiness * 1000 > 1)
        	{
        		epsilon = 0.00001;
        	}else if(last_avg_children_ftiness * 10000 >1)
        	{
        		epsilon = 0.000001;
        	}else if(last_avg_children_ftiness * 100000 >1)
        	{
        		epsilon = 0.0000001;
        	}
        	if(shock)
        	{
        		if(afterShock == 0)
        		{
        			antishock(pop);
        		}
        		
        		//offsprings = (int)(pop_size *0.3);
        		afterShock++;
        		if(afterShock > 30)
        		{
        			shock = false;
        			afterShock = 0;
        		}
        	}
        	
        	double avg_children_fitness = childrenSumFitness/offsprings;
    		if((Math.abs(last_avg_children_ftiness - avg_children_fitness) < epsilon) && !shock)
    		{
    			if(last_avg_children_ftiness > 1)
    			{
    				initPop(pop);
    				System.out.println("Pop Init -------------------------------------------------------------------------------------------------");
    				evals = 0;
			        while(evals < pop_size)
			        {
			            Double fitness = (double) evaluation_.evaluate(pop[evals][0]);
			           
			            FitnessIndex fitnessesIndexes = new FitnessIndex(evals, fitness);
			            fitnesses_pq.add(fitnessesIndexes);
			            
			            fitnesses[evals] = fitness;
			            ranks[evals] = fitness;
			            sumFitness += fitness;
			            evals++;
			            // Select survivors
			        }
			        
    			}else
    			{
    				shock = true;
        			//offsprings = pop_size;
        			shock(pop);
        			
        			System.out.println(" POPULATION HAVE BEEN SHOOOCKED !!! ");
    			}
    			
    			
    		}
        	
//        	if(runs % 32545252 ==0)
//        	{
//        		System.out.println();
//            	System.out.println("fitness of best indiv"+fi[0].fitness);
//            	
//            	for(int dim = 0; dim < 10; dim++)
//            	{
//            		int index = fi[0].index;
//            		if(index >= pop_size)
//            		{
//            			System.out.print("dim "+dim+" "+pop[index - pop_size][0][dim]);
//            		}
//            		else
//            		{
//            			System.out.print("dim "+dim+" "+pop[index][0][dim]);
//            		}
//            		
//            	}
//            	System.out.println();
//        	}
        	runs++;
        }  
	}
	
	private void shockIndiv(double[][] indiv)
	{
		for(int dim = 0; dim < 10; dim++)
    	{
			double chanse = random.nextDouble() * 2 - 1;
    		double sign = Math.signum(chanse);
    		double shock = random.nextDouble();
        	
    		indiv[0][dim] = indiv[0][dim] + sign * shock;
    	}
		
		for(int mutate_rate = 0; mutate_rate < 10; mutate_rate++)
    	{
    		//from 0.25 to 0.75
    		double mu_rate = random.nextDouble() * 0.5 + 0.25;
    		indiv[1][mutate_rate] = mu_rate;
    	}
    	
    	for(int mutate_size = 0; mutate_size < 10; mutate_size++)
    	{
    		//from -0.5 to 0.5
    		double chanse = random.nextDouble() * 1 - 0.5;
    		double sign = Math.signum(chanse);
    		double mu_size = random.nextDouble() * 0.5 + 0.1;
    		indiv[2][mutate_size] = sign * mu_size;
    	}
	}
	
	private void antishock(double[][][] pop) {
		// TODO Auto-generated method stub
		for(int individual = 0; individual < pop_size; individual++)
        {
//			
			for(int mutate_rate = 0; mutate_rate < 10; mutate_rate++)
        	{
        		//from 0.25 to 0.75
        		double mu_rate = random.nextDouble() * 0.2 + 0.1;
        		pop[individual][1][mutate_rate] = mu_rate;
        	}
        	
        	for(int mutate_size = 0; mutate_size < 10; mutate_size++)
        	{
        		//from -0.5 to 0.5
        		double chanse = random.nextDouble() * 1 - 0.5;
        		double sign = Math.signum(chanse);
        		double mu_size = random.nextDouble() * 0.2 + 0.1;
        		pop[individual][2][mutate_size] = sign * mu_size;
        	}
        }
	}
	
	private void shock(double[][][] pop) {
		// TODO Auto-generated method stub
		for(int individual = 0; individual < pop_size; individual++)
        {
			if(random.nextDouble() < 0.2) {
				for(int dim = 0; dim < 10; dim++)
	        	{
					
	        		double sign = getSign();
	        		double shock = random.nextDouble()*3;
	            	
	        		if(sign <0)
	        		{
	        			shock = Math.max(-5, shock);
	        		}
	        		else {
	        			shock = Math.min(5, shock);
	        		}
	        		pop[individual][0][dim] = pop[individual][0][dim] + sign * shock;
	        	}
			}
			
//			
			for(int mutate_rate = 0; mutate_rate < 10; mutate_rate++)
        	{
        		//from 0.25 to 0.75
        		double mu_rate = random.nextDouble() * 0.75 + 0.35;
        		pop[individual][1][mutate_rate] = mu_rate;
        	}
        	
        	for(int mutate_size = 0; mutate_size < 10; mutate_size++)
        	{
        		//from -0.5 to 0.5
        		double sign = getSign();
        		double mu_size = random.nextDouble() * 2 + 0.3;
        		pop[individual][2][mutate_size] = sign * mu_size;
        	}
        }
	}

	private void initPop(double[][][] pop)
	{
		for(int individual = 0; individual < pop_size; individual++)
        {
        	pop[individual] = new double[3][dimensions];
        	
        	for(int dim = 0; dim < 10; dim++)
        	{
        		double rangeMax = 10;
        		double rangeMin = -5;

            	double randomValue = (rangeMax * random.nextDouble()) + rangeMin;
            	
        		pop[individual][0][dim] = randomValue;
        	}
        	
        	for(int mutate_rate = 0; mutate_rate < 10; mutate_rate++)
        	{
        		//from 0.1 to 0.5
        		double mu_rate = random.nextDouble() * 0.4 + 0.1;
        		//double mu_rate = 0.2;
        		pop[individual][1][mutate_rate] = mu_rate;
        	}
        	
        	for(int mutate_size = 0; mutate_size < 10; mutate_size++)
        	{
        		//from -0.5 to 0.5
        		
        		double sign = getSign();
        		
        		double mu_size = random.nextDouble() * 0.5 + 0.05;
        		//double mu_size = 0.15;
        		pop[individual][2][mutate_size] = sign * mu_size;
        	}
        }
	}
	private double getSign()
	{
		double chanse = random.nextDouble() * 1 - 0.5;
		double sign = Math.signum(chanse);
		
		return sign;
	}
	
	private void print(double[][][] pop) {
		// TODO Auto-generated method stub
		for(int i =0; i< pop.length; i++)
		{
		    System.out.println();
		    System.out.println("indiv "+i+" dims: ");
			for(int dim = 0; dim < dimensions; dim++) {
				System.out.print(pop[i][0][dim] + " , ");
			}
			System.out.println();
		}
	}

	public void stochastic_transfer_gen(double[][][] pop, double[][][] next_gen,
								FitnessIndex[] fi, double[] fitnesses, double[] ranks)
	{	
		double[][][] sorted = new double[pop_size][3][dimensions];
		for(int i = 0; i < pop_size; i++)
    	{
			int indiv = 0;
			// elitism
			int elites = 30;
			if(i < elites) {
				indiv = i;
			}
			else if(random.nextDouble() > 0.98)
			{
				//5% of indivs are choosen by chanse
				indiv = random.nextInt(fi.length-i);
			}
			else {
				indiv = rouletteWheel(fi, fi.length-i);
			}
			ranks[i] =  (int) ((Math.abs(fi[indiv].getIndex() - fi.length + offsprings)+100)/100);
			fitnesses[i] = fi[indiv].fitness;
			int index = fi[indiv].getIndex();
			fi[indiv] = fi[fi.length-i-1];
			
			for(int attr = 0; attr < attributes; attr++)
			{
				for(int j = 0; j < dimensions; j++)
	    		{
	    			if(index >= pop_size)
	    			{
	    				sorted[i][attr][j] = next_gen[index - pop_size][attr][j];
	    			}
	    			else
	    			{
	    				sorted[i][attr][j] = pop[index][attr][j];
	    			}	
	    		}	
			}
    	}
		
		for(int ind = 0; ind < pop_size; ind++)
    	{
    		for(int attr = 0; attr < attributes; attr++)
    		{
    			for(int j = 0; j < dimensions; j++) {
        			pop[ind][attr][j] = sorted[ind][attr][j];
        		}
    		}
    		
    	}
	}
	
	public void transfer_gen(double[][][] pop, double[][][] next_gen,
			FitnessIndex[] fi, double[] fitnesses, double[] ranks)
	{
		double[][][] sorted = new double[pop_size][3][dimensions];
		for(int i = 0; i < pop_size; i++)
    	{
			fitnesses[i] = fi[i].fitness;
			ranks[i] = (int) ((Math.abs(fi[i].getIndex() - fi.length + offsprings)+100)/100);
			int index = fi[i].getIndex();
			
			for(int attr = 0; attr < attributes; attr++)
			{
				for(int j = 0; j < dimensions; j++)
	    		{
	    			if(index >= pop_size)
	    			{
	    				sorted[i][attr][j] = next_gen[index - pop_size][attr][j];
	    			}
	    			else
	    			{
	    				sorted[i][attr][j] = pop[index][attr][j];
	    			}	
	    		}	
			}
    	}
		
		for(int ind = 0; ind < pop_size; ind++)
    	{
    		for(int attr = 0; attr < attributes; attr++)
    		{
    			for(int j = 0; j < dimensions; j++) {
        			pop[ind][attr][j] = sorted[ind][attr][j];
        		}
    		}
    		
    	}
	}
	
	public FitnessIndex[] pq_to_array(PriorityQueue<FitnessIndex> pq, int offsprings)
	{
		//Iterator<FitnessIndex> fitness_iterator = pq.iterator();
		FitnessIndex[] fi = new FitnessIndex[pop_size+offsprings];
		int counter = 0;
		
		while(counter<(pop_size+offsprings))
		{
			
			fi[counter] = pq.poll();
//			if(counter < 100) {
//				//System.out.println("index "+ fi[counter].index +" fitness "+fi[counter].fitness);
//			}
				
			counter++;
		}
		
		return fi;
	}
	
	public void mutateSomeDims(double[] indiv, double chanse2)
	{
		for(int dim = 0; dim < dimensions; dim++)
		{
			double chanse_to_mutate = random.nextDouble();
			if(chanse2 > chanse_to_mutate) {
				double chanse = random.nextDouble();
				
				double change = random.nextDouble() * 0.2 + 0.05;
				
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
				}
				else 
				{
					indiv[dim] = indiv[dim] + change;
				}
			}
			
		}
	}
	
	public void mutateAllDims(double[][] indiv)
	{
		for(int dim = 0; dim < dimensions; dim++)
		{
			double chanse = random.nextDouble();
			double change = 0.15;
			
			if(chanse > 0.5)
			{
				if(indiv[0][dim]+change <5)
				{
					indiv[0][dim] = indiv[0][dim] + change;
				}
				else {
					indiv[0][dim] = indiv[0][dim] - change;
				}
			}
			else if(((indiv[0][dim]-change) > -5))
			{
				indiv[0][dim] = indiv[0][dim] - change;
			}else 
			{
				indiv[0][dim] = indiv[0][dim] + change;
			}
		}
		
	}
	
	public void mutate(double[] indiv)
	{
		int dim = random.nextInt(dimensions);
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
	
	public void arithmeticCrossOver(double[][] parent_a, double[][] parent_b, double[][] child)
	{
		// cross over
		for(int attribute = 0; attribute < attributes ; attribute++)
		{
			
			for(int i = 0; i < dimensions; i++)
			{
				double a = 0;
				double b = 0;
				double increment = 0;
				
				double max = Math.max(parent_a[attribute][i], parent_b[attribute][i]);
				double min = Math.min(parent_a[attribute][i], parent_b[attribute][i]);
				
				if(attribute == 0)
				{
					a = Math.min(5, max + increment);
					b = Math.max(-5, min - increment);
				}
				if(attribute == 1)
				{
					increment = 0.2;
					a = Math.min(1, max + increment);
					b = Math.max(0, min - increment);
				}
				else if(attribute == 2)
				{
					//increment = 0.2;
					a = max;
					b = min;
				}
				
				//double distance = random.nextDouble();
				double distance = 0.5;
				child[attribute][i] = (distance * a + (1 - distance) * b);
			}
		}
		
		//mutate forcebly at list 0.05 or -0.05
		//mutateEmbeded(child);
		double shock = random.nextDouble();
		if(shock > 0.99) {
			shockIndiv(child);
		}
	}
	
	public void uniformCrossOver(double[][] parent_a, double[][] parent_b, double[][] child)
	{
		for(int attribute = 0; attribute < attributes ; attribute++)
		{
			for(int i = 0; i < dimensions; i++)
			{
				double chanse = random.nextDouble();
				if(chanse > 0.5 )
				{
					child[attribute][i] = parent_a[attribute][i];
				}
				else {
					child[attribute][i] = parent_b[attribute][i];
				}
			}		
		}
		
		//mutateEmbeded(child);
		double shock = random.nextDouble();
		if(shock > 0.99) {
			shockIndiv(child);
		}
	}
	
	public void mutateEmbeded(double[][] child)
	{
		//mutate forcebly at list 0.05 or -0.05
		for(int i = 0; i < dimensions; i++)
		{
			double chanse = random.nextDouble();
			
			if(chanse < child[1][i])
			{
				double mutation_size = child[2][i];
				double val = child[0][i];
				
				if(mutation_size < 0)
				{
					double min_size = -0.05;
					mutation_size = Math.min(min_size, mutation_size);
					child[0][i] = Math.max(val + mutation_size, -5);
				}
				else 
				{
					double min_size = 0.05;
					mutation_size = Math.max(min_size, mutation_size);
					child[0][i] = Math.min(val + mutation_size, 5);
				}
				
			}
		}
	}
	
	public double[] crossOver(double[] parent_a, double[] parent_b)
	{
		double child[] = new double[dimensions];
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
	
//	public int island_wheel(double[] fi, int start, int interval)
//	{
//		
//		double[] fitnesses_interval = new double[interval];
//		//gather interval consecutive fitnesses
//		for(int i = 0; i < interval; i++)
//		{
//			fitnesses_interval[i] = fi[i+start];
//		}
//		
//		double sumFitness_interval = sumFitness(fitnesses_interval);
//		//int choosen_indiv = rouletteWheel(sumFitness_interval, fitnesses_interval);
//		int choosen_indiv = partial_rouletteWheel(fitnesses_interval, 2);
//		int parent = choosen_indiv+start;
//		return parent;
//	}
	
	
	public int island_wheel(FitnessIndex[] fi, int start, int interval)
	{
		double[] fitnesses_interval = new double[interval];
		//gather interval consecutive fitnesses
		for(int i = 0; i < interval; i++)
		{
			fitnesses_interval[i] = fi[i+start].fitness;
		}
		
		int choosen_indiv = rouletteWheel(fitnesses_interval);
		
		int parent = fi[choosen_indiv+start].index;
		return parent;
	}
	/*
	 * It will return an individual from start to start+interval proportionally to his fitness.
	 * The fitest individual in that interval has higher chanses to be picked.
	 */
	public int partial_rouletteWheel(double[] fitnesses, int interval, int participants)
	{
		double[] fitnesses_interval = new double[interval];
		int[] indexes_interval = new int[interval];
		//gather interval consecutive fitnesses
		for(int i = 0; i < interval; i++)
		{
			int choosen = random.nextInt(participants);
			
			indexes_interval[i] = choosen;
			fitnesses_interval[i] = fitnesses[choosen];
		}
		
		int choosen_indiv = rouletteWheel(fitnesses_interval);
	
		int parent = indexes_interval[choosen_indiv];
		return parent;
	}
	
	public int rouletteWheel(FitnessIndex[] fi, int participants)
	{
		
		double sum = 0;
		
		for(int j = 0; j < participants; j++)
		{
			sum += fi[j].fitness;
		}
		
		double rand = randomNumber(sum);
		
		for(int i = 0; i < participants; i++)
		{
			rand -= fi[i].fitness;
			
			if(rand < 0)
			{
				return i;
			}
		}
		System.out.println("error");
		return random.nextInt(fi.length);
	}
	
	public int rouletteWheel(double[] fitnesses)
	{
		double sumFitness = sumFitness(fitnesses);
		
		double rand = randomNumber(sumFitness);
		
		for(int i = 0; i < fitnesses.length; i++)
		{
			rand -= fitnesses[i];
			
			if(rand < 0)
			{
				return i;
			}
		}
		System.out.println("error----------------------------------------------------------------------------------------------------------------------------------");
		return random.nextInt(fitnesses.length);
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


