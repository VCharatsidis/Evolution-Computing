public class Individual {
	
	public int dimensions = 10;
	public double[] genome = new double[dimensions];
	public double[] mutation_steps = new double[dimensions]; 
	public double fitness;
	public int rank;
	public int index;
	
	public Individual(Individual to_copy)
	{
		genome = new double[dimensions];
		mutation_steps = new double[dimensions];
		
		for(int i = 0; i < dimensions; i++)
		{
			genome[i] = to_copy.genome[i];
			mutation_steps[i] = to_copy.mutation_steps[i];	
		}
		
		fitness = to_copy.fitness;
		rank = to_copy.rank;
		index = to_copy.index;
	}
	
	public Individual()
	{
		genome = new double[dimensions];
		mutation_steps = new double[dimensions];
		fitness = 0;
		rank = 0;
		index = 0;
	}
	
	/*
	 * Constructs an individual with random values as genome and mutation_steps inside the bounds.
	 */
	public Individual(double genome_upper_value, double genome_lower_value, double upper_mutation_step, double lower_mutation_step)
	{
		dimensions = 10;
		genome = new double[dimensions];
		mutation_steps = new double[dimensions];
		
		for(int dim = 0; dim < dimensions; dim++)
    	{
			genome[dim] = Utils.getSign() * Utils.double_in_range(genome_upper_value, genome_lower_value);
			mutation_steps[dim] = Utils.double_in_range(upper_mutation_step, lower_mutation_step);
    	}
		
		fitness = 0;
		rank = 0;
		index = 0;
	}
	
}
