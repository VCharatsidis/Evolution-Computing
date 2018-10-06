
public class ArithmeticCrossOver implements CrossOver {

	public Individual cross_over(Individual a, Individual b) {

		int dimensions = player72.dimensions;	
		Individual child = new Individual();
	
		// The distance can be fixed or random.
		double distance = 0.5;
		// Increment can be used with distance to be able to make genomes that are not always in the middle between parents genome.
		double genome_increment = 0;
		// Likewise for steps.
		double mutation_increment = 0.15;
		for(int dim = 0; dim < dimensions; dim++)
		{
			child.genome[dim] = fill_genome(a.genome[dim], b.genome[dim], child, distance, genome_increment);
			child.mutation_steps[dim] = fill_mutation_steps(a.mutation_steps[dim], b.mutation_steps[dim], child, mutation_increment);
		}	
		
		return child;
	}
	
	/*
	 * Combine parents genome to make the genome of the child.
	 */
	public double fill_genome(double parent_a_genome, double parent_b_genome, Individual child, double distance, double increment)
	{
		// We need to know the max and min values so that we add the increment to the max and we extract increment to the min.
		double max_genome_value = Math.max(parent_a_genome, parent_b_genome);
		double min_genome_value = Math.min(parent_a_genome, parent_b_genome);
		
		// This assures that they dont exceed 5 or -5 , the fixed bounds.
		parent_a_genome = Math.min(5, max_genome_value + increment);
		parent_b_genome = Math.max(-5, min_genome_value - increment);
		
		return (distance * parent_a_genome + (1 - distance) * parent_b_genome);
	}
	
	/*
	 * Combines the mutation step of the parents to make the mutation step of the child.
	 */
	public double fill_mutation_steps(double parent_a_steps, double parent_b_steps, Individual child, double increment)
	{
		// We need to know the max and min values so that we add the increment to the max and we extract increment to the min.
		double max_step = Math.max(parent_a_steps, parent_b_steps);
		double min_step = Math.min(parent_a_steps, parent_b_steps);
		
		// This reassures that the steps are between the bounds.
		double high_bound = 0.25;
		double low_bound  = 0;
		
		parent_a_steps = Math.min(high_bound, max_step + increment);
		parent_b_steps = Math.max(low_bound, min_step - increment);
		
		// The distance is fixed.
		double distance = player72.rnd_.nextDouble();
		return (distance * parent_a_steps + (1 - distance) * parent_b_steps);
	}
	
}
