import java.util.Random;

public class UniformCrossOver {
	
	private int dimensions = player72.dimensions;	
	
	public Individual cross_over(Individual a, Individual b)
	{
		Individual child = new Individual();
		
		for(int dim = 0; dim < dimensions; dim++)
		{
			double chanse = player72.rnd_.nextDouble();
			
			if(chanse > 0.5 ) 
			{
				child.genome[dim] = a.genome[dim];
				child.mutation_steps[dim] = a.mutation_steps[dim];
			}
			else
			{
				child.genome[dim] = b.genome[dim];
				child.mutation_steps[dim] = b.mutation_steps[dim];
			}
			
		}	
		
		return child;
	}
}
