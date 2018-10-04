import java.util.Comparator;

public class IndividualComparator implements Comparator<Individual> {

	@Override
	public int compare(Individual individual_1, Individual individual_2) {

		if(individual_1.fitness < individual_2.fitness)
		{
			return 1;
		}
		else if(individual_1.fitness > individual_2.fitness)
		{
			return -1;
		}
		
		return 0;
	}
}
