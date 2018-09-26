import java.util.Comparator;

public class FitnessComparator implements Comparator<FitnessIndex> {

	@Override
	public int compare(FitnessIndex individual_1, FitnessIndex individual_2) {
		// TODO Auto-generated method stub
		if(individual_1.fitness < individual_2.fitness)
		{
			return 1;
		}else if(individual_1.fitness > individual_2.fitness)
		{
			return -1;
		}
		return 0;
	}
	
}
