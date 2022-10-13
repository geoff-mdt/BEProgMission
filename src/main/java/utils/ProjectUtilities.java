package utils;

import fr.cnes.sirius.patrius.events.Phenomenon;
import fr.cnes.sirius.patrius.events.postprocessing.Timeline;
import reader.Site;

/**
 * [DO NOT MODIFY THIS CLASS] Utility class.
 * 
 * @author herberl
 */
public class ProjectUtilities {

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * 
	 * Print the input {@link Timeline} in the console.
	 * 
	 * @param timeline Input {@link Timeline} to print.
	 */
	public static void printTimeline(final Timeline timeline) {
		System.out.println("____ Printing Timeline ____");
		for (final Phenomenon phenom : timeline.getPhenomenaList()) {
			System.out.println(phenom);
		}
		System.out.println("___________________________");
	}

	// MODIFY FOR PERSONAL USE
	public static void printTimeline(final Timeline timeline, final Site site) {
		System.out.println("____ Printing Timeline " + site.getName() + " ____");
		for (final Phenomenon phenom : timeline.getPhenomenaList()) {
			System.out.println(phenom);
		}
		System.out.println("___________________________");
	}
}
