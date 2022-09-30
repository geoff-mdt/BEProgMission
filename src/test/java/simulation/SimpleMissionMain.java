package simulation;

import fr.cnes.sirius.patrius.utils.exception.PatriusException;
import progmission.SimpleMission;

/**
 * Main class for launching the simulation.
 *
 * @author herberl
 */
public class SimpleMissionMain {

	/**
	 * Main method.
	 * 
	 * @param args Nothing expected
	 * @throws PatriusException if a PatriusException occurs
	 */
	public static void main(final String[] args) throws PatriusException {

		System.out.println("##################################################");

		// Instantiating our mission using the SimpleMission object.
		final SimpleMission mission = new SimpleMission("BE Supaero mission", 20);
		System.out.println("Simple simulation starting ...");
		System.out.println(mission);

		// Creating simple VTS visualization by propagating our satellite with a simple
		// nadir pointing law
		mission.createSimpleVTSVisualization();

		System.out.println("\n\nSimulation done");
		System.out.println("##################################################");

	}
}
