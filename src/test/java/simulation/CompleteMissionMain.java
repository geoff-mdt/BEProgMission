package simulation;

import java.util.Map;

import fr.cnes.sirius.patrius.attitudes.AttitudeLawLeg;
import fr.cnes.sirius.patrius.attitudes.AttitudeLeg;
import fr.cnes.sirius.patrius.attitudes.StrictAttitudeLegsSequence;
import fr.cnes.sirius.patrius.events.postprocessing.Timeline;
import fr.cnes.sirius.patrius.utils.exception.PatriusException;
import progmission.CompleteMission;
import reader.Site;

public class CompleteMissionMain {

	/**
	 * Main method to run your simulation.
	 * 
	 * @param args No arg
	 * @throws PatriusException If a {@link PatriusException} occurs.
	 */
	public static void main(String[] args) throws PatriusException {
		System.out.println("##################################################");
		double t0 = System.currentTimeMillis();

		// Instantiating our mission using the CompleteMission object.
		final CompleteMission mission = new CompleteMission("BE Supaero mission", 100);
		System.out.println("Complete simulation starting ...");
		System.out.println(mission);

		// First step is to compute when the satellite can access the targets. Each
		// access is an observation opportunity to be considered in the latter scheduling
		// process.
		Map<Site, Timeline> accessPlan = mission.computeAccessPlan();
		// System.out.println(accessPlan.toString());

		/// Then we compute the observation plan, that is to say we fill a plan with
		// Observation objects that can be achieved one after each other by the
		// satellite without breaking the cinematic constraints imposed by the
		// satellite agility.
		Map<Site, AttitudeLawLeg> observationPlan = mission.computeObservationPlan();
		System.out.println(observationPlan.toString());

		// Then, we compute the cinematic plan, which is the whole cinematic sequence of
		// attitude law legs for our satellite during the mission horizon
		StrictAttitudeLegsSequence<AttitudeLeg> cinematicPlan = mission.computeCinematicPlan();
		System.out.println(cinematicPlan.toPrettyString());

		// Checking our cinematic plan
		boolean validity = mission.checkCinematicPlan();
		System.out.println("Plan validity : " + validity);

		// Only if the cinematic plan is valid : we compute the score of our
		// observationPlan
		System.out.println(mission.computeFinalScore());

		//Finally, we write the VTS outputs to visualize and validate our plan
		mission.generateVTSVisualization(cinematicPlan);


		System.out.println("\n\nSimulation done");

		// Computing the time of execution
		double t1 = System.currentTimeMillis();
		System.out.println("Total duration : " + 0.001 * (t1 - t0));

		System.out.println("##################################################");

		// To output files for AcessTimeline Visualization
		// mission.createSimpleVTSVisualization();
	}

}
