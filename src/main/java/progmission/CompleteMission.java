package progmission;

import java.io.File;
import java.util.*;
import java.util.Map.Entry;

import fr.cnes.sirius.patrius.assembly.models.SensorModel;
import fr.cnes.sirius.patrius.attitudes.*;
import fr.cnes.sirius.patrius.events.CodedEvent;
import fr.cnes.sirius.patrius.events.CodedEventsLogger;
import fr.cnes.sirius.patrius.events.GenericCodingEventDetector;
import fr.cnes.sirius.patrius.events.Phenomenon;
import fr.cnes.sirius.patrius.events.postprocessing.AndCriterion;
import fr.cnes.sirius.patrius.events.postprocessing.ElementTypeFilter;
import fr.cnes.sirius.patrius.events.postprocessing.PhenomenonDurationFilter;
import fr.cnes.sirius.patrius.events.postprocessing.Timeline;
import fr.cnes.sirius.patrius.events.sensor.SensorVisibilityDetector;
import fr.cnes.sirius.patrius.frames.FramesFactory;
import fr.cnes.sirius.patrius.frames.TopocentricFrame;
import fr.cnes.sirius.patrius.math.geometry.euclidean.threed.Vector3D;
import fr.cnes.sirius.patrius.math.util.FastMath;
import fr.cnes.sirius.patrius.orbits.Orbit;
import fr.cnes.sirius.patrius.orbits.pvcoordinates.PVCoordinatesProvider;
import fr.cnes.sirius.patrius.propagation.BoundedPropagator;
import fr.cnes.sirius.patrius.propagation.Propagator;
import fr.cnes.sirius.patrius.propagation.analytical.KeplerianPropagator;
import fr.cnes.sirius.patrius.propagation.events.ConstantRadiusProvider;
import fr.cnes.sirius.patrius.propagation.events.EventDetector;
import fr.cnes.sirius.patrius.propagation.events.ThreeBodiesAngleDetector;
import fr.cnes.sirius.patrius.time.AbsoluteDate;
import fr.cnes.sirius.patrius.time.AbsoluteDateInterval;
import fr.cnes.sirius.patrius.utils.exception.PatriusException;
import fr.cnes.sirius.patrius.utils.exception.PropagationException;
import reader.Site;
import utils.ConstantsBE;
import utils.ProjectUtilities;
import utils.VTSTools;

/**
 * This class encapsulates the context of an Earth Observation mission.
 *
 * @author herberl
 */
public class CompleteMission extends SimpleMission {

	/**
	 * Maximum checking interval (s) for the event detection during the orbit
	 * propagation.
	 */
	public static final double MAXCHECK_EVENTS = 120.0;

	/**
	 * Default convergence threshold (s) for the event computation during the orbit
	 * propagation.
	 */
	public static final double TRESHOLD_EVENTS = 1.e-4;

	/**
	 * This {@link Map} will be used to enumerate each site access {@link Timeline},
	 * that is to say a {@link Timeline} with access windows respecting all
	 * observation conditions. This object corresponds to the access plan, which
	 * will be computed in the computeAccessPlan() method.
	 */
	private final Map<Site, Timeline> accessPlan;

	/**
	 * This {@link Map} will be used to enumerate each site's programmed
	 * observation. We suggest to use an {@link AttitudeLawLeg} to encapsulate the
	 * guidance law of each observation. This object corresponds to the observation
	 * plan, which will be computed in the computeObservationPlan() method.
	 */
	private final Map<Site, AttitudeLawLeg> observationPlan;

	/**
	 * {@link StrictAttitudeLegsSequence} representing the cinematic plan during the
	 * whole mission horizon. Each {@link AttitudeLeg} corresponds to a diffrent
	 * attitude law : either nadir pointing, target pointing or a slew between two
	 * laws. This object corresponds to the cinematic plan, which will be computed
	 * in the computeCinematicPlan() method.
	 */
	private final StrictAttitudeLegsSequence<AttitudeLeg> cinematicPlan;

	/**
	 * Constructor for the {@link CompleteMission} class.
	 *
	 * @param missionName   Name of the mission
	 * @param numberOfSites Number of target {@link Site} to consider, please give a
	 *                      number between 1 and 100.
	 * @throws PatriusException      Can be raised by Patrius when building
	 *                               particular objects. Here it comes from
	 *                               {@link FramesFactory}
	 * @throws IllegalStateException if the mission horizon is too short
	 */
	public CompleteMission(final String missionName, int numberOfSites) throws PatriusException {

		// Since this class extends the SimpleMission class, we need to use the super
		// constructor to instantiate our instance of CompleteMission. All the
		// attributes of the super class will be instantiated during this step.
		super(missionName, numberOfSites);

		// Initialize the mission plans with empty maps. You will fill those HashMaps in
		// the "compute****Plan()" methods.
		this.accessPlan = new HashMap<>();
		this.observationPlan = new HashMap<>();
		this.cinematicPlan = new StrictAttitudeLegsSequence<>();

	}

	/**
	 * Compute the access plan.
	 * <p>
	 * Reminder : the access plan corresponds to the object gathering all the
	 * opportunities of access for all the sites of interest during the mission
	 * horizon. One opportunity of access is defined by an access window (an
	 * interval of time during which the satellite can observe the target and during
	 * which all the observation conditions are achieved : visibility, incidence
	 * angle, illumination of the scene,etc.). Here, we suggest you use the Patrius
	 * class {@link Timeline} to encapsulate all the access windows of each site of
	 * interest. One access window will then be described by the {@link Phenomenon}
	 * object, itself defined by two {@link CodedEvent} objects giving the start and
	 * end of the access window. Please find more tips and help in the submethods of
	 * this method.
	 *
	 * @return the sites access plan with one {@link Timeline} per {@link Site}
	 * @throws PatriusException If a {@link PatriusException} occurs during the
	 *                          computations
	 */
	public Map<Site, Timeline> computeAccessPlan() throws PatriusException {
		/**
		 * Here you need to compute one access Timeline per target Site. You can start
		 * with only one site and then try to compute all of them.
		 *
		 * Note : when computing all the sites, try to make sure you don't decrease the
		 * performance of the code too much. You might have some modifications to do in
		 * order to ensure a reasonable time of execution.
		 */

		// We are creating a site access timeline for each element in the mission siteList
		// and add each to the <Site,SiteAccessTimeline> HashMap
		for (Site targetSite : this.getSiteList()) {
			Timeline siteAccessTimeline = createSiteAccessTimeline(targetSite);
			this.accessPlan.put(targetSite, siteAccessTimeline);
			ProjectUtilities.printTimeline(siteAccessTimeline, targetSite);
		}

		/*Site targetSite = this.getSiteList().get(0);
		Timeline siteAccessTimeline = createSiteAccessTimeline(targetSite);
		accessPlan.put(targetSite, siteAccessTimeline);
		ProjectUtilities.printTimeline(siteAccessTimeline);*/

		return this.accessPlan;
	}

	/**
	 * Compute the observation plan.
	 * <p>
	 * Reminder : the observation plan corresponds to the sequence of observations
	 * programmed for the satellite during the mission horizon. Each observation is
	 * defined by an observation window (start date; end date defining an
	 * {@link AbsoluteDateInterval}), a target (target {@link Site}) and an
	 * {@link AttitudeLawLeg} giving the attitude guidance to observe the target.
	 *
	 * @return the sites observation plan with one {@link AttitudeLawLeg} per
	 * {@link Site}
	 * @throws PatriusException If a {@link PatriusException} occurs during the
	 *                          computations
	 */
	public Map<Site, AttitudeLawLeg> computeObservationPlan() throws PatriusException {

		List<Reservation> listResas = new ArrayList<>();

		// Sorting sites by score to maximize the total score
		this.getSiteList().sort(Comparator.comparing(Site::getScore).reversed());

		// Iterating through each site and getting its timeline
		for(Site target: this.getSiteList()){

			// Getting its access Timeline
			final Timeline timeline = this.accessPlan.get(target);
			boolean targetObserved = false;

			//
			for (final Phenomenon accessWindow : timeline.getPhenomenaList()) {

				AbsoluteDate access_start = accessWindow.getStartingEvent().getDate();
				AbsoluteDate access_end   = accessWindow.getEndingEvent().getDate();

				List<Reservation> listParallelsResas = new ArrayList<>();

				if(listResas.size() == 0) {
					listResas.add(new Reservation(access_start, target, getSatellite().getMaxSlewDuration()));
					break;
				} else{
					// Building an chronologically sorted list of reservation of interest, for a given access
					// Theses parallels reservations are those who overlap with the access
					for(Reservation resa: listResas){
						// If there is a reservation already booked for a target
						if(resa.getSite().equals(target)){
							targetObserved = true;
							break;
						}

						if((resa.getStartDate().compareTo(access_start) > 0
							&& resa.getStartDate().compareTo(access_end) < 0)
							||
							(resa.getEndDate().compareTo(access_end) < 0)
							&& (resa.getEndDate().compareTo(access_start) > 0)
							|| (resa.getStartDate().compareTo(access_start) < 0)
							&& (resa.getEndDate().compareTo(access_end) > 0)){
							listParallelsResas.add(resa);
						}
					}
					// We break the for a second time if we observed the target to escape
					if(targetObserved) { break; }
					listParallelsResas.sort(Comparator.comparing(Reservation::getStartDate));

					// Base case, if there's none, the beginnig of the access is set for the start of observation time
					if(listParallelsResas.size() == 0){
						listResas.add(new Reservation(access_start, target, getSatellite().getMaxSlewDuration()));
						break;
					} else {
						for (int i = 0; i < listParallelsResas.size(); i++) {

							// For 1 reservation overlapping the access duration, we try to put before or immediately after the next reservation
							if(listParallelsResas.size() == 1){
								// Before
								if(listParallelsResas.get(0).getStartDate().compareTo(access_start.shiftedBy(ConstantsBE.INTEGRATION_TIME+ getSatellite().getMaxSlewDuration())) > 0){
									listResas.add(new Reservation(access_start, target, getSatellite().getMaxSlewDuration()));
									break;
									// After
								} else if(listParallelsResas.get(0).getEndDate().shiftedBy(ConstantsBE.INTEGRATION_TIME+ getSatellite().getMaxSlewDuration()).compareTo(access_end) < 0 ){
									listResas.add(new Reservation(listParallelsResas.get(0).getEndDate(), target, getSatellite().getMaxSlewDuration()));
									break;
								}
							} else {
								// For 2 or more, we check BETWEEN successive reservations, thus needing the sort at the beginning.
								// Initial reservation case
								if(i == 0){
									if(listParallelsResas.get(0).getStartDate().compareTo(access_start.shiftedBy(ConstantsBE.INTEGRATION_TIME+ getSatellite().getMaxSlewDuration())) > 0){
										listResas.add(new Reservation(access_start, target, getSatellite().getMaxSlewDuration()));
										break;
									}
									// Final case
								} else if(i == listParallelsResas.size() -1){
									// N-2 to N-1
									if(listParallelsResas.get(i).getStartDate().compareTo(listParallelsResas.get(i-1).getEndDate().shiftedBy(ConstantsBE.INTEGRATION_TIME+ getSatellite().getMaxSlewDuration())) > 0){
										listResas.add(new Reservation(listParallelsResas.get(i-1).getEndDate(), target, getSatellite().getMaxSlewDuration()));
										break;
									}
									// N to the end of the access
									else if(listParallelsResas.get(i).getEndDate().shiftedBy(ConstantsBE.INTEGRATION_TIME+ getSatellite().getMaxSlewDuration()).compareTo(access_end) < 0 ){
										listResas.add(new Reservation(listParallelsResas.get(i).getEndDate(), target, getSatellite().getMaxSlewDuration()));
										break;
									}
									// Middle Case, between two reservations
								} else {
									// i-1 to i
									if(listParallelsResas.get(i).getStartDate().compareTo(listParallelsResas.get(i-1).getEndDate().shiftedBy(ConstantsBE.INTEGRATION_TIME+ getSatellite().getMaxSlewDuration())) > 0){
										listResas.add(new Reservation(listParallelsResas.get(i-1).getEndDate(), target, getSatellite().getMaxSlewDuration()));
										break;
									}
								}
							}
						}
					}
				}
			}
		}
		// Sorting ensure an easier debugging and is a negligible step in terms of computational time in our project
		listResas.sort(Comparator.comparing(Reservation::getStartDate));

		// We build the observationPlan with AttitudeLawLeg
		for(Reservation resa: listResas){
			AttitudeLaw observationLaw = createObservationLaw(resa.getSite());
			String legName = "OBS_" + resa.getSite().getName();
			AttitudeLawLeg obsLeg = new AttitudeLawLeg(observationLaw, new AbsoluteDateInterval(resa.getStartDate(), resa.getStartDate().shiftedBy(ConstantsBE.INTEGRATION_TIME)), legName);
			this.observationPlan.put(resa.getSite(), obsLeg);
		}
		return this.observationPlan;
	}


	/**
	 * Computes the cinematic plan...
	 * <p>
	 * Here you need to compute the cinematic plan, which is the cinematic chain of
	 * attitude law legs (observation, default law and slews) needed to perform the
	 * mission. Usually, we start and end the mission in default law and during the
	 * horizon, we alternate between default law, observation legs and slew legs.
	 *
	 * @return a {@link StrictAttitudeLegsSequence} that gives all the cinematic
	 * plan of the {@link Satellite}. It is a chronological sequence of all
	 * the {@link AttitudeLawLeg} that are necessary to define the
	 * {@link Attitude} of the {@link Satellite} during all the mission
	 * horizon. Those legs can have 3 natures : pointing a target site,
	 * pointing nadir and performing a slew between one of the two previous
	 * kind of legs.
	 * @throws PatriusException
	 */
	public StrictAttitudeLegsSequence<AttitudeLeg> computeCinematicPlan() throws PatriusException {

		/**
		 * Now we want to assemble a continuous attitude law which is valid during all
		 * the mission horizon. For that, we will use to object
		 * StrictAttitudeLegsSequence<AttitudeLeg> which is a chronological sequence of
		 * AttitudeLeg. In our case, each AttitudeLeg will be an AttitudeLawLeg, either
		 * a leg of site observation, a slew, or the nadir pointing attitude law (see
		 * the Satellite constructor and the BodyCenterGroundPointing class, it's the
		 * Satellite default attitude law). For more help about the Attitude handling,
		 * use the module 11 of the patrius formation.
		 *
		 * Tip 1 : Please give names to the different AttitudeLawLeg you build so that
		 * you can visualize them with VTS later on. For example "OBS_Paris" when
		 * observing Paris or "SlEW_Paris_Lyon" when adding a slew from Paris
		 * observation AttitudeLawLeg to Lyon observation AttitudeLawLeg.
		 *
		 * Tip 2 : the sequence you want to obtain should look like this :
		 * [nadir-slew-obs1-slew-obs2-slew-obs3-slew-nadir] for the simple version where
		 * you don't try to fit nadir laws between observations or
		 * [nadir-slew-obs1-slew-nadir-selw-obs2-slew-obs3-slew-nadir] for the more
		 * complexe version with nadir laws if the slew during two observation is long
		 * enough.
		 *
		 * Tip 3 : You can use the class ConstantSpinSlew(initialAttitude,
		 * finalAttitude, slewName) for the slews. This an AttitudeLeg so you will be
		 * able to add it to the StrictAttitudeLegsSequence as every other leg.
		 */

		// We get the start and end dates, which will be useful respectively for the first and last rounds of nadir laws
		AbsoluteDate start = this.getStartDate();
		AbsoluteDate end = this.getEndDate();

		// We create a list of all sites to be observed, to be processed in chronological order of observation
		List<Site> sortedSites = new ArrayList<>(this.observationPlan.keySet());

		//We create a custom comparator for two attitudes: att1>att2 if att2 starts before att1
		Collections.sort(sortedSites, Comparator.comparing((Site site) -> this.observationPlan.get(site).getTimeInterval()));
		// sortedSites is now actually sorted

		// Getting our nadir law behavior
		AttitudeLaw nadirLaw = this.getSatellite().getDefaultAttitudeLaw();

		// We calculate the maximum time it takes to get from a nadir law (0° relative to earth)
		// to a groundPointingLaw (maximum inclination of ConstantsBE.POINTING_CAPACITY)
		// We add 1s to make up for the non-null rotation speed at the end/beginning of the nadir law
		double MAX_TIME_TO_NADIR = this.getSatellite().computeSlewDuration(ConstantsBE.POINTING_CAPACITY)+1.0;

		// We want to know for a given observation if it's the first/last in the queue
		boolean isFirstObservation = true;
		boolean isLastObservation = false;

		// We store the previous last attitude and sites at the end of each iteration,
		// to be used to calculate the following slew
		Attitude endPreviousAttitude = null;
		Site previousSite = null;

		// We define a currentIndex from 1 to length(listKeys), which is iterated at each observation and
		// helps us calculate if an observation is the last in the queue
		int currentIndex = 0;
		for (final Site currentSite : sortedSites) {
			currentIndex += 1;
			if (currentIndex > 1) {
				isFirstObservation = false;
			}
			if (currentIndex == this.observationPlan.size()) {
				isLastObservation = true;
			}

			// Getting all the dates and observationLeg of the current reservation to compute our slews
			AttitudeLeg currentObsLeg = this.observationPlan.get(currentSite);
			AbsoluteDate obsStart = currentObsLeg.getDate();
			AbsoluteDate obsEnd = currentObsLeg.getEnd();

			// The propagator will be used to compute Attitudes
			KeplerianPropagator propagator = this.createDefaultPropagator();

			// Computing the Attitudes used to compute the slews
			Attitude startObsAttitude = currentObsLeg.getAttitude(propagator, obsStart, getEme2000());
			Attitude endObsAttitude = currentObsLeg.getAttitude(propagator, obsEnd, getEme2000());

			// if we are processing the first observation, we need to end the nadir law in time to point the site:
			// nadirLaw->slew->currentSite
			if (isFirstObservation) {
				AbsoluteDate endNadirLaw1 = obsStart.shiftedBy(-MAX_TIME_TO_NADIR);

				// We create the leg from the start to endNadirLaw
				Attitude endNadir1Attitude = nadirLaw.getAttitude(propagator, endNadirLaw1, getEme2000());
				AttitudeLawLeg nadir1 = new AttitudeLawLeg(nadirLaw, start, endNadirLaw1, "Nadir_Law_1");

				// We calculate the Nadir->firstSite slew
				ConstantSpinSlew slew1 = new ConstantSpinSlew(endNadir1Attitude, startObsAttitude, "Slew_Nadir_to_" + currentSite.getName());

				// We add the computed leg to the cinematicPlan
				cinematicPlan.add(nadir1);
				cinematicPlan.add(slew1);
			}else{
				// if it's not the first observation: we are currently following a GroundPointing law toward a site
				// and need check what to do at the end of the previous attitude

				// if we have the time to insert a nadir leg before the next observation law:
				// previousSite->slewInter1->nadirInter->slewInter2->currentSite
				if(startObsAttitude.getDate().durationFrom(endPreviousAttitude.getDate()) > 2*MAX_TIME_TO_NADIR){
					// We need to calculate the closest end of slewInter1
					AbsoluteDate endNadirSlewInter1 = endPreviousAttitude.getDate().shiftedBy(MAX_TIME_TO_NADIR);
					// and the latest start of slewInter2
					AbsoluteDate beginNadirSlewInter2 = obsStart.shiftedBy(-MAX_TIME_TO_NADIR);

					// We get the nadirAttitude at each of these dates
					Attitude beginNadirIntAttitude = nadirLaw.getAttitude(propagator, endNadirSlewInter1.getDate(), getEme2000());
					Attitude endNadirIntAttitude = nadirLaw.getAttitude(propagator, beginNadirSlewInter2.getDate(), getEme2000());

					// We compute the slewInter1, nadirInter and slewInter2 legs
					AttitudeLawLeg nadirInter = new AttitudeLawLeg(nadirLaw, endNadirSlewInter1, beginNadirSlewInter2, "Nadir_Law_Inter");
					ConstantSpinSlew slewInter1 = new ConstantSpinSlew(endPreviousAttitude, beginNadirIntAttitude, "Slew_"+previousSite.getName()+"_to_NadirInter");
					ConstantSpinSlew slewInter2 = new ConstantSpinSlew(endNadirIntAttitude, startObsAttitude, "Slew_NadirInter_to_"+currentSite.getName());

					// We add those to the cinematicPlan
					cinematicPlan.add(slewInter1);
					cinematicPlan.add(nadirInter);
					cinematicPlan.add(slewInter2);

				}else{
					// if we don't have to insert a nadir law, we go straight to the next observation
					// previousSite->slew->currentSite
					ConstantSpinSlew slew1 = new ConstantSpinSlew(endPreviousAttitude, startObsAttitude, "Slew_"+previousSite.getName()+"_to_"+currentSite.getName());
					cinematicPlan.add(slew1);
				}
			}

			// We add the leg of observation of the considered site to cinematicPlan
			cinematicPlan.add(currentObsLeg);

			// If it's the last observation, we have to insert a nadir law lasting until the end of the cinematicPlan
			// currentSite->slewNadirLaw2->nadir2
			if (isLastObservation) {
				// We calculate the earliest start of the final nadir law after a slew
				AbsoluteDate startNadirLaw2 = obsEnd.shiftedBy(MAX_TIME_TO_NADIR);
				Attitude startNadir2Attitude = nadirLaw.getAttitude(propagator, startNadirLaw2, getEme2000());

				// We compute the corresponding slew and nadir2 legs
				ConstantSpinSlew slew2 = new ConstantSpinSlew(endObsAttitude, startNadir2Attitude, "Slew_" + currentSite.getName() + "_to_Nadir");
				AttitudeLawLeg nadir2 = new AttitudeLawLeg(nadirLaw, startNadirLaw2, end, "Nadir_Law_2");

				// And add these to cinematicPlan
				cinematicPlan.add(slew2);
				cinematicPlan.add(nadir2);
			}

			// At the end of a leg, we store the last attitude and sites for later use
			endPreviousAttitude = endObsAttitude;
			previousSite = currentSite;
		}

		/**
		 * Now your job is finished, the two following methods will finish the job for
		 * you : checkCinematicPlan() will check that each slew's duration is longer
		 * than the theoretical duration it takes to perform the same slew. Then, if the
		 * cinematic plan is valid, computeFinalScore() will compute the score of your
		 * observation plan. Finally, generateVTSVisualization will write all the
		 * ephemeris (Position/Velocity + Attitude) and generate a VTS simulation that
		 * you will be able to play to visualize and validate your plans.
		 */
		return this.cinematicPlan;
	}

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * <p>
	 * Checks the cinematic plan and prints if it's ok or not.
	 * <p>
	 * We provide this method so that you can check if your cinematic plan doesn't
	 * violate a cinematic constraint. It returns a boolean saying if the plan is
	 * valid or not.
	 *
	 * @return A boolean indicating the cinematic validity of the plan.
	 * @throws PatriusException if an error occurs during propagation
	 */
	public boolean checkCinematicPlan() throws PatriusException {

		final KeplerianPropagator propagator = createDefaultPropagator();

		// Checking the cinematic plan validity
		boolean valid = true;
		for (final AttitudeLeg slew : this.cinematicPlan) {
			System.out.println(slew.getNature());

			final Attitude endAtt = slew.getAttitude(propagator, slew.getEnd(), this.getEme2000());
			final Attitude startAtt = slew.getAttitude(propagator, slew.getDate(), this.getEme2000());

			final boolean condition = slew.getTimeInterval().getDuration() > this.getSatellite()
					.computeSlewDuration(startAtt, endAtt);

			if (condition) {
				System.out.println("Cinematic is ok");
			} else {
				valid = false;
				System.out.println("WARNING : cinematic is not realist for this slew");
				System.out.println("Slew actual duration : " + slew.getTimeInterval().getDuration());
				System.out
						.println("Slew duration theory : " + this.getSatellite().computeSlewDuration(startAtt, endAtt));
			}
		}
		System.out.println("==== Is the cinematic plan valid ? => " + valid + " ==== ");
		return valid;
	}

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * <p>
	 * Compute the final score from the observation plan.
	 *
	 * <p>
	 * Note : the observation plan should have unique sites.<br>
	 * The duplicated elements won't be considered for the score (target with no
	 * value = lost opportunity)
	 * </p>
	 *
	 * @return the final score
	 */
	public double computeFinalScore() {

		// Convert the observation plan into a Set of sites to make sure a site will not
		// be evaluated twice
		// The Set makes sure there won't be duplicated sites
		final Set<Site> sitesSet = new HashSet<>(this.observationPlan.size());
		for (final Entry<Site, AttitudeLawLeg> entry : this.observationPlan.entrySet()) {
			sitesSet.add(entry.getKey());
		}

		// Loop over each site and sum its score
		double finalScore = 0.;
		for (final Site site : sitesSet) {
			finalScore += site.getScore();
		}

		return finalScore;
	}

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * <p>
	 * Writes the VTS output files : one CIC-POI file to print the sites of
	 * interest, one CIC-OEM file giving the position and velocity ephemeris of the
	 * satellite, one CIC-AEM file giving the attitude ephemeris of the satellite
	 * pointing Nadir only (to help visualize the access field of view of the
	 * satellite) and one CIC-AEM file giving the attitude ephemeris of the
	 * satellite cinematic plan. Also writes the cinematic plan as a sequence of
	 * pointing modes for the satellite in a CIC-MEM file.
	 *
	 * @param cinematicPlan Input cinematic plan.
	 * @throws PropagationException if an error happens during the {@link Orbit}
	 *                              propagation
	 */
	public void generateVTSVisualization(StrictAttitudeLegsSequence<AttitudeLeg> cinematicPlan)
			throws PropagationException {
		// First, create the propagator for the satellite's pointing capacity view
		// (nadir law)
		final KeplerianPropagator vtsPropagatorNadir = createDefaultPropagator();
		vtsPropagatorNadir.setEphemerisMode();
		vtsPropagatorNadir.propagate(this.getEndDate());

		// Get generated ephemeris
		final BoundedPropagator ephemerisNadir = vtsPropagatorNadir.getGeneratedEphemeris();

		// Then, we create the propagator for the cinematic plan visualization
		final KeplerianPropagator vtsPropagator = new KeplerianPropagator(this.getSatellite().getInitialOrbit(),
				cinematicPlan);

		vtsPropagator.setEphemerisMode();
		vtsPropagator.propagate(this.getEndDate());

		// Get generated ephemeris
		final BoundedPropagator ephemeris = vtsPropagator.getGeneratedEphemeris();

		// Writing the outputs
		final String pathPOI = ConstantsBE.PATH_VTS_DIRECTORY + File.separator + "BE_Supaero_Target_Sites_POI.txt";
		final String pathOEM = ConstantsBE.PATH_VTS_DIRECTORY + File.separator
				+ "BE_Supaero_Satellite_Trajectory_OEM.txt";
		final String pathAEMNadir = ConstantsBE.PATH_VTS_DIRECTORY + File.separator
				+ "BE_Supaero_Nadir_Pointing_AEM.txt";
		final String pathAEMCinematicPlan = ConstantsBE.PATH_VTS_DIRECTORY + File.separator
				+ "BE_Supaero_Cinematic_Plan_AEM.txt";
		final String pathMEMCinematicPlan = ConstantsBE.PATH_VTS_DIRECTORY + File.separator
				+ "BE_Supaero_Cinematic_Plan_Events_MEM.txt";

		System.out.println("\n\nWriting VTS outputs, please wait...");
		VTSTools.generatePOIFile(pathPOI, this.getSiteList());
		VTSTools.generateOEMFile(pathOEM, this.getStartDate(), this.getEndDate(), ephemeris);
		VTSTools.generateAEMFile(pathAEMNadir, this.getStartDate(), this.getEndDate(), ephemerisNadir);
		VTSTools.generateAEMFile(pathAEMCinematicPlan, this.getStartDate(), this.getEndDate(), ephemeris);
		VTSTools.generateLegSequenceMEMFile(pathMEMCinematicPlan, cinematicPlan);
		System.out.println("VTS outputs written");
	}

	/**
	 * This method should compute the input {@link Site}'s access {@link Timeline}.
	 * That is to say the {@link Timeline} which contains all the {@link Phenomenon}
	 * respecting the access conditions for this site : good visibility + corrrect
	 * illumination of the {@link Site}.
	 * <p>
	 * For that, we suggest you create as many {@link Timeline} as you need and
	 * combine them with logical gates to filter only the access windows phenomenon.
	 *
	 * @param targetSite Input target {@link Site}
	 * @return The {@link Timeline} of all the access {@link Phenomenon} for the
	 * input {@link Site}.
	 * @throws PatriusException If a {@link PatriusException} occurs.
	 */
	private Timeline createSiteAccessTimeline(Site targetSite) throws PatriusException {
		/**
		 * Step 1 :
		 *
		 * Create one Timeline per phenomenon you want to monitor.
		 */
		/*
		 * Use the createSiteXTimeline method to create a custom Timeline. More
		 * indication are given inside the method. Note that you will have to code one
		 * method per constraint, for example the method createSiteVisibilityTimeline
		 * for visibility constraint and createSiteIlluminationTimeline for illumination
		 * constraint. All the methods you code can be coded using the given
		 * createSiteXTimeline method as a basis.
		 */

		// We propagate to have a Timeline for each detector thanks to our loggers
		List<Timeline> listTimeline = propagateTimelines(targetSite);

		// Set a global Timeline to serve as final access timeline for our targetSite
		final Timeline siteAccessTimeline = new Timeline(
				new AbsoluteDateInterval(this.getStartDate(), this.getEndDate()));

		// Adding the phenomena of all the considered timelines
		for (Timeline timeline : listTimeline) {
			// ProjectUtilities.printTimeline(timeline, targetSite);
			for (final Phenomenon phenom : timeline.getPhenomenaList()) {
				siteAccessTimeline.addPhenomenon(phenom);
			}
		}

		/**
		 * Step 2 :
		 *
		 * Combine the timelines with logical gates and retrieve only the access
		 * conditions through a refined Timeline object.
		 *
		 * For that, you can use the classes in the events.postprocessing module : for
		 * example, the AndCriterion or the NotCriterion.
		 *
		 * Finally, you can filter only the Phenomenon matching a certain condition
		 * using the ElementTypeFilter
		 */

		// We define the criterion of simultaneous SunIncidence and Visibility
		AndCriterion andCriterionA = new AndCriterion("Visibility", "SunIncidence",
				"Visibility AND SunIncidence", "");

		// We define the criterion of simultaneous (SunIncidence and Visibility) and NonGlare
		AndCriterion andCriterionB = new AndCriterion("Visibility AND SunIncidence",
				"NonGlare",
				"Visibility AND SunIncidence AND NonGlare", "Obs conditions checked");

		// We add those criteria to the actual timeline
		andCriterionA.applyTo(siteAccessTimeline);
		andCriterionB.applyTo(siteAccessTimeline);

		// We filter out elements not respecting our ((SunIncidence and Visibility) and NonGlare) criterion
		final ElementTypeFilter obsThirdConditionFilter = new ElementTypeFilter("Visibility AND SunIncidence AND NonGlare", false);
		obsThirdConditionFilter.applyTo(siteAccessTimeline);

		// We filter out phenomenons lasting less than the minimum observation time INTEGRATION_TIME
		final PhenomenonDurationFilter integrationTimeFilter = new PhenomenonDurationFilter("Visibility AND SunIncidence AND NonGlare", ConstantsBE.INTEGRATION_TIME, true);
		integrationTimeFilter.applyTo(siteAccessTimeline);

		// Log the final access timeline associated to the current target
		// System.out.println("\n" + targetSite.getName());
		// ProjectUtilities.printTimeline(siteAccessTimeline);

		return siteAccessTimeline;
	}

	/**
	 * @param targetSite Input target {@link Site}
	 * @return The {@link List<Timeline>} containing a list of the propagated timeline of access,
	 * respectively phenomenonVisibilityTimeline, phenomenonSunIncidenceTimeline and phenomenonNonGlareTimeline.
	 * @throws PatriusException If a {@link PatriusException} occurs when creating
	 *                          the {@link Timeline}.
	 */
	private List<Timeline> propagateTimelines(Site targetSite) throws PatriusException {

		// We define a new propagator
		KeplerianPropagator propagator = this.createDefaultPropagator();

		// We create a CodedEventsLogger for each phenomenon
		CodedEventsLogger eventVisibilityLogger = createSiteVisibilityLogger(targetSite, propagator);
		CodedEventsLogger eventSunIncidenceLogger = createSiteSunIncidenceLogger(targetSite, propagator);
		CodedEventsLogger eventNonGlareLogger = createSiteNonGlareLogger(targetSite, propagator);

		// We run our propagator
		propagator.propagate(this.getEndDate());

		// We get the corresponding timelines for each type of event
		final Timeline phenomenonVisibilityTimeline = new Timeline(eventVisibilityLogger,
				new AbsoluteDateInterval(this.getStartDate(), this.getEndDate()), null);

		final Timeline phenomenonSunIncidenceTimeline = new Timeline(eventSunIncidenceLogger,
				new AbsoluteDateInterval(this.getStartDate(), this.getEndDate()), null);

		final Timeline phenomenonNonGlareTimeline = new Timeline(eventNonGlareLogger,
				new AbsoluteDateInterval(this.getStartDate(), this.getEndDate()), null);


		return Arrays.asList(phenomenonVisibilityTimeline, phenomenonSunIncidenceTimeline, phenomenonNonGlareTimeline);
	}

	/**
	 * @param targetSite Input target {@link Site}
	 * @param propagator Propagator to use
	 * @return The {@link CodedEventsLogger} setting the conditions for the Visibility phenomenon to monitor.
	 * @throws PatriusException If a {@link PatriusException} occurs when creating
	 *                          the {@link Timeline}.
	 */
	private CodedEventsLogger createSiteVisibilityLogger(Site targetSite, Propagator propagator) throws PatriusException {

		// Creating the relevant EventDetector and add it to our propagator
		EventDetector constraintVisibilityDetector = createConstraintVisibilityDetector(targetSite);
		propagator.addEventDetector(constraintVisibilityDetector);

		// Define the codes for our event
		GenericCodingEventDetector codingEventVisibilityDetector = new GenericCodingEventDetector(constraintVisibilityDetector,
				"Start of visibility", "End of visibility", true, "Visibility");
		CodedEventsLogger eventVisibilityLogger = new CodedEventsLogger();
		EventDetector eventVisibilityDetector = eventVisibilityLogger.monitorDetector(codingEventVisibilityDetector);

		propagator.addEventDetector(eventVisibilityDetector);

		return eventVisibilityLogger;
	}

	/**
	 * @param targetSite Input target {@link Site}
	 * @param propagator Propagator to use
	 * @return The {@link CodedEventsLogger} setting the conditions for the SunIncidence phenomenon to monitor.
	 * @throws PatriusException If a {@link PatriusException} occurs when creating
	 *                          the {@link Timeline}.
	 */
	private CodedEventsLogger createSiteSunIncidenceLogger(Site targetSite, Propagator propagator) throws PatriusException {
		// Creating the relevant EventDetector and add it to our propagator
		EventDetector constraintSunIncidenceDetector = createConstraintSunIncidenceDetector(targetSite);
		propagator.addEventDetector(constraintSunIncidenceDetector);

		// Define the codes for our event
		GenericCodingEventDetector codingEventSunIncidenceDetector = new GenericCodingEventDetector(constraintSunIncidenceDetector,
				"Start of illumination", "End of illumination", true, "SunIncidence");
		CodedEventsLogger eventSunIncidenceLogger = new CodedEventsLogger();
		EventDetector eventSunIncidenceDetector = eventSunIncidenceLogger.monitorDetector(codingEventSunIncidenceDetector);

		propagator.addEventDetector(eventSunIncidenceDetector);

		return eventSunIncidenceLogger;
	}

	/**
	 * @param targetSite Input target {@link Site}
	 * @param propagator Propagator to use
	 * @return The {@link CodedEventsLogger} setting the conditions for the NonGlare phenomenon to monitor.
	 * @throws PatriusException If a {@link PatriusException} occurs when creating
	 *                          the {@link Timeline}.
	 */
	private CodedEventsLogger createSiteNonGlareLogger(Site targetSite, Propagator propagator) throws PatriusException {
		// Creating the relevant EventDetector and add it to our propagator
		EventDetector constraintNonGlareDetector = createConstraintNonGlareDetector(targetSite);
		propagator.addEventDetector(constraintNonGlareDetector);

		// Define the codes for our event
		GenericCodingEventDetector codingEventNonGlareDetector = new GenericCodingEventDetector(constraintNonGlareDetector,
				"Start of correct phase angle", "End of correct phase angle", true, "NonGlare");
		CodedEventsLogger eventNonGlareLogger = new CodedEventsLogger();
		EventDetector eventNonGlareDetector = eventNonGlareLogger.monitorDetector(codingEventNonGlareDetector);

		propagator.addEventDetector(eventNonGlareDetector);

		return eventNonGlareLogger;
	}

	/**
	 * @param targetSite Input target {@link Site}
	 * @return The {@link EventDetector} containing all the {@link Phenomenon} relative
	 * to the Visibility phenomenon to monitor.
	 * @throws PatriusException If a {@link PatriusException} occurs when creating
	 *                          the {@link Timeline}.
	 */
	private EventDetector createConstraintVisibilityDetector(Site targetSite) {
		// Create the rel
		SensorModel sensorModel = new SensorModel(this.getSatellite().getAssembly(), Satellite.SENSOR_NAME);
		sensorModel.addMaskingCelestialBody(this.getEarth());

		PVCoordinatesProvider sitePVCoordinates = new TopocentricFrame(
				this.getEarth(),
				targetSite.getPoint(),
				targetSite.getName()
		);

		sensorModel.setMainTarget(sitePVCoordinates, new ConstantRadiusProvider(0.0)); // Radius = 0: point target

		return new SensorVisibilityDetector(sensorModel,
				MAXCHECK_EVENTS, TRESHOLD_EVENTS, EventDetector.Action.CONTINUE, EventDetector.Action.CONTINUE);
	}

	/**
	 * @param targetSite Input target {@link Site}
	 * @return The {@link EventDetector} containing all the {@link Phenomenon} relative
	 * to the SunIncidence phenomenon to monitor.
	 * @throws PatriusException If a {@link PatriusException} occurs when creating
	 *                          the {@link Timeline}.
	 */
	private EventDetector createConstraintSunIncidenceDetector(Site targetSite) {

		PVCoordinatesProvider sitePVCoordinates = new TopocentricFrame(
				this.getEarth(),
				targetSite.getPoint(),
				targetSite.getName()
		);

		return new ThreeBodiesAngleDetector(
				this.getEarth(),
				sitePVCoordinates,
				this.getSun(),
				FastMath.toRadians(180 - ConstantsBE.MAX_SUN_INCIDENCE_ANGLE),
				MAXCHECK_EVENTS, TRESHOLD_EVENTS,
				EventDetector.Action.CONTINUE
		);
	}

	/**
	 * @param targetSite Input target {@link Site}
	 * @return The {@link EventDetector} containing all the {@link Phenomenon} relative
	 * to the NonGlare phenomenon to monitor.
	 * @throws PatriusException If a {@link PatriusException} occurs when creating
	 *                          the {@link Timeline}.
	 */
	private EventDetector createConstraintNonGlareDetector(Site targetSite) {

		PVCoordinatesProvider sitePVCoordinates = new TopocentricFrame(
				this.getEarth(),
				targetSite.getPoint(),
				targetSite.getName()
		);

		return new ThreeBodiesAngleDetector(
				this.getSun(),
				sitePVCoordinates,
				ThreeBodiesAngleDetector.BodyOrder.SECOND,
				FastMath.toRadians(ConstantsBE.MAX_SUN_PHASE_ANGLE),
				MAXCHECK_EVENTS, TRESHOLD_EVENTS,
				EventDetector.Action.CONTINUE
		);
	}


	/**
	 * Create an observation leg, that is to say an {@link AttitudeLaw} that give
	 * the {@link Attitude} (pointing direction) of the {@link Satellite} in order
	 * to perform the observation of the input target {@link Site}.
	 * <p>
	 * An {@link AttitudeLaw} is an {@link AttitudeProvider} providing the method
	 * {@link AttitudeProvider #getAttitude()} which can be used to compute the
	 * {@link Attitude} of the {@link Satellite} at any given {@link AbsoluteDate}
	 * (instant) during the mission horizon.
	 * <p>
	 * An {@link AttitudeLaw} is valid at anu time in theory.
	 *
	 * @param target Input target {@link Site}
	 * @return An {@link AttitudeLawLeg} adapted to the observation.
	 */
	private AttitudeLaw createObservationLaw(Site target) {
		/**
		 * To perform an observation, the satellite needs to point the target for a
		 * fixed duration.
		 *
		 * Here, you will use the {@link TargetGroundPointing}. This law provides the
		 * Attitude of a Satellite that only points one target at the surface of a
		 * BodyShape. The earth object from the SimpleMission is a BodyShape and we
		 * remind you that the Site object has an attribute which is a GeodeticPoint.
		 * Use those information to your advantage to build a TargetGroundPointing.
		 */

		// We return a TargetGroundPointing, taking in account the 3D vector parameters to ensure the right orientation in space
		TargetGroundPointing targetGroundPointing = new TargetGroundPointing(this.getEarth(), target.getPoint(), Vector3D.MINUS_K, Vector3D.PLUS_I);

		return targetGroundPointing;
	}

	@Override
	public String toString() {
		return "CompleteMission [name=" + this.getName() + ", startDate=" + this.getStartDate() + ", endDate="
				+ this.getEndDate() + ", satellite=" + this.getSatellite() + "]";
	}
}
