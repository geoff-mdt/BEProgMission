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
import fr.cnes.sirius.patrius.time.AbsoluteDateIntervalsList;
import fr.cnes.sirius.patrius.utils.exception.PatriusException;
import fr.cnes.sirius.patrius.utils.exception.PropagationException;
import org.apache.commons.collections.comparators.ReverseComparator;
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
	 * [COMPLETE THIS METHOD TO ACHIEVE YOUR PROJECT]
	 * <p>
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
		/*
		 * We give a very basic example of incomplete code computing the first target
		 * site access Timeline and adding it to the accessPlan.
		 *
		 * Please complete the code below.
		 */
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
	 * [COMPLETE THIS METHOD TO ACHIEVE YOUR PROJECT]
	 * <p>
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

	public boolean checkIntervalAvailability(AbsoluteDateInterval obsInterval) {
		for (AttitudeLawLeg obsLeg : this.observationPlan.values()) {
			if (obsInterval.overlaps(obsLeg.getTimeInterval())) {
				return false;
			}
		}
		return true;
	}

	/**
	 * [COMPLETE THIS METHOD TO ACHIEVE YOUR PROJECT]
	 * <p>
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
		 * finalAttitude, slewName) for the slews. This an AtittudeLeg so you will be
		 * able to add it to the StrictAttitudeLegsSequence as every other leg.
		 */

		/*
		 * Example of code using our observation plan, let's say we only have one obs
		 * pointing Paris.
		 *
		 * Then we are going to create a very basic cinematic plan : nadir law => slew
		 * => obsParis => slew => nadir law
		 *
		 * To do that, we need to compute the slew duration from the end of nadir law to
		 * the beginning of Paris obs and then from the end of Paris obs to the beginning
		 * of nadir law. For that, we use the Satellite#computeSlewDurationMethod() as
		 * before. We know we have to the time to perform the slew thanks to the
		 * cinematic checks we already did during the observation plan computation.
		 */

		AbsoluteDate start = this.getStartDate();
		AbsoluteDate end = this.getEndDate();

		LinkedHashMap<Site, AttitudeLawLeg> sortedPlan = new LinkedHashMap<>();
		ArrayList<AttitudeLawLeg> list = new ArrayList<>();
		for (Map.Entry<Site, AttitudeLawLeg> entry : this.observationPlan.entrySet()) {
			list.add(entry.getValue());
		}
		Collections.sort(list, new Comparator<AttitudeLeg>() {
			@Override
			public int compare(AttitudeLeg att1, AttitudeLeg att2) {
				return (att1).getTimeInterval().compareTo(att2.getTimeInterval());
			}
		});
		for (AttitudeLawLeg att : list) {
			for (Entry<Site, AttitudeLawLeg> entry : this.observationPlan.entrySet()) {
				if (entry.getValue().equals(att)) {
					sortedPlan.put(entry.getKey(), att);
				}
			}
		}

		Set<Site> keySet = sortedPlan.keySet();
		List<Site> listKeys = new ArrayList<Site>(keySet);

		// Getting our nadir law
		AttitudeLaw nadirLaw = this.getSatellite().getDefaultAttitudeLaw();

		double MAX_TIME_TO_NADIR = this.getSatellite().computeSlewDuration(ConstantsBE.POINTING_CAPACITY)+1.0;
		boolean isFirstObservation = true;
		boolean isLastObservation = false;
		Attitude endPreviousAttitude = null;
		//AbsoluteDate endPreviousAttitudeLaw = null;
		Site previousSite = null;

		int currentIndex = 0;
		for (final Site currentSite : listKeys) {
			currentIndex += 1;
			if (currentIndex > 1) {
				isFirstObservation = false;
			}
			if (currentIndex == sortedPlan.size()) {
				isLastObservation = true;
			}

			AttitudeLeg currentObsLeg = sortedPlan.get(currentSite);
			//AttitudeLeg parisObsLeg = observationPlan.get(currentSite);

			// Getting all the dates we need to compute our slews
			AbsoluteDate obsStart = currentObsLeg.getDate();
			AbsoluteDate obsEnd = currentObsLeg.getEnd();

			// The propagator will be used to compute Attitudes
			KeplerianPropagator propagator = this.createDefaultPropagator();

			// For the slew nadir => paris and paris => nadir, we will use the maximum
			// duration because we have a lot of time here. In practice, you will use either
			// the maximum possible time if you have nothing else planned around or the
			// available time coming from the duration until next observation programmed.

			// Computing the Attitudes used to compute the slews
			Attitude startObsAttitude = currentObsLeg.getAttitude(propagator, obsStart, getEme2000());
			Attitude endObsAttitude = currentObsLeg.getAttitude(propagator, obsEnd, getEme2000());


			if (isFirstObservation) {
				AbsoluteDate endNadirLaw1 = obsStart.shiftedBy(-MAX_TIME_TO_NADIR);
				// We create our two Nadir legs using the dates we computed
				AttitudeLawLeg nadir1 = new AttitudeLawLeg(nadirLaw, start, endNadirLaw1, "Nadir_Law_1");
				// From nadir law 1 to current observation
				Attitude endNadir1Attitude = nadirLaw.getAttitude(propagator, endNadirLaw1, getEme2000());

				ConstantSpinSlew slew1 = new ConstantSpinSlew(endNadir1Attitude, startObsAttitude, "Slew_Nadir_to_" + currentSite.getName());

				cinematicPlan.add(nadir1);
				cinematicPlan.add(slew1);
			}else{
				// Si on a le temps de repasser au nadir, on le fait
				if(startObsAttitude.getDate().durationFrom(endPreviousAttitude.getDate()) > 2*MAX_TIME_TO_NADIR){
					AbsoluteDate endNadirSlewInter1 = endPreviousAttitude.getDate().shiftedBy(MAX_TIME_TO_NADIR);
					AbsoluteDate beginNadirSlewInter2 = obsStart.shiftedBy(-MAX_TIME_TO_NADIR);

					Attitude beginNadirIntAttitude = nadirLaw.getAttitude(propagator, endNadirSlewInter1.getDate(), getEme2000());
					Attitude endNadirIntAttitude = nadirLaw.getAttitude(propagator, beginNadirSlewInter2.getDate(), getEme2000());

					AttitudeLawLeg nadirInter = new AttitudeLawLeg(nadirLaw, endNadirSlewInter1, beginNadirSlewInter2, "Nadir_Law_Inter");
					ConstantSpinSlew slewInter1 = new ConstantSpinSlew(endPreviousAttitude, beginNadirIntAttitude, "Slew_"+previousSite.getName()+"_to_NadirInter");
					ConstantSpinSlew slewInter2 = new ConstantSpinSlew(endNadirIntAttitude, startObsAttitude, "Slew_NadirInter_to_"+currentSite.getName());

					cinematicPlan.add(slewInter1);
					cinematicPlan.add(nadirInter);
					cinematicPlan.add(slewInter2);

				}else{
					ConstantSpinSlew slew1 = new ConstantSpinSlew(endPreviousAttitude, startObsAttitude, "Slew_"+previousSite.getName()+"_to_"+currentSite.getName());
					cinematicPlan.add(slew1);
				}
			}


			cinematicPlan.add(currentObsLeg);

			// Finally computing the slews

			if (isLastObservation) {
				AbsoluteDate startNadirLaw2 = obsEnd.shiftedBy(+getSatellite().getMaxSlewDuration());
				Attitude startNadir2Attitude = nadirLaw.getAttitude(propagator, startNadirLaw2, getEme2000());
				// From currentObservation observation to nadir law 2
				ConstantSpinSlew slew2 = new ConstantSpinSlew(endObsAttitude, startNadir2Attitude, "Slew_" + currentSite.getName() + "_to_Nadir");
				AttitudeLawLeg nadir2 = new AttitudeLawLeg(nadirLaw, startNadirLaw2, end, "Nadir_Law_2");
				cinematicPlan.add(slew2);
				cinematicPlan.add(nadir2);
			}

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
	 * [COMPLETE THIS METHOD TO ACHIEVE YOUR PROJECT]
	 * <p>
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
		/*
		 * Code your logical operations on Timeline objects and filter only the access
		 * Phenomenon (gathering all constraints you need to define an access condition)
		 * below.
		 */
		// Combining all Timelines
		// Creating a global Timeline containing all phenomena, this Timeline will be
		// filtered and processed to that only the access Phenomenon remain, this is
		// our siteAccessTimeline

		AndCriterion andCriterionA = new AndCriterion("Visibility", "SunIncidence",
				"Visibility AND SunIncidence", "");
		// NotCriterion notNonGlareCriterion = new NotCriterion("NonGlare", "NotNonGlare", "");
		AndCriterion andCriterionB = new AndCriterion("Visibility AND SunIncidence",
				"NonGlare",
				"Visibility AND SunIncidence AND NonGlare", "Obs conditions checked");

		andCriterionA.applyTo(siteAccessTimeline);
		andCriterionB.applyTo(siteAccessTimeline);

		final ElementTypeFilter obsThirdConditionFilter = new ElementTypeFilter("Visibility AND SunIncidence AND NonGlare", false);
		obsThirdConditionFilter.applyTo(siteAccessTimeline);
		final PhenomenonDurationFilter integrationTimeFilter = new PhenomenonDurationFilter("Visibility AND SunIncidence AND NonGlare", ConstantsBE.INTEGRATION_TIME, true);
		integrationTimeFilter.applyTo(siteAccessTimeline);


		// Log the final access timeline associated to the current target
		// System.out.println("\n" + targetSite.getName());
		// ProjectUtilities.printTimeline(siteAccessTimeline);

		return siteAccessTimeline;
	}

	private List<Timeline> propagateTimelines(Site targetSite) throws PatriusException {

		KeplerianPropagator propagator = this.createDefaultPropagator();

		CodedEventsLogger eventVisibilityLogger = createSiteVisibilityLogger(targetSite, propagator);
		CodedEventsLogger eventSunIncidenceLogger = createSiteSunIncidenceLogger(targetSite, propagator);
		CodedEventsLogger eventNonGlareLogger = createSiteNonGlareLogger(targetSite, propagator);

		propagator.propagate(this.getEndDate());

		final Timeline phenomenonVisibilityTimeline = new Timeline(eventVisibilityLogger,
				new AbsoluteDateInterval(this.getStartDate(), this.getEndDate()), null);

		final Timeline phenomenonSunIncidenceTimeline = new Timeline(eventSunIncidenceLogger,
				new AbsoluteDateInterval(this.getStartDate(), this.getEndDate()), null);

		final Timeline phenomenonNonGlareTimeline = new Timeline(eventNonGlareLogger,
				new AbsoluteDateInterval(this.getStartDate(), this.getEndDate()), null);


		return Arrays.asList(phenomenonVisibilityTimeline, phenomenonSunIncidenceTimeline, phenomenonNonGlareTimeline);
	}

	/**
	 * [COPY-PASTE AND COMPLETE THIS METHOD TO ACHIEVE YOUR PROJECT]
	 * <p>
	 * This method should compute a {@link Timeline} object which encapsulates all
	 * the {@link Phenomenon} corresponding to an orbital phenomenon X relative to
	 * the input target {@link Site}. For example, X can be the {@link Site}
	 * visibility phenomenon.
	 * <p>
	 * You can copy-paste this method and adapt it for every X {@link Phenomenon}
	 * and {@link Timeline} you need to implement. The global process described here
	 * stays the same.
	 *
	 * @param targetSite Input target {@link Site}
	 * @return The {@link Timeline} containing all the {@link Phenomenon} relative
	 * to the X phenomenon to monitor.
	 * @throws PatriusException If a {@link PatriusException} occurs when creating
	 *                          the {@link Timeline}.
	 */
	private Timeline createSiteXTimeline(Site targetSite) throws PatriusException {
		/**
		 * Here is a quick idea of how to compute a Timeline. A Timeline contains a
		 * PhenomenaList, which is list of Phenomenon objects. Each Phenomenon object
		 * represents a phenomenon in orbit which is defined between two AbsoluteDate
		 * objects and their associated CodedEvent which define the begin and the end of
		 * the Phenomenon. For example, the Sun visibility can be defined as a
		 * phenomenon beginning with the start of visibility and ending with the end of
		 * visibility, itself defined using geometrical rules.
		 *
		 * Now, how to create a Phenomenon object matching the requirement of a given
		 * orbital phenomenon.
		 *
		 * For that, you can use Patrius possibilities with the
		 * "fr.cnes.sirius.patrius.propagation.events", "fr.cnes.sirius.patrius.events",
		 * "fr.cnes.sirius.patrius.events.sensor" and the
		 * "fr.cnes.sirius.patrius.events.postprocessing" modules. See the modules 05
		 * and 09 of the Patrius formation for those aspects, you have examples of codes
		 * using those modules and how to build a Timeline derived from other objects in
		 * a representative case.
		 *
		 * Below are some basic steps and tips to help you search for the right
		 * information in Javadoc and in the Patrius formation in order to compute your
		 * Timeline.
		 *
		 */

		/**
		 * Step 1 :
		 *
		 * Here we deal with event detection. As explain in the module 05, this is done
		 * with EventDetector objects. If you look at the Javadoc, you'll find you all
		 * sorts of detectors. You need to translate the X input constraint (for example
		 * an incidence angle between the sensor and the target, sun incidence angle,
		 * masking of the target by the Earth, etc.) into an EventDetector object.
		 * Scroll through the event detection modules to find the one adapted to your
		 * problem (represented by the X constraint which describe the X phenomenon you
		 * want to detect) and then look at the inputs you need to build it.
		 *
		 * Please note that in order to facilitate the task for you, we provide the
		 * object Satellite. If you look how the constructor build this object, you will
		 * find that our Satellite already has an Assembly filed with a lot of
		 * properties. Among those properties, there is a SensorProperty that you can
		 * use to your advantage when trying to build you detector (for example when
		 * trying to build a visibility detector). See the module 7 of the formation to
		 * learn more about the Assembly object. You can use the SensorProperty via the
		 * Assembly of the Satellite and its name to define appropriate detectors.
		 *
		 */
		/*
		 * Complete the method below to build your detector. More indications are given
		 * in the method.
		 */
		EventDetector constraintXDetector = createConstraintXDetector();

		/**
		 * Step 2 :
		 *
		 * When you have your detector, you can add it on an Orbit Propagator such as
		 * the KeplerianPropagator of your Satellite. If you give the detector the right
		 * parameters, you can then propagate the orbit (see the SimpleMission code and
		 * the module 03 from the Patrius formation) and the detector will automatically
		 * perform actions when a particular orbital event happens (you need to
		 * configure the right detector to detect the event you want to monitor).
		 *
		 * You can add several detectors to the propagator (one per constraint per Site
		 * for example).
		 */
		/*
		 * This is how you add a detector to a propagator, feel free to add several
		 * detectors to the satellite propagator !
		 */
		this.getSatellite().getPropagator().addEventDetector(constraintXDetector);

		/**
		 * Step 3 :
		 *
		 * Now you need to use the detector's ability to create CodedEvent objects to
		 * actually detect the events and visualize them. You can obtain CodedEvents
		 * with a CodedEventsLogger that you plug on an EventDetector with the
		 * CodedEventsLogger.monitorDetector() method. For that, you will need the
		 * GenericCodingEventDetector class. See the module 09 to understand how to use
		 * those objects in order to detect events.
		 */
		/*
		 * Develop the code in which you create your GenericCodingEventDetector and use
		 * it to create a CodedEventsLogger here. You have some example code to help.
		 */
		GenericCodingEventDetector codingEventXDetector = new GenericCodingEventDetector(constraintXDetector,
				"Event starting the X phenomenon", "Event ending the X phenomenon", true, "Name of the X phenomenon");
		CodedEventsLogger eventXLogger = new CodedEventsLogger();
		EventDetector eventXDetector = eventXLogger.monitorDetector(codingEventXDetector);
		// Then you add your logger to the propagator, it will monitor the event coded
		// by the codingEventDetector
		this.getSatellite().getPropagator().addEventDetector(eventXDetector);

		/**
		 * Step 4 :
		 *
		 * Now you can propagate your orbit and the propagator will use the added
		 * detectors and loggers the way you defined them, detecting all events you
		 * wanted to monitor.
		 */
		// Finally propagating the orbit
		this.getSatellite().getPropagator().propagate(this.getEndDate());
		/**
		 * Remark : since you can add as many EventDetectors as you want to a
		 * propagator, you might want to delay this step afterwards to propagate the orbit
		 * with all your detectors at once. Here we do it here to provide an example but
		 * feel free to code your own more performant version of it.
		 */

		/**
		 * Step 5 : WARNING : this can only be done after the propagation !
		 *
		 * Now, you have to post process all your events. That's when you actually
		 * create your Timeline object which contains the Phenomenon you want to
		 * monitor.
		 *
		 * Since you have propagated your orbit, the events that have been detected are
		 * stored inside the detector and logger. This mechanic is used to create a
		 * Timeline.
		 */
		/*
		 * See code below and create your own code to have your X Timeline describing
		 * all X phenomenon you want to detect.
		 */
		// Creating a Timeline to process the events : we are going to define one
		// visibility Phenomenon by couple of events "start -> end" (linked to the
		// increase and decrease of the g function of the visibility detector)
		final Timeline phenomenonXTimeline = new Timeline(eventXLogger,
				new AbsoluteDateInterval(this.getStartDate(), this.getEndDate()), null);

		return phenomenonXTimeline;
	}


	/**
	 * @param targetSite Input target {@link Site}
	 * @return The {@link Timeline} containing all the {@link Phenomenon} relative
	 * to the X phenomenon to monitor.
	 * @throws PatriusException If a {@link PatriusException} occurs when creating
	 *                          the {@link Timeline}.
	 */
	private CodedEventsLogger createSiteVisibilityLogger(Site targetSite, Propagator propagator) throws PatriusException {

		EventDetector constraintVisibilityDetector = createConstraintVisibilityDetector(targetSite);


		propagator.addEventDetector(constraintVisibilityDetector);


		GenericCodingEventDetector codingEventVisibilityDetector = new GenericCodingEventDetector(constraintVisibilityDetector,
				"Start of visibility", "End of visibility", true, "Visibility");
		CodedEventsLogger eventVisibilityLogger = new CodedEventsLogger();
		EventDetector eventVisibilityDetector = eventVisibilityLogger.monitorDetector(codingEventVisibilityDetector);

		propagator.addEventDetector(eventVisibilityDetector);


		// Finally propagating the orbit
		//this.getSatellite().getPropagator().propagate(this.getEndDate());
		/**
		 * Remark : since you can add as many EventDetectors as you want to a
		 * propagator, you might want to delay this step afterwards to propagate the orbit
		 * with all your detectors at once. Here we do it here to provide an example but
		 * feel free to code your own more performant version of it.
		 */
		//final Timeline phenomenonVisibilityTimeline = new Timeline(eventVisibilityLogger,
		//		new AbsoluteDateInterval(this.getStartDate(), this.getEndDate()), null);

		return eventVisibilityLogger;
	}

	/**
	 * @param targetSite Input target {@link Site}
	 * @return The {@link Timeline} containing all the {@link Phenomenon} relative
	 * to the X phenomenon to monitor.
	 * @throws PatriusException If a {@link PatriusException} occurs when creating
	 *                          the {@link Timeline}.
	 */
	private CodedEventsLogger createSiteSunIncidenceLogger(Site targetSite, Propagator propagator) throws PatriusException {
		EventDetector constraintSunIncidenceDetector = createConstraintSunIncidenceDetector(targetSite);


		propagator.addEventDetector(constraintSunIncidenceDetector);


		GenericCodingEventDetector codingEventSunIncidenceDetector = new GenericCodingEventDetector(constraintSunIncidenceDetector,
				"Start of illumination", "End of illumination", true, "SunIncidence");
		CodedEventsLogger eventSunIncidenceLogger = new CodedEventsLogger();
		EventDetector eventSunIncidenceDetector = eventSunIncidenceLogger.monitorDetector(codingEventSunIncidenceDetector);

		propagator.addEventDetector(eventSunIncidenceDetector);


		// Finally propagating the orbit
		// this.getSatellite().getPropagator().propagate(this.getEndDate());
		/**
		 * Remark : since you can add as many EventDetectors as you want to a
		 * propagator, you might want to delay this step afterwards to propagate the orbit
		 * with all your detectors at once. Here we do it here to provide an example but
		 * feel free to code your own more performant version of it.
		 */
		//final Timeline phenomenonSunIncidenceTimeline = new Timeline(eventSunIncidenceLogger,
		//		new AbsoluteDateInterval(this.getStartDate(), this.getEndDate()), null);

		return eventSunIncidenceLogger;
	}

	/**
	 * @param targetSite Input target {@link Site}
	 * @return The {@link Timeline} containing all the {@link Phenomenon} relative
	 * to the X phenomenon to monitor.
	 * @throws PatriusException If a {@link PatriusException} occurs when creating
	 *                          the {@link Timeline}.
	 */
	private CodedEventsLogger createSiteNonGlareLogger(Site targetSite, Propagator propagator) throws PatriusException {
		EventDetector constraintNonGlareDetector = createConstraintNonGlareDetector(targetSite);


		propagator.addEventDetector(constraintNonGlareDetector);


		GenericCodingEventDetector codingEventNonGlareDetector = new GenericCodingEventDetector(constraintNonGlareDetector,
				"Start of correct phase angle", "End of correct phase angle", true, "NonGlare");
		CodedEventsLogger eventNonGlareLogger = new CodedEventsLogger();
		EventDetector eventNonGlareDetector = eventNonGlareLogger.monitorDetector(codingEventNonGlareDetector);

		propagator.addEventDetector(eventNonGlareDetector);

		return eventNonGlareLogger;
	}


	/**
	 * [COPY-PASTE AND COMPLETE THIS METHOD TO ACHIEVE YOUR PROJECT]
	 * <p>
	 * Create an adapted instance of {@link EventDetector} matching the input need
	 * for monitoring the events defined by the X constraint. (X can be a lot of
	 * things).
	 * <p>
	 * You can copy-paste this method to adapt it to the {@link EventDetector} X
	 * that you want to create.
	 * <p>
	 * Note : this can have different inputs that we don't define here
	 *
	 * @return An {@link EventDetector} answering the constraint (for example a
	 * {@link SensorVisibilityDetector} for a visibility constraint).
	 */
	private EventDetector createConstraintXDetector() {
		/**
		 * Here you build an EventDetector object that correspond to the constraint X :
		 * visibility of the target from the satellite, target is in day time, whatever.
		 *
		 * Note that when you create a detector, you choose the actions that it will
		 * perform when the target event is detected. See the module 5 for more
		 * information about this.
		 *
		 * Tip 1 : For the visibility detector, you can use a SensorModel. You will have
		 * to add the Earth as a masking body with the method addMaskingCelestialBody
		 * and to set the main target of the SensorModel with the method setMainTarget.
		 * Then, you can use the class SensorVisibilityDetector with your SensorModel.
		 *
		 * Tip 2 : For the sun incidence angle detector (illumination condition), you
		 * can use the class ThreeBodiesAngleDetector, the three bodies being the ground
		 * target, the Earth and the Sun. See the inputs of this class to build the
		 * object properly.
		 *
		 * Tip 3 : When you create the detectors listed above, you can use the two
		 * public final static fields MAXCHECK_EVENTS and TRESHOLD_EVENTS to configure
		 * the detector (those values are often asked in input of the EventDectector
		 * classes). You will also indicate the Action to perform when the detection
		 * occurs, which is Action.CONTINUE.
		 */
		/*
		 * Create your detector and return it.
		 */
		return null;
	}

	private EventDetector createConstraintVisibilityDetector(Site targetSite) {

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
	 * [COMPLETE THIS METHOD TO ACHIEVE YOUR PROJECT]
	 * <p>
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


		/*
		 * Complete the code below to create your observation law and return it
		 */

		TargetGroundPointing targetGroundPointing = new TargetGroundPointing(this.getEarth(), target.getPoint(), Vector3D.MINUS_K, Vector3D.PLUS_I);

		return targetGroundPointing;
	}

	@Override
	public String toString() {
		return "CompleteMission [name=" + this.getName() + ", startDate=" + this.getStartDate() + ", endDate="
				+ this.getEndDate() + ", satellite=" + this.getSatellite() + "]";
	}
}
