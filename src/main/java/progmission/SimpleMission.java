package progmission;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.opencsv.exceptions.CsvValidationException;

import fr.cnes.sirius.addons.patriusdataset.PatriusDataset;
import fr.cnes.sirius.patrius.attitudes.AttitudeLawLeg;
import fr.cnes.sirius.patrius.attitudes.AttitudeLeg;
import fr.cnes.sirius.patrius.attitudes.StrictAttitudeLegsSequence;
import fr.cnes.sirius.patrius.bodies.CelestialBody;
import fr.cnes.sirius.patrius.bodies.CelestialBodyFactory;
import fr.cnes.sirius.patrius.bodies.ExtendedOneAxisEllipsoid;
import fr.cnes.sirius.patrius.frames.Frame;
import fr.cnes.sirius.patrius.frames.FramesFactory;
import fr.cnes.sirius.patrius.math.util.FastMath;
import fr.cnes.sirius.patrius.orbits.CircularOrbit;
import fr.cnes.sirius.patrius.orbits.PositionAngle;
import fr.cnes.sirius.patrius.propagation.analytical.KeplerianPropagator;
import fr.cnes.sirius.patrius.time.AbsoluteDate;
import fr.cnes.sirius.patrius.time.AbsoluteDateInterval;
import fr.cnes.sirius.patrius.time.TimeScale;
import fr.cnes.sirius.patrius.time.TimeScalesFactory;
import fr.cnes.sirius.patrius.utils.Constants;
import fr.cnes.sirius.patrius.utils.exception.PatriusException;
import fr.cnes.sirius.patrius.utils.exception.PropagationException;
import reader.Site;
import reader.SitesReader;
import utils.ConstantsBE;
import utils.VTSTools;

/**
 * Simple mission class to simulate an Earth observation mission with one
 * satellite.
 * 
 * @author herberl
 *
 */
public class SimpleMission {

	/** Name of the mission. */
	private final String name;

	/**
	 * EME2000 {@link Frame} will be the reference frame for our
	 * {@link SimpleMission} context.
	 */
	private final Frame eme2000;

	/** Reference TimeScale for the {@link SimpleMission} context. */
	private final TimeScale tai;

	/**
	 * Earth model. See {@link ExtendedOneAxisEllipsoid} for more details about this
	 * model.
	 */
	private final ExtendedOneAxisEllipsoid earth;

	/** Default Patrius Sun model. */
	private final CelestialBody sun;

	/** Initial date of the mission. */
	private final AbsoluteDate startDate;

	/** Final date of the mission. */
	private final AbsoluteDate endDate;

	/**
	 * Mission satellite object. The {@link Satellite} is a custom class developped
	 * in the context of the BEProgrammationMission project.
	 */
	private final Satellite satellite;

	/**
	 * {@link List} of {@link Site} objects, representing the list of observation
	 * targets for our Earth observation {@link Satellite}. The {@link Site} object
	 * is a custom class developped in the context of the BEProgrammationMission
	 * project.
	 */
	private final List<Site> siteList;

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * 
	 * Basic constructor for the {@link SimpleMission} object.
	 * 
	 * @param missionName Name of the mission
	 * @param numberOfSites Number of target {@link Site} to consider, please give a
	 *                    number between 1 and 100.
	 * @throws PatriusException If a {@link PatriusException} occurs
	 */
	public SimpleMission(final String missionName, int numberOfSites) throws PatriusException {
		// Naming our mission
		this.name = missionName;

		// Loading Patrius Dataset resources
		// Note that the PatriusDataset is necessary to provide all the physical context
		// (offsets bewteen Timescales, definition of the Frames, basic ephemeris of the
		// Celesital bodies, etc.
		PatriusDataset.addResourcesFromPatriusDataset();

		// Instantiating the space-time physical context
		this.eme2000 = FramesFactory.getEME2000();
		this.tai = TimeScalesFactory.getTAI();

		// Instantiating the Earth model as a one axis ellipsoid
		// Equatorial radius
		final double ae = Constants.WGS84_EARTH_EQUATORIAL_RADIUS;
		// Flattening : here we use a Sherical Earth (f=0)
		final double f = 0.;
		// ITRF is the Frame attached to the Earth (International Terrestrial Reference
		// Frame)
		this.earth = new ExtendedOneAxisEllipsoid(ae, f, FramesFactory.getITRF(), "Earth");

		// Instantiating the Sun, using the default Sun body of Patrius
		this.sun = CelestialBodyFactory.getSun();

		// Instantiating the mission horizon
		this.startDate = ConstantsBE.START_DATE;
		this.endDate = ConstantsBE.END_DATE;

		// Keplerian orbital parameters for a sun-synchronous orbit like Pleiades orbit
		final double a = Constants.WGS84_EARTH_EQUATORIAL_RADIUS + ConstantsBE.ALTITUDE;
		final double i = FastMath.toRadians(ConstantsBE.INCLINATION);
		final double raan = FastMath.toRadians(ConstantsBE.ASCENDING_NODE_LONGITUDE);
		final double e = FastMath.toRadians(ConstantsBE.MEAN_ECCENTRICITY);

		// Initial orbit, we use a simple circular orbit
		final CircularOrbit initialOrbit = new CircularOrbit(a, e, e, i, raan, 0.0, PositionAngle.TRUE, this.eme2000,
				this.startDate, Constants.WGS84_EARTH_MU);

		// Instantiating the satellite (see the Satellite object for more details)
		this.satellite = new Satellite(this, "Pleiades 1A", initialOrbit);

		// Read the sites list and extract only the top N ranking elements
		// Number of site to evaluate [1-100] : a low number improve the performance,
		// might be interesting for testing
		this.siteList = extractTopRankingSites(numberOfSites);

	}

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * 
	 * Read the sites list and extract only the top N ranking elements.
	 * 
	 * @param n Number of elements to return
	 * @return top N ranking sites list
	 * @throws IllegalStateException if an error occurs during reading the sites
	 *                               list file if the original sites list hasn't
	 *                               enough elements to be resized
	 */
	private static List<Site> extractTopRankingSites(final int n) {
		// Read the original sites list
		List<Site> fullSiteList;
		try {
			fullSiteList = SitesReader.readSites(SitesReader.OBSERVATION_SITES_FILE);
		} catch (CsvValidationException | NumberFormatException | IOException e) {
			throw new IllegalStateException(e);
		}

		// Check if the original sites list has enough elements to be resized
		if (fullSiteList.size() < n) {
			throw new IllegalStateException("The original sites list hasn't enough elements to be resized");
		}

		// Sort the list: the top rank elements come first
		Collections.sort(fullSiteList);

		// Resize the list with N elements
		return fullSiteList.subList(0, n);
	}

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * 
	 * Write the AEM and OEM CIC files to visualize a simple simulation with
	 * VTS.Also write the VTS ROI file to indicate to VTS to plot only the sites
	 * that have been selected for the mission. In this simulation we only have our
	 * satellite pointing nadir and its agility cone under it.
	 * 
	 * @throws PropagationException if the propagator fails
	 */
	public void createSimpleVTSVisualization() throws PropagationException {
		// Creating VTS outputs for basic trajectory with Nadir pointing with a
		// default propagator (pointing Nadir)
		final KeplerianPropagator vtsPropagator = createDefaultPropagator();

		// This object is used for the cinematic plan, don't bother understanding it at
		// this point, we use it as an AtittudeLaw object to write the orbital events
		// file in VTS simulation. In the Nadir case, this file is almost empty since we
		// have only on AtittudeLaw pointing only Nadir.
		final StrictAttitudeLegsSequence<AttitudeLeg> simpleSequence = new StrictAttitudeLegsSequence<>();
		simpleSequence.add(new AttitudeLawLeg(this.getSatellite().getDefaultAttitudeLaw(),
				new AbsoluteDateInterval(startDate, endDate)));

		// Writing the outputs using the VTSTools class

		// Declaring the paths of the output files
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

		// Calling the writing methods
		VTSTools.generatePOIFile(pathPOI, this.getSiteList());
		VTSTools.generateOEMFile(pathOEM, this.getStartDate(), this.getEndDate(), vtsPropagator);
		VTSTools.generateAEMFile(pathAEMNadir, this.getStartDate(), this.getEndDate(), vtsPropagator);
		VTSTools.generateAEMFile(pathAEMCinematicPlan, this.getStartDate(), this.getEndDate(), vtsPropagator);
		VTSTools.generateLegSequenceMEMFile(pathMEMCinematicPlan, simpleSequence);
		System.out.println("VTS outputs written");
	}

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * 
	 * Creating a new {@link KeplerianPropagator} based on the satellite
	 * propagator's configuration.
	 * 
	 * We use this method because the satellite propagator might be saturated with
	 * all its events detectors and loggers which slows the computations.
	 * 
	 * @return the new {@link KeplerianPropagator}
	 * @throws PropagationException if initial attitude cannot be computed
	 */
	public KeplerianPropagator createDefaultPropagator() throws PropagationException {
		return new KeplerianPropagator(this.getSatellite().getInitialOrbit(),
				this.getSatellite().getDefaultAttitudeLaw());
	}

	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * @return the eme2000
	 */
	public Frame getEme2000() {
		return eme2000;
	}

	/**
	 * @return the tai
	 */
	public TimeScale getTai() {
		return tai;
	}

	/**
	 * @return the earth
	 */
	public ExtendedOneAxisEllipsoid getEarth() {
		return earth;
	}

	/**
	 * @return the sun
	 */
	public CelestialBody getSun() {
		return sun;
	}

	/**
	 * @return the startDate
	 */
	public AbsoluteDate getStartDate() {
		return startDate;
	}

	/**
	 * @return the endDate
	 */
	public AbsoluteDate getEndDate() {
		return endDate;
	}

	/**
	 * @return the satellite
	 */
	public Satellite getSatellite() {
		return satellite;
	}

	/**
	 * @return the siteList
	 */
	public List<Site> getSiteList() {
		return siteList;
	}

	@Override
	public String toString() {
		return "CompleteMission [name=" + this.getName() + ", startDate=" + this.getStartDate() + ", endDate="
				+ this.getEndDate() + ", satellite=" + this.getSatellite() + "]";
	}

}
