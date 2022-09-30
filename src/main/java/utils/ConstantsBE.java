package utils;

import java.io.File;

import fr.cnes.sirius.patrius.time.AbsoluteDate;
import fr.cnes.sirius.patrius.time.TimeScale;
import fr.cnes.sirius.patrius.time.TimeScalesFactory;

/**
 * [DO NOT MODIFY THIS CLASS]
 * 
 * Set of useful physical constants and configuration constants.
 *
 * @author bonitt
 */
public interface ConstantsBE {

	// Configuration of the VTS path for the project

	/** Path of the VTS output directory. */
	public static final String PATH_VTS_DIRECTORY = ".." + File.separator + "Vts-WindowsNT-32bits-3.5.1"
			+ File.separator + "Data" + File.separator + "BESupaero" + File.separator + "Data";

	// Mission dates constants

	/** Reference time scale used for dates. */
	public static final TimeScale TIME_SCALE = TimeScalesFactory.getTAI();

	/** Initial analysis date. */
	public static final AbsoluteDate START_DATE = new AbsoluteDate("2000-09-01T00:00:00", TIME_SCALE);

	/** Final analysis date. */
	public static final AbsoluteDate END_DATE = new AbsoluteDate("2000-09-08T00:00:00", TIME_SCALE);

	// Orbital parameters

	/** Semi-major axis (radius in our case) of the orbit (m). */
	public static final double ALTITUDE = 694000;

	/** Inclination of the orbit (degrees). */
	public static final double INCLINATION = 98.2;

	/** Eccentricity (degrees). */
	public static final double MEAN_ECCENTRICITY = 0.0001153;

	/**
	 * RAAN of our Pleiades satellite's orbit. Source : http://celestrak.org/NORAD.
	 */
	public static final double ASCENDING_NODE_LONGITUDE = 267.6;

	// Algorithm configuration constants

	// Agility of the satellite
	/**
	 * Maximum pointing capacity of the satellite (degrees from nadir, cumulating
	 * roll and pitch).
	 *
	 * Justification : Each Pleiades satellite is capable of targeting images along
	 * any ground direction within 47° of vertical viewing position.
	 */
	public static final double POINTING_CAPACITY = 47.0;

	/**
	 * Agility durations to help compute the duration of a slew using interpolation.
	 * Agility - Roll or pitch: 5° in 8 seconds, 10° in 10 seconds, 60° in 25
	 * seconds. The last value has been added to increase the interval of definition
	 * of the interpolation function.
	 */
	public static final double[] POINTING_AGILITY_DURATIONS = { 0.0, 8.0, 10.0, 25.0, 150.0 };

	/**
	 * Agility rotations to help compute the duration of a slew using interpolation.
	 * Agility - Roll or pitch: 5° in 8 seconds, 10° in 10 seconds, 60° in 25
	 * seconds. The last value has been added to increase the interval of definition
	 * of the interpolation function.
	 */
	public static final double[] POINTING_AGILITY_ROTATIONS = { 0, 5, 10, 60, 360 };

	// Spacecraft properties

	/** Mass of the spacecraft (kg). */
	public static final double SPACECRAFT_MASS = 1000;

	// Observation constraints

	/**
	 * Maximum phase angle to prevent the detector's blinding with the Sun
	 * (degrees). Phase angle is the angle between the sun, the target and the
	 * detector.
	 */
	public static final double MAX_SUN_PHASE_ANGLE = 90;

	/**
	 * Maximum sun incidence angle in degrees. The solar incidence angle is the
	 * angle between the local normal and the sun. We prefer low solar incidence
	 * angles for illumination constraints (better SNR).
	 */
	public static final double MAX_SUN_INCIDENCE_ANGLE = 75;

	// Detector constants

	/**
	 * Here we make the assumption that we have a matrix CCD detector which takes
	 * observations by taking a photography of the target with a given time of
	 * integration. This assumption will simplify the attitude pointing law because
	 * we only need to point the location to image during that time and it's all.
	 * (No push broom pointing strategy and no surface imaging).
	 */
	public static final double INTEGRATION_TIME = 10.0;

}
