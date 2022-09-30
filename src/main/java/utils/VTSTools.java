package utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.TimeZone;

import fr.cnes.sirius.patrius.attitudes.Attitude;
import fr.cnes.sirius.patrius.attitudes.AttitudeLeg;
import fr.cnes.sirius.patrius.attitudes.StrictAttitudeLegsSequence;
import fr.cnes.sirius.patrius.frames.Frame;
import fr.cnes.sirius.patrius.frames.FramesFactory;
import fr.cnes.sirius.patrius.math.geometry.euclidean.threed.Rotation;
import fr.cnes.sirius.patrius.math.geometry.euclidean.threed.Vector3D;
import fr.cnes.sirius.patrius.math.util.FastMath;
import fr.cnes.sirius.patrius.propagation.Propagator;
import fr.cnes.sirius.patrius.time.AbsoluteDate;
import fr.cnes.sirius.patrius.time.TimeScale;
import fr.cnes.sirius.patrius.time.TimeScalesFactory;
import fr.cnes.sirius.patrius.utils.Constants;
import fr.cnes.sirius.patrius.utils.exception.PatriusException;
import reader.Site;

/**
 * [DO NOT MODIFY THIS CLASS]
 * 
 * This class develops method to write VTS output data.
 * 
 * @author lherbert
 */
public class VTSTools {

	/** Column separator. */
	private static final String SEP = " ";

	/** Useful strings. */
	private static final String NL = "\n";

	/** Unit conversion. */
	private static final double M_TO_KM = 1e-3;

	/** Time step in seconds for OEM files (PV ephemeris). */
	private static final double STEP_OEM = 60;

	/** Time step in seconds for AEM files (Attitude ephemeris). */
	private static final double STEP_AEM = 1;

	/** Reference date. */
	private static final AbsoluteDate REFERENCE_DATE = AbsoluteDate.MODIFIED_JULIAN_EPOCH;

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * 
	 * Creates a CIC MPM POI file which declare all the positions of the sites of
	 * interest for the VTS Simulation.
	 * 
	 * @param path     Path of the file
	 * @param siteList List of the {@link Site} objects to be written.
	 */
	public static void generatePOIFile(final String path, final List<Site> siteList) {

		// Create a new file
		final File file = new File(path);
		file.delete();

		try {

			// utc timescale
			final TimeScale utc = TimeScalesFactory.getUTC();

			// new file writer
			final FileWriter fw = new FileWriter(file);

			// header
			fw.append("CIC_MPM_VERS = 1.0");
			fw.append(NL);
			fw.append("CREATION_DATE = " + getCurrentDate(utc).toString(TimeScalesFactory.getUTC()));
			fw.append(NL);
			fw.append("ORIGINATOR = PATRIUS");
			fw.append(NL);
			fw.append("     ");
			fw.append(NL);
			fw.append("META_START");
			fw.append(NL);
			fw.append("");
			fw.append(NL);
			fw.append("USER_DEFINED_PROTOCOL = NONE");
			fw.append(NL);
			fw.append("USER_DEFINED_CONTENT = POINTS OF INTEREST");
			fw.append(NL);
			fw.append("USER_DEFINED_SIZE = 3");
			fw.append(NL);
			fw.append("USER_DEFINED_TYPE = STRING");
			fw.append(NL);
			fw.append("USER_DEFINED_UNIT = [n/a]");
			fw.append(NL);
			fw.append("META_STOP");
			fw.append(NL);

			// Sites writing process
			StringBuffer bf;
			for (Site site : siteList) {
				// write to file
				bf = new StringBuffer();
				bf.append(FastMath.toDegrees(site.getPoint().getLatitude()));
				bf.append(SEP);
				bf.append(FastMath.toDegrees(site.getPoint().getLongitude()));
				bf.append(SEP);
				bf.append("\"" + site.getName() + "\"");

				fw.append(bf.toString());
				fw.append(NL);

			}

			// close file writer
			fw.close();

		} catch (final IOException e) {
			e.printStackTrace();
		} catch (final PatriusException e) {
			e.printStackTrace();
		}
		System.out.println("VTS POI file created : " + path);
	}

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * 
	 * Create a position ephemeris file.
	 * 
	 * @param path       Path of the file
	 * @param start      Start date
	 * @param end        End date
	 * @param propagator Propagator
	 */
	public static void generateOEMFile(final String path, final AbsoluteDate start, final AbsoluteDate end,
			final Propagator propagator) {

		// Make sure the DEFAULT_DIR points to the classes
		// If you are running test classes, make sure it points to test-classes!
		// See static { ... }

		// create a new file
		final File file = new File(path);
		file.delete();

		try {

			// utc timescale
			final TimeScale utc = TimeScalesFactory.getUTC();
			final Frame eme2000 = FramesFactory.getEME2000();

			// new file writer
			final FileWriter fw = new FileWriter(file);

			// header
			fw.append("CIC_OEM_VERS = 2.0");
			fw.append(NL);
			fw.append("CREATION_DATE = " + getCurrentDate(utc).toString(TimeScalesFactory.getUTC()));
			fw.append(NL);
			fw.append("ORIGINATOR = PATRIUS");
			fw.append(NL);
			fw.append("     ");
			fw.append(NL);
			fw.append("META_START");
			fw.append(NL);
			fw.append("");
			fw.append(NL);
			fw.append("OBJECT_NAME = PATRIUS SATELLITE");
			fw.append(NL);
			fw.append("OBJECT_ID = PAT001");
			fw.append(NL);
			fw.append("CENTER_NAME = EARTH");
			fw.append(NL);
			fw.append("REF_FRAME = EME2000");
			fw.append(NL);
			fw.append("TIME_SYSTEM = UTC");
			fw.append(NL);
			fw.append("");
			fw.append(NL);
			fw.append("META_STOP");
			fw.append(NL);

			// ephemeris
			AbsoluteDate current = start.shiftedBy(0.1);
			Vector3D position;
			double seconds;
			int days;
			StringBuffer bf;
			while (current.offsetFrom(end, utc) < 0) {
				// elapsed time
				seconds = current.offsetFrom(REFERENCE_DATE, utc) + utc.offsetFromTAI(current);
				days = (int) (seconds / Constants.JULIAN_DAY);
				seconds = seconds - days * Constants.JULIAN_DAY;

				// position at time in EM2000
				position = propagator.getPVCoordinates(current, eme2000).getPosition();

				// write to file
				bf = new StringBuffer();
				bf.append(days);
				bf.append(SEP);
				bf.append(seconds);
				bf.append(SEP);
				bf.append(position.getX() * M_TO_KM);
				bf.append(SEP);
				bf.append(position.getY() * M_TO_KM);
				bf.append(SEP);
				bf.append(position.getZ() * M_TO_KM);

				fw.append(bf.toString());
				fw.append(NL);

				// move forward
				current = current.shiftedBy(STEP_OEM);
			}

			// close file writer
			fw.close();

		} catch (final IOException e) {
			e.printStackTrace();
		} catch (final PatriusException e) {
			e.printStackTrace();
		}
		System.out.println("VTS OEM file created : " + path);
	}

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * 
	 * Create an attitude ephemeris file.
	 * 
	 * @param path       Path of the file
	 * @param start      Start date
	 * @param end        End date
	 * @param propagator Propagator
	 */
	public static void generateAEMFile(final String path, final AbsoluteDate start, final AbsoluteDate end,
			final Propagator propagator) {

		// create a new file
		final File file = new File(path);

		try {

			// utc timescale
			final TimeScale utc = TimeScalesFactory.getUTC();

			// new file writer
			final FileWriter fw = new FileWriter(file);

			// header
			fw.append("CIC_AEM_VERS = 1.0");
			fw.append(NL);
			fw.append("CREATION_DATE = " + getCurrentDate(utc).toString(TimeScalesFactory.getUTC()));
			fw.append(NL);
			fw.append("ORIGINATOR = PATRIUS");
			fw.append(NL);
			fw.append("     ");
			fw.append(NL);
			fw.append("META_START");
			fw.append(NL);
			fw.append("");
			fw.append(NL);
			fw.append("OBJECT_NAME = PATRIUS SATELLITE");
			fw.append(NL);
			fw.append("OBJECT_ID = PAT001");
			fw.append(NL);
			fw.append("CENTER_NAME = EARTH");
			fw.append(NL);
			fw.append("REF_FRAME_A = EME2000");
			fw.append(NL);
			fw.append("REF_FRAME_B = SC_BODY_1");
			fw.append(NL);
			fw.append("ATTITUDE_DIR = A2B");
			fw.append(NL);
			fw.append("TIME_SYSTEM = UTC");
			fw.append(NL);
			fw.append("ATTITUDE_TYPE = QUATERNION");
			fw.append(NL);
			fw.append("");
			fw.append(NL);
			fw.append("META_STOP");
			fw.append(NL);

			// ephemeris
			AbsoluteDate current = start.shiftedBy(0.1);
			Attitude attitude;
			Rotation rot;
			double seconds;
			int days;
			StringBuffer bf;
			while (current.offsetFrom(end, utc) < 0) {
				// elapsed time
				seconds = current.offsetFrom(REFERENCE_DATE, utc) + utc.offsetFromTAI(current);
				days = (int) (seconds / Constants.JULIAN_DAY);
				seconds = seconds - days * Constants.JULIAN_DAY;

				// position at time in EM2000
				attitude = propagator.propagate(current).getAttitude();
				rot = attitude.getRotation().revert();

				// write to file
				bf = new StringBuffer();
				bf.append(days);
				bf.append(SEP);
				bf.append(seconds);
				bf.append(SEP);
				bf.append(-rot.getQuaternion().getQ0());
				bf.append(SEP);
				bf.append(rot.getQuaternion().getQ1());
				bf.append(SEP);
				bf.append(rot.getQuaternion().getQ2());
				bf.append(SEP);
				bf.append(rot.getQuaternion().getQ3());

				fw.append(bf.toString());
				fw.append(NL);

				// move forward
				// current = current.shiftedBy(STEP);
				current = current.shiftedBy(STEP_AEM);
			}

			// close file writer
			fw.close();

		} catch (final IOException e) {
			e.printStackTrace();
		} catch (final PatriusException e) {
			e.printStackTrace();
		}
		System.out.println("VTS AEM file created : " + path);
	}

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * 
	 * Writes a CIC MEM file to visualize the input
	 * {@link StrictAttitudeLegsSequence} in VTS.<br>
	 * Only the start and end of each leg and the leg's nature will be displayed.
	 * 
	 * @param path     Path of the output file
	 * @param sequence Input {@link StrictAttitudeLegsSequence} to be written
	 */
	public static void generateLegSequenceMEMFile(final String path,
			final StrictAttitudeLegsSequence<AttitudeLeg> sequence) {

		// Make sure the DEFAULT_DIR points to the classes
		// If you are running test classes, make sure it points to test-classes!
		// See static { ... }

		// create a new file
		final File file = new File(path);
		file.delete();

		try {

			// utc timescale
			final TimeScale utc = TimeScalesFactory.getUTC();

			// new file writer
			final FileWriter fw = new FileWriter(file);

			// header
			fw.append("CIC_MEM_VERS = 1.0");
			fw.append(NL);
			fw.append("CREATION_DATE = " + getCurrentDate(utc).toString(utc));
			fw.append(NL);
			fw.append("ORIGINATOR = PATRIUS");
			fw.append(NL);
			fw.append("     ");
			fw.append(NL);
			fw.append("META_START");
			fw.append(NL);
			fw.append("COMMENT SATELLITE POINTING MODES");
			fw.append(NL);
			fw.append("OBJECT_NAME = PATRIUS SATELLITE");
			fw.append(NL);
			fw.append("OBJECT_ID = PAT001");
			fw.append(NL);
			fw.append("USER_DEFINED_PROTOCOL	=	NONE");
			fw.append(NL);
			fw.append("USER_DEFINED_CONTENT	=	OEF");
			fw.append(NL);
			fw.append("USER_DEFINED_SIZE = 1");
			fw.append(NL);
			fw.append("USER_DEFINED_TYPE = STRING");
			fw.append(NL);
			fw.append("USER_DEFINED_UNIT = [n/a]");
			fw.append(NL);
			fw.append("TIME_SYSTEM = UTC");
			fw.append(NL);
			fw.append("");
			fw.append(NL);
			fw.append("META_STOP");
			fw.append(NL);

			double seconds;
			int days;
			StringBuffer bf;
			for (final AttitudeLeg leg : sequence) {
				// Writing the event

				// Getting the start date
				final AbsoluteDate start = leg.getDate();
				seconds = start.offsetFrom(REFERENCE_DATE, utc) + utc.offsetFromTAI(start);
				days = (int) (seconds / Constants.JULIAN_DAY);
				seconds = seconds - days * Constants.JULIAN_DAY;

				// Building the line for the start if the leg
				bf = new StringBuffer();
				bf.append(days);
				bf.append(SEP);
				bf.append(seconds);
				bf.append(SEP);
				bf.append("\"Start of " + leg.getNature() + "\"");
				bf.append(SEP);
				fw.append(bf.toString());
				fw.append(NL);

				// Getting the end date
				final AbsoluteDate end = leg.getEnd();
				seconds = end.offsetFrom(REFERENCE_DATE, utc) + utc.offsetFromTAI(end);
				days = (int) (seconds / Constants.JULIAN_DAY);
				seconds = seconds - days * Constants.JULIAN_DAY;

				// Building the line for the end if the leg
				bf = new StringBuffer();
				bf.append(days);
				bf.append(SEP);
				bf.append(seconds);
				bf.append(SEP);
				bf.append("\"End of " + leg.getNature() + "\"");
				bf.append(SEP);
				fw.append(bf.toString());
				fw.append(NL);

			}

			// close file writer
			fw.close();

		} catch (final IOException e) {
			e.printStackTrace();
		} catch (final PatriusException e) {
			e.printStackTrace();
		}
		System.out.println("VTS MEM file created : " + path);
	}

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * 
	 * Getter for the current date.
	 * 
	 * @return the current date
	 */
	@SuppressWarnings("deprecation")
	private static AbsoluteDate getCurrentDate(final TimeScale scale) {
		// LocalDateTime
		final GregorianCalendar cal = new GregorianCalendar(TimeZone.getTimeZone("GMT+00"));
		return new AbsoluteDate(cal.getTime(), scale);
	}
}
