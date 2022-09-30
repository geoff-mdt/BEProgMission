package progmission;

import fr.cnes.sirius.patrius.assembly.Assembly;
import fr.cnes.sirius.patrius.assembly.AssemblyBuilder;
import fr.cnes.sirius.patrius.assembly.properties.GeometricProperty;
import fr.cnes.sirius.patrius.assembly.properties.InertiaParallelepipedProperty;
import fr.cnes.sirius.patrius.assembly.properties.MassProperty;
import fr.cnes.sirius.patrius.assembly.properties.SensorProperty;
import fr.cnes.sirius.patrius.attitudes.Attitude;
import fr.cnes.sirius.patrius.attitudes.AttitudeLaw;
import fr.cnes.sirius.patrius.attitudes.BodyCenterGroundPointing;
import fr.cnes.sirius.patrius.fieldsofview.CircularField;
import fr.cnes.sirius.patrius.math.analysis.interpolation.LinearInterpolator;
import fr.cnes.sirius.patrius.math.analysis.polynomials.PolynomialSplineFunction;
import fr.cnes.sirius.patrius.math.complex.Quaternion;
import fr.cnes.sirius.patrius.math.geometry.euclidean.threed.Parallelepiped;
import fr.cnes.sirius.patrius.math.geometry.euclidean.threed.Rotation;
import fr.cnes.sirius.patrius.math.geometry.euclidean.threed.Vector3D;
import fr.cnes.sirius.patrius.math.util.FastMath;
import fr.cnes.sirius.patrius.orbits.Orbit;
import fr.cnes.sirius.patrius.propagation.SpacecraftState;
import fr.cnes.sirius.patrius.propagation.analytical.KeplerianPropagator;
import fr.cnes.sirius.patrius.utils.exception.PatriusException;
import utils.ConstantsBE;

/**
 * [DO NOT MODIFY THIS CLASS]
 * 
 * Satellite representation.
 * 
 * @author lherbert
 */
public class Satellite {

	/** Main body part name for the {@link Assembly} of the satellite. */
	public static final String MAIN_BODY_NAME = "main body";

	/** Sensor part name for the {@link Assembly} of the satellite. */
	public static final String SENSOR_NAME = "sensor";

	/** Solar panels part name for the {@link Assembly} of the satellite. */
	public static final String SOLAR_PANNELS_NAME = "solar pannels";

	/** Satellite name. */
	private final String name;

	/** Mission of the satellite. */
	private final SimpleMission mission;

	/** Satellite {@link Orbit}. */
	private final Orbit initialOrbit;

	/** Default {@link AttitudeLaw} for our satellite. */
	private final AttitudeLaw defaultAttitudeLaw;

	/** {@link SpacecraftState} of our satellite. */
	private final SpacecraftState spacecraftState;

	/**
	 * Propagator for our {@link Satellite}. Here we use a simple
	 * {@link KeplerianPropagator}
	 */
	private final KeplerianPropagator propagator;

	/** Origin point in the local frame of the satellite. */
	private final Vector3D frameOrigin;

	/** X axis in the local frame of the satellite. */
	private final Vector3D frameXAxis;

	/** Y axis in the local frame of the satellite. */
	private final Vector3D frameYAxis;

	/** Sensor axis in the local frame of the satellite. */
	private final Vector3D sensorAxis;

	/** {@link Assembly} object to encapsulate the satellite properties. */
	private final Assembly assembly;

	/**
	 * {@link PolynomialSplineFunction} used to compute the slews durations for our
	 * satellite, using interpolation. Here we assume that the initial and final
	 * rotation speeds are null.
	 */
	private final PolynomialSplineFunction slewCalculator;

	/**
	 * Maximum duration of a slew between two observations for our satellite in an
	 * observation context (the slew correspond to the rotation from the final
	 * Attitude of the first obs to the initial Attitude of the second obs). This
	 * field is calculated using the slewCalculator.
	 */
	private final double maxSlewDuration;

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * 
	 * Build the {@link Satellite} instance.
	 *
	 * @param mission Mission of the satellite
	 * @param name    Name of the satellite
	 * @param orbit   Orbit of the satellite
	 * @throws PatriusException if the Attitude computing fails
	 */
	public Satellite(final SimpleMission mission, final String name, final Orbit orbit) throws PatriusException {

		// Instantiating the context for our satellite : mission, name, orbit
		this.mission = mission;
		this.name = name;
		this.initialOrbit = orbit;

		// Instantiating the slew calculator for our satellite
		// This calculator simply interpolates the duration of the slew as a function of
		// the rotation angle that is needed to be performed. For that, we use the
		// Pleiades agility table in the ConstantsBE class.
		final LinearInterpolator interpolator = new LinearInterpolator();
		this.slewCalculator = interpolator.interpolate(ConstantsBE.POINTING_AGILITY_ROTATIONS,
				ConstantsBE.POINTING_AGILITY_DURATIONS);

		// Computing the maximum slew duration = duration the slew for a rotation of two
		// times the maximum pointing capacity. The satellite will never have to perform
		// a slew bigger than this one.
		this.maxSlewDuration = this.computeMaxSlewDuration();

		// Building an Assembly for our satellite (Patrius basic code)

		// Naming different parts
		final String mainBody = MAIN_BODY_NAME;
		final String sensor = SENSOR_NAME;

		// Creating an AssemblyBuilder for our satellite
		final AssemblyBuilder builder = new AssemblyBuilder();
		builder.addMainPart(mainBody);

		// Main vectors for our satellite shape and inertia
		this.frameOrigin = Vector3D.ZERO;
		this.frameXAxis = Vector3D.PLUS_I;
		this.frameYAxis = Vector3D.PLUS_J;
		// ZAxis is the vector product of YAxis and XAxis

		// Cubic shape
		final Parallelepiped cube = new Parallelepiped(this.frameOrigin, this.frameXAxis, this.frameYAxis, 1, 1, 1);

		// Main part properties
		final GeometricProperty geom = new GeometricProperty(cube);
		builder.addProperty(geom, mainBody);
		final MassProperty massProp = new MassProperty(ConstantsBE.SPACECRAFT_MASS);
		final InertiaParallelepipedProperty inertia = new InertiaParallelepipedProperty(1., 1., 1., massProp);

		// Adding the inertia property to the main part
		builder.addProperty(inertia, mainBody);

		// Creating the sensor
		// Here we specify the name, the position in the main body Frame, and the
		// rotation with respect to the main body Frame (satellite attached local Frame)
		builder.addPart(sensor, mainBody, new Vector3D(0, 0, -1), Rotation.IDENTITY);

		// Adding the sensor main axis = main direction of the sensor
		// Axis of the sensor line of sight = -z in the satellite local Frame
		this.sensorAxis = Vector3D.MINUS_K;
		final SensorProperty sensorProperty = new SensorProperty(this.sensorAxis);

		// Setting the satellite main field of view as a circular field of view
		// Here this fov is the global visibility of the satellite including its
		// agility, not the actual sensor fov which is smaller
		sensorProperty.setMainFieldOfView(new CircularField("SatelliteGlobalFoV",
				FastMath.toRadians(ConstantsBE.POINTING_CAPACITY), sensorProperty.getInSightAxis()));
		builder.addProperty(sensorProperty, sensor);

		// Default attitude law will be pointing the mass center of the Earth with the
		// sensor line of sight. It's the Nadir pointing law.
		this.defaultAttitudeLaw = new BodyCenterGroundPointing(mission.getEarth(), this.getSensorAxis(),
				this.getFrameXAxis());

		// Creating the propagator = the object which computes the coordinates
		// (position, velocity and attitude) of our satellite using its orbit and
		// attitude law.
		this.propagator = new KeplerianPropagator(this.initialOrbit, this.defaultAttitudeLaw);

		// Initial state of our satellite (position and velocity + attitude)
		this.spacecraftState = new SpacecraftState(orbit, this.defaultAttitudeLaw.getAttitude(orbit));

		// Assigning the main part frame to the SpacecraftState
		builder.initMainPartFrame(this.spacecraftState);

		// Finally building the satellite assembly
		this.assembly = builder.returnAssembly();
	}

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * 
	 * Returns the maximum slew duration of our {@link Satellite} : it corresponds
	 * by assumption to the duration to rotate the satellite with a rotation of 2
	 * times the maximum pointing capacity.
	 * 
	 * @return the maximum slew duration
	 */
	public double computeMaxSlewDuration() {
		return this.computeSlewDuration(2 * ConstantsBE.POINTING_CAPACITY);
	}

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * 
	 * Computes the duration of a slew given the input rotation to be achieved in
	 * degrees.
	 * <p>
	 * Here we assume that the initial and final rotation speeds are null.
	 * </p>
	 *
	 * @param rotationAngle Rotation angle (degrees).
	 * @return the duration of the slew (seconds)
	 */
	public double computeSlewDuration(final double rotationAngle) {
		return this.slewCalculator.value(rotationAngle);
	}

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * 
	 * Compute the duration of a slew between att1 and att2.
	 * <p>
	 * Here we assume that the initial and final rotation speeds are null.
	 * </p>
	 * 
	 * @param att1 Initial {@link Attitude}
	 * @param att2 Final {@link Attitude}
	 * @return the duration of the slew (seconds)
	 * @throws PatriusException if withReferenceFrame fails
	 */
	public double computeSlewDuration(final Attitude att1, final Attitude att2) throws PatriusException {
		// We will call computeSlewDuration so we need the angle between our two
		// attitudes
		// Computation of the angle between our two attitudes
		final Quaternion Q1 = att1.withReferenceFrame(this.mission.getEme2000()).getRotation().getQuaternion();
		final Rotation rot1 = new Rotation(true, Q1);
		final Quaternion Q2 = att2.withReferenceFrame(this.mission.getEme2000()).getRotation().getQuaternion();
		final Rotation rot2 = new Rotation(true, Q2);
		final double rotationAngle = Rotation.distance(rot1, rot2);
		return this.computeSlewDuration(FastMath.toDegrees(rotationAngle));
	}

	/**
	 * @return the name
	 */
	public String getName() {
		return this.name;
	}

	/**
	 * @return the mission
	 */
	public SimpleMission getMission() {
		return this.mission;
	}

	/**
	 * @return the orbit
	 */
	public Orbit getInitialOrbit() {
		return this.initialOrbit;
	}

	/**
	 * @return the defaultAttitudeLaw
	 */
	public AttitudeLaw getDefaultAttitudeLaw() {
		return this.defaultAttitudeLaw;
	}

	/**
	 * @return the propagator
	 */
	public KeplerianPropagator getPropagator() {
		return this.propagator;
	}

	/**
	 * @return the spacecraftState
	 */
	public SpacecraftState getSpacecraftState() {
		return this.spacecraftState;
	}

	/**
	 * @return the assembly
	 */
	public Assembly getAssembly() {
		return this.assembly;
	}

	/**
	 * @return the frameOrigin
	 */
	public Vector3D getFrameOrigin() {
		return this.frameOrigin;
	}

	/**
	 * @return the frameXAxis
	 */
	public Vector3D getFrameXAxis() {
		return this.frameXAxis;
	}

	/**
	 * @return the frameYAxis
	 */
	public Vector3D getFrameYAxis() {
		return this.frameYAxis;
	}

	/**
	 * @return the sensorAxis
	 */
	public Vector3D getSensorAxis() {
		return this.sensorAxis;
	}

	/**
	 * @return the slewCalculator
	 */
	public PolynomialSplineFunction getSlewCalculator() {
		return this.slewCalculator;
	}

	/**
	 * @return the maxSlewDuration
	 */
	public double getMaxSlewDuration() {
		return this.maxSlewDuration;
	}

	/**
	 * Returns a string representation of the satellite.
	 *
	 * @return a string representation of the satellite
	 */
	@Override
	public String toString() {
		return this.name;
	}
}
