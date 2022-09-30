package reader;

import fr.cnes.sirius.patrius.bodies.GeodeticPoint;
import fr.cnes.sirius.patrius.math.util.FastMath;

/**
 * Site description.
 *
 * @author bonitt
 */
public class Site implements Comparable<Site> {

    /** Site name. */
    private final String name;

    /** Site score. */
    private final double score;

    /** Geodetic point. */
    private final GeodeticPoint point;

    /**
     * [DO NOT MODIFY THIS METHOD]
     * 
     * Site constructor.
     * 
     * @param name
     *        Site name
     * @param score
     *        Site score
     * @param point
     *        Geodetic point
     */
    public Site(final String name, final double score, final GeodeticPoint point) {
        this.name = name;
        this.score = score;
        this.point = point;
    }

    /**
     * Getter for the site name.
     *
     * @return the site name
     */
    public String getName() {
        return this.name;
    }

    /**
     * Getter for the site score.
     *
     * @return the site score
     */
    public double getScore() {
        return this.score;
    }

    /**
     * Getter for the geodetic point.
     *
     * @return the geodetic point
     */
    public GeodeticPoint getPoint() {
        return this.point;
    }

    /**
     * Returns a string representation of the site.
     *
     * @return a string representation of the site
     */
    @Override
    public String toString() {
        // return "{" + this.name + " => " + this.score + ", location : " + this.getPoint() + "}";
        return this.name;
    }

    /**
     * Compares this site with the specified site for order.
     * 
     * @param site
     *        the site to be compared
     * @return a negative integer, zero, or a positive integer as this site
     *         has a lower, equal, or higher score than the specified site.
     */
    @Override
    public int compareTo(final Site site) {
        // We return the opposite comparison of the site score so that in a sorted
        // Collection of Sites, the highest scores will be first, property that is used
        // in the computation of the observation plan
        return (int) FastMath.signum(site.score - this.score);
    }
}
