package reader;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.opencsv.CSVReader;
import com.opencsv.exceptions.CsvValidationException;

import fr.cnes.sirius.patrius.bodies.GeodeticPoint;
import fr.cnes.sirius.patrius.math.util.MathLib;

/**
 * Read the reference sites.
 *
 * @author bonitt
 */
public class SitesReader {

    /** Path of the file containing our sites to read. */
	public static final String OBSERVATION_SITES_FILE = "src/test/resources/sites.csv";

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * 
     * Reads the input file containing the target sites for observation. The input
     * file must be a csv respecting particular format constraints.
     *
     * @param filename Relative path of the input file to read
     * @return sites data as a {@link List} of {@link Site}
     * @throws CsvValidationException
     *         if the csv is not valid
     * @throws NumberFormatException
     *         if an NumberFormatException occurs
     * @throws IOException
     *         if an IOException occurs
     */
    public static List<Site> readSites(final String filename)
			throws CsvValidationException, NumberFormatException, IOException {

        final List<Site> siteList = new ArrayList<>(100);
		final CSVReader reader = new CSVReader(new FileReader(OBSERVATION_SITES_FILE));
		String[] lineInArray;
		while ((lineInArray = reader.readNext()) != null) {
			if (!lineInArray[0].contains("ID")) {
				final String[] tab = lineInArray[0].trim().split(";");

				final String name = tab[2];
				final double score = Double.parseDouble(tab[1]);
				final double latitude = MathLib.toRadians(Double.parseDouble(tab[7]));
				final double longitude = MathLib.toRadians(Double.parseDouble(tab[6]));
				final double altitude = Double.parseDouble(tab[8]);
				final GeodeticPoint point = new GeodeticPoint(latitude, longitude, altitude);
				siteList.add(new Site(name, score, point));
			}
		}
		return siteList;
	}

	/**
	 * [DO NOT MODIFY THIS METHOD]
	 * 
     * Read the reference sites (test application).
     *
     * @param args
     *        Nothing expected
     * @throws CsvValidationException
     *         if the csv is not valid
     * @throws NumberFormatException
     *         if an NumberFormatException occurs
     * @throws IOException
     *         if an IOException occurs
     */
    public static void main(final String[] args) throws CsvValidationException, NumberFormatException, IOException {

        final List<Site> siteList = readSites(OBSERVATION_SITES_FILE);
		System.out.println("Loaded sites: " + siteList.size());
	}
}
