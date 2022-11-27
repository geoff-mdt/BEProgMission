package progmission;

import fr.cnes.sirius.patrius.time.AbsoluteDate;
import reader.Site;
import progmission.Satellite;

import java.util.Comparator;

public class Reservation implements Comparator<Reservation> {
	private AbsoluteDate startDate;
	private AbsoluteDate endDate;
	private Site site;

	public Reservation(AbsoluteDate startDate, Site site, double maxSlew) {
		this.startDate = startDate;
		this.endDate = startDate.shiftedBy(10+maxSlew);
		this.site = site;
	}

	public AbsoluteDate getStartDate() {
		return startDate;
	}

	public AbsoluteDate getEndDate() {
		return endDate;
	}

	public Site getSite() {
		return site;
	}

	@Override
	public int compare(Reservation r1, Reservation r2) {
		return r1.getStartDate().compareTo(r2.getStartDate());
	}
}
