package progmission;

import fr.cnes.sirius.patrius.time.AbsoluteDate;
import reader.Site;

public class Reservation {
	private AbsoluteDate startDate;
	private AbsoluteDate endDate;
	private Site site;

	public Reservation(AbsoluteDate startDate, Site site) {
		this.startDate = startDate;
		this.endDate = startDate.shiftedBy(10);
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
}
