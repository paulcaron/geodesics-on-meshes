import Jcg.polyhedron.*;
import Jcg.geometry.*;

public class Window implements Comparable<Window> {
	private double
			begin,	/* begin of window */
			end,	/* end of window */
			distBegin, /* distance from begin to (pseudo)source */
			distEnd, /* distance from end to (pseudo)source */
			distSource; /* distance from source to pseudosource */

	private Halfedge<Point_3> halfedge; /* halfedge corresponding to window */

 	/*
 	 * Is the path from this window towards
 	 * the (pseudo)source passing through the
 	 * face containing this halfedge?
 	 */
	private boolean side;

	public Window(
			double begin,
			double end,
			double distBegin,
			double distEnd,
			double distSource,
			Halfedge<Point_3> halfedge,
			boolean side
	) {
		this.begin = begin;
		this.end = end;
		this.distBegin = distBegin;
		this.distEnd = distEnd;
		this.distSource = distSource;
		this.halfedge = halfedge;
		this.side = side;
	}

	public double minimumDistance() {
		if (isAcuteTriangle())
			return distSource + height();
		else
			return distSource + Math.min(distBegin, distEnd);
	}

	private boolean isAcuteTriangle() {
		return true;
	}

	private double height() {
		return 0;
	}

	public int compareTo(Window other) {
		double thisDist = this.minimumDistance(),
			   otherDist = other.minimumDistance();
		if (thisDist < otherDist)
			return -1;
		else if (thisDist > otherDist)
			return 1;
		else
			return 0;
	}
}
