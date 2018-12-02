import Jcg.polyhedron.*;
import Jcg.geometry.*;

public class Window implements Comparable<Window> {
	private double
			begin,	/* begining of window (b0 in the article) */
			end,	/* end of window (b1 in the article) */
			distBegin, /* distance from begin to (pseudo)source (d0 in the article) */
			distEnd, /* distance from end to (pseudo)source (d1 in the article) */
			distSource; /* distance from source to pseudosource (sigma in the article) */

	private Halfedge<Point_3> halfedge; /* halfedge corresponding to window */

 	/*
 	 * Is the path from this window towards
 	 * the (pseudo)source passing through the
 	 * face containing this halfedge?
 	 */
	private boolean side; /* (tau in the article) */

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
	
	public double getBegin() {
		return this.begin;
	}

	public double getEnd() {
		return this.end;
	}

	public double getDistBegin() {
		return this.distBegin;
	}

	public double getDistEnd() {
		return this.distEnd;
	}

	public double getDistSource() {
		return this.distSource;
	}

	public boolean getSide() {
		return this.side;
	}
	
	public Halfedge<Point_3> getHalfedge() {
		return this.halfedge;
	}

	public double minimumDistance() {
		if (baseAnglesAcute())
			return distSource + getTriangleHeight();
		else
			return distSource + Math.min(distBegin, distEnd);
	}

	private boolean baseAnglesAcute() {
		return !((end - begin) * (end - begin) + distBegin * distBegin > distSource * distSource ||
				(end - begin) * (end - begin) + distSource * distSource > distBegin * distBegin);
	}

	private double getTriangleHeight() {
		double p = (end - begin + distBegin + distEnd) / 2;
		return 2 * Math.sqrt(p * (p - (end - begin)) * (p - distBegin) * (p - distEnd)) / (end - begin);
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
	
	public Point_3 getSource() {
		return new Point_3(); // to complete
	}
}
