import Jcg.polyhedron.*;
import Jcg.geometry.*;

public class Window implements Comparable<Window> {
	private double
			beginning,	/* begining of window (b0 in the article) */
			end,	/* end of window (b1 in the article) */
			distBeginning, /* distance from beginning to (pseudo)source (d0 in the article) */
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
			double beginning,
			double end,
			double distBeginning,
			double distEnd,
			double distSource,
			Halfedge<Point_3> halfedge,
			boolean side
	) {
		this.beginning = beginning;
		this.end = end;
		this.distBeginning = distBeginning;
		this.distEnd = distEnd;
		this.distSource = distSource;
		this.halfedge = halfedge;
		this.side = side;
	}
	
	public double getBeginning() {
		return this.beginning;
	}

	public double getEnd() {
		return this.end;
	}

	public double getDistBeginning() {
		return this.distBeginning;
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
		if (areBaseAnglesAcute())
			return distSource + GeoUtils.getTriangleHeight(end - beginning, distBeginning, distEnd);
		else
			return distSource + Math.min(distBeginning, distEnd);
	}

	private boolean areBaseAnglesAcute() {
		return !((end - beginning) * (end - beginning) + distBeginning * distBeginning > distSource * distSource ||
				(end - beginning) * (end - beginning) + distSource * distSource > distBeginning * distBeginning);
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
