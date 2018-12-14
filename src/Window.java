import Jcg.polyhedron.*;
import Jcg.geometry.*;

public class Window implements Comparable<Window> {
	private final double
			start,	/* begining of window (b0 in the article) */
			end,	/* end of window (b1 in the article) */
			distStart, /* distance from start to (pseudo)source (d0 in the article) */
			distEnd, /* distance from end to (pseudo)source (d1 in the article) */
			distSource; /* distance from source to pseudosource (sigma in the article) */

	private final Halfedge<Point_3> halfedge; /* halfedge corresponding to window */

 	/*
 	 * Is the path from this window towards
 	 * the (pseudo)source passing through the
 	 * face containing this halfedge?
 	 */
	private final boolean side; /* (tau in the article) */

	public Window(
			double start,
			double end,
			double distStart,
			double distEnd,
			double distSource,
			Halfedge<Point_3> halfedge,
			boolean side
	) {
		this.start = start;
		this.end = end;
		this.distStart = distStart;
		this.distEnd = distEnd;
		this.distSource = distSource;
		this.halfedge = halfedge;
		this.side = side;
	}

	public Window(Window other) {
		this(other.start,
			other.end,
			other.distStart,
			other.distEnd,
			other.distSource,
			other.halfedge,
			other.side);
	}
	
	public double getStart() {
		return this.start;
	}

	public double getEnd() {
		return this.end;
	}

	public double getDistStart() {
		return this.distStart;
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

	public double getLength() {
		return end - start;
	}
	
	public Halfedge<Point_3> getHalfedge() {
		return this.halfedge;
	}

	public double getMinimumDistanceFromEndpointToSource() {
		if (areBaseAnglesAcute())
			return distSource + GeoUtils.getTriangleHeight(end - start, distStart, distEnd);
		else
			return distSource + Math.min(distStart, distEnd);
	}

	/* 
	 * Calculates and returns the abscissa of the source point
	 * considering a cartesian plan having origin at start,
	 * and x-axis passing through the start and end points of
	 * the window. The result might be negative.
	 */
	public double getAbscisseSourcePlaneProjection() {
		return GeoUtils.getTriangleVertexPlaneProjection(end - start, distStart, distEnd).getX().doubleValue();
	}

	private boolean areBaseAnglesAcute() {
		return !((end - start) * (end - start) + distStart * distStart > distSource * distSource ||
				(end - start) * (end - start) + distSource * distSource > distStart * distStart);
	}

	public int compareTo(Window other) {
		double thisDist = this.getMinimumDistanceFromEndpointToSource(),
			   otherDist = other.getMinimumDistanceFromEndpointToSource();
		if (thisDist < otherDist)
			return -1;
		else if (thisDist > otherDist)
			return 1;
		else
			return 0;
	}
	
	public Point_3 getSource() {
		throw new IllegalArgumentException("Not implemented");
	}

	public Window setStart(double newStart) {
		if (newStart >= end)
			throw new IllegalArgumentException("New start value should be smaller than end");
		else if (newStart < 0)
			throw new IllegalArgumentException("New start value should be non-negative");

		double newDistStart = GeoUtils.getCevianLength(distSource, distEnd, newStart - start, end - newStart);
		return new Window(newStart, end, newDistStart, distEnd, distSource, halfedge, side);
	}

	public Window setEnd(double newEnd) {
		if (newEnd <= start)
			throw new IllegalArgumentException("New end value should be greater than start");
		else if (newEnd > GeoUtils.getHalfedgeLength(this.halfedge))
			throw new IllegalArgumentException("New end value should be non-negative");

		double newDistEnd = GeoUtils.getCevianLength(distSource, distEnd, newEnd - start, end - newEnd);
		return new Window(start, newEnd, distStart, newDistEnd, distSource, halfedge, side);
	}

	public Window getOpposite() {
		double halfedgeLength = GeoUtils.getHalfedgeLength(getHalfedge());
		return new Window(
			halfedgeLength - getStart(),
			halfedgeLength - getEnd(),
			getDistEnd(),
			getDistStart(),
			getDistSource(),
			getHalfedge().getOpposite(),
			!getSide());
	}

}
