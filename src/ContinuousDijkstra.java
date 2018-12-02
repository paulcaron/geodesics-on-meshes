import java.util.ArrayList;
import java.util.PriorityQueue;
import java.util.HashMap;
import java.lang.NullPointerException;

import Jcg.polyhedron.*;
import Jcg.geometry.*;
import Jama.*;

public class ContinuousDijkstra {
	private final double THRESHOLD = 1e-9;
	private final SurfaceMesh mesh;
	private HashMap<Halfedge, ArrayList<Window>> halfedgeToWindowsList;

	public ContinuousDijkstra(SurfaceMesh mesh) {
		if (mesh == null)
			throw new NullPointerException("Mesh is null");

		this.mesh = mesh;
	}

	public void buildDistances(Point_3 source) {
		if (source == null)
			throw new NullPointerException("Source point is null");

		PriorityQueue<Window> pq = new PriorityQueue<>();

		initializePriorityQueue(source, pq);

		while (!pq.isEmpty()) {
			Window window = pq.poll();

			propagateWindow(window, pq);
		}
	}

	private void initializePriorityQueue(
			Point_3 source,
			PriorityQueue<Window> pq
	) {
		Face<Point_3> face = getFaceContainingPoint(source);

		return;
	}

	private void propagateWindow(Window window, PriorityQueue<Window> pq) {
	}

	private double getLengthHalfedge(Halfedge<Point_3> halfedge) {
		assert halfedge != null;
		Vertex<Point_3> destination = halfedge.getVertex();
		Vertex<Point_3> origin = halfedge.getOpposite().getVertex();

		return (double) destination.getPoint().distanceFrom(origin.getPoint());
	}

	public double getMinDistance(Point_3 destination) {
		if (destination == null)
			throw new NullPointerException("Destination point is null");

		Face<Point_3> faceContainingPoint = getFaceContainingPoint(destination);

		if (isVertex(destination, faceContainingPoint)) {
			return getMinDistanceFromVertex(destination, faceContainingPoint);
		} else if (isHalfedgePoint(destination, faceContainingPoint)) {
			return getMinDistanceFromHalfedge(destination, faceContainingPoint);
		} else {
			return getMinDistanceFromFace(destination, faceContainingPoint);
		}
	}

	/* The destination point is a vertex of face */
	private double getMinDistanceFromVertex(Point_3 destination, Face<Point_3> face) {
		return 0;
	}

	/* The destination point lies on an edge of face */
	private double getMinDistanceFromHalfedge(Point_3 destination, Face<Point_3> face) {
		return 0;
	}

	/* The destination point lies in the interior of face */
	private double getMinDistanceFromFace(Point_3 destination, Face<Point_3> face) {
		double minDistance = Double.MAX_VALUE;

		Halfedge<Point_3> halfedge = face.getEdge();
		for (int i = 0; i < 3; i++) {
			double minDistanceThroughHalfedge = getMinDistance(destination, halfedge);

			if (minDistanceThroughHalfedge < minDistance)
				minDistance = minDistanceThroughHalfedge;

			halfedge = halfedge.getNext();
		}
		assert halfedge == face.getEdge();

		assert positive(minDistance);

		return minDistance;
	}

	private double getMinDistance(Point_3 destination, Halfedge<Point_3> halfedge) {
		double minDistance = Double.MAX_VALUE;
		ArrayList<Window> windowsList = halfedgeToWindowsList.get(halfedge);
		assert windowsList.size() > 0;

		for (int i = 0; i < windowsList.size(); i++) {
			double minDistanceForWindow = getMinDistance(destination, windowsList.get(i));

			if (minDistanceForWindow < minDistance)
				minDistance = minDistanceForWindow;
		}

		assert positive(minDistance);

		return minDistance;
	}

	private double getMinDistance(Point_3 destination, Window window) {
		double minDistance = Double.MAX_VALUE;

		Halfedge<Point_3> halfedge = window.getHalfedge();
		double lengthHalfedge = getLengthHalfedge(halfedge);
		Point_3 halfedgeOrigin = halfedge.getOpposite().getVertex().getPoint();
		Point_3 halfedgeEnd = halfedge.getVertex().getPoint();
		Vector_3 halfedgeUnaryVector = getUnaryVector(new Vector_3(halfedgeOrigin, halfedgeEnd));

		double dx = lengthHalfedge / 100;
		for (double x = window.getBegin(); x <= window.getEnd(); x += dx) {
			Point_3 halfedgePoint = new Point_3(halfedgeOrigin);
			halfedgePoint.translateOf(halfedgeUnaryVector.multiplyByScalar(x));

			double distance = (double) destination.distanceFrom(halfedgePoint) + window.getDistSource() +
				stewart(window.getDistBegin(), window.getDistEnd(), x - window.getBegin(), window.getEnd() - x);

			if (distance < minDistance)
				minDistance = distance;
		}
		assert positive(minDistance);

		return minDistance;
	}

	private Vector_3 getUnaryVector(Vector_3 vector) {
		Vector_3 unaryVector = vector.divisionByScalar(Math.sqrt(vector.squaredLength().doubleValue()));
		assert zero(unaryVector.squaredLength().doubleValue() - 1);
		return unaryVector;
	}

	/*
	 * Following convention employed by Wikipedia to name triangle sides
	 * https://en.wikipedia.org/wiki/Stewart's_theorem
	 */
	private double stewart(double b, double c, double n, double m) {
		double a = n + m;
		double d = Math.sqrt((b * b * m + c * c * n) / a - m * n);
		assert d >= 0;
		return d;
	}

	private boolean isVertex(Point_3 point, Face<Point_3> face) {
		Halfedge<Point_3> halfedge = face.getEdge();
		for (int i = 0; i < 3; i++) {
			if (zero(point.distanceFrom(halfedge.getVertex().getPoint()).doubleValue()))
				return true;

			halfedge = halfedge.getNext();
		}
		assert halfedge == face.getEdge();

		return false;
	}

	private boolean isHalfedgePoint(Point_3 point, Face<Point_3> face) {
		Halfedge<Point_3> halfedge = face.getEdge();
		for (int i = 0; i < 3; i++) {
			if (zero(distancePointHalfedge(point, halfedge)))
				return true;

			halfedge = halfedge.getNext();
		}
		assert halfedge == halfedge.getNext();

		return false;
	}

	private double distancePointHalfedge(Point_3 point, Halfedge<Point_3> halfedge) {
		Point_3 a = halfedge.getVertex().getPoint();
		Point_3 b = halfedge.getOpposite().getVertex().getPoint();

		Vector_3 ab = new Vector_3(a, b);
		Vector_3 ap = new Vector_3(a, point);

		Vector_3 crossProduct = ab.crossProduct(ap);

		return Math.sqrt(crossProduct.squaredLength().doubleValue() / ab.squaredLength().doubleValue());
	}

	private Face<Point_3> getFaceContainingPoint(Point_3 point) {
		for (Face<Point_3> polyhedronFace : mesh.getFaces()) {
			if (insideFace(point, polyhedronFace))
				return polyhedronFace;
		}

		throw new IllegalArgumentException("Point " + point.toString() + " not contained in any face of the mesh");
	}

	private boolean insideFace(Point_3 point, Face<Point_3> face) {
		Halfedge<Point_3> halfedge = face.getEdge();
		Point_3 a = halfedge.getVertex().getPoint();

		halfedge = halfedge.getNext();
		Point_3 b = halfedge.getVertex().getPoint();

		halfedge = halfedge.getNext();
		Point_3 c = halfedge.getVertex().getPoint();

		Vector_3 ab = new Vector_3(a, b);
		Vector_3 ac = new Vector_3(a, c);
		Vector_3 ap = new Vector_3(a, point);

		Matrix matrix = new Matrix(new double[][] {
			{ab.getX().doubleValue(), ab.getY().doubleValue(), ab.getZ().doubleValue()},
			{ac.getX().doubleValue(), ac.getY().doubleValue(), ac.getZ().doubleValue()},
			{ap.getX().doubleValue(), ap.getY().doubleValue(), ap.getZ().doubleValue()}});

		/* Point point is not in the plane containing face */
		if (!zero(matrix.det()))
			return false;

		/*
		 * Tries to solve the system alpha * ab + beta * ac = ap
		 * with 0 <= alpha, beta <= 1 and alpha + beta <= 1
		 */
		matrix = new Matrix(new double[][] {
			{ab.getX().doubleValue(), ac.getX().doubleValue()},
			{ab.getY().doubleValue(), ac.getY().doubleValue()}});
		if (!zero(matrix.det())) {
			Matrix coeffs = matrix.solve(new Matrix (new double[][] {{ap.getX().doubleValue()}, {ap.getY().doubleValue()}}));
			double alpha = coeffs.get(0, 0),
				   beta = coeffs.get(1, 0);
			return positive(alpha) && positive(beta) && negative(alpha + beta - 1);
		}

		matrix = new Matrix(new double[][] {
			{ab.getX().doubleValue(), ac.getX().doubleValue()},
			{ab.getZ().doubleValue(), ac.getZ().doubleValue()}});
		if (!zero(matrix.det())) {
			Matrix coeffs = matrix.solve(new Matrix (new double[][] {{ap.getX().doubleValue()}, {ap.getZ().doubleValue()}}));
			double alpha = coeffs.get(0, 0),
				   beta = coeffs.get(1, 0);
			return positive(alpha) && positive(beta) && negative(alpha + beta - 1);
		}

		matrix = new Matrix(new double[][] {
			{ab.getZ().doubleValue(), ac.getZ().doubleValue()},
			{ab.getY().doubleValue(), ac.getY().doubleValue()}});
		if (!zero(matrix.det())) {
			Matrix coeffs = matrix.solve(new Matrix (new double[][] {{ap.getZ().doubleValue()}, {ap.getY().doubleValue()}}));
			double alpha = coeffs.get(0, 0),
				   beta = coeffs.get(1, 0);
			return positive(alpha) && positive(beta) && negative(alpha + beta - 1);
		}

		return false;
	}

	private boolean zero(double x) {
		return Math.abs(x) <= THRESHOLD;
	}

	private boolean positive(double x) {
		return x > THRESHOLD;
	}

	private boolean negative(double x) {
		return x < -THRESHOLD;
	}

	
	public ArrayList<Double> computePoints(Window window) {
		ArrayList<Double> arr = new ArrayList<Double>(3);
		Point_3 p0, p1, p2, b0, b1;
		Point_3 source = window.computeSource();
		p0 = window.getHalfedge().getVertex().getPoint();
		p1 = window.getHalfedge().getOpposite().getVertex().getPoint();
		p2 = window.getHalfedge().getNext().getVertex().getPoint();
		double lengthHalfedge = getLengthHalfedge(window.getHalfedge());
		boolean p20, p21; //p20=true if p2 is on the left side of the window : 
		Number[] coefficients0 = {1-window.getBegin() / lengthHalfedge, window.getBegin() / lengthHalfedge };
		Number[] coefficients1 = {1-window.getEnd()/ lengthHalfedge, window.getEnd()/ lengthHalfedge};
		Point_3[] points = {p0, p1};
		b0 = Point_3.linearCombination(points, coefficients0);
		b1 = Point_3.linearCombination(points, coefficients1);
		Vector_3 b0s = new Vector_3(b0, source);
		Vector_3 b1s = new Vector_3(b1, source);
		Vector_3 b0p2 = new Vector_3(b0, p2);
		Vector_3 b1p2 = new Vector_3(b1, p2);
		Vector_3 p0p2 = new Vector_3(p0, p2);
		Vector_3 p2p1 = new Vector_3(p2, p1);
		Vector_3 p0s = new Vector_3(p0, source);
		Vector_3 p1s = new Vector_3(p1, source);
		Vector_3 p2s = new Vector_3(p2, source);
		Vector_3 n0 = new Vector_3(-b0s.getY().doubleValue(), b0s.getX(), 0.);
		Vector_3 n1 = new Vector_3(-b1s.getY().doubleValue(), b1s.getX(), 0.);
		p20 = b0p2.innerProduct(n0).doubleValue() > 0;
		p21 = b1p2.innerProduct(n1).doubleValue() > 0;
		if(p20) {
			arr.set(0, -1.);
			arr.set(1, Math.sqrt(p2p1.squaredLength().doubleValue())*p2s.innerProduct(n0).doubleValue() / p2p1.innerProduct(n0).doubleValue());
			arr.set(2, Math.sqrt(p2p1.squaredLength().doubleValue())*p2s.innerProduct(n1).doubleValue() / p2p1.innerProduct(n1).doubleValue());
		}
		else if (p21){
			arr.set(0, Math.sqrt(p0p2.squaredLength().doubleValue())*p0s.innerProduct(n0).doubleValue() / p0p2.innerProduct(n0).doubleValue());
			arr.set(1, -1.);
			arr.set(2, Math.sqrt(p2p1.squaredLength().doubleValue())*p0s.innerProduct(n1).doubleValue() / p2p1.innerProduct(n1).doubleValue());
		}
		else {
			arr.set(0, Math.sqrt(p0p2.squaredLength().doubleValue())*p0s.innerProduct(n0).doubleValue() / p0p2.innerProduct(n0).doubleValue());
			arr.set(1, Math.sqrt(p0p2.squaredLength().doubleValue())*p0s.innerProduct(n1).doubleValue() / p0p2.innerProduct(n1).doubleValue());
			arr.set(2, -1.);
		}
		return arr;
	}
}
