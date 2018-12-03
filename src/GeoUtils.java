import Jcg.polyhedron.*;
import Jcg.geometry.*;
import Jama.*;

public final class GeoUtils {
	private static final double THRESHOLD = 1e-9;

	private GeoUtils() { }

	public static boolean zero(double x) {
		return Math.abs(x) <= THRESHOLD;
	}

	public static boolean positive(double x) {
		return x > THRESHOLD;
	}

	public static boolean negative(double x) {
		return x < -THRESHOLD;
	}

	public static boolean equals(double x, double y) {
		return zero(x - y);
	}

	public static double distancePointHalfedge(Point_3 point, Halfedge<Point_3> halfedge) {
		Point_3 a = halfedge.getVertex().getPoint();
		Point_3 b = halfedge.getOpposite().getVertex().getPoint();

		Vector_3 ab = new Vector_3(a, b);
		Vector_3 ap = new Vector_3(a, point);

		Vector_3 crossProduct = ab.crossProduct(ap);

		return Math.sqrt(crossProduct.squaredLength().doubleValue() / ab.squaredLength().doubleValue());
	}

	public static boolean insideFace(Point_3 point, Face<Point_3> face) {
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

	public static Halfedge<Point_3> closestHalfedgeFromPoint(Point_3 point, Face<Point_3> face) {
		Halfedge<Point_3> halfedge = face.getEdge();
		for (int i = 0; i < 3; i++) {
			if (zero(distancePointHalfedge(point, halfedge)))
				return halfedge;

			halfedge = halfedge.getNext();
		}
		assert halfedge == halfedge.getNext();

		return null;
	}

	public static boolean isHalfedgePoint(Point_3 point, Face<Point_3> face) {
		return zero(distancePointHalfedge(point, closestHalfedgeFromPoint(point, face)));
	}

	public static boolean isFaceVertex(Point_3 point, Face<Point_3> face) {
		Halfedge<Point_3> halfedge = face.getEdge();
		for (int i = 0; i < 3; i++) {
			if (zero(point.distanceFrom(halfedge.getVertex().getPoint()).doubleValue()))
				return true;

			halfedge = halfedge.getNext();
		}
		assert halfedge == face.getEdge();

		return false;
	}

	/*
	 * Following convention employed by Wikipedia to name triangle sides
	 * https://en.wikipedia.org/wiki/Stewart's_theorem
	 */
	public static double stewart(double b, double c, double n, double m) {
		double a = n + m;
		double d = Math.sqrt((b * b * m + c * c * n) / a - m * n);
		assert d >= 0;
		return d;
	}

	public static Vector_3 unaryVector(Vector_3 vector) {
		Vector_3 unaryVector = vector.divisionByScalar(Math.sqrt(vector.squaredLength().doubleValue()));
		assert equals(unaryVector.squaredLength().doubleValue(), 1);
		return unaryVector;
	}

	public static double lengthHalfedge(Halfedge<Point_3> halfedge) {
		Vertex<Point_3> destination = halfedge.getVertex();
		Vertex<Point_3> origin = halfedge.getOpposite().getVertex();

		return (double) destination.getPoint().distanceFrom(origin.getPoint());
	}

	public static double triangleHeight(double base, double firstSide, double secondSide) {
		double p = (base + firstSide + secondSide) / 2;
		return 2 * Math.sqrt(p * (p - base) * (p - firstSide) * (p - secondSide)) / base;
	}

	public static Point_3 sumPointVector(Point_3 origin, Vector_3 direction, double norm) {
		Vector_3 unaryVector = unaryVector(direction);
		Point_3 result = new Point_3(origin);
		result.translateOf(unaryVector.multiplyByScalar(norm));

		return result;
	}

	public static Point_3 halfedgeMidpoint(Halfedge<Point_3> halfedge, double norm) {
		Point_3 vectorOrigin = halfedge.getOpposite().getVertex().getPoint();
		Point_3 vectorEnd = halfedge.getVertex().getPoint();
		Vector_3 vector = new Vector_3(vectorOrigin, vectorEnd);

		assert equals(vector.squaredLength().doubleValue(), norm * norm);

		return sumPointVector(vectorOrigin, vector, norm);
	}

	public static Point_3 halfedgeOrigin(Halfedge<Point_3> halfedge) {
		return halfedge.getOpposite().getVertex().getPoint();
	}

	public static Point_3 halfedgeEnd(Halfedge<Point_3> halfedge) {
		return halfedge.getVertex().getPoint();
	}
}
