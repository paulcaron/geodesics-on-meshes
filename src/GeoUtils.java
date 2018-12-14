import Jcg.polyhedron.*;
import Jcg.geometry.*;
import Jama.*;

public final class GeoUtils {
	private static final double THRESHOLD = 1e-9;

	private GeoUtils() { }

	public static boolean isZero(double x) {
		return Math.abs(x) <= THRESHOLD;
	}

	public static boolean isPositive(double x) {
		return x > THRESHOLD;
	}

	public static boolean isNegative(double x) {
		return x < -THRESHOLD;
	}

	public static boolean isEqual(double x, double y) {
		return isZero(x - y);
	}

	public static double getDistanceFromPointToHalfedge(Point_3 point, Halfedge<Point_3> halfedge) {
		Point_3 a = halfedge.getVertex().getPoint();
		Point_3 b = halfedge.getOpposite().getVertex().getPoint();

		Vector_3 ab = new Vector_3(a, b);
		Vector_3 ap = new Vector_3(a, point);

		Vector_3 crossProduct = ab.crossProduct(ap);

		return Math.sqrt(crossProduct.squaredLength().doubleValue() / ab.squaredLength().doubleValue());
	}

	public static boolean isPointOnFace(Point_3 point, Face<Point_3> face) {
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
		if (!isZero(matrix.det()))
			return false;

		/*
		 * Tries to solve the system alpha * ab + beta * ac = ap
		 * with 0 <= alpha, beta <= 1 and alpha + beta <= 1
		 */
		matrix = new Matrix(new double[][] {
			{ab.getX().doubleValue(), ac.getX().doubleValue()},
			{ab.getY().doubleValue(), ac.getY().doubleValue()}});
		if (!isZero(matrix.det())) {
			Matrix coeffs = matrix.solve(new Matrix (new double[][] {{ap.getX().doubleValue()}, {ap.getY().doubleValue()}}));
			double alpha = coeffs.get(0, 0),
				   beta = coeffs.get(1, 0);
			return !isNegative(alpha) && !isNegative(beta) && !isPositive(alpha + beta - 1);
		}

		matrix = new Matrix(new double[][] {
			{ab.getX().doubleValue(), ac.getX().doubleValue()},
			{ab.getZ().doubleValue(), ac.getZ().doubleValue()}});
		if (!isZero(matrix.det())) {
			Matrix coeffs = matrix.solve(new Matrix (new double[][] {{ap.getX().doubleValue()}, {ap.getZ().doubleValue()}}));
			double alpha = coeffs.get(0, 0),
				   beta = coeffs.get(1, 0);
			return !isNegative(alpha) && !isNegative(beta) && !isPositive(alpha + beta - 1);
		}

		matrix = new Matrix(new double[][] {
			{ab.getZ().doubleValue(), ac.getZ().doubleValue()},
			{ab.getY().doubleValue(), ac.getY().doubleValue()}});
		if (!isZero(matrix.det())) {
			Matrix coeffs = matrix.solve(new Matrix (new double[][] {{ap.getZ().doubleValue()}, {ap.getY().doubleValue()}}));
			double alpha = coeffs.get(0, 0),
				   beta = coeffs.get(1, 0);
			return !isNegative(alpha) && !isNegative(beta) && !isPositive(alpha + beta - 1);
		}

		return false;
	}

	public static Halfedge<Point_3> getClosestHalfedgeOfFaceToPoint(Point_3 point, Face<Point_3> face) {
		Halfedge<Point_3> halfedge = face.getEdge();
		Halfedge<Point_3> closestHalfedge = null;
		double distFromClosestHalfedge = Double.MAX_VALUE;

		for (int i = 0; i < 3; i++) {
			double distFromHalfdedge = getDistanceFromPointToHalfedge(point, halfedge);
			if (distFromHalfdedge < distFromClosestHalfedge) {
				distFromClosestHalfedge = distFromHalfdedge;
				closestHalfedge = halfedge;
			}

			halfedge = halfedge.getNext();
		}
		assert closestHalfedge != null;

		return closestHalfedge;
	}

	public static boolean isPointOnHalfedge(Point_3 point, Face<Point_3> face) {
		return isZero(getDistanceFromPointToHalfedge(point, getClosestHalfedgeOfFaceToPoint(point, face)));
	}

	public static boolean isPointVertexOfFace(Point_3 point, Face<Point_3> face) {
		Halfedge<Point_3> halfedge = face.getEdge();
		for (int i = 0; i < 3; i++) {
			if (areSamePoints(point, halfedge.getVertex().getPoint()))
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
	public static double getCevianLength(double b, double c, double n, double m) {
		double a = n + m;
		double d = Math.sqrt((b * b * m + c * c * n) / a - m * n);
		assert d >= 0;
		return d;
	}

	public static Vector_3 getUnaryVector(Vector_3 vector) {
		Vector_3 unaryVector = vector.divisionByScalar(Math.sqrt(vector.squaredLength().doubleValue()));
		assert isEqual(unaryVector.squaredLength().doubleValue(), 1);
		return unaryVector;
	}

	public static double getHalfedgeLength(Halfedge<Point_3> halfedge) {
		Vertex<Point_3> destination = halfedge.getVertex();
		Vertex<Point_3> origin = halfedge.getOpposite().getVertex();

		return (double) destination.getPoint().distanceFrom(origin.getPoint());
	}

	public static double getTriangleHeight(double base, double firstSide, double secondSide) {
		double p = (base + firstSide + secondSide) / 2;
		return 2 * Math.sqrt(p * (p - base) * (p - firstSide) * (p - secondSide)) / base;
	}

	public static Point_3 sumPointVector(Point_3 origin, Vector_3 direction, double norm) {
		Vector_3 unaryVector = getUnaryVector(direction);
		Point_3 result = new Point_3(origin);
		result.translateOf(unaryVector.multiplyByScalar(norm));

		return result;
	}

	public static Point_3 getHalfedgePoint(Halfedge<Point_3> halfedge, double norm) {
		Point_3 vectorOrigin = halfedge.getOpposite().getVertex().getPoint();
		Point_3 vectorEnd = halfedge.getVertex().getPoint();
		Vector_3 vector = new Vector_3(vectorOrigin, vectorEnd);

		assert isEqual(vector.squaredLength().doubleValue(), norm * norm);

		return sumPointVector(vectorOrigin, vector, norm);
	}

	public static Point_3 getHalfedgeOrigin(Halfedge<Point_3> halfedge) {
		return halfedge.getOpposite().getVertex().getPoint();
	}

	public static Point_3 getHalfedgeEnd(Halfedge<Point_3> halfedge) {
		return halfedge.getVertex().getPoint();
	}

	public static Pair<Double, Double> solveSecondDegreeEquation(double B, double C) {
		double delta = B * B - 4 * C;
		if (delta < 0)
			return new Pair<Double, Double> (Double.MAX_VALUE, Double.MAX_VALUE);
		double average = - B / 2;
		return new Pair<Double, Double> (average - 0.5 * Math.sqrt(delta), average + 0.5 * Math.sqrt(delta));
	}

	public static boolean areSamePoints(Point_3 first, Point_3 second) {
		return isZero(first.distanceFrom(second).doubleValue());
	}

	public static Vertex<Point_3> identifyVertex(Point_3 point, Face<Point_3> face) {
		Halfedge<Point_3> halfedge = face.getEdge();
		Vertex<Point_3> vertex = null;
		for (int i = 0; i < 3; i++) {
			Point_3 vertexPoint = halfedge.getVertex().getPoint();
			if (areSamePoints(point, vertexPoint))
				vertex = halfedge.getVertex();
			else halfedge = halfedge.getNext();
		}
		assert vertex != null;

		return vertex;
	}

	/*
	 * Returns the vertex point of a triangle opposite
	 * to a given triangle side. The origin is supposed to be
	 * the vertex adjacent to the base and the left side of the
	 * triangle, with axis Ox-axis supporting the base, and oriented
	 * from the origin to the other base vertex. The desired vertex is
	 * admitted to be on the positive half-plane y > 0.
	 *
	 * @param base length of side of triangle opposite to desired vertex
	 * @param sideLeft length of left side of the triangle
	 * @param sideRight length of right side of the triangle
	 */
	public Point_2 getTriangleVertexPlaneProjection(double base, double lengthLeft, double lengthRight) {
		double cos = Math.pow(lengthLeft, 2) - Math.pow(lengthRight, 2) + Math.pow(base, 2) / \
					 (2 * lengthLeft * base);
		double sin = Math.sqrt(1 - Math.pow(cos, 2));

		double x = lengthLeft * cos;
		double y = lengthLeft * sin;

		return new Point_2(x, y);
	}
}
