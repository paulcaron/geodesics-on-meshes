import java.util.ArrayList;
import java.util.PriorityQueue;
import java.util.HashMap;
import java.lang.NullPointerException;

import Jcg.polyhedron.*;
import Jcg.geometry.*;
import Jama.*;

public class ContinuousDijkstra {
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
		Face<Point_3> face = faceContainingPoint(source);

		return;
	}

	public double minDistance(Point_3 destination) {
		if (destination == null)
			throw new NullPointerException("Destination point is null");

		Face<Point_3> faceContainingPoint = faceContainingPoint(destination);

		if (GeoUtils.faceVertex(destination, faceContainingPoint)) {
			return minDistanceFromVertex(destination, faceContainingPoint);
		} else if (GeoUtils.halfedgePoint(destination, faceContainingPoint)) {
			return minDistanceFromHalfedge(destination, faceContainingPoint);
		} else {
			return minDistanceFromFace(destination, faceContainingPoint);
		}
	}

	/* The destination point is a vertex of face */
	private double minDistanceFromVertex(Point_3 destination, Face<Point_3> face) {
		Halfedge<Point_3> halfedge = face.getEdge();
		for (int i = 0; i < 3; i++) {
			if (GeoUtils.zero(destination.distanceFrom(halfedge.getVertex().getPoint()).doubleValue()))
				return minDistance(destination, halfedge);

			halfedge = halfedge.getNext();
		}
		throw new AssertionError();
	}

	/* The destination point lies on an edge of face */
	private double minDistanceFromHalfedge(Point_3 destination, Face<Point_3> face) {
		Halfedge<Point_3> halfedge = GeoUtils.closestHalfedgeFromPoint(destination, face);
		return minDistance(destination, halfedge);
	}

	/* The destination point lies in the interior of face */
	private double minDistanceFromFace(Point_3 destination, Face<Point_3> face) {
		double minDistance = Double.MAX_VALUE;

		Halfedge<Point_3> halfedge = face.getEdge();
		for (int i = 0; i < 3; i++) {
			double minDistanceThroughHalfedge = minDistance(destination, halfedge);

			if (minDistanceThroughHalfedge < minDistance)
				minDistance = minDistanceThroughHalfedge;

			halfedge = halfedge.getNext();
		}
		assert halfedge == face.getEdge();

		assert GeoUtils.positive(minDistance);

		return minDistance;
	}

	private double minDistance(Point_3 destination, Halfedge<Point_3> halfedge) {
		double minDistance = Double.MAX_VALUE;
		ArrayList<Window> windowsList = halfedgeToWindowsList.get(halfedge);
		assert windowsList.size() > 0;

		for (int i = 0; i < windowsList.size(); i++) {
			double minDistanceForWindow = minDistance(destination, windowsList.get(i));

			if (minDistanceForWindow < minDistance)
				minDistance = minDistanceForWindow;
		}

		assert GeoUtils.positive(minDistance);

		return minDistance;
	}

	private double minDistance(Point_3 destination, Window window) {
		double minDistance = Double.MAX_VALUE;

		Halfedge<Point_3> halfedge = window.getHalfedge();
		double lengthHalfedge = GeoUtils.lengthHalfedge(halfedge);
		Point_3 halfedgeOrigin = halfedge.getOpposite().getVertex().getPoint();

		double dx = lengthHalfedge / 100;
		for (double x = window.getBeginning(); x <= window.getEnd(); x += dx) {
			Point_3 halfedgePoint = GeoUtils.halfedgeMidpoint(halfedge, x);

			double distance = destination.distanceFrom(halfedgePoint).doubleValue() + window.getDistSource() +
				GeoUtils.stewart(window.getDistBeginning(), window.getDistEnd(), x - window.getBeginning(), window.getEnd() - x);

			if (distance < minDistance)
				minDistance = distance;
		}
		assert GeoUtils.positive(minDistance);

		return minDistance;
	}

	private Face<Point_3> faceContainingPoint(Point_3 point) {
		for (Face<Point_3> polyhedronFace : mesh.getFaces()) {
			if (GeoUtils.insideFace(point, polyhedronFace))
				return polyhedronFace;
		}

		throw new IllegalArgumentException("Point " + point.toString() + " not contained in any face of the mesh");
	}

	private void propagateWindow(Window window, PriorityQueue<Window> pq) {
		ArrayList<Double> propagatedExtremities = propagatedExtremities(window);
		Halfedge<Point_3> propagatingHalfedge = window.getHalfedge();
		assert propagatedExtremities.size() == 3;
		double firstExtremity = propagatedExtremities.get(0);
		double secondExtremity = propagatedExtremities.get(1);
		double thirdExtremity = propagatedExtremities.get(2);

		if (GeoUtils.negative(firstExtremity)) {
			/* The window propagates to only one window */

			Halfedge<Point_3> halfedge = propagatingHalfedge.getNext().getNext();

			assert GeoUtils.positive(secondExtremity);
			assert secondExtremity < thirdExtremity;
			assert thirdExtremity <= GeoUtils.lengthHalfedge(halfedge);

			Window propagatedWindow = new Window(
					secondExtremity,
					thirdExtremity,
					GeoUtils.halfedgeMidpoint(propagatingHalfedge, window.getEnd()).distanceFrom(
						GeoUtils.halfedgeMidpoint(halfedge, secondExtremity)).doubleValue() + window.getDistEnd(),
					GeoUtils.halfedgeMidpoint(propagatingHalfedge, window.getBeginning()).distanceFrom(
						GeoUtils.halfedgeMidpoint(halfedge, thirdExtremity)).doubleValue() + window.getDistBeginning(),
					window.getDistSource(),
					halfedge,
					true);
			mergeWindow(propagatedWindow, pq);

			/*
			 * If the window is adjacent to a saddle/boundary vertex
			 * we must propagate through its adjacent edge
			 */
			if (GeoUtils.equals(window.getEnd(), GeoUtils.lengthHalfedge(propagatingHalfedge)) &&
				specialVertex(propagatingHalfedge.getVertex())) {
				halfedge = propagatingHalfedge.getNext();
				Window firstExtraWindow = new Window(
						0,
						GeoUtils.lengthHalfedge(halfedge),
						0, /* the new pseudosource is the saddle/boundary vertex itself */
						GeoUtils.lengthHalfedge(halfedge),
						window.getDistEnd() + window.getDistSource(),
						halfedge,
						true);
				mergeWindow(firstExtraWindow, pq);

				halfedge = halfedge.getNext();
				Window secondExtraWindow = new Window(
						0,
						secondExtremity,
						GeoUtils.lengthHalfedge(propagatingHalfedge.getNext()),
						GeoUtils.halfedgeMidpoint(halfedge, secondExtremity).distanceFrom(
							GeoUtils.halfedgeEnd(propagatingHalfedge)).doubleValue(),
						window.getDistEnd() + window.getDistSource(),
						halfedge,
						true);
				mergeWindow(secondExtraWindow, pq);
			}

			assert GeoUtils.positive(window.getBeginning()) ||
				!specialVertex(propagatingHalfedge.getOpposite().getVertex());
		} else if (GeoUtils.negative(thirdExtremity)) {
			/* The window propagates to only one window */

			Halfedge<Point_3> halfedge = propagatingHalfedge.getNext();

			assert GeoUtils.positive(firstExtremity); 
			assert firstExtremity < secondExtremity;
			assert secondExtremity <= GeoUtils.lengthHalfedge(halfedge);

			Window propagatedWindow = new Window(
					firstExtremity,
					secondExtremity,
					GeoUtils.halfedgeMidpoint(propagatingHalfedge, window.getEnd()).distanceFrom(
						GeoUtils.halfedgeMidpoint(halfedge, firstExtremity)).doubleValue() + window.getDistEnd(),
					GeoUtils.halfedgeMidpoint(propagatingHalfedge, window.getBeginning()).distanceFrom(
						GeoUtils.halfedgeMidpoint(halfedge, secondExtremity)).doubleValue() + window.getDistBeginning(),
					window.getDistSource(),
					halfedge,
					true);
			mergeWindow(propagatedWindow, pq);

			/*
			 * If the window is adjacent to a saddle/boundary vertex
			 * we must propagate through its adjacent edge
			 */
			if (GeoUtils.zero(window.getBeginning()) ||
				specialVertex(propagatingHalfedge.getOpposite().getVertex())) {
				halfedge = propagatingHalfedge.getNext().getNext();
				Window firstExtraWindow = new Window(
						0,
						GeoUtils.lengthHalfedge(halfedge),
						GeoUtils.lengthHalfedge(halfedge), /* the new pseudosource is the saddle/boundary vertex itself */
						0,
						window.getDistBeginning() + window.getDistSource(),
						halfedge,
						true);
				mergeWindow(firstExtraWindow, pq);

				halfedge = propagatingHalfedge.getNext();
				Window secondExtraWindow = new Window(
						GeoUtils.lengthHalfedge(halfedge) - secondExtremity,
						GeoUtils.lengthHalfedge(halfedge),
						GeoUtils.halfedgeMidpoint(halfedge, secondExtremity).distanceFrom(
							GeoUtils.halfedgeOrigin(propagatingHalfedge)).doubleValue(),
						GeoUtils.lengthHalfedge(halfedge.getNext()),
						window.getDistBeginning() + window.getDistSource(),
						halfedge,
						true);
				mergeWindow(secondExtraWindow, pq);
			}

			assert GeoUtils.negative(window.getEnd() - GeoUtils.lengthHalfedge(propagatingHalfedge)) ||
				!specialVertex(propagatingHalfedge.getVertex());
		} else {
			/* 
			 * The window propagates to two windows,
			 * each in a different edge
			 */
			assert GeoUtils.negative(secondExtremity);

			Halfedge<Point_3> firstHalfedge = propagatingHalfedge.getNext();

			assert GeoUtils.positive(firstExtremity); 
			assert firstExtremity < GeoUtils.lengthHalfedge(firstHalfedge);

			Window firstPropagatedWindow = new Window(
					firstExtremity,
					GeoUtils.lengthHalfedge(firstHalfedge),
					GeoUtils.halfedgeMidpoint(firstHalfedge, firstExtremity).distanceFrom(
						GeoUtils.halfedgeMidpoint(propagatingHalfedge, window.getEnd())).doubleValue() + window.getDistEnd(),
					minDistance(GeoUtils.halfedgeEnd(firstHalfedge), window), // is it really correct?
					window.getDistSource(),
					firstHalfedge,
					true);
			mergeWindow(firstPropagatedWindow, pq);

			Halfedge<Point_3> secondHalfedge = propagatingHalfedge.getNext().getNext();

			assert GeoUtils.positive(thirdExtremity); 
			assert thirdExtremity < GeoUtils.lengthHalfedge(secondHalfedge);

			Window secondPropagatedWindow = new Window(
					thirdExtremity,
					GeoUtils.lengthHalfedge(secondHalfedge),
					minDistance(GeoUtils.halfedgeOrigin(secondHalfedge), window), // is it really correct?
					GeoUtils.halfedgeMidpoint(propagatingHalfedge, window.getBeginning()).distanceFrom(
						GeoUtils.halfedgeMidpoint(secondHalfedge, thirdExtremity)).doubleValue() + window.getDistBeginning(),
					window.getDistSource(),
					secondHalfedge,
					true);
			mergeWindow(secondPropagatedWindow, pq);

			assert !specialVertex(propagatingHalfedge.getVertex()) ||
				GeoUtils.negative(window.getEnd() - GeoUtils.lengthHalfedge(propagatingHalfedge));

			assert !specialVertex(propagatingHalfedge.getOpposite().getVertex()) ||
				GeoUtils.positive(window.getBeginning());
		}
	}

	private void mergeWindow(Window window, PriorityQueue<Window> pq) {
		return;
	}

	private boolean specialVertex(Vertex<Point_3> vertex) {
		return false;
	}

	public ArrayList<Double> propagatedExtremities(Window window) {
		ArrayList<Double> arr = new ArrayList<Double>(3);
		Point_3 p0, p1, p2, b0, b1;
		Point_3 source = window.getSource();
		p0 = window.getHalfedge().getVertex().getPoint();
		p1 = window.getHalfedge().getOpposite().getVertex().getPoint();
		p2 = window.getHalfedge().getNext().getVertex().getPoint();
		double lengthHalfedge = GeoUtils.lengthHalfedge(window.getHalfedge());
		boolean p20, p21; //p20=true if p2 is on the left side of the window : 
		Number[] coefficients0 = {1-window.getBeginning() / lengthHalfedge, window.getBeginning() / lengthHalfedge };
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
