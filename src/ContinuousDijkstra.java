import java.util.ArrayList;
import java.util.PriorityQueue;


import java.util.HashMap;
import java.lang.NullPointerException;

import Jcg.polyhedron.*;
import Jcg.geometry.*;

public class ContinuousDijkstra {
	private final SurfaceMesh mesh;
	private HashMap<Halfedge<Point_3>, ArrayList<Window>> halfedgeToWindowsList;

	public ContinuousDijkstra(SurfaceMesh mesh) {
		if (mesh == null)
			throw new NullPointerException("Mesh is null");

		this.mesh = mesh;
	}

	public void buildDistances(Point_3 source) {
		if (source == null)
			throw new NullPointerException("Source point is null");

		PriorityQueue<Window> pq = new PriorityQueue<>();

		initializeHashMap();
		initializePriorityQueue(source, pq);

		int i = 0;
		while (!pq.isEmpty()) {
			Window window = pq.poll();
			if (isValidWindow(window)) {
				System.out.println(i++);
				propagateWindow(window, pq);
			}		
		}
	}
	
	private boolean isValidWindow(Window window) {
		ArrayList<Window> windowsList;
		if (halfedgeToWindowsList.containsKey(window.getHalfedge()))
			windowsList = halfedgeToWindowsList.get(window.getHalfedge());
		else
			windowsList = halfedgeToWindowsList.get(window.getHalfedge().getOpposite());

		for (Window listWindow : windowsList) 
			if (listWindow == window)
				return true;

		return false;
	}

	private void initializeHashMap() {
		halfedgeToWindowsList = new HashMap<>();

		ArrayList<Halfedge<Point_3>> halfedges = mesh.getHalfedges();
		for (Halfedge<Point_3> halfedge : halfedges) {
			if (!halfedgeToWindowsList.containsKey(halfedge.getOpposite()))
				halfedgeToWindowsList.put(halfedge, new ArrayList<Window>());
		}
		
		assert allHalfedgesIncluded();
	}

	private void initializePriorityQueue(
			Point_3 source,
			PriorityQueue<Window> pq
	) {
		if (source == null)
			throw new NullPointerException("Source point is null");
			
		Face<Point_3> faceContainingPoint = getFaceContainingPoint(source);

		if (GeoUtils.isPointVertexOfFace(source, faceContainingPoint)) {
			initializeWithVertex(source, faceContainingPoint, pq);
		} else if (GeoUtils.isPointOnHalfedge(source, faceContainingPoint)) {
			initializeWithHalfedgePoint(source, faceContainingPoint, pq);
		} else {
			initializeWithFacePoint(source, faceContainingPoint, pq);
		}
	}
	
	private boolean allHalfedgesIncluded() {
		for (Halfedge<Point_3> hf : mesh.getHalfedges())
			if (halfedgeToWindowsList.containsKey(hf) || halfedgeToWindowsList.containsKey(hf.getOpposite()))
				return true;
		return false;
	}

	private void initializeWithVertex(Point_3 source, Face<Point_3> face, PriorityQueue<Window> pq) {
		Vertex<Point_3> vertex = GeoUtils.identifyVertex(source, face);
		Halfedge<Point_3> halfedge = vertex.getHalfedge();
		
		assert allHalfedgesIncluded();
		do {
			Window window = new Window(
					0,
					GeoUtils.getHalfedgeLength(halfedge),
					GeoUtils.getHalfedgeLength(halfedge),
					0,
					0,
					halfedge,
					true);
			ArrayList<Window> windowList;
			if (halfedgeToWindowsList.containsKey(halfedge)) {
				windowList = halfedgeToWindowsList.get(halfedge);
			}
			else {
				windowList = halfedgeToWindowsList.get(halfedge.getOpposite());
				window = window.getOpposite();
			}
			pq.add(window);
			windowList.add(window);
			halfedge = halfedge.getNext().getOpposite();
		}while(halfedge != vertex.getHalfedge());
	}

	private void initializeWithHalfedgePoint(Point_3 source, Face<Point_3> face, PriorityQueue<Window> pq) {
		// Find the halfedge the source is on
		Halfedge<Point_3> halfedge = GeoUtils.getClosestHalfedgeOfFaceToPoint(source, face);
		
		//Initialze pq with 3 windows
		
		//window 1 (initialHalfedge)		
		Vector_3 sp0 = new Vector_3(source, halfedge.getOpposite().getVertex().getPoint());
		Vector_3 sp1 = new Vector_3(source, halfedge.getVertex().getPoint());
		Window window = new Window(
				0, 
				GeoUtils.getHalfedgeLength(halfedge),
				Math.sqrt(sp0.squaredLength().doubleValue()),
				Math.sqrt(sp1.squaredLength().doubleValue()),
				0,
				halfedge,
				true);
		ArrayList<Window> windowList;
		if (halfedgeToWindowsList.containsKey(halfedge)) {
			windowList = halfedgeToWindowsList.get(halfedge);
		}
		else {
			windowList = halfedgeToWindowsList.get(halfedge.getOpposite());
			window = window.getOpposite();
		}
		pq.add(window);
		windowList.add(window);
		
		//window 2 (initialHalfedge.opposite.next)		
		halfedge = halfedge.getOpposite().getNext();
		sp0 = new Vector_3(source, halfedge.getOpposite().getVertex().getPoint());
		sp1 = new Vector_3(source, halfedge.getVertex().getPoint());
		window = new Window(
				0, 
				GeoUtils.getHalfedgeLength(halfedge),
				Math.sqrt(sp0.squaredLength().doubleValue()),
				Math.sqrt(sp1.squaredLength().doubleValue()),
				0,
				halfedge,
				true);
		if (halfedgeToWindowsList.containsKey(halfedge)) {
			windowList = halfedgeToWindowsList.get(halfedge);
		}
		else {
			windowList = halfedgeToWindowsList.get(halfedge.getOpposite());
			window = window.getOpposite();
		}
		pq.add(window);
		windowList.add(window);
		
		
		//window 3 (initialHalfedge.opposite.next.next)	
		halfedge = halfedge.getNext();
		sp0 = new Vector_3(source, halfedge.getOpposite().getVertex().getPoint());
		sp1 = new Vector_3(source, halfedge.getVertex().getPoint());
		window = new Window(
				0, 
				GeoUtils.getHalfedgeLength(halfedge),
				Math.sqrt(sp0.squaredLength().doubleValue()),
				Math.sqrt(sp1.squaredLength().doubleValue()),
				0,
				halfedge,
				true);
		if (halfedgeToWindowsList.containsKey(halfedge)) {
			windowList = halfedgeToWindowsList.get(halfedge);
		}
		else {
			windowList = halfedgeToWindowsList.get(halfedge.getOpposite());
			window = window.getOpposite();
		}
		pq.add(window);
		windowList.add(window);
	}

	private void initializeWithFacePoint(Point_3 source, Face<Point_3> face, PriorityQueue<Window> pq) {
		Halfedge<Point_3> halfedge = face.getEdge();
		Halfedge<Point_3> pEdge = halfedge;
		do {
			Vector_3 sp0 = new Vector_3(source, halfedge.getOpposite().getVertex().getPoint());
			Vector_3 sp1 = new Vector_3(source, halfedge.getVertex().getPoint());
			Window window = new Window(
					0,
					GeoUtils.getHalfedgeLength(halfedge),
					Math.sqrt(sp0.squaredLength().doubleValue()),
					Math.sqrt(sp1.squaredLength().doubleValue()),
					0,
					halfedge,
					true);
			ArrayList<Window> windowList;
			if (halfedgeToWindowsList.containsKey(halfedge)) {
				windowList = halfedgeToWindowsList.get(halfedge);
			}
			else {
				windowList = halfedgeToWindowsList.get(halfedge.getOpposite());
				window = window.getOpposite();
			}
			pq.add(window);
			windowList.add(window);
			
			halfedge=halfedge.getNext();
		}while(!halfedge.equals(pEdge));
	}

	public double getDistanceToSource(Point_3 destination) {
		if (destination == null)
			throw new NullPointerException("Destination point is null");

		Face<Point_3> faceContainingPoint = getFaceContainingPoint(destination);

		if (GeoUtils.isPointVertexOfFace(destination, faceContainingPoint)) {
			return getMinDistanceFromVertexPointToSource(destination, faceContainingPoint);
		} else if (GeoUtils.isPointOnHalfedge(destination, faceContainingPoint)) {
			return getMinDistanceFromHalfedgePointToSource(destination, faceContainingPoint);
		} else {
			return getMinDistanceToSource(destination, faceContainingPoint);
		}
	}

	/* The destination point is a vertex of face */
	private double getMinDistanceFromVertexPointToSource(Point_3 destination, Face<Point_3> face) {
		Vertex<Point_3> vertex = GeoUtils.identifyVertex(destination, face);
		return getDistanceToSourcePassingOverHalfedge(destination, vertex.getHalfedge());
	}

	/* The destination point lies on an edge of face */
	private double getMinDistanceFromHalfedgePointToSource(Point_3 destination, Face<Point_3> face) {
		Halfedge<Point_3> halfedge = GeoUtils.getClosestHalfedgeOfFaceToPoint(destination, face);
		return getDistanceToSourcePassingOverHalfedge(destination, halfedge);
	}

	/* The destination point lies in the interior of face */
	private double getMinDistanceToSource(Point_3 destination, Face<Point_3> face) {
		double distanceToSource = Double.MAX_VALUE;

		Halfedge<Point_3> halfedge = face.getEdge();
		for (int i = 0; i < 3; i++) {
			double minDistanceThroughHalfedge = getDistanceToSourcePassingOverHalfedge(destination, halfedge);

			if (minDistanceThroughHalfedge < distanceToSource)
				distanceToSource = minDistanceThroughHalfedge;

			halfedge = halfedge.getNext();
		}
		assert halfedge == face.getEdge();

		assert GeoUtils.isPositive(distanceToSource);

		return distanceToSource;
	}

	private double getDistanceToSourcePassingOverHalfedge(Point_3 destination, Halfedge<Point_3> halfedge) {
		double distanceToSource = Double.MAX_VALUE;
		ArrayList<Window> windowsList;
		if (halfedgeToWindowsList.containsKey(halfedge))
			windowsList = halfedgeToWindowsList.get(halfedge);
		windowsList = halfedgeToWindowsList.get(halfedge.getOpposite());
		
		assert windowsList.size() > 0;

		for (int i = 0; i < windowsList.size(); i++) {
			double minDistanceForWindow = getDistanceToSourcePassingOverWindow(destination, windowsList.get(i));

			if (minDistanceForWindow < distanceToSource)
				distanceToSource = minDistanceForWindow;
		}

		assert GeoUtils.isPositive(distanceToSource);

		return distanceToSource;
	}

	private double getDistanceToSourcePassingOverWindow(Point_3 destination, Window window) {
		double distanceToSource = Double.MAX_VALUE;

		Halfedge<Point_3> halfedge = window.getHalfedge();
		double halfedgeLength = GeoUtils.getHalfedgeLength(halfedge);

		double dx = halfedgeLength / 100;
		for (double x = window.getStart(); x <= window.getEnd(); x += dx) {
			Point_3 halfedgePoint = GeoUtils.getHalfedgePoint(halfedge, x);

			double distance = destination.distanceFrom(halfedgePoint).doubleValue() + window.getDistSource();
			if (GeoUtils.isEqual(x, window.getStart()))
				distance += window.getDistStart();
			else if (GeoUtils.isEqual(x, window.getEnd()))
				distance += window.getDistEnd();
			else
				GeoUtils.getCevianLength(window.getDistStart(), window.getDistEnd(), x - window.getStart(), window.getEnd() - x);

			if (distance < distanceToSource)
				distanceToSource = distance;
		}
		assert GeoUtils.isPositive(distanceToSource);

		return distanceToSource;
	}

	private Face<Point_3> getFaceContainingPoint(Point_3 point) {
		for (Face<Point_3> polyhedronFace : mesh.getFaces()) {
			if (GeoUtils.isPointOnFace(point, polyhedronFace))
				return polyhedronFace;
		}
		
		throw new IllegalArgumentException("Point " + point.toString() + " not contained in any face of the mesh");
	}

	private void propagateWindow(Window window, PriorityQueue<Window> pq) {
		if (!window.getSide())
			window = window.getOpposite();
		assert window.getSide();

		ArrayList<Double> propagatedExtremities = getPropagatedExtremities(window);
		Halfedge<Point_3> propagatingHalfedge = window.getHalfedge();

		assert propagatedExtremities.size() == 3;

		double firstExtremity = propagatedExtremities.get(0);
		double secondExtremity = propagatedExtremities.get(1);
		double thirdExtremity = propagatedExtremities.get(2);

		if (GeoUtils.isNegative(firstExtremity)) {
			/* The window propagates to only one window */

			Halfedge<Point_3> halfedge = propagatingHalfedge.getNext().getNext();

			assert !GeoUtils.isNegative(secondExtremity);
			assert secondExtremity < thirdExtremity;
			assert thirdExtremity <= GeoUtils.getHalfedgeLength(halfedge);

			Window propagatedWindow = new Window(
					secondExtremity,
					thirdExtremity,
					GeoUtils.getHalfedgePoint(propagatingHalfedge, window.getEnd()).distanceFrom(
						GeoUtils.getHalfedgePoint(halfedge, secondExtremity)).doubleValue() + window.getDistEnd(),
					GeoUtils.getHalfedgePoint(propagatingHalfedge, window.getStart()).distanceFrom(
						GeoUtils.getHalfedgePoint(halfedge, thirdExtremity)).doubleValue() + window.getDistStart(),
					window.getDistSource(),
					halfedge,
					true);
			insertWindow(propagatedWindow, pq);

			/*
			 * If the window is adjacent to a saddle/boundary vertex
			 * we must propagate through its adjacent edge
			 */
			if (GeoUtils.isEqual(window.getEnd(), GeoUtils.getHalfedgeLength(propagatingHalfedge)) &&
				isSpecialVertex(propagatingHalfedge.getVertex())) {
				halfedge = propagatingHalfedge.getNext();
				Window firstExtraWindow = new Window(
						0,
						GeoUtils.getHalfedgeLength(halfedge),
						0, /* the new pseudosource is the saddle/boundary vertex itself */
						GeoUtils.getHalfedgeLength(halfedge),
						window.getDistEnd() + window.getDistSource(),
						halfedge,
						true);
				insertWindow(firstExtraWindow, pq);

				halfedge = halfedge.getNext();
				Window secondExtraWindow = new Window(
						0,
						secondExtremity,
						GeoUtils.getHalfedgeLength(propagatingHalfedge.getNext()),
						GeoUtils.getHalfedgePoint(halfedge, secondExtremity).distanceFrom(
							GeoUtils.getHalfedgeEnd(propagatingHalfedge)).doubleValue(),
						window.getDistEnd() + window.getDistSource(),
						halfedge,
						true);
				insertWindow(secondExtraWindow, pq);
			}

		} else if (GeoUtils.isNegative(thirdExtremity)) {
			/* The window propagates to only one window */

			Halfedge<Point_3> halfedge = propagatingHalfedge.getNext();

			assert !GeoUtils.isNegative(firstExtremity); 
			assert firstExtremity < secondExtremity;
			assert secondExtremity <= GeoUtils.getHalfedgeLength(halfedge);

			Window propagatedWindow = new Window(
					firstExtremity,
					secondExtremity,
					GeoUtils.getHalfedgePoint(propagatingHalfedge, window.getEnd()).distanceFrom(
						GeoUtils.getHalfedgePoint(halfedge, firstExtremity)).doubleValue() + window.getDistEnd(),
					GeoUtils.getHalfedgePoint(propagatingHalfedge, window.getStart()).distanceFrom(
						GeoUtils.getHalfedgePoint(halfedge, secondExtremity)).doubleValue() + window.getDistStart(),
					window.getDistSource(),
					halfedge,
					true);
			insertWindow(propagatedWindow, pq);

			/*
			 * If the window is adjacent to a saddle/boundary vertex
			 * we must propagate through its adjacent edge
			 */
			if (GeoUtils.isZero(window.getStart()) ||
				isSpecialVertex(propagatingHalfedge.getOpposite().getVertex())) {
				halfedge = propagatingHalfedge.getNext().getNext();
				Window firstExtraWindow = new Window(
						0,
						GeoUtils.getHalfedgeLength(halfedge),
						GeoUtils.getHalfedgeLength(halfedge), /* the new pseudosource is the saddle/boundary vertex itself */
						0,
						window.getDistStart() + window.getDistSource(),
						halfedge,
						true);
				insertWindow(firstExtraWindow, pq);

				halfedge = propagatingHalfedge.getNext();
				Window secondExtraWindow = new Window(
						GeoUtils.getHalfedgeLength(halfedge) - secondExtremity,
						GeoUtils.getHalfedgeLength(halfedge),
						GeoUtils.getHalfedgePoint(halfedge, secondExtremity).distanceFrom(
							GeoUtils.getHalfedgeOrigin(propagatingHalfedge)).doubleValue(),
						GeoUtils.getHalfedgeLength(halfedge.getNext()),
						window.getDistStart() + window.getDistSource(),
						halfedge,
						true);
				insertWindow(secondExtraWindow, pq);
			}
		} else {
			/* 
			 * The window propagates to two windows,
			 * each in a different edge
			 */
			assert GeoUtils.isNegative(secondExtremity);
			
			// comecou aqui
			double halfedgeLength = GeoUtils.getHalfedgeLength(propagatingHalfedge);
			Point_2 p0 = new Point_2(0, 0);
			Point_2 p1 = new Point_2(halfedgeLength, 0);
			Point_2 b0 = new Point_2(halfedgeLength - window.getEnd(), 0);
			Point_2 source;
			if (GeoUtils.isZero(window.getDistEnd())) 
				source = p0;
			else if (GeoUtils.isZero(window.getDistStart()))
				source = p1;
			else
				source = GeoUtils.getTriangleVertexPlaneProjection(window.getLength(), window.getDistEnd(), window.getDistStart());
			Point_2 p2 = GeoUtils.getTriangleVertexPlaneProjection(halfedgeLength, GeoUtils.getHalfedgeLength(propagatingHalfedge.getNext()), GeoUtils.getHalfedgeLength(propagatingHalfedge.getNext().getNext()));
			source.setX(source.getX().doubleValue() + b0.getX().doubleValue());
			p2.setY(-p2.getY().doubleValue());
			Vector_2 p2s = new Vector_2(p2, source);
			// terminou aqui

			Halfedge<Point_3> firstHalfedge = propagatingHalfedge.getNext();

			assert firstExtremity < GeoUtils.getHalfedgeLength(firstHalfedge);
			
			Window firstPropagatedWindow = new Window(
					firstExtremity,
					GeoUtils.getHalfedgeLength(firstHalfedge),
					GeoUtils.getHalfedgePoint(firstHalfedge, firstExtremity).distanceFrom(
						GeoUtils.getHalfedgePoint(propagatingHalfedge, window.getEnd())).doubleValue() + window.getDistEnd(),
					Math.sqrt(p2s.squaredLength().doubleValue()), // is it really correct?
					window.getDistSource(),
					firstHalfedge,
					true);
			insertWindow(firstPropagatedWindow, pq);

			Halfedge<Point_3> secondHalfedge = propagatingHalfedge.getNext().getNext();

			assert GeoUtils.isPositive(thirdExtremity); 
			assert !GeoUtils.isPositive(thirdExtremity - GeoUtils.getHalfedgeLength(secondHalfedge));
			
			Window secondPropagatedWindow = new Window(
					0,
					thirdExtremity,
					Math.sqrt(p2s.squaredLength().doubleValue()), // is it really correct?
					GeoUtils.getHalfedgePoint(propagatingHalfedge, window.getStart()).distanceFrom(
						GeoUtils.getHalfedgePoint(secondHalfedge, thirdExtremity)).doubleValue() + window.getDistStart(),
					window.getDistSource(),
					secondHalfedge,
					true);
			insertWindow(secondPropagatedWindow, pq);
		}
	}

	private void insertWindow(Window newWindow, PriorityQueue<Window> pq) {
		Halfedge<Point_3> halfedge = newWindow.getHalfedge();
		if (!halfedgeToWindowsList.containsKey(halfedge)) {
			halfedge = halfedge.getOpposite();
			newWindow = newWindow.getOpposite();
		}

		ArrayList<Window> oldWindows = halfedgeToWindowsList.get(halfedge);
		ArrayList<Window> newWindows = new ArrayList<>();
		//assert areWindowsInIncreasingOrder(oldWindows);
		/*
		 * This linear scan is too slow, we should optimize
		 * it using binary search to achieve the time complexity
		 * indicated in the paper.
		 */
		int i = 0;
		while (i < oldWindows.size()) {
			Window oldWindow = oldWindows.get(i);

			if (oldWindow.getEnd() <= newWindow.getStart()) {
				newWindows.add(oldWindow);
				i++;
			} else if (newWindow.getEnd() <= oldWindow.getStart()) {
				break;
			} else if (newWindow.getStart() < oldWindow.getStart()) {
				Window auxWindow = newWindow.setEnd(oldWindow.getStart());
				pq.add(auxWindow);
				newWindows.add(auxWindow);
				newWindow = newWindow.setStart(oldWindow.getStart());
			}  else {
				pq.remove(oldWindow);
				Window adjustedWindow = mergeWindows(oldWindow, newWindow);
				newWindows.add(adjustedWindow);
				
				if (adjustedWindow.compareTo(oldWindow) < 0)
					pq.add(adjustedWindow);

				if (GeoUtils.isNegative(adjustedWindow.getEnd() - oldWindow.getEnd()))
					oldWindows.set(i, oldWindow.setStart(adjustedWindow.getEnd()));
				else
					i++;
				if (GeoUtils.isNegative(adjustedWindow.getEnd() - newWindow.getEnd())) 
					newWindow = newWindow.setStart(adjustedWindow.getEnd());
				else {
					newWindow = null;
					break;
				}
			}
		}

		if (newWindow != null) {
			newWindows.add(newWindow);
			pq.add(newWindow);
		}
		while (i < oldWindows.size())
			newWindows.add(oldWindows.get(i++));

		halfedgeToWindowsList.put(halfedge, newWindows);
	}

	private Window mergeWindows(Window leftWindow, Window rightWindow) {
		assert leftWindow.getStart() <= rightWindow.getStart() &&
			!areWindowsDisjoint(leftWindow, rightWindow);

		double leftSourceProj = leftWindow.getAbscisseSourcePlaneProjection();
		double rightSourceProj = rightWindow.getAbscisseSourcePlaneProjection();

		double heightLeftTriangle = GeoUtils.getTriangleHeight(leftWindow.getLength(), leftWindow.getDistStart(), leftWindow.getDistEnd());
		double heightRightTriangle = GeoUtils.getTriangleHeight(rightWindow.getLength(), rightWindow.getDistStart(), rightWindow.getDistEnd());

		double leftSourceNormSquared = leftSourceProj * leftSourceProj + heightLeftTriangle * heightLeftTriangle;
		double rightSourceNormSquared = rightSourceProj * rightSourceProj + heightRightTriangle * heightRightTriangle;

		double alpha = rightSourceProj - leftSourceProj;
		double beta = rightWindow.getDistSource() - leftWindow.getDistSource();
		double gamma = leftSourceNormSquared - rightSourceNormSquared - beta * beta;

		double A = alpha * alpha - beta * beta;
		double B = gamma * alpha + 2 * rightSourceProj * beta * beta;
		double C = 0.25 * gamma * gamma - rightSourceNormSquared * beta * beta;

		double maxPossible = Math.min(leftWindow.getEnd(), rightWindow.getEnd());
		double solution = getSolutionInRange(GeoUtils.solveSecondDegreeEquation(B / A, C / A), maxPossible);

		if (leftWindow.getDistStart() + leftWindow.getDistSource() < rightWindow.getDistStart() + rightWindow.getDistSource()) {
			if (solution > 0)
				return leftWindow.setEnd(leftWindow.getStart() + solution);
			else
				return leftWindow;
		}
		else {
			if (solution > 0)
				return rightWindow.setEnd(rightWindow.getStart() + solution);
			else
				return rightWindow;
		}
	}

	private double getSolutionInRange(Pair<Double, Double> solutions, double maxPossible) {
		double first = solutions.first();
		double second = solutions.second();
		if (first > 0 && first <= maxPossible)
			return first;
		if (second > 0 && second <= maxPossible)
			return second;

		return -1;
	}

	private boolean areWindowsDisjoint(Window leftWindow, Window rightWindow) {
		return leftWindow.getEnd() <= rightWindow.getStart();
	}

	private boolean areWindowsInIncreasingOrder(ArrayList<Window> windows) {
		for (int i = 1; i < windows.size(); i++)
			if (!GeoUtils.isNegative(windows.get(i).getStart() - windows.get(i - 1).getEnd())) {
				System.out.printf("%f\n", windows.get(i).getStart() - windows.get(i - 1).getEnd());
				return false;
			}
		return true;
	}

	private boolean isSpecialVertex(Vertex<Point_3> vertex) {
		return true;
	}

	public ArrayList<Double> getPropagatedExtremities(Window window) {
		ArrayList<Double> arr = new ArrayList<Double>(3);
		
		if(GeoUtils.isZero(window.getDistStart()) || GeoUtils.isZero(window.getDistEnd())) { // s=p0 or s=p1
			arr.add(0.);
			arr.add(-1.);
			arr.add(GeoUtils.getHalfedgeLength(window.getHalfedge().getNext().getNext()));
			return arr;
		}

		Halfedge<Point_3> halfedge = window.getHalfedge();
		double halfedgeLength = GeoUtils.getHalfedgeLength(halfedge);
		Point_2 p0 = new Point_2(0, 0);
		Point_2 p1 = new Point_2(halfedgeLength, 0);
		Point_2 b0 = new Point_2(halfedgeLength - window.getEnd(), 0);
		Point_2 b1 = new Point_2(halfedgeLength - window.getStart(), 0);
		Point_2 source = GeoUtils.getTriangleVertexPlaneProjection(window.getLength(), window.getDistEnd(), window.getDistStart());
		Point_2 p2 = GeoUtils.getTriangleVertexPlaneProjection(halfedgeLength, GeoUtils.getHalfedgeLength(halfedge.getNext()), GeoUtils.getHalfedgeLength(halfedge.getNext().getNext()));
		source.setX(source.getX().doubleValue() + b0.getX().doubleValue());
		p2.setY(-p2.getY().doubleValue());
		
		Vector_2 b0s = new Vector_2(b0, source);
		Vector_2 b1s = new Vector_2(b1, source);
		Vector_2 b0p2 = new Vector_2(b0, p2);
		Vector_2 b1p2 = new Vector_2(b1, p2);
		Vector_2 p0p2 = new Vector_2(p0, p2);
		Vector_2 p2p1 = new Vector_2(p2, p1);
		Vector_2 p0s = new Vector_2(p0, source);
		Vector_2 p2s = new Vector_2(p2, source);
		Vector_2 n0 = new Vector_2(-b0s.getY().doubleValue(), b0s.getX());
		Vector_2 n1 = new Vector_2(-b1s.getY().doubleValue(), b1s.getX());
		boolean p20, p21; //p20=true if p2 is on the left side of the window
		p20 = b0p2.innerProduct(n0).doubleValue() > 0;
		p21 = b1p2.innerProduct(n1).doubleValue() > 0;
		if(p20) {
			arr.add(-1.);
			arr.add(Math.sqrt(p2p1.squaredLength().doubleValue())*p2s.innerProduct(n0).doubleValue() / p2p1.innerProduct(n0).doubleValue());
			arr.add(Math.sqrt(p2p1.squaredLength().doubleValue())*p2s.innerProduct(n1).doubleValue() / p2p1.innerProduct(n1).doubleValue());
		}
		else if (p21){
			arr.add(Math.sqrt(p0p2.squaredLength().doubleValue())*p0s.innerProduct(n0).doubleValue() / p0p2.innerProduct(n0).doubleValue());
			arr.add(-1.);
			arr.add(Math.sqrt(p2p1.squaredLength().doubleValue())*p0s.innerProduct(n1).doubleValue() / p2p1.innerProduct(n1).doubleValue());
		}
		else {
			arr.add(Math.sqrt(p0p2.squaredLength().doubleValue())*p0s.innerProduct(n0).doubleValue() / p0p2.innerProduct(n0).doubleValue());
			arr.add(Math.sqrt(p0p2.squaredLength().doubleValue())*p0s.innerProduct(n1).doubleValue() / p0p2.innerProduct(n1).doubleValue());
			arr.add(-1.);
		}
		return arr;
		
	}

	
}
