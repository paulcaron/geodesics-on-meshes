import java.util.ArrayList;
import java.util.PriorityQueue;

import java.util.HashMap;
import java.lang.NullPointerException;

import Jcg.polyhedron.*;
import Jcg.geometry.*;

public class ContinuousDijkstra {
	private final SurfaceMesh mesh;
	private HashMap<Halfedge<Point_3>, ArrayList<Window>> halfedgeToWindowsList;
	int var;

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

		var = 0;
		while (!pq.isEmpty()) {
			Window window = pq.poll();
			if (isValidWindow(window)) {
				System.out.println(var++);
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
			
			Halfedge<Point_3> oppositeHalfedge = halfedge.getOpposite().getNext();
			Window oppositeWindow = new Window(
					0,
					GeoUtils.getHalfedgeLength(oppositeHalfedge),
					GeoUtils.getHalfedgeLength(halfedge),
					GeoUtils.getHalfedgeLength(oppositeHalfedge.getNext()),
					0,
					oppositeHalfedge,
					true);
			if (halfedgeToWindowsList.containsKey(oppositeHalfedge)) {
				windowList = halfedgeToWindowsList.get(oppositeHalfedge);
			}
			else {
				windowList = halfedgeToWindowsList.get(oppositeHalfedge.getOpposite());
				oppositeWindow = oppositeWindow.getOpposite();
			}
			pq.add(oppositeWindow);
			windowList.add(oppositeWindow);
			
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
		else 
			windowsList = halfedgeToWindowsList.get(halfedge.getOpposite());
		
		assert windowsList.size() > 0;
		
		Pair<Double, Point_3> optimalPair = null;
		for (int i = 0; i < windowsList.size(); i++) {
			Pair<Double, Point_3> distancePoint = getDistanceToSourcePassingOverWindow(destination, windowsList.get(i));

			double minDistanceForWindow = distancePoint.first();
			if (minDistanceForWindow < distanceToSource) {
				distanceToSource = minDistanceForWindow;
				optimalPair = distancePoint;
			}
		}

		assert !GeoUtils.isNegative(distanceToSource);
		this.mesh.intermediate = optimalPair.second();
		return distanceToSource;
	}

	private Pair<Double, Point_3> getDistanceToSourcePassingOverWindow(Point_3 destination, Window window) {
		double distanceToSource = Double.MAX_VALUE;

		Halfedge<Point_3> halfedge = window.getHalfedge();
		double halfedgeLength = GeoUtils.getHalfedgeLength(halfedge);

		double dx = halfedgeLength / 100;
		double optimalX = window.getStart();
		for (double x = window.getStart(); x <= window.getEnd(); x += dx) {
			Point_3 halfedgePoint = GeoUtils.getHalfedgePoint(halfedge, x);

			double distance = destination.distanceFrom(halfedgePoint).doubleValue() + window.getDistSource();
			if (GeoUtils.isEqual(x, window.getStart()))
				distance += window.getDistStart();
			else if (GeoUtils.isEqual(x, window.getEnd()))
				distance += window.getDistEnd();
			else
				distance += GeoUtils.getCevianLength(window.getDistStart(), window.getDistEnd(), x - window.getStart(), window.getEnd() - x);

			if (distance < distanceToSource) {
				distanceToSource = distance;
				optimalX = x;
			}
		}
		assert !GeoUtils.isNegative(distanceToSource);

		return new Pair<Double, Point_3> (distanceToSource, GeoUtils.getHalfedgePoint(window.getHalfedge(), optimalX));
	}

	private Face<Point_3> getFaceContainingPoint(Point_3 point) {
		for (Face<Point_3> polyhedronFace : mesh.getFaces()) {
			if (GeoUtils.isPointOnFace(point, polyhedronFace))
				return polyhedronFace;
		}
		
		throw new IllegalArgumentException("Point " + point.toString() + " not contained in any face of the mesh");
	}

	private void propagateWindow(Window window, PriorityQueue<Window> pq) {
		if (GeoUtils.isZero(window.getLength()))
			return;
		
		if (window.getSide())
			window = window.getOpposite();
		assert !window.getSide();

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
			assert !GeoUtils.isPositive(thirdExtremity - GeoUtils.getHalfedgeLength(halfedge));

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
						secondExtremity,
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

			Halfedge<Point_3> firstHalfedge = propagatingHalfedge.getNext();

			assert firstExtremity < GeoUtils.getHalfedgeLength(firstHalfedge);
			
			Window firstPropagatedWindow = new Window(
					firstExtremity,
					GeoUtils.getHalfedgeLength(firstHalfedge),
					GeoUtils.getHalfedgePoint(firstHalfedge, firstExtremity).distanceFrom(
						GeoUtils.getHalfedgePoint(propagatingHalfedge, window.getEnd())).doubleValue() + window.getDistEnd(),
					Math.sqrt(p2s.squaredLength().doubleValue()), 
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
					Math.sqrt(p2s.squaredLength().doubleValue()),
					GeoUtils.getHalfedgePoint(propagatingHalfedge, window.getStart()).distanceFrom(
						GeoUtils.getHalfedgePoint(secondHalfedge, thirdExtremity)).doubleValue() + window.getDistStart(),
					window.getDistSource(),
					secondHalfedge,
					true);
			insertWindow(secondPropagatedWindow, pq);
		}
	}

	private void insertWindow(Window newWindow, PriorityQueue<Window> pq) {
		
		if(!newWindow.isValid()) return;
		
		System.out.println("in");
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
			if (!newWindow.isValid()) {
				newWindow = null;
				break;
			}
			Window oldWindow = oldWindows.get(i);
			boolean pqContainsOldWindow = pq.contains(oldWindow);
			if (oldWindow.getEnd() <= newWindow.getStart()) {
				newWindows.add(oldWindow);
				i++;
			} else if (newWindow.getEnd() <= oldWindow.getStart()) {
				break;
			} 
			else if(GeoUtils.isZero(newWindow.getStart() - oldWindow.getStart())) {
				
				if(GeoUtils.isZero(newWindow.getEnd() - oldWindow.getEnd()) &&
						GeoUtils.isZero(newWindow.getDistSource() - oldWindow.getDistSource()) &&
						GeoUtils.isZero(newWindow.getDistStart() - oldWindow.getDistStart()) &&
						GeoUtils.isZero(newWindow.getDistEnd() - oldWindow.getDistEnd())) {
					newWindow = null;
					break;
				}
				
				double p = mergeWindows(newWindow, oldWindow);
				if(p < 0) {
					if(newWindow.compareTo(oldWindow) < 0) {
						
						if(newWindow.getEnd() < oldWindow.getEnd()) {
							if(pqContainsOldWindow) pq.remove(oldWindow);
							oldWindow = oldWindow.setStart(newWindow.getEnd());
							oldWindows.set(i, oldWindow);
							if(pqContainsOldWindow) pq.add(oldWindow);
							
							break;
						}
						else {
							pq.add(newWindow);
							if(pqContainsOldWindow) pq.remove(oldWindow);
							i++;
						}
					}
					else {
						if(newWindow.getEnd() <= oldWindow.getEnd()) {
							newWindow = null;
							break;
						}
						else {
							newWindow = newWindow.setStart(oldWindow.getEnd());
							newWindows.add(oldWindow);
							i++;
						}
					}
				}
				else {
					if(newWindow.getDistSource() + newWindow.getDistStart() < oldWindow.getDistSource() + oldWindow.getDistStart()) {
						if(newWindow.getEnd() < oldWindow.getEnd()) {
							newWindow = newWindow.setEnd(p);
							if(pqContainsOldWindow) pq.remove(oldWindow);
							oldWindow = oldWindow.setStart(p);
							if(pqContainsOldWindow) pq.add(oldWindow);
							break;
						}
						else if(newWindow.getEnd() > oldWindow.getEnd()) {
							newWindows.add(newWindow.setEnd(p));
							pq.add(newWindow.setEnd(p));
							if(pqContainsOldWindow) pq.remove(oldWindow);
							oldWindow = oldWindow.setStart(p);
							if(pqContainsOldWindow) pq.add(oldWindow);
							newWindows.add(oldWindow);
							newWindow = newWindow.setStart(oldWindow.getEnd());
							i++;
						} 
						else {
							newWindow = newWindow.setEnd(p);
							if(pqContainsOldWindow) pq.remove(oldWindow);
							oldWindow = oldWindow.setStart(p);
							if(pqContainsOldWindow) pq.add(oldWindow);
							i++;
							break;
						}
					}
					else {
						if(newWindow.getEnd() < oldWindow.getEnd()) {
							newWindow.setStart(p);
							if(pqContainsOldWindow) pq.remove(oldWindow);
							Window newOldWindow = oldWindow;
							oldWindow = oldWindow.setEnd(p);
							newOldWindow.setStart(newWindow.getEnd());
							if(pqContainsOldWindow) {
								pq.add(oldWindow);
								pq.add(newOldWindow);
							}
							newWindows.add(oldWindow);
							newWindows.add(newWindow);
							pq.add(newWindow);
							newWindows.add(newOldWindow);
							i++;
							newWindow = null;
							break;
						}
						else {
							newWindow = newWindow.setStart(p);
							if(pqContainsOldWindow) pq.remove(oldWindow);
							oldWindow = oldWindow.setEnd(p);
							if(pqContainsOldWindow) pq.add(oldWindow);
							newWindows.add(oldWindow);
							i++;
						}
						

					}					
				}
				
			}
			else if (newWindow.getStart() <= oldWindow.getStart()) {
				Window auxWindow = newWindow.setEnd(oldWindow.getStart());
				pq.add(auxWindow);
				newWindows.add(auxWindow);
				newWindow = newWindow.setStart(oldWindow.getStart());
			}  
			else {
				if(pqContainsOldWindow) pq.remove(oldWindow);
				Window adjustedWindow = oldWindow.setEnd(newWindow.getStart());
				newWindows.add(adjustedWindow);
				oldWindow = oldWindow.setStart(newWindow.getStart());
				oldWindows.set(i, oldWindow);
				if(pqContainsOldWindow) {
					pq.add(adjustedWindow);
					pq.add(oldWindow);
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
		System.out.println("out");
	}
	
	private double mergeWindows(Window window1, Window window2){ // leftWindows and rightWindows have the same start
		
		assert GeoUtils.isZero(window1.getStart() - window2.getStart()) &&
				!areWindowsDisjoint(window1, window2);
		
		double s0x = window1.getAbscisseSourcePlaneProjection() + window1.getStart(); 
		double s1x = window2.getAbscisseSourcePlaneProjection() + window2.getStart();;
		
		double s0y = GeoUtils.getTriangleHeight(window1.getLength(), window1.getDistStart(), window1.getDistEnd());
		double s1y = GeoUtils.getTriangleHeight(window2.getLength(), window2.getDistStart(), window2.getDistEnd());

		double s0squared = s0x * s0x + s0y * s0y;
		double s1squared = s1x * s1x + s1y * s1y;
		
		double alpha = s1x - s0x;
		double beta = window2.getDistSource() - window1.getDistSource();
		double gamma = s0squared - s1squared - beta * beta;
		
		double A = alpha * alpha - beta * beta;
		double B = gamma * alpha + 2 * s1x * beta * beta;
		double C = 0.25 * gamma * gamma - s1squared * beta * beta;
		
		double minPossible = Math.max(window1.getStart(), window2.getStart());
		double maxPossible = Math.min(window1.getEnd(), window2.getEnd());
		double solution = getSolutionInRange(GeoUtils.solveSecondDegreeEquation(B / A, C / A), minPossible, maxPossible);
		
		return solution;
	}

	private double getSolutionInRange(Pair<Double, Double> solutions, double minPossible, double maxPossible) {
		double first = solutions.first();
		double second = solutions.second();
		if (first > minPossible && first < maxPossible)
			return first;
		if (second > minPossible && second < maxPossible)
			return second;

		return -1;
	}

	private boolean areWindowsDisjoint(Window leftWindow, Window rightWindow) {
		return leftWindow.getEnd() <= rightWindow.getStart();
	}

	private boolean isSpecialVertex(Vertex<Point_3> vertex) {
		return true;
	}

	public ArrayList<Double> getPropagatedExtremities(Window window) {
		assert !window.getSide();
		ArrayList<Double> arr = new ArrayList<Double>(3);
		
		

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
		Vector_2 p0p1 = new Vector_2(p0, p1);
		Vector_2 p2p1 = new Vector_2(p2, p1);
		Vector_2 p0s = new Vector_2(p0, source);
		Vector_2 p2s = new Vector_2(p2, source);
		Vector_2 n0 = new Vector_2(-b0s.getY().doubleValue(), b0s.getX());
		Vector_2 n1 = new Vector_2(-b1s.getY().doubleValue(), b1s.getX());
		boolean p20, p21; //p20=true if p2 is on the left side of the window
		p20 = b0p2.innerProduct(n0).doubleValue() > 0;
		p21 = b1p2.innerProduct(n1).doubleValue() > 0;
		
		if(GeoUtils.isZero(Math.abs(p0s.innerProduct(p0p1).doubleValue()) - Math.sqrt(p0s.squaredLength().doubleValue()*p0p1.squaredLength().doubleValue())) || GeoUtils.isZero(b0p2.innerProduct(n0).doubleValue()) || GeoUtils.isZero(b1p2.innerProduct(n1).doubleValue())  ){
			arr.add(0.);
			arr.add(-1.);
			arr.add(GeoUtils.getHalfedgeLength(window.getHalfedge().getNext().getNext()));
		}
	
		
		else if(p20) {
			arr.add(-1.);
			arr.add(Math.sqrt(p2p1.squaredLength().doubleValue()) * p2s.innerProduct(n0).doubleValue() / p2p1.innerProduct(n0).doubleValue());
			arr.add(Math.sqrt(p2p1.squaredLength().doubleValue()) * p2s.innerProduct(n1).doubleValue() / p2p1.innerProduct(n1).doubleValue());
			assert arr.get(1) < arr.get(2);
		}
		else if (p21){
			arr.add(Math.sqrt(p0p2.squaredLength().doubleValue()) * p0s.innerProduct(n0).doubleValue() / p0p2.innerProduct(n0).doubleValue());
			arr.add(-1.);
			arr.add(Math.sqrt(p2p1.squaredLength().doubleValue()) * p2s.innerProduct(n1).doubleValue() / p2p1.innerProduct(n1).doubleValue());
		}
		else {
			arr.add(Math.sqrt(p0p2.squaredLength().doubleValue()) * p0s.innerProduct(n0).doubleValue() / p0p2.innerProduct(n0).doubleValue());
			arr.add(Math.sqrt(p0p2.squaredLength().doubleValue()) * p0s.innerProduct(n1).doubleValue() / p0p2.innerProduct(n1).doubleValue());
			arr.add(-1.);
			assert arr.get(0) < arr.get(1);
		}
		
		assert !Double.isNaN(arr.get(0));
		assert !Double.isNaN(arr.get(1));
		assert !Double.isNaN(arr.get(2));

		return arr;
		
	}

	
}
