import java.util.PriorityQueue;
import Jcg.polyhedron.*;
import Jcg.geometry.*;

public class ContinuousDijkstra {
	private SurfaceMesh mesh;

	public ContinuousDijkstra(SurfaceMesh mesh) {
		this.mesh = mesh;
	}

	public void buildDistances(Vertex<Point_3> source) {
		PriorityQueue<Window> pq = new PriorityQueue<>();

		initializePriorityQueue(source, pq);

		while (!pq.isEmpty()) {
			Window window = pq.poll();

			propagateWindow(window, pq);
		}
	}

	private void initializePriorityQueue(
			Vertex<Point_3> source,
			PriorityQueue<Window> pq
	) {
		Halfedge<Point_3> firstHalfedge = source.getHalfedge();
		Halfedge<Point_3> halfedge = firstHalfedge;
		do {
			double lengthHalfedge = calculateLengthHalfedge(halfedge);
			Window fullEdgeWindow =
				new Window(0, lengthHalfedge, lengthHalfedge, 0, 0, halfedge, false);
			pq.add(fullEdgeWindow);

			halfedge = halfedge.getNext().getOpposite();
		} while (halfedge != firstHalfedge);
	}

	private void propagateWindow(Window window, PriorityQueue<Window> pq) {
	}

	private double calculateLengthHalfedge(Halfedge<Point_3> halfedge) {
		Vertex<Point_3> destination = halfedge.getVertex();
		Vertex<Point_3> origin = halfedge.getOpposite().getVertex();

		return (double) destination.getPoint().distanceFrom(origin.getPoint());
	}

	public double calculateDistance(Point_3 destination) {
		if (isVertex(destination)) {
			return 1;
		} else if (isInHalfedge(destination)) {
			return 1;
		} else {
			Halfedge<Point_3> halfedgeMinDistance = getFaceContainingPoint(destination).getEdge();
			Pair<Double, Point_3> distancePair = getMinimumDistanceThroughHalfedge(destination, halfedgeMinDistance);
			double minDistance = distancePair.first();
			Point_3 minDistanceHalfedgePoint = distancePair.second();

			Halfedge<Point_3> halfedge = halfedgeMinDistance;
			for (int i = 0; i < 2; i++) {
				halfedge = halfedge.getNext();
				distancePair = getMinimumDistanceThroughHalfedge(destination, halfedge);

				if (distancePair.first() < minDistance) {
					minDistance = distancePair.first();
					halfedgeMinDistance = halfedge;
				}
			}

			return 1;
		}
	}

	public boolean isVertex(Point_3 point) {
		return true;
	}

	public boolean isInHalfedge(Point_3 point) {
		return true;
	}

	public Face<Point_3> getFaceContainingPoint(Point_3 point) {
		return null;
	}

	public Pair<Double, Point_3> getMinimumDistanceThroughHalfedge(Point_3 point, Halfedge<Point_3> halfedge) {
		return null;
	}
}
