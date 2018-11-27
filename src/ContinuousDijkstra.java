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

		while (!pq.empty()) {
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
				new Window(0, lengthHalfedge, lengthHalfedge, 0, 0, halfedge, 0);
			pq.add(fullEdgeWindow);

			halfedge = halfedge.getNext().getOpposite();
		} while (halfedge != firstHalfedge);
	}

	private double calculateLengthHalfedge(Halfedge<Point_3> halfedge) {
		Vertex<Point_3> destination = halfedge.getVertex();
		Vertex<Point_3> origin = halfedge.getOpposite().getVertex();

		return (double) destination.getPoint().distanceFrom(origin.getPoint());
	}

	public double calculateDistance(Vertex<Point_3> destination) {
		return 0;
	}
}
