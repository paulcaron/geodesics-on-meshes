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

	public double calculateDistance(Vertex<Point_3> destination) {
		return 0;
	}
}
