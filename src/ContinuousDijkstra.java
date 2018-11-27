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
	}

	public double calculateDistance(Vertex<Point_3> destination) {
		return 0;
	}
}
