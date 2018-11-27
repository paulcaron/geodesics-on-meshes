import java.util.ArrayList;
import java.util.PriorityQueue;
import Jcg.polyhedron.*;
import Jcg.geometry.*;
import Jama.*;

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
				new Window(0., lengthHalfedge, lengthHalfedge, 0., 0., halfedge, false);
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
	
	public ArrayList<Double> computePoints(Window window) {
		ArrayList<Double> arr = new ArrayList<Double>(3);
		Point_3 p0, p1, p2, b0, b1;
		Point_3 source = window.computeSource();
		p0 = window.getHalfedge().getVertex().getPoint();
		p1 = window.getHalfedge().getOpposite().getVertex().getPoint();
		p2 = window.getHalfedge().getNext().getVertex().getPoint();
		boolean p20, p21; //p20=true if p2 is on the left side of the window : 
		Number[] coefficients0 = {1-window.getBegin(), window.getBegin()};
		Number[] coefficients1 = {window.getEnd(), 1-window.getEnd()};
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
	
	public void propagateWindow(Window window, PriorityQueue<Window>pq) {
		ArrayList<Double> arr = computePoints(window);
		
		
	}
}
