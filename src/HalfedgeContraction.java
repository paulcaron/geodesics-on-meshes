import Jcg.geometry.*;
import Jcg.polyhedron.*;

/**
 * @author Luca Castelli Aleardi (INF555, 2012)
 *
 */
public class HalfedgeContraction extends MeshSimplification {
	
	public HalfedgeContraction(Polyhedron_3<Point_3> polyhedron3D) {
		super(polyhedron3D);
	}
	
	/**
	 * Basic example of simplification based on edge contractions
	 * Simply select at random edges to be contracted
	 */
	public void simplify() {
		System.out.println("To be completed");		
	}
	
	
	/**
	 * Check whether a given halfedge can be contracted
	 */
	boolean isLegal(Halfedge<Point_3> h){
		throw new Error("To be completed");
	}
	
}
