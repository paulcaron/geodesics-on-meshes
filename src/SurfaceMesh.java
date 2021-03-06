import java.util.ArrayList;

import Jcg.geometry.*;
import Jcg.polyhedron.*;


/**
 * Class for rendering a surface triangle mesh
 */
public class SurfaceMesh {
	
	public double scaleFactor = 50;
	public Point_3 intermediate = null;
	private final MeshViewer view;
	private Polyhedron_3<Point_3> polyhedron3D;
	
	/**
	 * Create a surface mesh from an OFF file
	 */	
	public SurfaceMesh(MeshViewer view, String filename) {
		this.view = view;

		// shared vertex representation of the mesh
    	SharedVertexRepresentation sharedVertex = new SharedVertexRepresentation(filename);
    	LoadMesh<Point_3> load3D = new LoadMesh<Point_3>();
    	
    	polyhedron3D=load3D.createTriangleMesh(sharedVertex.points,sharedVertex.faceDegrees,
				sharedVertex.faces,sharedVertex.sizeHalfedges);

    	polyhedron3D.isValid(false);
    	    	
    	this.scaleFactor = this.computeScaleFactor();
	}
	
	public void drawSegment(Point_3 p, Point_3 q) {
		float s = (float) this.scaleFactor;
		this.view.line(p.getX().floatValue() * s, p.getY().floatValue() * s, 
				p.getZ().floatValue() * s, q.getX().floatValue() * s, 
				q.getY().floatValue() * s, q.getZ().floatValue() * s);
	}

	public void drawTriangle(Point_3 p, Point_3 q, Point_3 r) {
		float s = (float)this.scaleFactor;
		view.vertex((p.getX().floatValue() * s), (p.getY().floatValue() * s), (p.getZ().floatValue() * s));
		view.vertex((q.getX().floatValue() * s), (q.getY().floatValue() * s), (q.getZ().floatValue() * s));
		view.vertex((r.getX().floatValue() * s), (r.getY().floatValue() * s), (r.getZ().floatValue() * s));
	}
	
	public void drawVertex(Point_3 p) {
		float s = (float) this.scaleFactor;
		float x1=(float)p.getX().doubleValue()*s;
		float y1=(float)p.getY().doubleValue()*s;
		float z1=(float)p.getZ().doubleValue()*s;
		
		view.translate(x1, y1, z1);
		view.sphere(s/25f);
		view.translate(-x1, -y1, -z1);
	}

	
	public void draw() {
		this.drawAxis();
		
		view.beginShape(view.TRIANGLES);
		for(Face<Point_3> f : this.polyhedron3D.facets) {
			Halfedge<Point_3> e = f.getEdge();
			Point_3 p = e.vertex.getPoint();
			Point_3 q = e.getNext().vertex.getPoint();
			Point_3 r = e.getNext().getNext().vertex.getPoint();
			
			view.noStroke();
			view.fill(200, 200, 200, 255); // color of the triangle
			this.drawTriangle(p, q, r); // draw a triangle face
		}
		view.endShape();
		
		view.strokeWeight(2); // line width (for edges)
		view.stroke(20);
		for(Halfedge<Point_3> e : this.polyhedron3D.halfedges) {
			Point_3 p = e.vertex.getPoint();
			Point_3 q = e.opposite.vertex.getPoint();
			
			this.drawSegment(p, q); // draw edge (p,q)
		}
		view.strokeWeight(1);
	}
	
	public void drawAxis() {
		double s = 1;
		Point_3 p000 = new Point_3(0, 0, 0);
		Point_3 p100 = new Point_3(s, 0, 0);
		Point_3 p010 = new Point_3(0, s, 0);
		Point_3 p011 = new Point_3(0, 0, s);
		
		drawSegment(p000, p100);
		drawSegment(p000, p010);
		drawSegment(p000, p011);
	}


	public static double round(double x, int precision) {
		return ((int) (x * precision) / (double) precision);
	}
	
	/**
	 * Compute the scale factor (depending on the max distance of the point set)
	 */
	public double computeScaleFactor() {
		if(this.polyhedron3D == null || this.polyhedron3D.vertices.size() < 1)
			return 1;

		double maxDistance = 0;
		Point_3 origin = new Point_3(0, 0, 0);
		for(Vertex<Point_3> v : this.polyhedron3D.vertices) {
			double distance = Math.sqrt(v.getPoint().squareDistance(origin).doubleValue());
			maxDistance = Math.max(maxDistance, distance);
		}

		return Math.sqrt(3) / maxDistance * 150;
	}
	
	public Polyhedron_3<Point_3> getPolyhedron() {
		return polyhedron3D;
	}

	public ArrayList<Face<Point_3>> getFaces() {
		return polyhedron3D.facets;
	}

	public ArrayList<Vertex<Point_3>> getVertices() {
		return polyhedron3D.vertices;
	}

	public ArrayList<Halfedge<Point_3>> getHalfedges() {
		return polyhedron3D.halfedges;
	}
}
