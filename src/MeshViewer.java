import processing.core.*;

import Jcg.geometry.*;
import Jcg.polyhedron.*;

/**
 * A simple 3d viewer for visualizing surface meshes
 * 
 * @author Luca Castelli Aleardi (INF555, 2012)
 *
 */
public class MeshViewer extends PApplet {

	SurfaceMesh mesh;
	ContinuousDijkstra continuousDijsktra;
	MeshSimplification ms;
	//String filename="OFF/high_genus.off";
	//String filename="OFF/sphere.off";
	String filename="OFF/cube.off";
	//String filename="OFF/torus_33.off";
	//String filename="OFF/tore.off";
	//String filename="OFF/tri_round_cube.off";
	//String filename="OFF/tri_hedra.off";
	//String filename="OFF/tri_horse.off";
	
	int simplificationMethod=0;
	int nMethods=3; // number of simplification methods proposed
	Point_3 origin;
	Point_3 destination;
	int renderType = 0;
	
	public void setup() {
		size(800,600,P3D);
		ArcBall arcball = new ArcBall(this);
		
		this.mesh=new SurfaceMesh(this, filename);
		this.ms=new EdgeContraction(this.mesh.getPolyhedron());
		//this.ms=new HalfedgeContraction(this.mesh.polyhedron3D);
		//this.ms=new QuadricErrorMetrics(this.mesh.polyhedron3D);
		
		//ms.simplify();
		
		this.continuousDijsktra = new ContinuousDijkstra(this.mesh);

		origin = mesh.getFaces().get(0).getEdge().getVertex().getPoint();
		//origin.multiply(50);
		destination = mesh.getVertices().get(1).getPoint();
		//this.continuousDijsktra.buildDistances(origin);
		//double distance = continuousDijsktra.getDistanceToSource(destination);

		//System.out.printf("Calculated distance is " + distance + "\n");
	}
	
	public void updatedMethod() {
		if(this.simplificationMethod==0) {
		  this.ms=new EdgeContraction(this.mesh.getPolyhedron());
		  System.out.println("Simplification method changed: edge contraction");
		}
		else if(this.simplificationMethod==1) {
			System.out.println("Simplification method changed: halfedge contraction");
		  this.ms=new HalfedgeContraction(this.mesh.getPolyhedron());
		}
		else {
			System.out.println("Simplification method changed: quadric error metrics");
			this.ms=new QuadricErrorMetrics(this.mesh.getPolyhedron());
		}
	}

		 
		public void draw() {
		  background(0);
		  //this.lights();
		  directionalLight(101, 204, 255, -1, 0, 0);
		  directionalLight(51, 102, 126, 0, -1, 0);
		  directionalLight(51, 102, 126, 0, 0, -1);
		  directionalLight(102, 50, 126, 1, 0, 0);
		  directionalLight(51, 50, 102, 0, 1, 0);
		  directionalLight(51, 50, 102, 0, 0, 1);
		  		 
		  translate(width/2.f,height/2.f,-1*height/2.f);
		  this.strokeWeight(1);
		  stroke(150,150,150);
		  //translate(origin.getX().floatValue(), origin.getY().floatValue(), origin.getZ().floatValue());
		  //sphere(10);
		  float s = 50;
		  stroke(200);
		  
		  this.mesh.drawVertex(origin);
		  this.mesh.drawVertex(destination);
		  this.mesh.drawSegment(origin, destination);
		  
		  this.mesh.draw();
		}
		
		
		public void keyPressed(){
			  switch(key) {
			    case('s'):case('S'): this.simplify(); 
			    break;
			    
			    case('c'):this.simplificationMethod=(this.simplificationMethod+1)%this.nMethods; 
			    this.updatedMethod();
			    break;
			  }
		}
		
		public void simplify() {
			this.ms.simplify();
			//this.mesh.updateScaleFactor();
			//this.mesh.polyhedron3D.isValid(false);
		}

		public void getDistanceToSource(Point_3 source, Point_3 destination) {
			continuousDijsktra.buildDistances(source);
			double distance = continuousDijsktra.getDistanceToSource(destination);

			System.out.printf("Calculated distance is " + distance);
		}
		
		/**
		 * For running the PApplet as Java application
		 */
		public static void main(String args[]) {
			//PApplet pa=new MeshViewer();
			//pa.setSize(400, 400);
			PApplet.main(new String[] { "MeshViewer" });
		}
		
}
