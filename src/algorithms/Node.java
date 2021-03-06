package algorithms;
import java.awt.Graphics;


public class Node {
	public static final int SIZE = 6;
	
	private double x, y;
	private double scalar;
	
	public Node(double x, double y) {
		super();
		this.x = x;
		this.y = y;
		this.scalar = 1;
	}
	
	//calculates the distance from one node to another node
	public double distanceToOtherNode(Node other) {
		return Math.sqrt(Math.pow(this.y - other.y, 2) + Math.pow(this.x - other.x, 2));
	}

	//scales the node by a scalar	
	public void scale(double scalar) {
		this.scalar = scalar;
		x /= scalar;
		y /= scalar;
	}

	//returns x coordinate of node	
	public double getX() {
		return scalar*x;
	}
	
	//returns y coordinate of node
	public double getY() {
		return scalar*y;
	}
		
	public void setScalar(double scalar) {
		this.scalar=scalar;
	}
	
	//draws the node!
	public void draw(Graphics g) {
		g.fillRect((int) x-2, (int) y-2, 4, 4);
	}

	//copy node
	public Node copy() {
		Node n = new Node(x,y);
		return n;
	}
	
	@Override
	public String toString() {
		return "(" + x + "," + y + ")";
	}
}
