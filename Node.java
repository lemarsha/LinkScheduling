import java.awt.Color;
import java.awt.Graphics;
import java.awt.Point;


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
	
	//draws the node!
	public void draw(Graphics g) {
		g.setColor(Color.blue);
		g.fillOval((int)(scalar*x) - SIZE/2, MaxLinksDemo.HEIGHT - (int)(scalar*y) - SIZE/2 , SIZE, SIZE);
	}
	
	@Override
	public String toString() {
		return "(" + x + "," + y + ")";
	}
}
