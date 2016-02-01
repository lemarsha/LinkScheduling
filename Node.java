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
	
	public double distanceToOtherNode(Node other) {
		return Math.sqrt(Math.pow(this.y - other.y, 2) + Math.pow(this.x - other.x, 2));
	}
	
	public void scale(double scalar) {
		this.scalar = scalar;
		x /= scalar;
		y /= scalar;
	}
	
	public double getX() {
		return scalar*x;
	}
	
	public double getY() {
		return scalar*y;
	}
	
	public void draw(Graphics g) {
		g.setColor(Color.blue);
		g.fillOval((int)(scalar*x) - SIZE/2, MaxLinksDemo.HEIGHT - (int)(scalar*y) - SIZE/2 , SIZE, SIZE);
	}
	
	@Override
	public String toString() {
		return "(" + x + "," + y + ")";
	}
}
