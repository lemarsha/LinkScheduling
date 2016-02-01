import java.awt.Color;
import java.awt.Graphics;


public class Link implements Comparable {
	private Node sender, receiver;
	private double length;
	private boolean selected = false;
	
	public Link(Node sender, Node receiver) {
		this.sender = sender;
		this.receiver = receiver;
		this.length = sender.distanceToOtherNode(receiver);
	}
	
	public Link(Node sender, Node receiver, double length) {
		this.sender = sender;
		this.receiver = receiver;
		this.length = length;
	}
	
	public Node getSender() {
		return sender;
	}
	
	public Node getReceiver() {
		return receiver;
	}
	
	public double getLength() {
		return length;
	}
	
	public void scaleLengths(double scalar) {
		length /= scalar;
	}
	
	public void draw(Graphics g) {
		if(selected) g.setColor(Color.red);
		else g.setColor(Color.black);
		g.drawLine((int)sender.getX(), MaxLinksDemo.HEIGHT - (int)sender.getY(),
					(int)receiver.getX(), MaxLinksDemo.HEIGHT - (int)receiver.getY());
	}
	
	public void select() {
		selected = true;
	}
	
	public Link copy() {
		return new Link(sender, receiver, length);
	}

	@Override
	public int compareTo(Object o) {
		if(this.length < ((Link) o).length) return -1;
		else if(this.length > ((Link) o).length) return 1;
		return 0;
	}
	
	@Override
	public String toString() {
		return "(" + sender + "," + receiver + ")";
	}

}
