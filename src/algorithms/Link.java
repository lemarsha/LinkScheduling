package algorithms;
import java.awt.Color;
import java.awt.Graphics;


@SuppressWarnings("rawtypes")
public class Link implements Comparable {
	private Node sender, receiver;
	private double length;
	private boolean selected = false;
	private boolean isPrimaryLink = false;
	private double tolerance;
	
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

	//returns the node that is sending the signal	
	public Node getSender() {
		return sender;
	}
	
	//returns the node that is recieving the signal
	public Node getReceiver() {
		return receiver;
	}
	
	//returns the length of the node
	public double getLength() {
		return length;
	}

	public double getTolerance() {
		return tolerance;
	}
	
	public void setTolerance(double tolerance) {
		this.tolerance=tolerance;
	}
	
	public void decrementTolerance(double tolerance) {
		this.tolerance-=tolerance;
	}
	
	//draws the link!	
	public void draw(Graphics g) {
		if (isPrimaryLink) g.setColor(Color.RED);
		else if (selected) g.setColor(Color.GREEN);
		g.drawLine((int) sender.getX(), (int) sender.getY(), (int) receiver.getX(), (int) receiver.getY());
		g.setColor(Color.BLACK);
	}

	//indicates if link was selected in algo	
	public void select() {
		selected = true;
	}

	//makes the link the Primary Link
	public void setPrimaryLink() {
		isPrimaryLink=true;
	}

	//creates new link with same variables as link calling function	
	public Link copy() {
		return new Link(sender, receiver, length);
	}

	//compares the lenght of the links
	//returns 0 if the lenghts of two links are equal
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
