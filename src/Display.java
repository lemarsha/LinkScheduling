import java.awt.Graphics;
import java.util.ArrayList;
import java.util.Set;

import javax.swing.JPanel;


public class Display extends JPanel {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private ArrayList<Node> nodes;
	private ArrayList<Link> links;
	private Set<Link> selectedLinks;
	private Link primaryLink;

	public Display(ArrayList<Node> nodes, ArrayList<Link> links, Set<Link> selectedLinks, Link primaryLink) {
		super();
		this.nodes = nodes;
		this.links = links;
		this.selectedLinks = selectedLinks;
		this.primaryLink=primaryLink;
		if (primaryLink==null) System.out.println("hi");
	}

	public Display(ArrayList<Node> nodes, ArrayList<Link> links, Set<Link> selectedLinks) {
		super();
		this.nodes=nodes;
		this.links=links;
		this.selectedLinks=selectedLinks;

	}


	//draws the display
	//first draws all nodes
	//then draws all the links in links
	//then draws all the selected links	
	@Override
	public void paintComponent(Graphics g) {
		g.drawLine(MaxLinksDemo.WIDTH-1, 0, MaxLinksDemo.WIDTH-1, MaxLinksDemo.HEIGHT + 23);
		for(Node n : nodes) {
			n.draw(g);
		}
		for(Link l : links) {
			l.draw(g);
		}
		for(Link l : selectedLinks) {
			l.draw(g);
		}
		if (primaryLink!=null) {
			primaryLink.getSender().draw(g);
			primaryLink.getReceiver().draw(g);
			primaryLink.draw(g);
		} else {
			System.out.println("No PL");
		}
		primaryLink.draw(g);
	}

}
