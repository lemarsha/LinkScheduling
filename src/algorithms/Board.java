package algorithms;

import java.awt.Graphics;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JPanel;


public class Board extends JPanel {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private ArrayList<Link> links_s;
	private Set<Link> selectedLinks;
	private ArrayList<Link> links_p;
	
	public Board(ArrayList<Link> links_s, ArrayList<Link> links_p) {
		super();
		this.links_s = links_s;
		this.links_p = links_p;
		this.selectedLinks = new HashSet<Link>();
	}

	public Board(ArrayList<Link> links_s, ArrayList<Link> links_p, Set<Link> selectedLinks) {
		super();
		this.links_s = links_s;
		this.links_p=links_p;
		this.selectedLinks = selectedLinks;
	}
	
	public void setSelectedLinks(Set<Link> selectedLinks) {
		this.selectedLinks=selectedLinks;
	}
	
	public void clearSelectedLinks() {
		selectedLinks.clear();
	}

	//draws the display
	//first draws all nodes
	//then draws all the links in links
	//then draws all the selected links	
	@Override
	public void paintComponent(Graphics g) {
		g.drawLine(Simulation.WIDTH-1, 0, Simulation.WIDTH-1, Simulation.HEIGHT + 23);
		//draw secondary links
		for(Link l : links_s) {
			l.getSender().draw(g);
			l.getReceiver().draw(g);
			l.draw(g);
		}
		for(Link l : selectedLinks) {
			l.select();
			l.getSender().draw(g);
			l.getReceiver().draw(g);
			l.draw(g);
		}
		for (Link l: links_p) {
			l.setPrimaryLink();
			l.getSender().draw(g);
			l.getReceiver().draw(g);
			l.draw(g);
		}
	}
}
