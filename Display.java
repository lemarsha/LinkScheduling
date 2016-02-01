import java.awt.Graphics;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JPanel;


public class Display extends JPanel {
	private ArrayList<Node> nodes;
	private ArrayList<Link> links;
	private Set<Link> selectedLinks;
	
	public Display(ArrayList<Node> nodes, ArrayList<Link> links, Set<Link> selectedLinks) {
		super();
		this.nodes = nodes;
		this.links = links;
		this.selectedLinks = selectedLinks;
	}
	
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
	}

}
