import java.awt.BorderLayout;
import java.awt.Graphics;
import java.awt.GridLayout;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

import javax.swing.JFrame;
import javax.swing.JPanel;

public class MaxLinksDemo3 extends JFrame implements KeyListener {

	public static final int WIDTH = 600;
	public static final int HEIGHT = 600;
	public static final int NUM_NODES = 100;

	private Link primaryLink; //primary link
	private ArrayList<Link> links; //all links
	private ArrayList<Node> nodes; //all nodes
	private HashSet<Link> selected_links; //links selected in algorithm
	private HashSet<Node> selected_nodes; //selected sender nodes

	private double phi;
	private double kappa;
	private double R;
	private double sigma;
	
	private Display display;

	public MaxLinksDemo3() {
		super();
		links=new ArrayList<Link>();
		nodes = new ArrayList<Node>();
		selected_links = new HashSet<Link>();
		selected_nodes = new HashSet<Node>();

		phi = 0.5;
		kappa = 4;
		R = links.get(links.size()-1).getLength()+1;
		sigma = 16;
		
		setUp();

		setSize(WIDTH*2,HEIGHT+23);	
		display = new Display(nodes, links, selected_links, primaryLink);	
		setLayout(new GridLayout(1,2));
		add(display);

	}

	//initialize links
	private void setUp() {
		//clear old vals
		nodes.clear();
		links.clear();
		selected_links.clear();
		selected_nodes.clear();
		//generate new vals
		generatePrimaryLink();
		generateNodes();
		generateLinks();	
		Collections.sort(links);
		scale();
	}

	//generate primary link
	private void generatePrimaryLink(){
		Random r = new Random();
		Node sender = new Node(r.nextInt(WIDTH), r.nextInt(HEIGHT));
		Node receiver = new Node(r.nextInt(WIDTH), r.nextInt(HEIGHT));
		primaryLink = new Link(sender, receiver);
		primaryLink.setPrimaryLink();
	}

	//generate nodes
	private void generateNodes() {
		Random r = new Random();
		for (int i = 0; i<NUM_NODES; i++) {
			Node n = new Node(r.nextInt(WIDTH), r.nextInt(HEIGHT));
			nodes.add(n);
		}
	}

	//generate links
	private void generateLinks() {
		for(int i = 0; i<nodes.size(); i++) {
			for(int j = 0; j<nodes.size(); j++) {
				if (i!=j) {
					Link s = new Link(nodes.get(i), nodes.get(j));
					links.add(s);
				}
			}
		}
	}

	//scale the links
	public void scale() {
		double shortest = links.get(0).getLength();
		if (primaryLink.getLength()<shortest) shortest=primaryLink.getLength();
		for (Link l: links) {
			l.scaleLengths(shortest);
		}
		for (Node n: nodes) {
			n.scale(shortest);
		}	
		primaryLink.scaleLengths(shortest);
		primaryLink.getSender().scale(shortest);
		primaryLink.getReceiver().scale(shortest);
	}

	//removes links in first removal
	private void first_removal(Link a, double rho) {
		ArrayList<Link> S_1 = new ArrayList<Link>();
		for (Link l: links) {
			if (a.getSender().distanceToOtherNode(l.getSender())<rho*a.getLength()) S_1.add(l);
		}
		for (Link l: S_1) {
			links.remove(l);
		}
	}

	//removes links in second removal
	private void second_removal(Link a) {
		ArrayList<Link> S_2 = new ArrayList<Link>();
		double ri = 0;
		if (a.getLength()<primaryLink.getLength()) {
			for (Link l: links) {
				if((relativeInterference(selected_nodes,l)+relativeInterference(primaryLink.getSender(),l))>(1-phi)) S_2.add(l);
			}
		}else{
			for (Link l: links) {
				if(relativeInterference(selected_nodes,l)>(1-phi)) S_2.add(l);
			}
		}
		for (Link l: S_2) {
			links.remove(l);
		}
	}	

	//calculates the relative interference of a set of nodes on another node
	private double relativeInterference(HashSet<Node> W, Link a) {
		double totalInterference = 0;
		for(Node w: W) {
			totalInterference+=relativeInterference(w,a);
		}
		return totalInterference;
	}

	//calculates the relative inteference of a node on another node
	private double relativeInterference(Node w, Link a) {
		return sigma*Math.pow(a.getLength()/w.distanceToOtherNode(a.getReceiver()),kappa)
			/(1-Math.pow(a.getLength()/R,kappa));
	}

	@Override
	public void keyTyped(KeyEvent e) {

	}
	@Override
	public void keyPressed(KeyEvent e) {

	}
	@Override
	public void keyReleased(KeyEvent e) {

	}
	
	public static void main(String[] args) {
		MaxLinksDemo3 test = new MaxLinksDemo3();
		test.setDefaultCloseOperation(EXIT_ON_CLOSE);
		test.setVisible(true);
	}

}

