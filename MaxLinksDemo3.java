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

	private Link primaryLink, primaryLink2; //primary link
	private ArrayList<Link> links, links2; //all links
	private ArrayList<Node> nodes, nodes2; //all nodes
	private HashSet<Link> selected_links;  //links selected in algorithm
	private HashSet<Node> selected_nodes; //selected sender nodes

	private double phi;
	private double kappa;
	private double R;
	private double sigma;
	private double rho_0;
	private double shortest;
	
	private Display display1, display2;

	public MaxLinksDemo3() {
		super();
		links=new ArrayList<Link>();
		nodes = new ArrayList<Node>();
		selected_links = new HashSet<Link>();
		selected_nodes = new HashSet<Node>();
	
		setUp();
		setUpConstants();

		setSize(WIDTH*2,HEIGHT+23);
		display1 = new Display(nodes, links, selected_links, primaryLink);	
		display2 = new Display(nodes2, links2, selected_links, primaryLink);
		setLayout(new GridLayout(1,2));
		add(display1);
		add(display2);
		addKeyListener(this);

	}

	//initialize links
	private void setUp() {
		//clear old vals
		nodes.clear();
		links.clear();
		//generate new vals
		generatePrimaryLink();
		generateNodes();
		generateLinks();	
		Collections.sort(links);
		scale();
		reset();
	}

	//set up stuff that gets edited
	private void reset() {
		selected_nodes.clear();
		selected_links.clear();
		nodes2 = new ArrayList<Node>();
		for (Node n: nodes) {
			Node n2 = n.copy();
			n2.setScalar(shortest);
			nodes2.add(n2);
		}
		links2 = new ArrayList<Link>();
		for (Link l: links) {
			Link l2 = l.copy();
			links2.add(l2);
		}

	}

	//initialize constants
	private void setUpConstants() {
		phi = 0.5;
		kappa = 4;
		R = links.get(links.size()-1).getLength()+1;
		sigma = 16;
		rho_0=1+Math.pow((sigma*(16*1.202+8*Math.pow(Math.PI,4)/90-6)/phi),1/kappa);
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
		shortest = links.get(0).getLength();
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
	private void firstRemoval(Link a, double rho) {
		ArrayList<Link> S_1 = new ArrayList<Link>();
		for (Link l: links2) {
			if (a.getSender().distanceToOtherNode(l.getSender())<rho*a.getLength()) S_1.add(l);
		}
		for (Link l: S_1) {
			links2.remove(l);
		}
	}

	//removes links in second removal
	private void secondRemoval(Link a, boolean use_primary) {
		ArrayList<Link> S_2 = new ArrayList<Link>();
		double ri = 0;
		if (a.getLength()<primaryLink.getLength() && use_primary) {
			for (Link l: links2) {
				if((relativeInterference(selected_nodes,l)+relativeInterference(primaryLink.getSender(),l))>(1-phi)) S_2.add(l);
			}
		}else{
			for (Link l: links2) {
				if(relativeInterference(selected_nodes,l)>(1-phi)) S_2.add(l);
			}
		}
		for (Link l: S_2) {
			links2.remove(l);
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

	//original algorithm
	private void algorithmAAMISL() {
		Link a; 
		double rho;
		while (links2.size()>0) {
			a=links2.get(0);
			a.select();
			links2.remove(0);
			selected_links.add(a);
			selected_nodes.add(a.getSender());
			rho = 1+(rho_0-1)/Math.pow(1-Math.pow(a.getLength()/R,kappa),1/kappa);
			firstRemoval(a,rho);
			secondRemoval(a,false);
		}
	}

	//algorithm with primary link
	private void algorithmPLMISL() {
		Link a;
		double rho;
		boolean first=true;
		while(links2.size()>0) {
			a=links2.get(0);
			links2.remove(0);
			if (a.getLength()<primaryLink.getLength()) {
				if (relativeInterference(a.getSender(),primaryLink)+relativeInterference(selected_nodes,primaryLink)<=(1-phi)) {
					a.select();
					selected_links.add(a);
					selected_nodes.add(a.getSender());
					rho = 1+(rho_0-1)/Math.pow(1-Math.pow(a.getLength()/R,kappa),1/kappa);
					firstRemoval(a,rho);
					secondRemoval(a,true);
				}
			}
			else {
				if (first) {
					first = false;
					a=primaryLink;
					selected_links.add(a);
					selected_nodes.add(a.getSender());
					rho = 1+(rho_0-1)/Math.pow(1-Math.pow(a.getLength()/R,kappa),1/kappa);
					firstRemoval(a,rho);
					secondRemoval(a,true);
				}
				else {
					a.select();
					selected_links.add(a);
					selected_nodes.add(a.getSender());
					rho = 1+(rho_0-1)/Math.pow(1-Math.pow(a.getLength()/R,kappa),1/kappa);
					firstRemoval(a,rho);
					secondRemoval(a,true);

				}
	
			}
		}
	}

	//Circle algorithm
	public void algorithmCMISL() {
		Link a; 
		double rho;
		selected_links.add(primaryLink);
		selected_nodes.add(primaryLink.getSender());
		double rho_s = 1+(rho_0-1)/Math.pow(1-Math.pow(1/R,kappa),1/kappa);
		double c = 24*sigma*Math.pow(primaryLink.getLength(),kappa);
		c/=( (kappa-2)*(1-Math.pow(primaryLink.getLength()/R,kappa))*Math.pow(rho_s,kappa) );
		c = Math.ceil(Math.pow(c,1/(kappa-2))+0.5);
		ArrayList<Link> S_1 = new ArrayList<Link>();
		ArrayList<Link> S_2 = new ArrayList<Link>();
		for (Link l: links2) {
			if (primaryLink.getSender().distanceToOtherNode(l.getSender())<c*rho_s) S_1.add(l);
		}
		for(Link l: S_1) {
			links2.remove(l);
		}
		for (Link l: links2) {
			if (relativeInterference(primaryLink.getSender(),l)>1-phi) S_2.add(l);
		}
		for (Link l: S_2) {
			links2.remove(l);
		}
		
		while (links2.size()>0) {
			a=links2.get(0);
			a.select();
			links2.remove(0);
			selected_links.add(a);
			selected_nodes.add(a.getSender());
			rho = 1+(rho_0-1)/Math.pow(1-Math.pow(a.getLength()/R,kappa),1/kappa);
			firstRemoval(a,rho);
			secondRemoval(a,false);
		}
	}	

	//verifies independence
	public boolean checkIndependent() {
		HashSet<Node> all_selected_nodes = new HashSet<Node>();
		for (Link l: selected_links) {
			//check that all nodes unique
			if (all_selected_nodes.contains(l.getSender())) return false;
			all_selected_nodes.add(l.getSender());
			if (all_selected_nodes.contains(l.getReceiver())) return false;
			all_selected_nodes.add(l.getReceiver());
			//check interference
			selected_nodes.remove(l.getSender());
			double ri = relativeInterference(selected_nodes,l);
			selected_nodes.add(l.getSender());
			if (ri>1-phi) return false;
		}
		return true;
	}


	@Override
	public void keyPressed(KeyEvent e) {
		if (e.getKeyChar()=='a') {
			reset();
			algorithmAAMISL();
		}else if (e.getKeyChar()=='s') {
			reset();
			algorithmPLMISL();
		}else if (e.getKeyChar()=='d'){
			reset();
			algorithmCMISL();
		}else {
			reset();
		}
		System.out.println(checkIndependent());
		repaint();
	}
	@Override
	public void keyTyped(KeyEvent e) {

	}
	@Override
	public void keyReleased(KeyEvent e) {

	}
	
	public HashSet<Node> getSelectedNodes() {
		return selected_nodes;
	}

	public HashSet<Link> getSelectedLinks() {
		return selected_links;
	}

	public static void main(String[] args) {
		MaxLinksDemo3 test = new MaxLinksDemo3();
		test.setDefaultCloseOperation(EXIT_ON_CLOSE);
		test.setVisible(true);
	}

}

