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
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.JPanel;


public class MaxLinksDemo2 extends JFrame implements KeyListener {
	public static final int WIDTH = 600;
	public static final int HEIGHT = 600;
	public static final int NUM_NODES = 30;

	private Link primaryLink; //primary link
	private ArrayList<Link> S1; //S1 is list of all possible nodes
	private Set<Link> I1; //I1 is list of all selected nodes
	private Set<Node> U1; //U1 is the list of all the senders
	private ArrayList<Node> all; //holds all the nodes
	private Display display1, display2;
	//R is the length of the longest node + 1
	private double P, R;
	private double nu = 0.99;
	private double xi = 0.01;
	private double kappa = 4;
	private double sigma = 16;
	private double phi;
	private double rho_01, rho_02; //rho_0 constant

	public MaxLinksDemo2() {
		super();
		//initialize variables
		S1 = new ArrayList<Link>();
		I1 = new HashSet<Link>();
		U1 = new HashSet<Node>();
		all = new ArrayList<Node>();
		P = 1000;
		phi = calcPhi();
		System.out.println("phi = " + phi);
		rho_01 = 1 + Math.pow((sigma*(16*1.202 + 8*Math.pow(Math.PI, 4)/90 - 6)/phi), 1/kappa);
		rho_02 = rho_01 - 1;
		//generate nodes and add all possible links to original set
		addElements();
		//sort links from smallest to largest
		Collections.sort(S1);
		//scales the links by the lenght of the smallest link
		scale();
		//set R to lenght of longest node + 1
		R = S1.get(S1.size()-1).getLength() + 1;
		//set display parameters
		setSize(WIDTH*2, HEIGHT + 23);
		//Display is jpanel that shows all the links
		//S1=all links, I1=selected links
		display1 = new Display(all, S1, I1, primaryLink);
		display2 = new Display(all, S1, I1, primaryLink);
		setLayout(new GridLayout(1,2));
		add(display1);
		add(display2);
		addKeyListener(this);
	}

	//rho of 0 is a constant used in the rho calc
	private double rho_0(double phi) {
		return 1 + Math.pow((sigma*(16*1.202 + 8*Math.pow(Math.PI, 4)/90 - 6)/phi), 1/kappa);
	}

	private double f1(double phi) {
		double rho_0 = rho_0(phi);
		double denom = Math.pow(sigma, 1/kappa) - 1;
		return Math.PI*2/Math.pow(3, 0.5)*Math.pow(rho_0/denom, 2) + Math.PI*rho_0/denom + 1;
	}

	private double f2(double phi) {
		return 5*Math.pow(Math.pow(1-phi, -1/kappa) + 2*Math.pow(sigma, -1/kappa), kappa);
	}

	private double f(double phi) {
		return f1(phi) + f2(phi);
	}

	//calculates the optimal value for phi
	private double calcPhi() {
		double low = 0;
		double high = 1, highValue = f(high);
		double mid = 0.5, midValue = f(mid);
		while(midValue > highValue) {
			low = mid;
			mid = (high + low)/2;
			midValue = f(mid);
		}
		while(high - mid > 0.001) {
			double temp = (mid + high)/2;
			if(f(temp) > midValue) {
				high = temp;
				highValue = f(temp);
			}
			else {
				low = mid;
				mid = temp;
				midValue = f(temp);
			}
		}
		return mid;
	}

	//this function creates n nodes that each have random positions 
	//then it adds each possible combination of links (n^2) to S1
	private void addElements() {
		Random randomGenerator = new Random();
		for(int i=0; i<NUM_NODES; i++) {
			all.add(new Node(randomGenerator.nextInt(WIDTH), randomGenerator.nextInt(HEIGHT)));
		}
		for(int i=0; i<all.size(); i++) {
			for(int j=0; j<all.size(); j++) {
				if(i != j) S1.add(new Link(all.get(i), all.get(j)));
			}
		}
		generatePrimaryLink();
	}

	private void generatePrimaryLink() {
		Random randomGenerator = new Random();
		Node sender = new Node(randomGenerator.nextInt(WIDTH), randomGenerator.nextInt(HEIGHT));
		Node receiver = new Node(randomGenerator.nextInt(WIDTH), randomGenerator.nextInt(HEIGHT));
		primaryLink = new Link(sender,receiver);
		primaryLink.setPrimaryLink();
	}

	//This function scales all of the lengths of the links by the smallest one
	//also scales position of node but idk why
	public void scale() {
		double shortest = S1.get(0).getLength();
		if (primaryLink.getLength()<shortest) shortest=primaryLink.getLength();
		for(Link l : S1) {
			l.scaleLengths(shortest);
		}
		for(Node n : all) {
			n.scale(shortest);
		}
		primaryLink.scaleLengths(shortest);
		primaryLink.getSender().scale(shortest);
		primaryLink.getReceiver().scale(shortest);
		P /= shortest;
	}

	//prints all the links in S1
	public void printS() {
		for(Link l : S1) {
			System.out.println(l);
		}	
		System.out.println(primaryLink);
	}

	//selects all links in I1, sets color to red on screen
	public void showSelectedLinks() {
		for(Link l : I1) {
			l.select();
		}
		repaint();
	}

	//Approximation Algorithm for Maximum Independent Set of Links
	public void algorithmAAMISL() {
		Link a;
		double rho;
		//gets smallest link in S1
		//removes a from S1, adds it to I1, and adds sender to U1
		//calculates rho of the link
		//removes links that don't have mutual distance of rho||uv|| or RI(U;l)>1-phi
		//continues until S1 contains no more links
		while(S1.size() > 0) {
			a = S1.get(0);
			S1.remove(0);
			I1.add(a);
			U1.add(a.getSender());
			rho = 1 + (rho_01 - 1)/Math.pow(1 - Math.pow(a.getLength()/R, kappa), 1/kappa);
			ArrayList<Link> S_1 = new ArrayList<Link>(); //links removed in first itt
			ArrayList<Link> S_2 = new ArrayList<Link>(); //links removed in second itt
			//adds link to S_1 if distance from senders is less than rho||uv||
			for(Link l : S1) {
				if(a.getSender().distanceToOtherNode(l.getSender()) < rho * a.getLength()) {
					S_1.add(l);
				}
			}
			//removes links that are in S_1 from S1
			for(Link l : S_1) {
				S1.remove(l);
			}
			//adds link to S_2 if RI(U;l) is greater than 1 minus phi
			for(Link l : S1) {
				if(relativeInterference(U1, l) > 1-phi) S_2.add(l);
			}
			//removes links that are in S_2 from S1
			for(Link l : S_2) {
				S1.remove(l);
			}
		}
	}

	//Algorithm for MISL given existing Link
	private void algorithmELMISL() {
		Link a;
		double rho;
		boolean first=true;
		while(S1.size()>0) {
			a = S1.get(0);
			S1.remove(0);
			if (a.getLength()<primaryLink.getLength()) {
				if (relativeInterference(a.getSender(),primaryLink)+relativeInterference(U1,primaryLink)<=(1-phi)) {
					I1.add(a);
					U1.add(a.getSender());	
					rho = 1+(rho_01-1)/Math.pow(1-Math.pow(a.getLength()/R,kappa), 1/kappa);
					ArrayList<Link> S_1 = new ArrayList<Link>();
					ArrayList<Link> S_2 = new ArrayList<Link>();

					for (Link l: S1) {
						if(a.getSender().distanceToOtherNode(l.getSender())<rho*a.getLength()) S_1.add(l);
					}
					for (Link l: S_1) {
						S1.remove(l);
					}
					for (Link l: S1) {
						if((relativeInterference(U1,l)+relativeInterference(primaryLink.getSender(),l))>(1-phi)) S_2.add(l);
					}
					for (Link l: S_2) {
						S1.remove(l);
					}
				}
			}else {
				if (first) {
					first = false;
					a=primaryLink;
					I1.add(a);
					U1.add(a.getSender());	
					rho = 1+(rho_01-1)/Math.pow(1-Math.pow(a.getLength()/R,kappa), 1/kappa);
					ArrayList<Link> S_1 = new ArrayList<Link>();
					ArrayList<Link> S_2 = new ArrayList<Link>();

					for (Link l: S1) {
						if(a.getSender().distanceToOtherNode(l.getSender())<rho*a.getLength()) S_1.add(l);
					}
					for (Link l: S_1) {
						S1.remove(l);
					}
					for (Link l: S1) {
						if(relativeInterference(U1,l)>(1-phi)) S_2.add(l);
					}
					for (Link l: S_2) {
						S1.remove(l);
					}
				}else {
					I1.add(a);
					U1.add(a.getSender());	
					rho = 1+(rho_01-1)/Math.pow(1-Math.pow(a.getLength()/R,kappa), 1/kappa);
					ArrayList<Link> S_1 = new ArrayList<Link>();
					ArrayList<Link> S_2 = new ArrayList<Link>();

					for (Link l: S1) {
						if(a.getSender().distanceToOtherNode(l.getSender())<rho*a.getLength()) S_1.add(l);
					}
					for (Link l: S_1) {
						S1.remove(l);
					}
					for (Link l: S1) {
						if(relativeInterference(U1,l)>(1-phi)) S_2.add(l);
					}
					for (Link l: S_2) {
						S1.remove(l);
					}
				}
			}
		}
	}
	
	//Calculates the RI of a set of nodes, W on a link, a
	private double relativeInterference(Set<Node> W, Link a) {
		double totalInterference = 0;
		for(Node w : W) {
			totalInterference += relativeInterference(w, a);
		}
		return totalInterference;
	}
	
	//Calculates the RI of a node, w, on a link, a
	private double relativeInterference(Node w, Link a) {
		return sigma*Math.pow(a.getLength()/w.distanceToOtherNode(a.getReceiver()), kappa)
				/(1-Math.pow(a.getLength()/R, kappa));
	}

	//resets all of the variables
	private void reset() {
		S1.clear();
		I1.clear();
		U1.clear();
		all.clear();
		addElements();
		Collections.sort(S1);
		scale();
		R = S1.get(S1.size()-1).getLength() + 1;
	}

	@Override
	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub

	}

	//when a key is pressed this happens!
	@Override
	public void keyPressed(KeyEvent e) {
		algorithmELMISL();
		showSelectedLinks();
	}

	@Override
	public void keyReleased(KeyEvent e) {

	}

	public static void main(String[] args) {
		MaxLinksDemo2 test = new MaxLinksDemo2();
		test.setDefaultCloseOperation(EXIT_ON_CLOSE);
		test.setVisible(true);
	}
}
