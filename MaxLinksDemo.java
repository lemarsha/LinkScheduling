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


public class MaxLinksDemo extends JFrame implements KeyListener {
	public static final int WIDTH = 600;
	public static final int HEIGHT = 600;
	public static final int NUM_NODES = 5;

	private ArrayList<Link> S1, S2;
	private Set<Link> I1, I2;
	private Set<Node> U1, U2;
	private ArrayList<Node> all;
	private Display display1, display2;
	private double P, R;
	private double nu = 0.99;
	private double xi = 0.01;
	private double kappa = 4;
	private double sigma = 16;
	private double phi;
	private double rho_01, rho_02;

	public MaxLinksDemo() {
		super();
		S1 = new ArrayList<Link>();
		I1 = new HashSet<Link>();
		U1 = new HashSet<Node>();
		S2 = new ArrayList<Link>();
		I2 = new HashSet<Link>();
		U2 = new HashSet<Node>();
		all = new ArrayList<Node>();
		P = 1000;
		phi = calcPhi();
		System.out.println("phi = " + phi);
		rho_01 = 1 + Math.pow((sigma*(16*1.202 + 8*Math.pow(Math.PI, 4)/90 - 6)/phi), 1/kappa);
		rho_02 = rho_01 - 1;
		addElements();
		Collections.sort(S1);
		scale();
		for(Link l : S1) {
			S2.add(l.copy());
		}
		R = S1.get(S1.size()-1).getLength() + 1;
		setSize(WIDTH*2, HEIGHT + 23);
		display1 = new Display(all, S1, I1);
		display2 = new Display(all, S2, I2);
		setLayout(new GridLayout(1,2));
		add(display1);
		add(display2);
		addKeyListener(this);
	}

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
	}

	public void scale() {
		double shortest = S1.get(0).getLength();
		for(Link l : S1) {
			l.scaleLengths(shortest);
			//System.out.println(l.getLength());
		}
		for(Node n : all) {
			n.scale(shortest);
		}
		P /= shortest;
	}

	public void printS() {
		for(Link l : S1) {
			System.out.println(l);
		}
	}

	public void showSelectedLinks() {
		//I.add(S.get(S.size()-1));
		for(Link l : I1) {
			l.select();
		}
		for(Link l : I2) {
			l.select();
		}
		repaint();
	}

	public void pickLinksOld() {
		//double R = Math.pow((nu*P)/(sigma*xi), 1/kappa);
		//System.out.println(R);
		Link a;
		double rho;
		//System.out.println(S1.size());
		while(S1.size() > 0) {
			a = S1.get(0);
			S1.remove(0);
			I1.add(a);
			U1.add(a.getSender());
			rho = 1 + (rho_01 - 1)/Math.pow(1 - Math.pow(a.getLength()/R, kappa), 1/kappa);
			//System.out.println("rho = " + rho + ", length = " + a.getLength());
			ArrayList<Link> S_1 = new ArrayList<Link>();
			ArrayList<Link> S_2 = new ArrayList<Link>();
			for(Link l : S1) {
				if(a.getSender().distanceToOtherNode(l.getSender()) < rho * a.getLength()) {
					S_1.add(l);
				}
			}
			for(Link l : S_1) {
				S1.remove(l);
			}
			for(Link l : S1) {
				if(relativeInterference(U1, l) > 1-phi) S_2.add(l);
			}
			for(Link l : S_2) {
				S1.remove(l);
			}
			//System.out.println(S1.size());
		}
	}

	public void pickLinksNew() {
		Link a;
		double rho;
		while(S2.size() > 0) {
			a = S2.get(0);
			S2.remove(0);
			I2.add(a);
			U2.add(a.getSender());
			rho = rho_02/Math.pow(1 - Math.pow(a.getLength()/R, kappa), 1/kappa);
			ArrayList<Link> S_1 = new ArrayList<Link>();
			ArrayList<Link> S_2 = new ArrayList<Link>();
			for(Link l : S2) {
				if(a.getReceiver().distanceToOtherNode(l.getSender()) < rho * a.getLength()) {
					S_1.add(l);
				}
			}
			for(Link l : S_1) {
				S2.remove(l);
			}
			for(Link l : S2) {
				if(relativeInterference(U2, l) > 1-phi) S_2.add(l);
			}
			for(Link l : S_2) {
				S2.remove(l);
			}
		}
	}
	
	public void newAlgo() {
		int size = S2.size();
		List<List<Double>> mat = new ArrayList<List<Double>>();
		for(int row = 0; row < size; row++) {
			mat.add(new ArrayList<Double>());
			//mat.set(row, new ArrayList<Double>(size));
			for(int col = 0; col < size; col++) {
				if(col != row) {
					mat.get(row).add(relativeInterference(S2.get(row).getSender(), S2.get(col)));
				}
				else mat.get(row).add(0.0);
			}
		}
		/*for(List<Double> lt : mat) {
			for(Double d : lt) {
				System.out.print(d + " ");
			}
			System.out.print("\n");
		}*/
		
	}

	private double relativeInterference(Set<Node> W, Link a) {
		double totalInterference = 0;
		for(Node w : W) {
			totalInterference += relativeInterference(w, a);
		}
		return totalInterference;
	}
	
	private double relativeInterference(Node w, Link a) {
		return sigma*Math.pow(a.getLength()/w.distanceToOtherNode(a.getReceiver()), kappa)
				/(1-Math.pow(a.getLength()/R, kappa));
	}

	private void reset() {
		S1.clear();
		I1.clear();
		U1.clear();
		S2.clear();
		I2.clear();
		U2.clear();
		all.clear();
		addElements();
		Collections.sort(S1);
		scale();
		for(Link l : S1) {
			S2.add(l.copy());
		}
		R = S1.get(S1.size()-1).getLength() + 1;
	}

	@Override
	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void keyPressed(KeyEvent e) {
		pickLinksOld();
		//pickLinksNew();
		newAlgo();
		showSelectedLinks();
		if(e.getKeyChar() == 'g') {
			int iterations = 1000;
			double avg1 = 0, avg2 = 0;
			for(int i=0; i<iterations; i++) {
				reset();
				pickLinksOld();
				pickLinksNew();
				avg1 += I1.size();
				avg2 += I2.size();
			}
			showSelectedLinks();
			avg1 /= iterations;
			avg2 /= iterations;
			System.out.println("Average number of links chosen by method 1: " + avg1);
			System.out.println("Average number of links chosen by method 2: " + avg2);
		}
	}

	@Override
	public void keyReleased(KeyEvent e) {
		// TODO Auto-generated method stub

	}

	public static void main(String[] args) {
		MaxLinksDemo test = new MaxLinksDemo();
		//test.printS();
		test.setDefaultCloseOperation(EXIT_ON_CLOSE);
		test.setVisible(true);
	}
}
