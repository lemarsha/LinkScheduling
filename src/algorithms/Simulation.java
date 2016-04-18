package algorithms;


import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

public class Simulation {
	
	//dimensions of area
	public static final int WIDTH = 100;
	public static final int HEIGHT = 100;
	
	private ArrayList<Node> senders_s;
	private Node sender_p;
	private ArrayList<Link> links_s, links_p;
	
	//constants
	double phi, exp, max_rad, sigma, power;
	double noise = 0.000000001; //1.0*Math.pow(10, -9);
	double path_loss = 1;
	
	//results
	Set<Link> algorithmC, algorithmPLMISL, algorithmTolerance, algorithmMatrix;
	
	double rho_0;
	
	//constructor
	@SuppressWarnings("unchecked")
	public Simulation(int numPrimaryLinks, int numSecondaryLinks, int minLinkLength, int maxLinkLength) {
		//generate senders and links for secondary links
		senders_s = generateNodes(numSecondaryLinks);
		links_s = generateLinks(senders_s, minLinkLength, maxLinkLength);
		Collections.sort(links_s);
		//generate senders and links for primary links
		sender_p = generatePrimaryNode(); //in the middle
		links_p = generateLinks(sender_p, numPrimaryLinks, minLinkLength, maxLinkLength);
		Collections.sort(links_p);
	}
	
	//constructor
	@SuppressWarnings("unchecked")
	public Simulation(int numPrimaryLinks, int numSecondaryLinks, int minLinkLength, int maxLinkLength,
			double phi, double exp, double sigma, double power) {
		//generate senders and links for secondary links
		senders_s = generateNodes(numSecondaryLinks);
		links_s = generateLinks(senders_s, minLinkLength, maxLinkLength);
		Collections.sort(links_s);
		//generate senders and links for primary links
		sender_p = generatePrimaryNode(); //in the middle
		links_p = generateLinks(sender_p, numPrimaryLinks, minLinkLength, maxLinkLength);
		Collections.sort(links_p);
		
		this.phi=phi;
		this.exp = exp;
		this.max_rad = Math.pow((path_loss*power)/(sigma*noise),1/exp);
		this.sigma = sigma;
		this.power=power;
		rho_0=1+Math.pow((sigma*(16*1.202+8*Math.pow(Math.PI,4)/90-6)/phi),1/exp);
	}
	
	void run() {
		ArrayList<Link> links_s2 = new ArrayList<Link>(links_s);
		ArrayList<Link> links_p2 = new ArrayList<Link>(links_p);
		algorithmC = new HashSet<Link>(algorithmC(links_p2, links_s2));
		links_s2 = new ArrayList<Link>(links_s);
		links_p2 = new ArrayList<Link>(links_p);
		algorithmPLMISL = new HashSet<Link>(algorithmPLMISL(links_p2, links_s2));
		links_s2 = new ArrayList<Link>(links_s);
		links_p2 = new ArrayList<Link>(links_p);
		algorithmTolerance = new HashSet<Link>(algorithmTolerance(links_p2, links_s2));
		links_s2 = new ArrayList<Link>(links_s);
		links_p2 = new ArrayList<Link>(links_p);
		algorithmMatrix = new HashSet<Link>(algorithmMatrix(links_p2, links_s2));
	}
	
	//generate nodes
	public ArrayList<Node> generateNodes(int numNodes){
		ArrayList<Node> nodes = new ArrayList<Node>();
		Random r = new Random();
		//loop over number of  nodes and generate random senders
		for (int i = 0; i<numNodes; i++) {
			nodes.add(new Node(r.nextInt(WIDTH), r.nextInt(HEIGHT)));
		}
		return nodes;
	}
	
	//generate primary node
	public Node generatePrimaryNode() {
		return new Node(WIDTH/2, HEIGHT/2);
	}
	
	//generate secondary links
	public ArrayList<Link> generateLinks(ArrayList<Node> nodes, int minLinkLength, int maxLinkLength) {
		ArrayList<Link> links = new ArrayList<Link>();
		Random r = new Random();
		double x=0, y=0;
		for (Node n: nodes) {
			boolean possible=false;
			int length = r.nextInt(maxLinkLength-minLinkLength)+minLinkLength;
			while (!possible) {
				int angle = r.nextInt(360);
				x = n.getX()+length*Math.cos(angle*Math.PI/180.0);
				y = n.getY()-length*Math.sin(angle*Math.PI/180.0);
				if (x>0 && x<WIDTH && y>0 && y<HEIGHT) possible=true;
			}
			Node n2 = new Node(x,y);
			links.add(new Link(n, n2, length));
		}
		return links;
	}
	
	//generate primary links
	public ArrayList<Link> generateLinks(Node n, int numNodes, int minLinkLength, int maxLinkLength) {
		ArrayList<Link> links = new ArrayList<Link>();
		Random r = new Random();
		for (int i = 0; i<numNodes; i++) {
			int angle = r.nextInt(360);
			int length = r.nextInt(maxLinkLength-minLinkLength)+minLinkLength;
			double x = n.getX()+length*Math.cos(angle*Math.PI/180.0);
			double y = n.getY()-length*Math.sin(angle*Math.PI/180.0);
			Node n2 = new Node(x,y);
			Link l = new Link(n,n2,length);
			l.setPrimaryLink();
			links.add(l);
		}
		return links;
	}
	
	//calculates the relative interference of a set of nodes on another node
	public double relativeInterference(Set<Link> selected_links, Link a) {
		double totalInterference = 0;
		for(Link l: selected_links) {
			Node w = l.getSender();
			totalInterference+=relativeInterference(w,a);
		}
		return totalInterference;
	}

	//calculates the relative interference of a node on another node
	public double relativeInterference(Node w, Link a) {
		return sigma*Math.pow(a.getLength()/w.distanceToOtherNode(a.getReceiver()),exp)
			/(1-Math.pow(a.getLength()/max_rad,exp));
	}
	
	//removes links in first removal
	public ArrayList<Link> firstRemoval(ArrayList<Link> links_s, Link a, double rho) {
		ArrayList<Link> S_1 = new ArrayList<Link>();
		for (Link l: links_s) {
			if (a.getSender().distanceToOtherNode(l.getSender())<rho*a.getLength()) S_1.add(l);
		}
		for (Link l: S_1) {
			links_s.remove(l);
		}
		return links_s;
	}

	//removes links in second removal
	public ArrayList<Link> secondRemoval(ArrayList<Link> links_s, Set<Link> selected_links, Link primaryLink,
			Link a, boolean use_primary) {
		
		ArrayList<Link> S_2 = new ArrayList<Link>();
		if (a.getLength()<primaryLink.getLength() && use_primary) {
			for (Link l: links_s) {
				if((relativeInterference(selected_links,l)+relativeInterference(primaryLink.getSender(),l))>(1-phi)) S_2.add(l);
			}
		}else{
			for (Link l: links_s) {
				if(relativeInterference(selected_links,l)>(1-phi)) S_2.add(l);
			}
		}
		for (Link l: S_2) {
			links_s.remove(l);
		}
		return links_s;
	}
	
	
	//Circle algorithm
	public Set<Link> algorithmC(ArrayList<Link> links_p, ArrayList<Link> links_s) {
		Set<Link> selected_links = new HashSet<Link>();
		//add all primary links to set
		selected_links.addAll(links_p);
		//set up variables
		Link a;
		double rho;
		double rho_s = 1+(rho_0-1)/Math.pow(1-Math.pow(1/max_rad,exp),1/exp);
		for (Link primaryLink: links_p) {
			//calculate c
			double c = 24*sigma*Math.pow(primaryLink.getLength(),exp);
			c/=((exp-2)*(1-Math.pow(primaryLink.getLength()/max_rad,exp))*Math.pow(rho_s,exp));
			c = Math.ceil(Math.pow(c,1/(exp-2))+0.5);
			//removal of links in c*rho
			ArrayList<Link> S_1 = new ArrayList<Link>();
			for (Link l: links_s) {
				if (primaryLink.getSender().distanceToOtherNode(l.getSender())<c*rho_s) S_1.add(l);
			}
			for(Link l: S_1) {
				links_s.remove(l);
			}

		}
		//removals with interference too high
		ArrayList<Link> S_2 = new ArrayList<Link>();
		for (Link l: links_s) {
			if (relativeInterference(selected_links,l)>1-phi) S_2.add(l);
		}
		for (Link l: S_2) {
			links_s.remove(l);
		}
		//rest of algorithm
		while (links_s.size()>0) {
			a=links_s.get(0);
			links_s.remove(0);
			selected_links.add(a);
			rho = 1+(rho_0-1)/Math.pow(1-Math.pow(a.getLength()/max_rad,exp),1/exp);
			links_s = firstRemoval(links_s, a, rho);
			links_s = secondRemoval(links_s, selected_links, links_p.get(0), a, false);
		}
		
		return selected_links;
	}
	
	//calculates rho
	double calcrho(Link a) {
		return 1+(rho_0-1)/Math.pow(1-Math.pow(a.getLength()/max_rad,exp),1/exp);
	}
	
	//algorithm with primary link
	@SuppressWarnings("unchecked")
	public Set<Link> algorithmPLMISL(ArrayList<Link> links_p, ArrayList<Link> links_s) {
		Set<Link> selected_links = new HashSet<Link>();
		Set<Link> selected_primary_links = new HashSet<Link>();
		//set up variables
		Link a;
		double rho;
		boolean first=true;
		//sort primary links and get the smallest
		Collections.sort(links_p);
		Link primaryLink = links_p.get(0);
		//run algorithm
		while(links_s.size()>0) {
			a=links_s.get(0);
			links_s.remove(0);
			if (a.getLength()<primaryLink.getLength()) {
				if (relativeInterference(a.getSender(),primaryLink)+relativeInterference(selected_links,primaryLink)<=(1-phi)) {
					selected_links.add(a);
					rho = calcrho(a);
					links_s = firstRemoval(links_s, a,rho);
					links_s = secondRemoval(links_s, selected_links, primaryLink, a,true);
				}
			}
			else {
				if (first) {
					a=primaryLink;
					selected_links.add(a);
					selected_primary_links.add(a);
					rho = calcrho(a);
					links_s = firstRemoval(links_s, a,rho);
					links_s = secondRemoval(links_s, selected_links, primaryLink, a,true);
					//check if more primary links
					links_p.remove(0);
					if (links_p.size()==0) {
						first = false;
					} else {
						primaryLink = links_p.get(0);
					}
				}
				else {
					selected_links.add(a);
					rho = calcrho(a);
					links_s = firstRemoval(links_s, a,rho);
					links_s = secondRemoval(links_s, selected_links, primaryLink, a,true);

				}
			}
		}
		return selected_links;
	}
	
	public double powerOfLink(Node sender, Node reciever) {
		double multiplier = Math.pow(sender.distanceToOtherNode(reciever), exp);
		if (multiplier<1) multiplier =1;
		return power/multiplier;
	}
	
	public Set<Link> algorithmTolerance(ArrayList<Link> links_p, ArrayList<Link> links_s) {
		Set<Link> selected_links = new HashSet<Link>();
		//add primary links to selected links and initialize their tolerance
		for (Link l: links_p) {
			l.setTolerance(powerOfLink(l.getSender(), l.getReceiver())/sigma-noise);
			selected_links.add(l);
		}
		ArrayList<Link> removals = new ArrayList<Link>();
		//set the tolerances of all the secondary links
		for (Link l: links_s) {
			l.setTolerance(powerOfLink(l.getSender(), l.getReceiver())/sigma-noise);
			//decrement tolerance by interference of primary links
			for (Link p: links_p) {
				l.decrementTolerance(powerOfLink(p.getSender(), l.getReceiver()));
			}
			//remove link if tolerance is less than 0
			if (l.getTolerance()<0) {
				removals.add(l);
			}
		}
		//remove links whose tolerance is less than 0
		for (Link l: removals) links_s.remove(l);
		
		//pick links while there are still feasible links
		while (links_s.size()!=0) {
			removals.clear();
			//remove links that aren't feasible
			for (Link l: links_s) {
				for (Link s: selected_links) {
					double interferenceOfLink = powerOfLink(l.getSender(),s.getReceiver());
					if (interferenceOfLink>s.getTolerance()) {
						removals.add(l);
						continue;
					}
				}
			}
			//remove links that aren't feasible
			for (Link l: removals) links_s.remove(l);
			removals.clear();
			if (links_s.size()<=0) break;
			//pick link with highest max tolerance/max interference ratio
			Link current_link = links_s.get(0);
			double max_ratio = -1;
			for (Link l: links_s) {
				//find maximum interference
				double max_interference = -1;
				for (Link s: selected_links) {
					double interference = powerOfLink(l.getSender(),s.getReceiver());
					if (interference>max_interference) max_interference=interference;
				}
				//calculate the ratio and see if it is the maximum
				double ratio = l.getTolerance()/max_interference;
				if (ratio>max_ratio) {
					current_link = l;
					max_ratio=ratio;
				}
			}
			//remove link from list
			links_s.remove(current_link);
			//decrement tolerance of secondary links
			for (Link l: links_s) {
				l.decrementTolerance(powerOfLink(current_link.getSender(),l.getReceiver()));
				if (l.getTolerance()<0) removals.add(l);
			}
			for (Link l: removals) links_s.remove(l);
			//decrement tolerance of selected links by current one
			double interference=0;
			for (Link l: selected_links) {
				interference+=powerOfLink(l.getSender(),current_link.getReceiver());
				l.decrementTolerance(powerOfLink(current_link.getSender(),l.getReceiver()));
				if (l.getTolerance()<0) System.out.println("WTF");
			}
			if (interference>powerOfLink(current_link.getSender(),current_link.getReceiver())) System.out.println("WTF2");
			selected_links.add(current_link);

		}
		
		return selected_links;
	}
	
	//finds the link that should be removed based on row+col
	public int findRemoval_rc(double[][] matrix, int psize) {
		double max_interference = 0;
		int link=psize;
		for (int i = psize; i<matrix.length-1; i++) {
			double interference = matrix[i][matrix.length-1]+matrix[matrix.length-1][i];
			if (interference>max_interference) {
				max_interference=interference;
				link=i;
			}
		}
		return link;
	}
	
	//finds the link that should be removed based on row
	//interference that the link has on other links
	public int findRemoval_r(double[][] matrix, int psize) {
		double max_interference = 0;
		int link=psize;
		for (int i = psize; i<matrix.length-1; i++) {
			double interference = matrix[i][matrix.length-1];
			if (interference>max_interference) {
				max_interference=interference;
				link=i;
			}
		}
		return link;
	}
	
	//finds the link that should be removed based on col
	//interference of other links on the link
	public int findRemoval_c(double[][] matrix, int psize) {
		double max_interference = 0;
		int link=0;
		for (int i = psize; i<matrix.length-1; i++) {
			double interference = matrix[matrix.length-1][i];
			if (interference>max_interference) {
				max_interference=interference;
				link=i;
			}
		}
		return link;
	}
	
	public boolean independent(double[][] matrix) {
		for (int i = 0; i<matrix.length-1; i++) {
			if (matrix[matrix.length-1][i]>1) return false;
		}
		return true;
	}
	
	public Set<Link> algorithmMatrix(ArrayList<Link> links_p, ArrayList<Link> links_s) {
		Set<Link> selected_links = new HashSet<Link>();
		ArrayList<Link> links = new ArrayList<Link>();
		links.addAll(links_p);
		links.addAll(links_s);
		selected_links.addAll(links);
		int psize = links_p.size(), ssize = links_s.size();
		int size = psize+ssize+1;
		double[][] matrix = new double[size][size];
		//initialize the matrix
		for (int i = 0; i<size-1; i++) {
			for (int j = 0; j<size-1; j++) {
				if (i==j) {
					matrix[i][j]=0;
				}else if (i<psize && j<psize) {
					matrix[i][j]=0;
				}else {
					matrix[i][j] = relativeInterference(links.get(i).getSender(), links.get(j));
					matrix[i][size-1]+=matrix[i][j];
					matrix[size-1][j]+=matrix[i][j];
				}
			}
		}
		matrix[size-1][size-1]=0;
		//while the links aren't independent
		
		while (!independent(matrix)) {
			//find a link to remove
			int link_id = findRemoval_rc(matrix,psize);
			selected_links.remove(links.get(link_id));
			//subtract interference of that link from the other links
			for (int i = 0; i<size-1; i++) {
				//get the interference of the link on the link i
				double interference = matrix[link_id][i];
				//subtract the interference that link had on i from the total interference on i
				matrix[size-1][i]-=interference;
				//set that value = 0
				matrix[link_id][i] = 0;
			}
			matrix[size-1][link_id] = 0;
			matrix[link_id][size-1]=0;
		}
		return selected_links;
	}
	
	public ArrayList<Link> getSecondaryLinks() {
		return links_s;
	}
	
	public ArrayList<Link> getPrimaryLinks() {
		return links_p;
	}
	
	public Set<Link> getAlgorithmC() {
		return algorithmC;
	}

	public Set<Link> getAlgorithmPLMISL() {
		return algorithmPLMISL;
	}

	public Set<Link> getAlgorithmTolerance() {
		return algorithmTolerance;
	}

	public Set<Link> getAlgorithmMatrix() {
		return algorithmMatrix;
	}

	public static void main(String[] args) {
		//FIX R EQUATION
		Simulation s = new Simulation(10, 100, 1, 30, 0.5, 4, 16, 0.2);
		s.run();
		System.out.println("Algorithm C size: " + s.getAlgorithmC().size());
		System.out.println("Algorithm PLMISL size: " + s.getAlgorithmPLMISL().size());
		System.out.println("Algorithm Tolerance size: " + s.getAlgorithmTolerance().size());
		System.out.println("Algorithm Matrix size: " + s.getAlgorithmMatrix().size());
	}

}
