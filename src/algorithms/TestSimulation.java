package algorithms;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Set;

import org.junit.Test;

public class TestSimulation {

	private int minlen=1, maxlen=30;
	private int WIDTH=Simulation.WIDTH, HEIGHT=Simulation.HEIGHT;
	
	//tests link creation for secondary links
	@Test
	public void testSecondaryLinkGeneration() {
		Simulation s = new Simulation(1, 500, minlen, maxlen);
		ArrayList<Node> nodes = s.generateNodes(500); //500 links
		ArrayList<Link> links = s.generateLinks(nodes, minlen, maxlen); // b/w 100 and 300 in size
		
		assertEquals(500,nodes.size());
		assertEquals(500,links.size());
		
		for (Link l: links) {
			Node n1 = l.getSender();
			Node n2 = l.getReceiver();
			//assert link length and distance between nodes is equal
			assertEquals(n1.distanceToOtherNode(n2), l.getLength(), .01);
			//assert link length is greater than min and less than max
			assertTrue(l.getLength()<=maxlen);
			assertTrue(l.getLength()>=minlen);
			//assert all nodes are in bounds
			assertTrue(n1.getX()>=0);
			assertTrue(n1.getX()<=WIDTH);
			assertTrue(n1.getY()>=0);
			assertTrue(n1.getY()<=HEIGHT);
			assertTrue(n2.getX()>=0);
			assertTrue(n2.getX()<=WIDTH);
			assertTrue(n2.getY()>=0);
			assertTrue(n2.getY()<=HEIGHT);
		}
	}
	
	//tests link creation for primary link
	@Test
	public void testPrimaryLinkGeneration() {
		//1 link
		Simulation s = new Simulation(1, 500, minlen, maxlen);
		Node sender = s.generatePrimaryNode();
		ArrayList<Link> links = s.generateLinks(sender, 1, minlen, maxlen);
		
		//assert only one link
		assertEquals(links.size(),1);
		//assert link length reasonable and equal to distance b/w nodes
		assertEquals(links.get(0).getSender().distanceToOtherNode(links.get(0).getReceiver()),
				links.get(0).getLength(), 0.01);
		assertTrue(links.get(0).getLength()>=minlen);
		assertTrue(links.get(0).getLength()<=maxlen);
		//assert all nodes in bounds
		assertTrue(links.get(0).getSender().getX()==WIDTH/2);
		assertTrue(links.get(0).getSender().getY()==HEIGHT/2);
		assertTrue(links.get(0).getReceiver().getX()>0);
		assertTrue(links.get(0).getReceiver().getY()>0);
		assertTrue(links.get(0).getReceiver().getX()<=WIDTH);
		assertTrue(links.get(0).getReceiver().getY()<=HEIGHT);
		
		//more than 1 link
		s = new Simulation(10, 500, minlen, maxlen);
		sender = s.generatePrimaryNode();
		links = s.generateLinks(sender, 10, minlen, maxlen);
		
		//assert 10 links
		assertEquals(links.size(),10);
		Node n = links.get(0).getSender();
		//assert sender in middle
		assertTrue(n.getX()==WIDTH/2);
		assertTrue(n.getY()==HEIGHT/2);
		
		for (Link l: links) {
			Node n2 = l.getReceiver();
			//assert link length reasonable and equal to distance b/w nodes
			assertEquals(n.distanceToOtherNode(n2), l.getLength(), 0.01);
			assertTrue(l.getLength()>=minlen);
			assertTrue(l.getLength()<=maxlen);
			//assert all senders are the same
			assertTrue(l.getSender().getX()==n.getX());
			assertTrue(l.getSender().getY()==n.getY());
			//assert all nodes in bounds
			assertTrue(n2.getX()>=0);
			assertTrue(n2.getX()<=WIDTH);
			assertTrue(n2.getY()>=0);
			assertTrue(n2.getY()<=HEIGHT);
		}
		
	}

	
	@Test
	public void testIndependent_C() {
		Simulation s = new Simulation(1, 100, minlen, maxlen, 0.5, 4, 16, 0.2);
		s.run();
		Set<Link> selected_links = s.getAlgorithmC();
		System.out.println("Algorithm C size: " + selected_links.size());
		for (Link l1: selected_links) {
			double totalInterference=0.0;
			//check that all nodes unique
			for (Link l2: selected_links) {
				if (l1.equals(l2)) continue;
				totalInterference+=s.relativeInterference(l2.getSender(), l1);
			}
			assertTrue(totalInterference<=1);
		}
	}
	
	@Test
	public void testIndependent_PLMISL() {
		Simulation s = new Simulation(1, 100, minlen, maxlen, 0.5, 4, 16, 0.2);
		s.run();
		Set<Link> selected_links = s.getAlgorithmPLMISL();
		System.out.println("Algorithm PLMISL size: " + selected_links.size());
		for (Link l1: selected_links) {
			double totalInterference=0.0;
			//check that all nodes unique
			for (Link l2: selected_links) {
				if (l1.equals(l2)) continue;
				totalInterference+=s.relativeInterference(l2.getSender(), l1);
			}
			assertTrue(totalInterference<=1);
		}
	}
	
	@Test
	public void testIndependent_Tolerance() {
		Simulation s = new Simulation(1, 100, minlen, maxlen, 0.5, 4, 16, 0.2);
		s.run();
		Set<Link> selected_links = s.getAlgorithmTolerance();
		System.out.println("Algorithm Tolerance size: " + selected_links.size());
		for (Link l1: selected_links) {
			double totalInterference=0.0;
			//check that all nodes unique
			for (Link l2: selected_links) {
				if (l1.equals(l2)) continue;
				totalInterference+=s.relativeInterference(l2.getSender(), l1);
			}
			assertTrue(totalInterference<=1);
		}
	}
	
	@Test
	public void testIndependent_Matrix() {
		Simulation s = new Simulation(1, 100, minlen, maxlen, 0.5, 4, 16, 0.2);
		s.run();
		Set<Link> selected_links = s.getAlgorithmMatrix();
		System.out.println("Algorithm Matrix size: " + selected_links.size());
		for (Link l1: selected_links) {
			double totalInterference = 0.0;
			for (Link l2: selected_links) {
				if (l1.equals(l2)) continue;
				totalInterference+=s.relativeInterference(l2.getSender(),l1);
			}
			assertTrue(totalInterference<=1);
		}
	}
	
}
