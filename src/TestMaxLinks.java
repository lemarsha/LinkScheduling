import static org.junit.Assert.*;

import java.util.ArrayList;


import org.junit.Test;


public class TestMaxLinks {
	
	MaxLinksDemo3 test = new MaxLinksDemo3();
	ArrayList<Node> nodes = new ArrayList<Node>();
	ArrayList<Link> links = new ArrayList<Link>();

	
	//tests that the constant calculations are correct
	@Test
	public void testSetUpConstants() {
		assertEquals(6.1446, test.getRho0(),0.0001);
	}
	


}
