import static org.junit.Assert.*;

import java.util.Collections;

import org.junit.Before;
import org.junit.Test;


public class TestMaxLinkCalculations {
	
	MaxLinksDemo3 test;

	@Before
	public void setUp() throws Exception {
		test = new MaxLinksDemo3();
		/*
		 * 	methods called so far: 
		 * 	setUp()
		 * 		generatePrimaryLink();
		 * 		generateNodes();
		 * 		generateLinks();	
		 * 		Collections.sort(links);
		 * 		scale();
		 * 		reset();
		 * 	setUpConstants()
		 */
	}
	
	//tests that the constant calculations are correct
	@Test
	public void testSetUpConstants() {
		assertEquals(6.1446, test.getRho0(),0.0001);
	}
	
	//tests that the first removal is correct
	@Test
	public void testFirstRemocal() {
		
	}
	
	//tests that the second removal is correct
	@Test
	public void testSecondRemoval() {
		
	}
	}

}
