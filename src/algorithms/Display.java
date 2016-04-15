package algorithms;

import java.awt.GridLayout;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;

import javax.swing.JFrame;

public class Display extends JFrame implements KeyListener{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	Board board1;
	Simulation s;
	
	public Display(int numPrimaryLinks, int numSecondaryLinks, int minLinkLength, int maxLinkLength,
			double phi, double exp, double sigma, double power) {
		super();
		s = new Simulation(numPrimaryLinks, numSecondaryLinks, minLinkLength, maxLinkLength, phi, exp, sigma, power);
		s.run();
		setSize(Simulation.WIDTH,Simulation.HEIGHT);
		board1 = new Board(s.getSecondaryLinks(), s.getPrimaryLinks());	
		setLayout(new GridLayout(1,1));
		add(board1);
		addKeyListener(this);
	}
	
	
	@Override
	public void keyPressed(KeyEvent e) {
		if (e.getKeyChar()=='a') {
			System.out.println("AlgorithmC size: "+s.getAlgorithmC().size());
			board1.setSelectedLinks(s.getAlgorithmC());
		}else if (e.getKeyChar()=='s') {
			System.out.println("AlgorithmPLMISL size: "+s.getAlgorithmPLMISL().size());
			board1.setSelectedLinks(s.getAlgorithmPLMISL());
		}else if (e.getKeyChar()=='d'){
			System.out.println("AlgorithmTolerance size: "+s.getAlgorithmTolerance().size());
			board1.setSelectedLinks(s.getAlgorithmTolerance());
		}else if (e.getKeyChar()=='f'){
			System.out.println("AlgorithmMatrix size: "+s.getAlgorithmMatrix().size());
			board1.setSelectedLinks(s.getAlgorithmMatrix());
		}else if (e.getKeyChar()=='g') {
			board1.clearSelectedLinks();
		}
		repaint();
	}

	@Override
	public void keyReleased(KeyEvent e) {
	}

	@Override
	public void keyTyped(KeyEvent e) {	
	}
	
	public static void main(String[] args) {
		Display test = new Display(10, 100, 1, 30, 0.5, 4, 16, 0.2);
		test.setDefaultCloseOperation(EXIT_ON_CLOSE);
		test.setVisible(true);
	}

}
