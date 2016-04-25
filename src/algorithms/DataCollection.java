package algorithms;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;

public class DataCollection {
	
	/*
	 * Region size = 100,000 x 100,000
	 * Path loss exp = 4
	 * Threshold = 16
	 * Power = 20
	 * Noise = 10 E -16
	 * Links [1,000 to 10,000]
	 * Secondary links: 100, 150, 200, ... 500
	 * Primary links: 1, 2, 5, 10, 20, 30
	 * Average of 1000 instances
	 */

	
	//array list will hold selected link averages
	//inner array list holds values for number of secondary links
	//outer array list holds values for number of primary links
	ArrayList<ArrayList<Double>> selected_c, selected_t, selected_m, selected_p; 	
	//array list will hold selected time averages
	ArrayList<ArrayList<Long>> time_c, time_t, time_m, time_p;
	
	
	public DataCollection() {
		selected_c = new ArrayList<ArrayList<Double>>(); 
		selected_t = new ArrayList<ArrayList<Double>>();
		selected_m = new ArrayList<ArrayList<Double>>();
		selected_p = new ArrayList<ArrayList<Double>>();
		time_c = new ArrayList<ArrayList<Long>>();
		time_t = new ArrayList<ArrayList<Long>>();
		time_m = new ArrayList<ArrayList<Long>>();
		time_p = new ArrayList<ArrayList<Long>>();
	}
	
	public void CollectData_Base() {
		int minLinkLength = 1000, maxLinkLength = 10000;
		double phi=0.5, exp=4, sigma=16, power=20;
		
		//open the file to write to
		PrintWriter writer;
		try {
			writer = new PrintWriter("results.txt", "UTF-8");
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
		
		//print file thing to writer
		writer.println("Results of simulation");
		writer.println("Algorithm, NumPrimaryLinks, NumSecondaryLinks, NumLinksSelected, TimeToExecute(ns)");
		
		//1, 2, 5, 10, 20, 30 primary links
		for (int k = 0; k<6; k++) {
			int numPrimaryLinks = (k-2)*10;
			if (k==0) numPrimaryLinks =1;
			if (k==1) numPrimaryLinks = 2;
			if (k==2) numPrimaryLinks = 5;
			
			//loops over secondary links, inner array
			ArrayList<Double> selected_c_inner = new ArrayList<Double>();
			ArrayList<Double> selected_t_inner = new ArrayList<Double>();
			ArrayList<Double> selected_m_inner = new ArrayList<Double>();
			ArrayList<Double> selected_p_inner = new ArrayList<Double>();
			ArrayList<Long> time_c_inner = new ArrayList<Long>();
			ArrayList<Long> time_t_inner = new ArrayList<Long>();
			ArrayList<Long> time_m_inner = new ArrayList<Long>();
			ArrayList<Long> time_p_inner = new ArrayList<Long>();
			
			//100, 150, ... 300 secondary links
			for (int i = 0; i<9; i++) {
				int numSecondaryLinks = (i+1)*50+50;
				
				//average over the calculations
				int s_c=0;
				int s_t=0;
				int s_m=0;
				int s_p=0;
				long t_c=0;
				long t_t=0;
				long t_m=0;
				long t_p=0;
				
				//make sure stuff is correct
				System.out.println("Primary Links: "+numPrimaryLinks+", Secondary Links: "+numSecondaryLinks);
				
				//average over ten calculations
				int num_calcs = 1000;
				for (int j = 0; j<num_calcs; j++) {
					Simulation s = new Simulation(numPrimaryLinks, numSecondaryLinks, minLinkLength, maxLinkLength, phi, exp, sigma, power);
					
					// C algorithm
					long startTime = System.nanoTime();
					s.runC();
					long endTime = System.nanoTime();
					s_c+=s.getAlgorithmC().size();
					t_c+=(endTime-startTime);
					
					//Tolerance Algorithm
					startTime = System.nanoTime();
					s.runTolerance();
					endTime = System.nanoTime();
					s_t+=s.getAlgorithmTolerance().size();
					t_t+=(endTime-startTime);
					
					//Matrix Algorithm
					startTime = System.nanoTime();
					s.runMatrix();
					endTime = System.nanoTime();
					s_m+=s.getAlgorithmMatrix().size();
					t_m+=(endTime-startTime);
					
					//PLMISL algorithm
					startTime = System.nanoTime();
					s.runPLMISL();
					endTime = System.nanoTime();
					s_p+=s.getAlgorithmPLMISL().size();
					t_p+=(endTime-startTime);
					
					
				} //end loop for multiple trials
				
				
				//write the results to the files
				//System.out.println("C Algorithm, "+numPrimaryLinks+", "+numSecondaryLinks+", "+1.0*s_c/num_calcs+", "+t_c/num_calcs);
				writer.println("C Algorithm, "+numPrimaryLinks+", "+numSecondaryLinks+", "+1.0*s_c/num_calcs+", "+t_c/num_calcs);
				writer.println("T Algorithm, "+numPrimaryLinks+", "+numSecondaryLinks+", "+1.0*s_t/num_calcs+", "+t_t/num_calcs);
				writer.println("M Algorithm, "+numPrimaryLinks+", "+numSecondaryLinks+", "+1.0*s_m/num_calcs+", "+t_m/num_calcs);
				writer.println("P Algorithm, "+numPrimaryLinks+", "+numSecondaryLinks+", "+1.0*s_p/num_calcs+", "+t_p/num_calcs);
				
				//average the times and push them to array
				selected_c_inner.add(1.0*s_c/num_calcs);
				selected_t_inner.add(1.0*s_t/num_calcs);
				selected_m_inner.add(1.0*s_m/num_calcs);
				selected_p_inner.add(1.0*s_p/num_calcs);
				time_c_inner.add(t_c/num_calcs);
				time_t_inner.add(t_t/num_calcs);
				time_m_inner.add(t_m/num_calcs);
				time_p_inner.add(t_p/num_calcs);
				
				
				
			} //end loop for all the secondary links
			
			
			//push the arrays to the big arrays
			selected_c.add(selected_c_inner);
			selected_t.add(selected_t_inner);
			selected_m.add(selected_m_inner);
			selected_p.add(selected_p_inner);
			time_c.add(time_c_inner);
			time_t.add(time_t_inner);
			time_m.add(time_m_inner);
			time_p.add(time_p_inner);
			
			
		} //end loop for all the primary links
		
		//close the file
		writer.close();
	}
	
	public String printResults(int numPrimaryLinks, int numSecondaryLinks, double selected_links, String algo) {
		String results = "Algo, NumPrimary, NumSecondary, numSelected: ";
		results += (algo+", "+numPrimaryLinks+", "+numSecondaryLinks+", "+selected_links);
		return results;
	}
	
	public static void main(String[] args) {
		DataCollection dc = new DataCollection();
		dc.CollectData_Base();
		
	}
	

}

//public Simulation(int numPrimaryLinks, int numSecondaryLinks, int minLinkLength, int maxLinkLength, double phi, double exp, double sigma, double power) 
