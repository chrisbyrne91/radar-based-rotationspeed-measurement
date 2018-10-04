/*
 * 	
 * 
 * 	© Christopher Byrne   <c13712341@mydit.ie> 
 *	http://ahfr.dit.ie/node/540
 *
 *	Dublin Institute of Technology 
 *	School of Electronic & Electrical Engineering
 * 	Antenna and High Frequency Research Centre
 * 
 */

package main.java;

import krause.vna.data.calibrated.VNACalibratedSampleBlock;
import krause.vna.library.VNALibrary;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.*;

import org.knowm.xchart.QuickChart;
import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XChartPanel;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.Styler.LegendPosition;
import org.knowm.xchart.style.colors.XChartSeriesColors;
import org.knowm.xchart.style.markers.SeriesMarkers;

public class LibrarySampleRunner {
	
	static ArrayList<Complex> z;					
	static ArrayList<Complex> dz;					
	static ArrayList<Complex> dzDt;					

	static ArrayList<Double> timeElapsed;			
	static ArrayList<Double> dt;					
	static ArrayList<Double> s11Phase;				
	static ArrayList<Double> s11Magnitude;			
	static ArrayList<Double> angleDzDt;				
	static ArrayList<Double> bCoefficients;			
	static ArrayList<Double> signalStrength;		
	static ArrayList<Double> filteredSignalStrength;
	static ArrayList<Double> unwrappedPhaseAngle;	
	static ArrayList<Double> rotationSpeed;			
	static ArrayList<Double> filteredRotationSpeed;	
	static ArrayList<Double> acceleration;	
	static ArrayList<Double> jerk;
	static ArrayList<Double> dt2;
	
	public static void main(String[] args) {
		
		VNALibrary lib = null;
		VNACalibratedSampleBlock rc = null;
			
		try {
			
			lib = new VNALibrary();		
			int frequency = 0;
			float sweepTime = 0;
			int windowSize = 0;
			int filterItterations = 0; 
			double singnalThreshold = 0.0;
		
			
			//lib.loadDriverByName("miniVNA-pro-extender", "COM6");								
			//lib.loadCalibrationFile("D:\\miniVNA Calibration Files\\1metreCopperCable.cal");	
			
			System.out.print("Please enter the MiniVNA's associated COM port: ");
			String comPort = userInput.next( );
			System.out.print("Please enter the classpath of the calibration file: ");
			String calibrationFileDirectory = userInput.next( );
			lib.loadDriverByName("miniVNA-pro-extender", comPort);				
			System.out.println("MiniVNA-pro-extender successfully connected to: " + comPort); 
			lib.loadCalibrationFile(calibrationFileDirectory);									
			System.out.println("Calibration file successfully installed from: " + calibrationFileDirectory);
			
			String userInputCondition = "Y";
			while(userInputCondition.equals("Y")) {
				
				System.out.print("\nPlease enter scan frequency in Hz: ");
				String frequencyString = userInput.next( );
				frequency = Integer.parseInt(frequencyString);
				System.out.print("Please enter sweep time in ms: ");
				String sweepTimeString = userInput.next( );
				sweepTime = Float.parseFloat(sweepTimeString);
				System.out.print("Please enter window size for digtial filter: ");
				String windowSizeString = userInput.next( );
				windowSize = Integer.parseInt(windowSizeString);
				System.out.print("Please enter number of filter iterations: ");
				String filterItterationsString = userInput.next( );
				filterItterations = Integer.parseInt(filterItterationsString);
				System.out.print("Please enter signal threshold: ");
				String filterItterationsString = userInput.next( );
				filterItterations = Double.parseDouble(singnalThreshold);
				
				
				System.out.print("\nYou have entered the following parameters:\n");
				System.out.println("Scan Frequency:\t\t" + frequency + "Hz");
				System.out.println("Sweep Time:\t\t" + sweepTime + "ms");
				System.out.println("FIR Window Size:\t" + windowSize);
				System.out.println("Filter Iterations:\t" + filterItterations);
				System.out.println("Signal Threshold:\t" + singnalThreshold);
				System.out.print("Would you like to change these parameteres? Y/N:\t");
				userInputCondition = userInput.next();
			}
			
			// Modify if required 
			//int frequency = 1218200000; 	
			//float sweepTime = 5000;	
			//int windowSize = 20;		
			//int filterItterations = 40;
			
			String results = "Time (s)" + "," + "Signal Strength" + "," + "Angular Velocity (RPM)" + "," 
					+ "Acceleration (Revs.s^-2)" + "," + "Jerk (Revs.s^-3\n";	
			
			double[] xData = new double[1];	// X-axis of Graphical Chart for Plotting Results
			double[] yData = new double[1];	// Y-axis of Graphical Chart for Plotting Results 
			  
			List <XYChart> charts = new ArrayList<XYChart>();
			final XYChart chartSignal = QuickChart.getChart("Strength of Rotation Detection", "Time [s]",
				"Magnitude", "Signal", xData, yData); 
			chartSignal.getStyler().setLegendPosition(LegendPosition.InsideNW);
			charts.add(chartSignal);
			
			final XYChart chartRPM = QuickChart.getChart("Angular Velocity", "Time [s]", "Speed [RPM]",
				"Angular Velocity", xData, yData); 
			chartRPM.getStyler().setLegendPosition(LegendPosition.InsideNW);
			charts.add(chartRPM);
			
			final XYChart chartAcceleration = QuickChart.getChart("Angular Acceleration", "Time [s]", 
				"Acceleration [Revs/min2]", "Angular Acceleration", xData, yData); 
			chartAcceleration.getStyler().setLegendPosition(LegendPosition.InsideNW);
			charts.add(chartAcceleration);
			
			final XYChart chartJerk = QuickChart.getChart("Angular Jerk", "Time [s]", "Jerk [Revs/min3]", 
				"Angular Jerk", xData, yData); 
			chartJerk.getStyler().setLegendPosition(LegendPosition.InsideNW);
			charts.add(chartJerk);
			
			SwingWrapper<XYChart> sw = new SwingWrapper<XYChart>(charts); 
			sw.displayChartMatrix();
			
			ArrayList<Double> accumulativeTime = new ArrayList<Double>();
			ArrayList<Double> accumulativeSignalStrength = new ArrayList<Double>();
			ArrayList<Double> accumulativeRotationSpeed = new ArrayList<Double>();
			ArrayList<Double> accumulativeAcceleration = new ArrayList<Double>();
			ArrayList<Double> accumulativeJerk = new ArrayList<Double>();
			
			boolean subsequentLoop = false;
			double lastSampleTime = 0;		
			int globalLoopIndex = 0;		
			while(true){
				 
				initialiseVariables();								
			
				getSParameters(lib, theFrequency, theSweepTime);	
				
				double nf = computeIQNormalisationFactor(s11Phase, s11Magnitude);	
				z = computeComplexVectorArray(nf, s11Phase, s11Magnitude);			
				dt = computeChangeInVector(timeElapsed);
				dzDt = complexDifferenitation(timeElapsed, z);					
				angleDzDt = computeVectorPhaseAngle(dzDt);		
				
				signalStrength = computeVectorMagnitude(dzDt);
				bCoefficients = computeBCoefficients(theWindowSize);				
				filteredSignalStrength = filter(theWindowSize, theFilterIterations, signalStrength, bCoefficients);	
				double averageSignalStrength =  round(getArrayAverage(filteredSignalStrength),2);			
				accumulativeSignalStrength.add(averageSignalStrength);
				
				if (averageSignalStrength < singnalThreshold) {
					for (int i = 0; i < dt.size(); i++) {
						filteredRotationSpeed.add(0.0);
						acceleration.add(0.0);
						jerk.add(0.0);
					}
				}else {
				unwrappedPhaseAngle = unwrapPhaseAngleArray(angleDzDt);				
				rotationSpeed = calculateRotationSpeed(unwrappedPhaseAngle, dt);	
				filteredRotationSpeed = filter(theWindowSize, theFilterIterations, rotationSpeed, bCoefficients); 
				}
				double averageRotationSpeed = round(getArrayAverage(filteredRotationSpeed),2);	
				accumulativeRotationSpeed.add(averageRotationSpeed); 
				
				if (subsequentLoop) {							
					lastSampleTime = (accumulativeTime.get(accumulativeTime.size()-1));
					accumulativeTime.add( round( (timeElapsed.get(timeElapsed.size()-1)+lastSampleTime),3) );
					
					double deltaLoopTime = ( (accumulativeTime.get(globalLoopIndex)) 
						- (accumulativeTime.get(globalLoopIndex-1)) );
					double deltaRotationSpeed = ( (accumulativeRotationSpeed.get(globalLoopIndex)) - 
						(accumulativeRotationSpeed.get(globalLoopIndex-1)) );
					accumulativeAcceleration.add( round(deltaRotationSpeed/deltaLoopTime,2) );
					double deltaAcceleration = ( (accumulativeAcceleration.get(globalLoopIndex)) - 
						(accumulativeAcceleration.get(globalLoopIndex-1)) );
					accumulativeJerk.add( round(deltaAcceleration/deltaLoopTime,2) );
				} 
				else{
					accumulativeTime.add(round(timeElapsed.get(timeElapsed.size()-1),3));
					accumulativeAcceleration.add(0.0);
					accumulativeJerk.add(0.0);
				}
								
				plotData(accumulativeTime, accumulativeSignalStrength, charts.get(0),
					sw.getXChartPanel(0), "Signal");
				plotData(accumulativeTime, accumulativeRotationSpeed, charts.get(1), 
					sw.getXChartPanel(1), "Angular Velocity");
				plotData(accumulativeTime, accumulativeAcceleration, charts.get(2), 
					sw.getXChartPanel(2), "Angular Acceleration");
				plotData(accumulativeTime, accumulativeJerk, charts.get(3), 
					sw.getXChartPanel(3), "Angular Jerk");
				
				System.out.println("\nAt time: " + lastSampleTime + "s " + "with signal 
					amplitude of " + averageSignalStrength);
				System.out.println("The rotation speed is " + averageRotationSpeed + "RPM" + 
					" with acceleration of " + accumulativeAcceleration.get(globalLoopIndex) + 
						"Revs/min2 and jerk of " + accumulativeJerk.get(globalLoopIndex) + "Revs/min3");
				
				results += lastSampleTime + "," + averageSignalStrength + "," + averageRotationSpeed 
					+ "," + accumulativeAcceleration.get(globalLoopIndex) + "," + 
						accumulativeJerk.get(globalLoopIndex) + "\n";	
				outputResults(results);
				
				subsequentLoop = true;
				globalLoopIndex++;
			}
		} catch (Exception e) {
				
				System.out.println("failed with " + e.getMessage());
				
		} finally {
			if (lib != null) {
				lib.shutdown();
			}
		}
	}
	
	private static void initialiseVariables() {

		z = new ArrayList<Complex>();					
		dz = new ArrayList<Complex>();
		dzDt = new ArrayList<Complex>();
		
		s11Phase = new ArrayList<Double>();			
		s11Magnitude = new ArrayList<Double>();			
		angleDzDt = new ArrayList<Double>();			
		signalStrength = new ArrayList<Double>();		
		filteredSignalStrength = new ArrayList<Double>();	
		unwrappedPhaseAngle = new ArrayList<Double>();
		rotationSpeed = new ArrayList<Double>();		
		bCoefficients = new ArrayList<Double>();		
		filteredRotationSpeed = new ArrayList<Double>();
		timeElapsed = new ArrayList<Double>();			
		dt = new ArrayList<Double>();					
		dt2 = new ArrayList<Double>();	
		acceleration = new ArrayList<Double>();
		jerk = new ArrayList<Double>();
	}
	
	private static void getSParameters(VNALibrary lib, int frequency, double sweepTime) {
		VNACalibratedSampleBlock rc = null ;
		double currentTime = System.currentTimeMillis();  	
		double sweepTimeStart = currentTime;				
		double sweepTimeEnd = sweepTimeStart + sweepTime;
		
		while(System.currentTimeMillis() < sweepTimeEnd) {
			
			try { 
				rc = lib.scan(frequency, frequency+1088, 1, "REFL");				
			} catch (Exception e) {
				e.printStackTrace();
			}
			currentTime = System.currentTimeMillis();
			double sampleTime = (float) ((currentTime - sweepTimeStart)/1000f); 	
			
			// Store the Reflection Loss, Phase Loss and Sample time into arrays
			s11Magnitude.add(rc.getCalibratedSamples()[0].getReflectionLoss());		
			s11Phase.add(rc.getCalibratedSamples()[0].getReflectionPhase());	  	
			timeElapsed.add(sampleTime);									      
		}
	}

	private static double computeIQNormalisationFactor(ArrayList<Double> phase, 
		ArrayList<Double> magnitude) {
		double minRL = Collections.min(magnitude);		
		double maxRL = Collections.max(magnitude);		
		double minPL = Collections.min(phase); 			
		double maxPL = Collections.max(phase);			
		double normilisationFactor = (Math.abs(maxRL-minRL))/(maxPL-minPL);	
		return normilisationFactor;
	}
	
	private static ArrayList<Complex> computeComplexVectorArray(double normilisationFactor, 
		ArrayList<Double> ImaginaryPart, ArrayList<Double> realPart ) {
		ArrayList<Complex> z = new ArrayList<Complex>(); 
		for (int i = 0; i <= realPart.size()-1; i++) {
			z.add( new Complex( realPart.get(i), (ImaginaryPart.get(i))*normilisationFactor ) );			
		}
		return (z);
	}
	
	private static ArrayList<Double> computeChangeInVector(ArrayList<Double> vector) {
		ArrayList<Double> result = new ArrayList<Double>();
		for (int i = 0; i < vector.size()-1; i++) {
			result.add( (vector.get(i+1)) - (vector.get(i)) );
			}
		return (result);
	}
	
	@SuppressWarnings("unused")
	private static ArrayList<Double> differentiate(ArrayList<Double> X, ArrayList<Double> Y) {
		ArrayList<Double> dX = new ArrayList<Double>();
		ArrayList<Double> dY = new ArrayList<Double>();
		ArrayList<Double> result = new ArrayList<Double>();
		dX = computeChangeInVector(X);
		dY = computeChangeInVector(Y);
		for (int i = 0; i < Y.size()-1; i++) {
			result.add( (dX.get(i)) / (dY.get(i)) );
			}
		return (result);
	}
	private static ArrayList<Complex> complexDifferenitation(ArrayList<Double> 
		timeElapsed, ArrayList<Complex> z) {
		ArrayList<Complex>complexDt = new ArrayList<Complex>();
		ArrayList<Complex>dz = new ArrayList<Complex>();
		ArrayList<Complex>result = new ArrayList<Complex>();
		
		for (int i = 0; i < timeElapsed.size()-1; i++) {
			// complex cannot be divided by float - set imaginary component to zero
			complexDt.add(new Complex((timeElapsed.get(i+1)) - (timeElapsed.get(i)), 0));
			dz.add((z.get(i+1).minus(z.get(i))));		
			result.add((dz.get(i)).divides(complexDt.get(i)));	
		}
		return (result);
	}
	
	private static ArrayList<Double> computeVectorPhaseAngle(ArrayList<Complex> complexVector) {
		ArrayList<Double> phaseAngle = new ArrayList<Double>();
		for (int i = 0; i < complexVector.size()-1; i++) {
			phaseAngle.add(Math.atan2(((complexVector.get(i)).im()), ((complexVector.get(i)).re())));
		}
		return(phaseAngle);
	}
	
	private static ArrayList<Double> computeVectorMagnitude(ArrayList<Complex> complexVector) {
		ArrayList<Double> magnitude = new ArrayList<Double>();
		for (int i = 0; i < complexVector.size()-1; i++) {
			magnitude.add((complexVector.get(i)).abs());
		}
		return(magnitude);
	}
	
	private static ArrayList<Double> unwrapPhaseAngleArray(ArrayList<Double> phaseAngle) {
		
		ArrayList<Double> unwrappedPhase = phaseAngle;
		for (int i = 2; i < phaseAngle.size(); i++) {
			double difference = (phaseAngle.get(i)) - (phaseAngle.get(i-1));
			if (difference > (Math.PI)) {
				for (int j = i; j < phaseAngle.size(); j++) {
					unwrappedPhase.set(j, (unwrappedPhase.get(j) - (2*(Math.PI))));
				}
			}		
			if (difference < -(Math.PI)) {
				for (int j = i; j < phaseAngle.size(); j++) {
					unwrappedPhase.set(j, (unwrappedPhase.get(j) + (2*(Math.PI))));
				}
			}	
		}
		return (unwrappedPhase);
	}
	
	private static ArrayList<Double> calculateRotationSpeed(ArrayList<Double> 
		unwrappedPhaseAngle, ArrayList<Double> deltaTime) {
		ArrayList<Double> dPhase = new ArrayList<Double>();
		ArrayList<Double> dPhaseDt= new ArrayList<Double>();
		ArrayList<Double> result = new ArrayList<Double>();
		for (int i = 0; i < unwrappedPhaseAngle.size()-1; i++) {
			dPhase.add( (unwrappedPhaseAngle.get(i+1)) - (unwrappedPhaseAngle.get(i)) );
			dPhaseDt.add( (dPhase.get(i)) / (deltaTime.get(i)) );
			result.add( (60*(1/(4*Math.PI)))*(dPhaseDt.get(i)) );
		}
		return (result);
	}

	private static ArrayList<Double> computeBCoefficients(int windowSize) {
		ArrayList<Double> coefficients = new ArrayList<Double>();
		for (int i = 0; i < windowSize; i++) {
			Double temp = (1.0/windowSize);
			coefficients.add(temp);		
		}
		return(coefficients);
	}

	private static ArrayList<Double> filter(int windowSize, int filterItterations, 
		ArrayList<Double> unfilteredArray,ArrayList<Double> coefficients) {
		
		ArrayList<Double> filteredArray = new ArrayList<Double>();
		for (int k = 1; k <= filterItterations; k++) {
			if (k == 1){
				for (int n = windowSize; n < unfilteredArray.size(); n++) {
					int m = 0;
					double tempYSum = 0;
					for (int j = 0; j < windowSize; j++) {
						double tempY = (coefficients.get(j))*(unfilteredArray.get(n-m));
						tempYSum += tempY;
						m++;
					}
					filteredArray.add(tempYSum);
				}
			}
			else{
				int i = 1;
				for (int n = windowSize; n < filteredArray.size(); n++) {
					int m = 0;
					double tempYSum = 0;
					for (int j = 0; j < windowSize; j++) {
						double tempY = (coefficients.get(j))*(filteredArray.get(n-m));
						tempYSum += tempY;
						m++;
					}
					filteredArray.set(i, tempYSum);
					i++;
				}
			}
		}
		return(filteredArray);
	}
	
	private static void plotData(ArrayList<Double> xData, ArrayList<Double> yData,
			final XYChart theChart, final XChartPanel<XYChart> xChartPanel, String theSeries) {
		
		javax.swing.SwingUtilities.invokeLater(new Runnable() {
			  
			public void run() {
				
				if(yData.size()>xData.size()) {
					double[] xAxis = new double[xData.size()];
					double[] yAxis = new double[xData.size()];
				
					for (int i = 0; i < xData.size(); i++) {
						xAxis[i] = xData.get(i);
						yAxis[i] = yData.get(i);
						
						theChart.updateXYSeries(theSeries, xAxis, yAxis, null);
						
					}
				}else {
					double[] xAxis = new double[yData.size()];
					double[] yAxis = new double[yData.size()];
				
					for (int i = 0; i < yData.size(); i++) {
						xAxis[i] = xData.get(i);
						yAxis[i] = yData.get(i);
						
						theChart.updateXYSeries(theSeries, xAxis, yAxis, null);
						
						}
					}
				xChartPanel.repaint();
			    }
		});
	}

	private static void outputResults(String output) throws FileNotFoundException, 
		UnsupportedEncodingException {
		PrintWriter writer1 = new PrintWriter("results.csv", "UTF-8");
		writer1.println(output);
		writer1.close();
	}
	
	public static double round(double value, int places) {
	    if (places < 0) throw new IllegalArgumentException();
	
	    long factor = (long) Math.pow(10, places);
	    value = value * factor;
	    long tmp = Math.round(value);
	    return (double) tmp / factor;
	}
	
	public static double getArrayAverage(ArrayList<Double> arrayValues) {
	double sum = 0;
	for (int i = 0; i < arrayValues.size(); i++) {
		double temp = arrayValues.get(i);
		sum = temp + sum; 
	}
	return (sum/arrayValues.size());
	}
}
