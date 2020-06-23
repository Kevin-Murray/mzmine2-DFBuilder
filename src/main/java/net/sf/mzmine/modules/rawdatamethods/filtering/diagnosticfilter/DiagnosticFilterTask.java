/*
 * Copyright 2006-2020 The MZmine Development Team
 * 
 * This file is part of MZmine.
 * 
 * MZmine is free software; you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 * 
 * MZmine is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with MZmine; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
 * USA
 */

package net.sf.mzmine.modules.rawdatamethods.filtering.diagnosticfilter;

import java.util.logging.Level;
import java.util.logging.Logger;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.Ostermiller.util.CSVParser;
import com.google.common.collect.Range;

import net.sf.mzmine.datamodel.DataPoint;
import net.sf.mzmine.datamodel.RawDataFile;
import net.sf.mzmine.datamodel.Scan;
import net.sf.mzmine.parameters.ParameterSet;
import net.sf.mzmine.parameters.parametertypes.tolerances.MZTolerance;
import net.sf.mzmine.parameters.parametertypes.selectors.ScanSelection;
import net.sf.mzmine.taskcontrol.AbstractTask;
import net.sf.mzmine.taskcontrol.TaskPriority;
import net.sf.mzmine.taskcontrol.TaskStatus;

public class DiagnosticFilterTask extends AbstractTask {
	
  private Logger logger = Logger.getLogger(this.getClass().getName());

  private RawDataFile rawDataFile;
  private int totalScans, processedScans;
  private Scan[] scans;
  
  // User Parameters
  private ScanSelection scanSelection;
  private File diagnosticFile;
  private Boolean useExclusion;
  private File exclusionFile;
  private Range<Double> mzRange;
  private MZTolerance mzDifference;
  private File fileName;
  private Double basePeakPercent;
  private Double minIntensity;
  private int finishedLines = 0;

  DiagnosticFilterTask(final RawDataFile rawDataFile, final ParameterSet parameters) {

    this.rawDataFile = rawDataFile;
    
    this.scanSelection = parameters.getParameter(DiagnosticFilterParameters.scanSelection).getValue();
    this.diagnosticFile = parameters.getParameter(DiagnosticFilterParameters.diagnosticFile).getValue();
    this.useExclusion = parameters.getParameter(DiagnosticFilterParameters.exclusionFile).getValue();
    this.exclusionFile = !useExclusion ? null
            : parameters.getParameter(DiagnosticFilterParameters.exclusionFile).getEmbeddedParameter()
            .getValue();
    this.mzRange = parameters.getParameter(DiagnosticFilterParameters.mzRange).getValue();
    this.mzDifference = parameters.getParameter(DiagnosticFilterParameters.mzDifference).getValue();
    this.basePeakPercent = parameters.getParameter(DiagnosticFilterParameters.basePeakPercent).getValue() / 100;
    this.minIntensity = parameters.getParameter(DiagnosticFilterParameters.minIntensity).getValue();
    this.fileName = parameters.getParameter(DiagnosticFilterParameters.fileName).getValue();
  }

  public void run() {

    setStatus(TaskStatus.PROCESSING);
    
    logger.info("Started diagnostic filter on " + rawDataFile);
    
    scans = scanSelection.getMatchingScans(rawDataFile);
    totalScans = scans.length;
    
    List<DiagnosticInformation> targetList = this.readDiagnostic();
    
    
    List<ExclusionInformation> exclusionList = !useExclusion ? null 
            : this.readExclusion();
    
    // dataList that will contain output m/z values, RT, and scan number for
    // ID, for use in targeted feature detection
    List<String> dataList = new ArrayList<String>();

    for (Scan scan : scans) {

      // Cancel?
      if (isCanceled())
        return;
      
      if (scan.getMSLevel() != 2) {
        processedScans++;
    	continue;
      }
      
      // check parent m/z
      if (!mzRange.contains(scan.getPrecursorMZ())) {
        processedScans++;
        continue;
      }
      
      // check precursor-rt in exclusion list
      if(useExclusion){
          
          boolean exclude = false;
          
          for(ExclusionInformation target : exclusionList){
              Double mz = target.getMZ();
              Range<Double> rtRange = target.getRTRange();
              if(mzDifference.checkWithinTolerance(mz, scan.getPrecursorMZ())){
                  if(rtRange.contains(scan.getRetentionTime())){
                      exclude = true;
                  }
              }
          }
          
          if(exclude){
              processedScans++;
              continue;
          }
      }
      
      
      // get m/z and intensity values
      DataPoint scanDataPoints[] = scan.getDataPoints();

      // skip empty scans
      if (scan.getHighestDataPoint() == null) {
        processedScans++;
        continue;
      }

      // topPeaks will contain indexes to mzValues in scan above a
      // threshold defined as : 'scan
      // basePeak Intensity' * percent of base Peak to include
      List<Integer> topPeaksList = new ArrayList<Integer>();
      
      double highestIntensity = 0;
      if(basePeakPercent == 0){
          highestIntensity = minIntensity;
      } else if (minIntensity == 0) {
          highestIntensity = scan.getHighestDataPoint().getIntensity() * basePeakPercent;
      } else {
          highestIntensity = Double.min(scan.getHighestDataPoint().getIntensity() * basePeakPercent,
              minIntensity);
      }
            
      for (int i = 0; i < scanDataPoints.length; i++) {
        // Cancel?
        if (isCanceled())
          return;

        if ((scanDataPoints[i].getIntensity()) > highestIntensity) {
          // add the peaks
          topPeaksList.add(i);
        }
      }

      // Transfer topPeakList over to array
      Integer[] topPeaks = topPeaksList.toArray(new Integer[topPeaksList.size()]);

      HashMap<String, Boolean> targetMap = new HashMap<>();
      
      for (DiagnosticInformation target : targetList) {
    	  
    	  String targetName = target.getName();
    	  double[] targetedMZ = target.getMZList();
    	  double[] targetedNF = target.getNFList();
    	  
          // Default set to pass scan and not add to list
          targetMap.put(targetName, false);
    	  
          /**
           * Depending on filter conditions these if statements will filter based off of product m/z or
           * neutral loss or both within a scan. Pass becomes set to true if filter conditions are met
           * and scan is added to output file and visual plot
           */
          
          // Filter based off both m/z and neutral loss if both are not equal
          // to 0
          if (targetedMZ[0] != 0 && targetedNF[0] != 0) {
            boolean passA = false;
            boolean passB = false;
            boolean[] booleanValuesA = new boolean[targetedMZ.length];
            boolean[] booleanValuesB = new boolean[targetedNF.length];

            // scan through each m/z within scan m/z peaks
            for (int h = 0; h < topPeaks.length; h++) {
              // Cancel?
              if (isCanceled())
                return;

              int peakIndex = topPeaks[h];
              if (peakIndex < 0)
                break;
              double neutralLoss = scan.getPrecursorMZ() - scanDataPoints[peakIndex].getMZ();

              // scan for all m/z values if more than one, set pass to
              // true if all m/z values are found
              for (int j = 0; j < targetedMZ.length; j++) {
                // Cancel?
                if (isCanceled())
                  return;

                if (mzDifference.getToleranceRange(targetedMZ[j])
                    .contains(scanDataPoints[peakIndex].getMZ()) == true) {
                  booleanValuesA[j] = true;
                }
              }

              if (isAllTrue(booleanValuesA)) {
                passA = true;
              }

              // scan for all neutral loss values if more than one, set
              // pass to true if all neutral loss
              // values are found
              for (int j = 0; j < targetedNF.length; j++) {
                // Cancel?
                if (isCanceled())
                  return;

                if (mzDifference.getToleranceRange(targetedNF[j])
                    .contains(neutralLoss) == true) {
                  booleanValuesB[j] = true;
                }
              }
              if (isAllTrue(booleanValuesB)) {
                passB = true;
              }

            }
            // if both m/z and neutral loss pass, then total pass becomes
            // set to true, and scan is added
            if (passA && passB) {
            	targetMap.put(targetName, true);
            }

            // if only m/z requirements set, search for m/z and set to pass
            // if found in scan
          } else if (targetedMZ[0] != 0) {
            boolean[] booleanValues = new boolean[targetedMZ.length];
            for (int h = 0; h < topPeaks.length; h++) {
              int peakIndex = topPeaks[h];
              if (peakIndex < 0)
                break;
              for (int j = 0; j < targetedMZ.length; j++) {
                // Cancel?
                if (isCanceled())
                  return;

                if (mzDifference.getToleranceRange(targetedMZ[j])
                    .contains(scanDataPoints[peakIndex].getMZ()) == true) {
                  booleanValues[j] = true;
                }
              }
              if (isAllTrue(booleanValues)) {
            	  targetMap.put(targetName, true);
              }
            }

            // scan for n/f if both are not searched for and m/z is not
            // searched for
          } else if (targetedNF[0] != 0) {
            boolean[] booleanValues = new boolean[targetedMZ.length];
            for (int h = 0; h < topPeaks.length; h++) {
              // Cancel?
              if (isCanceled())
                return;

              int peakIndex = topPeaks[h];
              if (peakIndex < 0)
                break;
              double neutralLoss = scan.getPrecursorMZ() - scanDataPoints[peakIndex].getMZ();
              for (int j = 0; j < targetedNF.length; j++) {
                // Cancel?
                if (isCanceled())
                  return;

                if (mzDifference.getToleranceRange(targetedNF[j])
                    .contains(neutralLoss) == true) {
                  booleanValues[j] = true;
                }
              }
              if (isAllTrue(booleanValues)) {
            	  targetMap.put(targetName, true);
              }

            }

            // If no requirements set, continue
          } else {
            continue;
          }
          

      }
      
      // If target found, output precursor mz, rt, scan, and target info
      if (targetMap.containsValue(true)) {

        // add precursor m/z, retention time, and scan number to output
        // .csv file
        String dataMZ = Double.toString(scan.getPrecursorMZ());
        String dataRT = Double.toString(scan.getRetentionTime());
        
        String temp = dataMZ + "," + dataRT + ",";
        
        for(Map.Entry<String, Boolean> entry : targetMap.entrySet()) {
        	if(entry.getValue()) {
        		temp = temp + "target=" + entry.getKey() + ';';
        	}
        }
        temp = temp.substring(0, temp.length() - 1);

        dataList.add(temp);
      }
     
      processedScans++;
    }

    // Write output to csv file - for targeted feature detection module.
    try {
      // Cancel?
      if (isCanceled())
        return;

      String namePattern = "{}";
      File curFile = fileName;

      if (fileName.getPath().contains(namePattern)) {

        String cleanPlName = rawDataFile.getName().replaceAll("[^a-zA-Z0-9.-]", "_");
        // Substitute
        String newFilename = fileName.getPath().replaceAll(Pattern.quote(namePattern), cleanPlName);
        curFile = new File(newFilename);
      }

      FileWriter writer = new FileWriter(curFile, true);

      String collect = dataList.stream().collect(Collectors.joining("\n"));

      writer.write(collect);
      writer.write("\n");
      writer.close();

    } catch (IOException e) {
      System.out.print("Could not output to file");
      System.out.print(e.getStackTrace());

      setStatus(TaskStatus.FINISHED);
    }
    
    setStatus(TaskStatus.FINISHED);
  }

  /**
   * @see org.jfree.data.general.AbstractSeriesDataset#getSeriesKey(int)
   */
  public Comparable<Integer> getSeriesKey(int series) {
    return series;
  }

  public void cancel() {
    setStatus(TaskStatus.CANCELED);
  }

  public double getFinishedPercentage() {
    if (totalScans == 0)
      return 0;
    else
      return ((double) processedScans / totalScans);
  }

  public String getTaskDescription() {
    return "Updating fragment filter visualizer of " + rawDataFile;
  }

  public static boolean isAllTrue(boolean[] array) {
    for (boolean b : array)
      if (!b)
        return false;
    return true;
  }

  @Override
  public TaskPriority getTaskPriority() {
    return TaskPriority.NORMAL;
  }
  
  public List<DiagnosticInformation> readDiagnostic(){
	  
	  FileReader dbFileReader = null;
	  try {
		  List<DiagnosticInformation> list = new ArrayList<DiagnosticInformation>();
		  dbFileReader = new FileReader(diagnosticFile);
		  
		  String[][] diagnosticListValue = CSVParser.parse(dbFileReader, ',');
		  
		  for (; finishedLines < diagnosticListValue.length; finishedLines++) {
			  try {
				  String name = diagnosticListValue[finishedLines][0].trim();
				  
				  // Remove FEFF character from CSV
				  String mzString = diagnosticListValue[finishedLines][1].replace("\uFEFF", "").trim();
				  String nlString = diagnosticListValue[finishedLines][2].replace("\uFEFF", "").trim();
				  
				  double[] mz = {0};
				  double[] nl = {0};
				  
				  if(!mzString.isEmpty()) {
					  mz = Stream.of(mzString.split(";")).mapToDouble(Double::parseDouble).toArray();
				  }
				  
				  if(!nlString.isEmpty()) {
					  nl = Stream.of(nlString.split(";")).mapToDouble(Double::parseDouble).toArray();			
				  }
				  
				  list.add(new DiagnosticInformation(name, mz, nl));
			  } catch (Exception e) {
				  e.printStackTrace();
			  }
		  }
		  
		  dbFileReader.close();
		  return list;
	  } catch (Exception e) {
		  e.printStackTrace();
		  logger.log(Level.WARNING, "Could not read file " + diagnosticFile, e);
		  setStatus(TaskStatus.ERROR);
		  setErrorMessage(e.toString());
		  return null;
	  }
	  
  }
  
   public List<ExclusionInformation> readExclusion(){
	  
	  FileReader dbFileReader = null;
	  try {
		  List<ExclusionInformation> list = new ArrayList<ExclusionInformation>();
		  dbFileReader = new FileReader(exclusionFile);
		  
		  String[][] exclusionListValue = CSVParser.parse(dbFileReader, ',');
		  
		  for (; finishedLines < exclusionListValue.length; finishedLines++) {
			  try {
				  
				  // Remove FEFF character from CSV
				  String mzString = exclusionListValue[finishedLines][0].replace("\uFEFF", "").trim();
				  String rtStart = exclusionListValue[finishedLines][1].replace("\uFEFF", "").trim();
				  String rtEnd = exclusionListValue[finishedLines][2].replace("\uFEFF", "").trim();

				  Double mz = Double.parseDouble(mzString);
				  Range<Double> rtRange = Range.closed(Double.parseDouble(rtStart), Double.parseDouble(rtEnd));
                                                                    
				  list.add(new ExclusionInformation(mz,rtRange));
			  } catch (Exception e) {
				  e.printStackTrace();
			  }
		  }
		  
		  dbFileReader.close();
		  return list;
	  } catch (Exception e) {
		  e.printStackTrace();
		  logger.log(Level.WARNING, "Could not read file " + diagnosticFile, e);
		  setStatus(TaskStatus.ERROR);
		  setErrorMessage(e.toString());
		  return null;
	  }
	  
  }
   
   public void writeDiagnostic(PeakInformation peak){
       
   }
}
