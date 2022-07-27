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
package net.sf.mzmine.modules.rawdatamethods.peakpicking.diagnosticfilter;

import java.util.logging.Level;
import java.util.logging.Logger;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.Ostermiller.util.CSVParser;
import com.google.common.collect.Range;

import net.sf.mzmine.datamodel.MZmineProject;
import net.sf.mzmine.datamodel.DataPoint;
import net.sf.mzmine.datamodel.Feature;
import net.sf.mzmine.datamodel.PeakList;
import net.sf.mzmine.datamodel.PeakListRow;
import net.sf.mzmine.datamodel.RawDataFile;
import net.sf.mzmine.datamodel.Scan;
import net.sf.mzmine.datamodel.impl.SimpleDataPoint;
import net.sf.mzmine.datamodel.impl.SimpleFeature;
import net.sf.mzmine.datamodel.impl.SimplePeakIdentity;
import net.sf.mzmine.datamodel.impl.SimplePeakList;
import net.sf.mzmine.datamodel.impl.SimplePeakListAppliedMethod;
import net.sf.mzmine.datamodel.impl.SimplePeakListRow;
import net.sf.mzmine.modules.peaklistmethods.qualityparameters.QualityParameters;
import net.sf.mzmine.parameters.ParameterSet;
import net.sf.mzmine.parameters.parametertypes.tolerances.MZTolerance;
import net.sf.mzmine.parameters.parametertypes.tolerances.RTTolerance;
import net.sf.mzmine.parameters.parametertypes.selectors.ScanSelection;
import net.sf.mzmine.taskcontrol.AbstractTask;
import net.sf.mzmine.taskcontrol.TaskPriority;
import net.sf.mzmine.taskcontrol.TaskStatus;
import net.sf.mzmine.util.scans.ScanUtils;

public class DiagnosticFilterTask extends AbstractTask {

    private final MZmineProject project;
    private Logger logger = Logger.getLogger(this.getClass().getName());

    private RawDataFile rawDataFile;
    private int totalScans, processedScans;
    private Scan[] scans;
    private PeakList targetPeakList;


    // User Parameters
    private ParameterSet parameters;
    private ScanSelection scanSelection;
    private File diagnosticFile;
    private Boolean useExclusion;
    private File exclusionFile;
    private Range<Double> mzRange;
    private MZTolerance mzDifference;
    private Boolean exportFile;
    private File fileName;
    private Double basePeakPercent;
    private Double minIntensity;
    private RTTolerance rtTolerance;
    private int finishedLines = 0;

    DiagnosticFilterTask(MZmineProject project, final RawDataFile rawDataFile, final ParameterSet parameters) {

        this.project = project;
        this.parameters = parameters;
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
        this.rtTolerance = parameters.getParameter(DiagnosticFilterParameters.rtTolerance).getValue();

        this.exportFile = parameters.getParameter(DiagnosticFilterParameters.exportFile).getValue();
        this.fileName = !exportFile ? null
                : parameters.getParameter(DiagnosticFilterParameters.exportFile).getEmbeddedParameter()
                        .getValue();
    }

    public void run() {

        setStatus(TaskStatus.PROCESSING);

        logger.info("Started diagnostic filtering on " + rawDataFile);

        scans = scanSelection.getMatchingScans(rawDataFile);
        totalScans = scans.length;

        // Create new target feature list
        targetPeakList = new SimplePeakList(rawDataFile.getName() + " targetChromatograms", rawDataFile);

        List<DiagnosticInformation> targetList = this.readDiagnostic();

        List<ExclusionInformation> exclusionList = !useExclusion ? null
                : this.readExclusion();

        // exportList that will contain output m/z values, RT, and scan number for
        // ID for export to a CSV
        List<String> exportList = new ArrayList<>();

        for (Scan scan : scans) {

            // Cancel?
            if (isCanceled()) {
                return;
            }

            // Only use MS2
            if (scan.getMSLevel() != 2) {
                processedScans++;
                continue;
            }

            // check parent m/z in mass range
            if (!mzRange.contains(scan.getPrecursorMZ())) {
                processedScans++;
                continue;
            }

            // check precursor-rt in exclusion list
            if (useExclusion) {

                boolean exclude = false;

                for (ExclusionInformation target : exclusionList) {
                    Double mz = target.getMZ();
                    Range<Double> rtRange = target.getRTRange();
                    if (mzDifference.checkWithinTolerance(mz, scan.getPrecursorMZ())) {
                        if (rtRange.contains(scan.getRetentionTime())) {
                            exclude = true;
                        }
                    }
                }
                if (exclude) {
                    processedScans++;
                    continue;
                }
            }
            
            // Get intensity threshold - basepeak vs relative
            double highestIntensity = 0;
            if (basePeakPercent == 0) {
                highestIntensity = minIntensity;
            } else if (minIntensity == 0) {
                highestIntensity = scan.getHighestDataPoint().getIntensity() * basePeakPercent;
            } else {
                highestIntensity = Double.max(scan.getHighestDataPoint().getIntensity() * basePeakPercent,
                        minIntensity);
            }
            
            DataPoint[] scanDataPoints = scan.getDataPointsOverIntensity(highestIntensity);

            // If no data points meet threshold, skip scan
            if (scanDataPoints.length == 0) {
                processedScans++;
                continue;
            }
            
            // Get Product Ions
            List<Double> fragmentIon = new ArrayList<>();
            for(DataPoint dataPoint : scanDataPoints){
                fragmentIon.add(dataPoint.getMZ());
            }

            // Target information
            HashMap<String, Boolean> targetMap = new HashMap<>();

            // Search Ions and Neutral Loss for targets
            for (DiagnosticInformation target : targetList) {
                
                String targetName = target.getName();
                double[] targetedMZ = target.getMZList();
                double[] targetedNF = target.getNFList();
                
                List<Boolean> foundMZ = new ArrayList<>();
                List<Boolean> foundNF = new ArrayList<>();
                
                // Check if fragment ion in scan
                // If no fragment ions to be searched, return true
                if (targetedMZ[0] != 0) {
                    for (double key : targetedMZ) {
                        List<Boolean> check = new ArrayList<>();
                        for (double MZ : fragmentIon) {
                            check.add(mzDifference.getToleranceRange(key).contains(MZ));
                        }
                        foundMZ.add(check.contains(true));
                    }
                } else {
                    foundMZ.add(true);
                }
                
                // Check if neutral loss in scan
                // If no neutral loss to be searched, return true
                if (targetedNF[0] != 0) {
                    for (double key : targetedNF) {
                        List<Boolean> check = new ArrayList<>();
                        double targetNL = scan.getPrecursorMZ() - key;
                        
                        for (double MZ : fragmentIon) {
                            check.add(mzDifference.getToleranceRange(targetNL).contains(MZ));
                        }
                        foundNF.add(check.contains(true));
                    }
                } else {
                    foundNF.add(true);
                }
                
                // If all fragment ions and neutral losses found, add target
                if (!foundMZ.contains(false) && !foundNF.contains(false)) {
                    targetMap.put(targetName, true);
                } else {
                    targetMap.put(targetName, false);
                }
            }

            // If target found, build chromatogram in RT range
            if (targetMap.containsValue(true)) {

                // Get target info
                Range<Double> mzRange = mzDifference.getToleranceRange(scan.getPrecursorMZ());
                Range<Double> rtRange = rtTolerance.getToleranceRange(scan.getRetentionTime());
                String id = "";

                // Format all diagnostic targets names detected
                for (Map.Entry<String, Boolean> entry : targetMap.entrySet()) {
                    if (entry.getValue()) {
                        id = id + "target=" + entry.getKey() + ';';
                    }
                }
                
                // Remove tailing semi-colon
                id = id.substring(0, id.length() - 1);
                
                // Initilaize new row in peaklist
                int rowID = targetPeakList.getNumberOfRows() + 1;
                PeakListRow newRow = new SimplePeakListRow(rowID);
                targetPeakList.addRow(newRow);
                
                // Build chromatogram for target precursor
                PeakListRow row = targetPeakList.findRowByID(rowID);
                row.addPeakIdentity(new SimplePeakIdentity(id), true);
                buildChromatogram(row, mzRange, rtRange, id);

                if (exportFile) {
                    // add precursor m/z, retention time, and scan number to output
                    // .csv file
                    String dataMZ = Double.toString(scan.getPrecursorMZ());
                    String dataRT = Double.toString(scan.getRetentionTime());

                    String temp = dataMZ + "," + dataRT + ",";

                    for (Map.Entry<String, Boolean> entry : targetMap.entrySet()) {
                        if (entry.getValue()) {
                            temp = temp + "target=" + entry.getKey() + ';';
                        }
                    }
                    temp = temp.substring(0, temp.length() - 1);

                    exportList.add(temp);
                }
            }

            processedScans++;
        }

        // Append processed feature list to the project
        project.addPeakList(targetPeakList);

        // Add quality parameters to peaks
        QualityParameters.calculateQualityParameters(targetPeakList);

        // Add task description to peakList
        targetPeakList.addDescriptionOfAppliedTask(
                new SimplePeakListAppliedMethod("DFBuilder Target chromatogram builder ", parameters));

        if (exportFile) {
            writeDiagnostic(exportList);
        }

        logger.log(Level.INFO, "Finished diagnostic fragmentation screnning on {0}", this.rawDataFile);
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
        if (totalScans == 0) {
            return 0;
        } else {
            return ((double) processedScans / totalScans);
        }
    }

    public String getTaskDescription() {
        return "Screening for fragment patterns in " + rawDataFile;
    }

    public static boolean isAllTrue(boolean[] array) {
        for (boolean b : array) {
            if (!b) {
                return false;
            }
        }
        return true;
    }

    @Override
    public TaskPriority getTaskPriority() {
        return TaskPriority.NORMAL;
    }

    public List<DiagnosticInformation> readDiagnostic() {

        FileReader dbFileReader = null;
        try {
            List<DiagnosticInformation> list = new ArrayList<DiagnosticInformation>();
            dbFileReader = new FileReader(diagnosticFile);

            String[][] diagnosticListValue = CSVParser.parse(dbFileReader, ',');

            for (; finishedLines < diagnosticListValue.length; finishedLines++) {
                try {
                    
                    String name = new String();
                    String mzString = new String();
                    String nlString = new String();
                    
                    name = diagnosticListValue[finishedLines][0].trim();

                    // Remove FEFF character from CSV
                    mzString = diagnosticListValue[finishedLines][1].replace("\uFEFF", "").trim();
                    
                    // If no neutral loss specifed, file no third column to parse
                    // TODO - better handling of this. 
                    try {
                        nlString = diagnosticListValue[finishedLines][2].replace("\uFEFF", "").trim();
                    } catch (Exception e) {
                        e.printStackTrace();
                    }

                    double[] mz = {0};
                    double[] nl = {0};

                    if (!mzString.isEmpty()) {
                        mz = Stream.of(mzString.split(";")).mapToDouble(Double::parseDouble).toArray();
                    }

                    if (!nlString.isEmpty()) {
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

    public List<ExclusionInformation> readExclusion() {

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

                    list.add(new ExclusionInformation(mz, rtRange));
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

    public void buildChromatogram(PeakListRow targetListRow, Range<Double> mzRange, Range<Double> rtRange, String ID) {

        int[] scanRange = rawDataFile.getScanNumbers(1, rtRange);

        List<ScanPoint> targetChromatogram = new ArrayList<>();
        for (int num : scanRange) {
            Scan scan = rawDataFile.getScan(num);
            double mz = (mzRange.lowerEndpoint() + mzRange.upperEndpoint()) / 2.0;
            double intensity = 0;

            DataPoint basePeak = ScanUtils.findBasePeak(scan, mzRange);

            if (basePeak != null) {
                mz = basePeak.getMZ();
                intensity = basePeak.getIntensity();
            }

            ScanPoint targetScan = new ScanPoint(scan.getScanNumber(), mz, scan.getRetentionTime(), intensity);

            targetChromatogram.add(targetScan);
        }

        double area = 0, height = 0, mz = 0, rt = 0;
        int scanNumbers[] = new int[targetChromatogram.size()];
        DataPoint finalDataPoint[] = new DataPoint[targetChromatogram.size()];
        Range<Double> finalRTRange = null, finalMZRange = null, finalIntensityRange = null;
        int representativeScan = 0;

        // Process all datapoints
        for (int i = 0; i < targetChromatogram.size(); i++) {

            ScanPoint dp = targetChromatogram.get(i);

            if (i == 0) {
                finalRTRange = Range.singleton(dp.getRT());
                finalMZRange = Range.singleton(dp.getMZ());
                finalIntensityRange = Range.singleton(dp.getIntensity());
            } else {
                assert finalRTRange != null && finalMZRange != null && finalIntensityRange != null;
                finalRTRange = finalRTRange.span(Range.singleton(dp.getRT()));
                finalMZRange = finalMZRange.span(Range.singleton(dp.getMZ()));
                finalIntensityRange = finalIntensityRange.span(Range.singleton(dp.getIntensity()));
            }

            scanNumbers[i] = targetChromatogram.get(i).getScanNumber();
            finalDataPoint[i] = new SimpleDataPoint(dp.getMZ(), dp.getIntensity());
            mz += targetChromatogram.get(i).getMZ();

            // Check height
            if (targetChromatogram.get(i).getIntensity() > height) {
                height = targetChromatogram.get(i).getIntensity();
                rt = targetChromatogram.get(i).getRT();
                representativeScan = targetChromatogram.get(i).getScanNumber();
            }

            // Skip last data point
            if (i == targetChromatogram.size() - 1) {
                break;
            }

            // X axis interval length
            double rtDifference
                    = targetChromatogram.get(i + 1).getRT() - targetChromatogram.get(i).getRT();

            // Convert the RT scale to seconds
            rtDifference *= 60d;

            // intensity at the beginning and end of the interval
            double intensityStart = targetChromatogram.get(i).getIntensity();
            double intensityEnd = targetChromatogram.get(i + 1).getIntensity();

            // calculate area of the interval
            area += (rtDifference * (intensityStart + intensityEnd) / 2);
        }

        // Calculate average m/z value
        mz /= targetChromatogram.size();

        // Find the best fragmentation scan, if available
        int fragmentScan = ScanUtils.findBestFragmentScan(rawDataFile, finalRTRange, finalMZRange);

        // Find all MS2 fragment scans, if available
        int[] allMS2fragmentScanNumbers
                = ScanUtils.findAllMS2FragmentScans(rawDataFile, finalRTRange, finalMZRange);

        SimpleFeature newPeak = new SimpleFeature(rawDataFile, mz, rt, height, area, scanNumbers,
                finalDataPoint, Feature.FeatureStatus.ESTIMATED, representativeScan, fragmentScan,
                allMS2fragmentScanNumbers, finalRTRange, finalMZRange, finalIntensityRange);

        // Fill the gap
        targetListRow.addPeak(rawDataFile, newPeak);
    }

    public void writeDiagnostic(List<String> export) {
        // Write output to csv file - for DFBuilder target chromatogram builder module.
        try {
            // Cancel?
            if (isCanceled()) {
                return;
            }

            String namePattern = "{}";
            File curFile = fileName;

            if (fileName.getPath().contains(namePattern)) {

                String cleanPlName = rawDataFile.getName().replaceAll("[^a-zA-Z0-9.-]", "_");
                // Substitute
                String newFilename = fileName.getPath().replaceAll(Pattern.quote(namePattern), cleanPlName);
                curFile = new File(newFilename);
            }

            FileWriter writer = new FileWriter(curFile, true);

            String collect = export.stream().collect(Collectors.joining("\n"));

            writer.write(collect);
            writer.write("\n");
            writer.close();

        } catch (IOException e) {
            System.out.print("Could not output to file");
            System.out.print(e.getStackTrace());

            setStatus(TaskStatus.FINISHED);
        }
    }
}
