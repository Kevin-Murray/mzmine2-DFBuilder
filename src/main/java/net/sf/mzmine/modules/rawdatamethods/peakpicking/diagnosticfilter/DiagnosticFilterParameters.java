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

import net.sf.mzmine.main.MZmineCore;
import net.sf.mzmine.parameters.Parameter;
import net.sf.mzmine.parameters.impl.SimpleParameterSet;
import net.sf.mzmine.parameters.parametertypes.DoubleParameter;
import net.sf.mzmine.parameters.parametertypes.OptionalParameter;
import net.sf.mzmine.parameters.parametertypes.filenames.FileNameParameter;
import net.sf.mzmine.parameters.parametertypes.ranges.MZRangeParameter;
import net.sf.mzmine.parameters.parametertypes.tolerances.RTToleranceParameter;
import net.sf.mzmine.parameters.parametertypes.selectors.RawDataFilesParameter;
import net.sf.mzmine.parameters.parametertypes.tolerances.MZToleranceParameter;
import net.sf.mzmine.parameters.parametertypes.selectors.ScanSelectionParameter;

public class DiagnosticFilterParameters extends SimpleParameterSet {

  public static final RawDataFilesParameter dataFiles = new RawDataFilesParameter();

  public static final ScanSelectionParameter scanSelection = new ScanSelectionParameter();
  
  public static final MZRangeParameter mzRange =
      new MZRangeParameter("Precursor m/z", "Range of precursor m/z values");
 
  public static final FileNameParameter diagnosticFile = new FileNameParameter("Diagnostic feature list file",
		  "CSV file containing diagnostic filter targets. See Help for more info.",
		  "csv");
  
  public static final OptionalParameter<FileNameParameter> exclusionFile =
      new OptionalParameter<>(new FileNameParameter("(Optional) Exclusion feature list file",
		  "Optional CSV file of mass-RT combinations to exclude. See Help for more info.",
		  "csv"));
  
  public static final MZToleranceParameter mzDifference = new MZToleranceParameter();

  public static final DoubleParameter basePeakPercent = new DoubleParameter(
      "Minimum ion intensity (% base peak)",
      "Minimum ion intesity for screening, scaled to base peak. Will choose non-zero maximum between %base peak and relative abundance.",
      MZmineCore.getConfiguration().getRTFormat(), 5.0, 0.0, 100.0);
    
  public static final DoubleParameter minIntensity = new DoubleParameter(
      "Minimum ion intensity (relative abundance)",
      "Minimum ion intesity for screening, relative abundance. Will choose non-zero maximum between %base peak and relative abundance.",
      MZmineCore.getConfiguration().getRTFormat(), 25000.0);
  
  public static final RTToleranceParameter rtTolerance = new RTToleranceParameter();

  public static final OptionalParameter<FileNameParameter> exportFile =
      new OptionalParameter<>(new FileNameParameter("(Optional) Export detected precursor list",
		  "Optional export of precursor hits of interest. See Help for more info.",
		  "csv"));

  public DiagnosticFilterParameters() {
    super(new Parameter[] {dataFiles, scanSelection, mzRange, mzDifference,
        diagnosticFile, basePeakPercent, minIntensity, rtTolerance,
        exclusionFile, exportFile});
  }

}
