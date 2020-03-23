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

import net.sf.mzmine.main.MZmineCore;
import net.sf.mzmine.parameters.Parameter;
import net.sf.mzmine.parameters.impl.SimpleParameterSet;
import net.sf.mzmine.parameters.parametertypes.DoubleParameter;
import net.sf.mzmine.parameters.parametertypes.OptionalParameter;
import net.sf.mzmine.parameters.parametertypes.filenames.FileNameParameter;
import net.sf.mzmine.parameters.parametertypes.ranges.MZRangeParameter;
import net.sf.mzmine.parameters.parametertypes.selectors.RawDataFilesParameter;
import net.sf.mzmine.parameters.parametertypes.tolerances.MZToleranceParameter;
import net.sf.mzmine.parameters.parametertypes.selectors.ScanSelectionParameter;

public class DiagnosticFilterParameters extends SimpleParameterSet {

  public static final RawDataFilesParameter dataFiles = new RawDataFilesParameter();

  public static final ScanSelectionParameter scanSelection = new ScanSelectionParameter();
  
  public static final MZRangeParameter mzRange =
      new MZRangeParameter("Precursor m/z", "Range of precursor m/z values");
 
  public static final FileNameParameter diagnosticFile = new FileNameParameter("Diagnostic feature list input file",
		  "Name of the input diagnostic features containing the name, neutral loss shift, and product ions of feature",
		  "csv");
  
  public static final OptionalParameter<FileNameParameter> exclusionFile =
      new OptionalParameter<>(new FileNameParameter("Optional exclusion feature list input file",
		  "Spectral massess and retentions windows to exclude from export, containing the mz, start RT, and end RT",
		  "csv"));

  public static final MZToleranceParameter mzDifference = new MZToleranceParameter();

  public static final DoubleParameter basePeakPercent = new DoubleParameter(
      "Minimum diagnostic ion intensity (% base peak)",
      "Percent of scan base peak of which ms/ms product ions must be above to be included in analysis",
      MZmineCore.getConfiguration().getRTFormat(), 5.0);

  public static final FileNameParameter fileName = new FileNameParameter("Peaklist output file",
      "Name of the output CSV file containing m/z and RT of selected precursor ions. "
          + "If the file already exists, it will be overwritten.",
      "csv");

  public DiagnosticFilterParameters() {
    super(new Parameter[] {dataFiles, scanSelection, diagnosticFile, exclusionFile,
        mzRange, mzDifference, basePeakPercent, fileName});
  }

}
