/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package net.sf.mzmine.modules.rawdatamethods.peakpicking.diagnosticfilter;

import com.google.common.collect.Range;

public class ExclusionInformation {
    	
	private Double mz;
	private Range<Double> rtRange;
	
	public ExclusionInformation(Double mz, Range<Double> rtRange){
		this.mz = mz;
		this.rtRange = rtRange;
	}
	
	public Double getMZ() {
		return mz;
	}
	
	public Range<Double> getRTRange(){
		return rtRange;
	}
}
