package umontreal.ssj.stat;


/*
 * Class:        TallyHistogram
 * Description:  Histogram of a tally
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        January 2011
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
// package umontreal.ssj.stat;
//import umontreal.ssj.stat.*;
//import java.util.logging.Level;
import java.util.logging.Logger;
import umontreal.ssj.util.PrintfFormat;

/**
 * This class is similar to @ref TallyHistogram, except that it does not maintain the min, max,
 * average, and variance of the observations. Only the counters for the histogram are maintained.
 * There are also no extra virtual bins on the left and on the right to count the observations that
 * fall outside [a,b].
 *
 * <div class="SSJ-bigskip"></div>
 */
public class HistogramOnly extends StatProbe implements Cloneable {
	protected int numObs;  // number of observations
	protected int numBins; // number of bins
	protected int[] count;   // count[j] = num obs in bin[j], j=0,...,numbins-1.
	protected double m_h; // width of 1 bin
	protected double m_a; // left boundary of first bin
	protected double m_b; // right boundary of last bin
	private Logger log = Logger.getLogger("umontreal.ssj.stat");

	/**
	 * Constructs a `TallyHistogram` statistical probe. Divide the interval
	 * 
	 * @f$[a,b]@f$ into @f$s@f$ bins of equal width and initializes a counter to 0 for each bin.
	 *             Whenever an observation falls into a bin, the bin counter is increased by 1.
	 *             There are two extra bins (and counters) that count the number of
	 *             observations @f$x@f$ that fall outside the interval @f$[a,b]@f$: one for
	 *             those @f$x< a@f$, and the other for those @f$x > b@f$.
	 * @param a
	 *            left boundary of interval
	 * @param b
	 *            right boundary of interval
	 * @param s
	 *            number of bins
	 */
	public HistogramOnly(double a, double b, int s) {
		super();
		init(a, b, s);
	}

	/**
	 * Constructs a new `TallyHistogram` statistical probe with name `name`.
	 * 
	 * @param name
	 *            the name of the tally.
	 * @param a
	 *            left boundary of interval
	 * @param b
	 *            right boundary of interval
	 * @param s
	 *            number of bins
	 */
	public HistogramOnly(String name, double a, double b, int s) {
		super();
		this.name = name;
		init(a, b, s);
	}

	/**
	 * Initializes this object. Divide the interval @f$[a,b]@f$ into
	 * 
	 * @f$s@f$ bins of equal width and initializes all counters to 0.
	 * @param s
	 *            number of bins
	 * @param a
	 *            left boundary of interval
	 * @param b
	 *            right boundary of interval
	 */
	public void init(double a, double b, int s) {
		/*
		 * The counters co[1] to co[s] contains the number of observations falling in the interval
		 * [a, b]. co[0] is the number of observations < a, and co[s+1] is the number of
		 * observations > b.
		 */
		if (b <= a)
			throw new IllegalArgumentException("   b <= a");
		count = new int[s];
		numBins = s;
		m_h = (b - a) / s;
		m_a = a;
		m_b = b;
		for (int i = 0; i < s; i++)
			count[i] = 0;
		numObs = 0;
	}

	public void init() {
		for (int i = 0; i < numBins; i++)
			count[i] = 0;
		numObs = 0;
	}

	/**
	 * Fills this object from the first n observations in array obs.
	 */
	public void fillFromArray(double[] obs, int numObs) {
		for (int i = 0; i < numBins; i++)
			count[i] = 0;
		for (int i = 0; i < numObs; i++)
			add(obs[i]);
	}

	/**
	 * Fills this object from the entire array obs.
	 */
	public void fillFromArray(double[] obs) {
		fillFromArray(obs, obs.length);
	}

	/**
	 * Fills this object from the observations in a TallyStore object.
	 */
	public void fillFromTallyStore(TallyStore ts) {
		fillFromArray(ts.getArray(), ts.numberObs());
	}

	/**
	 * Gives a new observation @f$x@f$ to the statistical collectors. Increases by 1 the bin counter
	 * in which value @f$x@f$ falls. Values that fall outside the interval @f$[a,b]@f$ are added in
	 * extra bin counter bin[0] if @f$x < a@f$, and in bin[@f$s+1@f$] if @f$x > b@f$.
	 * 
	 * @param x
	 *            observation value
	 */
	public void add(double x) {
		numObs++;
		if ((x >= m_a) & (x <= m_b))
			++count[(int) ((x - m_a) / m_h)];
	}

	/**
	 * Returns the number of observations given to histogram since its last initialization.
	 * 
	 * @return the number of collected observations
	 */
	public int numberObs() {
		return numObs;
	}

	/**
	 * Remove empty bins in the tails (left and right), without changing the bin size.
	 */
	public void trimHistogram() {
		int cpL = 0; // number of empty bins from left
		int cpR = 0; // number of empty bins from right
		while (count[cpL] == 0)
			cpL++;
		while (count[numBins - cpR - 1] == 0)
			cpR++;
		if (cpL + cpR > 0) {
			numBins -= cpL + cpR;
			int[] coco = new int[numBins];
			System.arraycopy(count, cpL, coco, 0, numBins);
			count = coco;
			m_a = m_a + (cpL * m_h);
			m_b = m_b - (cpR * m_h);
		}
	}

	/**
	 * Merges this histogram with the other histogram, by adding the bin counts of the two
	 * histograms. Returns the result in a new histogram.
	 * 
	 * @param other
	 *            the histogram to add
	 * 
	 *            Returns the merged histogram.
	 */
	public HistogramOnly addHistograms(HistogramOnly other) {
		if (this.numBins != other.numBins)
			throw new IllegalArgumentException(
			        "trying to add two histograms with different numbers of bin");
		HistogramOnly image;
		try {
			image = (HistogramOnly) super.clone();
		} catch (CloneNotSupportedException e) {
			throw new IllegalStateException("Tally can't clone");
		}
		int[] coco = new int[numBins];
		System.arraycopy(count, 0, coco, 0, numBins);
		int countOther[] = other.getCounters();
		for (int i = 0; i < numBins; i++)
			coco[i] = count[i] + countOther[i];
		image.count = coco;
		image.m_h = m_h;
		image.m_a = m_a;
		image.m_b = m_b;
		image.numBins = numBins;
		image.numObs = numObs + other.numObs;
		return image;
	}

	/**
	 * Merges bins by groups of size $g$. If there are $m$ bins initially, the new number of bins
	 * will be $\lceil m/g\rceil$. The last bin may regroup less than $g$ original bins if $m$ is
	 * not a multiple of $g$
	 **/
	public HistogramOnly aggregateBins(int g) {
		HistogramOnly image;
		try {
			image = (HistogramOnly) super.clone();
		} catch (CloneNotSupportedException e) {
			throw new IllegalStateException("Tally can't clone");
		}
		int numBinsNew = (int) Math.ceil((double) numBins / (double) g);
		// int[] coco = new int[numBins];
		// System.arraycopy(count, 0, coco, 0, numBins);
		int[] countNew = new int[numBinsNew];
		int b = 0;
		for (int j = 0; j < numBinsNew - 1; j++) {
			for (int i = b; i < b + g; i++)
				countNew[j] += count[i];
			b = b + g;
		}
		while (b < numBins) {
			countNew[numBinsNew] += count[b];
			b++;
		}
		image.count = countNew;
		image.m_h = m_h * g;
		image.m_a = m_a;
		image.m_b = m_h * numBinsNew;
		image.numBins = numBinsNew;
		return image;
	}

	/**
	 * Returns the bin counters. Each counter contains the number of observations that fell in its
	 * corresponding bin. The counters bin[@f$i@f$], @f$i=1, 2, …, s@f$ contain the number of
	 * observations that fell in each subinterval of @f$[a,b]@f$. Values that fell outside the
	 * interval @f$[a,b]@f$ were added in extra bin counter bin[0] if @f$x < a@f$, and in
	 * bin[@f$s+1@f$] if @f$x > b@f$. There are thus @f$s+2@f$ counters.
	 * 
	 * @return the array of counters
	 */
	public int[] getCounters() {
		return count;
	}

	/**
	 * Returns the number of bins @f$s@f$ dividing the interval
	 * 
	 * @f$[a,b]@f$. Does not count the two extra bins for the values of
	 * @f$x<a@f$ or @f$x>b@f$.
	 * @return the number of bins
	 */
	public int getNumBins() {
		return numBins;
	}

	/**
	 * Returns the left boundary @f$a@f$ of interval @f$[a,b]@f$.
	 * 
	 * @return left boundary of interval
	 */
	public double getA() {
		return m_a;
	}

	/**
	 * Returns the right boundary @f$b@f$ of interval @f$[a,b]@f$.
	 * 
	 * @return right boundary of interval
	 */
	public double getB() {
		return m_b;
	}

	/**
	 * Clones this object and the array which stores the counters.
	 */
	public HistogramOnly clone() {
		HistogramOnly image;
		try {
			image = (HistogramOnly) super.clone();
		} catch (CloneNotSupportedException e) {
			throw new IllegalStateException("Tally can't clone");
		}
		int[] countNew = new int[numBins];
		System.arraycopy(count, 0, countNew, 0, numBins);
		image.count = countNew;
		image.m_h = m_h;
		image.m_a = m_a;
		image.m_b = m_b;
		image.numBins = numBins;
		image.numObs = numObs;
		return image;
	}

	/**
	 * Returns the bin counters as a `String`.
	 */
	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append("---------------------------------------" + PrintfFormat.NEWLINE);
		sb.append(name + PrintfFormat.NEWLINE);
		sb.append("Interval = [ " + m_a + ", " + m_b + " ]" + PrintfFormat.NEWLINE);
		sb.append("Number of bins = " + numBins + PrintfFormat.NEWLINE);
		sb.append(PrintfFormat.NEWLINE + "Counters = {" + PrintfFormat.NEWLINE);
		for (int i = 0; i < numBins; i++) {
			double a = m_a + (i - 1) * m_h;
			double b = m_a + i * m_h;
			sb.append("   (" + PrintfFormat.f(6, 3, a) + ", " + PrintfFormat.f(6, 3, b) + ")    "
			        + count[i] + PrintfFormat.NEWLINE);
		}
		sb.append("}" + PrintfFormat.NEWLINE);
		return sb.toString();
	}

	@Override
	public double average() {
		// TODO Auto-generated method stub
		throw new IllegalStateException("HistogramOnly.average() is not supported.");
	}

	@Override
	public String report() {
		// TODO Auto-generated method stub
		throw new IllegalStateException("HistogramOnly.report() is not supported.");
	}

	@Override
	public String shortReport() {
		// TODO Auto-generated method stub
		throw new IllegalStateException("HistogramOnly.shortReport() is not supported.");
	}

	@Override
	public String shortReportHeader() {
		// TODO Auto-generated method stub
		throw new IllegalStateException("HistogramOnly.shortReportHeader() is not supported.");
	}

}