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
//package umontreal.ssj.stat;
//import umontreal.ssj.stat.*;
//import java.util.logging.Level;
import java.util.logging.Logger;
import umontreal.ssj.util.PrintfFormat;

/**
 * This class is an extension of @ref Tally which gives a more detailed view of the observations
 * statistics. The individual observations are assumed to fall into different bins (boxes) of equal
 * width on an interval. The total number of observations falling into the bins are kept in an array
 * of counters. This is useful, for example, if one wish to build a histogram from the observations.
 * One must access the array of bin counters to compute quantities not supported by the methods
 * in @ref Tally.
 *
 * *Never add or remove observations directly* on the array of bin counters because this would put
 * the @ref Tally counters in an inconsistent state.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class TallyHistogram extends Tally {
	protected int numBins; // number of bins
	protected int[] count; // counter: num of values in bin[i].
	// count[0] and count[numBins] count the values outside interval.
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
	 * @param numBins
	 *            number of bins
	 */
	public TallyHistogram(double a, double b, int numBins) {
		super();
		init(a, b, numBins);
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
	 * @param numBins
	 *            number of bins
	 */
	public TallyHistogram(String name, double a, double b, int numBins) {
		super(name);
		init(a, b, numBins);
	}

	/**
	 * Initializes this object. Divide the interval @f$[a,b]@f$ into
	 * 
	 * @f$s@f$ bins of equal width and initializes all counters to 0.
	 * @param numBins
	 *            number of bins
	 * @param a
	 *            left boundary of interval
	 * @param b
	 *            right boundary of interval
	 */
	public void init(double a, double b, int numBins) {
		/*
		 * The counters count[1] to count[s] contains the number of observations falling in the
		 * interval [a, b]. count[0] is the number of observations < a, and count[numBins+1] is the
		 * number of observations > b.
		 */
		super.init();
		if (b <= a)
			throw new IllegalArgumentException("   b <= a");
		count = new int[numBins + 2];
		this.numBins = numBins;
		m_h = (b - a) / numBins;
		m_a = a;
		m_b = b;
		for (int i = 0; i <= numBins + 1; i++)
			count[i] = 0;
	}

	/**
	 * Fills this object from the first n observations in array obs.
	 */
	public void fillFromArray(double[] obs, int numObs) {
		super.init();
		for (int i = 0; i <= numBins + 1; i++)
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
		super.add(x);
		if (x < m_a)
			++count[0];
		else if (x > m_b)
			++count[1 + numBins];
		else {
			int i = 1 + (int) ((x - m_a) / m_h);
			++count[i];
		}
	}

	/**
	 * Remove empty bins in the tails (left and right), without changing the bin size.
	 */
	public TallyHistogram trimHistogram() {
		TallyHistogram image = (TallyHistogram) super.clone();

		int i = 1;
		int j = numBins; // number of bins in the initial histogram
		int cpL = 0; // number of empty bin from left initialized to zero
		int cpR = 0; // number of empty bin from right initialized to zero
		while (count[i] == 0) {
			i++;
			cpL++;
		}
		while (count[j] == 0) {
			j--;
			cpR++;
		}
		int[] coco = new int[2 + numBins - cpL - cpR];
		System.arraycopy(count, i, coco, 1, j - i + 1);
		System.arraycopy(count, 0, coco, 0, 1);
		System.arraycopy(count, count.length - 1, coco, coco.length - 1, 1);
		image.count = coco;
		image.m_h = m_h;
		image.m_a = m_a + (cpL * m_h);
		image.m_b = m_b - (cpR * m_h);
		image.numBins = numBins - cpL - cpR;
		return image;
	}

	/**
	 * Merges this histogram with the other histogram, by adding the bin counts of the two
	 * histograms.
	 * 
	 * @param other
	 *            the histogram to add
	 * 
	 *            Returns the merged histogram.
	 */
	public TallyHistogram addHistograms(TallyHistogram other) {
		if (this.numBins != other.numBins)
			throw new IllegalArgumentException("different number of bin in two histogram to merge");
		TallyHistogram image = (TallyHistogram) super.clone();
		int[] countNew = new int[2 + numBins];
		System.arraycopy(count, 0, countNew, 0, 2 + numBins);
		int coOther[] = other.getCounters();
		for (int i = 0; i < countNew.length; i++)
			countNew[i] = countNew[i] + coOther[i];
		image.count = countNew;
		image.m_h = m_h;
		image.m_a = m_a;
		image.m_b = m_b;
		image.numBins = numBins;
		return image;

	}

	/**
	 * Merges bins by groups of size $g$. If there are $m$ bins initially, the new number of bins
	 * will be $\lceil m/g\rceil$. The last bin may regroup less than $g$ original bins if $m$ is
	 * not a multiple of $g$
	 **/
	public TallyHistogram aggregateBins(int g) {
		TallyHistogram image = (TallyHistogram) super.clone();
		double numB = Math.ceil((double) numBins / (double) g);
		int numBinsNew = (int) numB;
		int[] countNew = new int[2 + numBinsNew];
		int b = 1;
		for (int j = 1; j < numBinsNew; j++) {
			for (int i = b; i < b + g; i++)
				countNew[j] = countNew[j] + count[i];
			b = b + g;
		}
		while (b < numBins - 1) {
			countNew[numBinsNew] = countNew[numBinsNew] + count[b];
			b++;
		}
		countNew[0] = count[0];
		countNew[numBinsNew - 1] = count[numBins - 1];
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
	 * Returns the width @f$h@f$ of rectangles.
	 * 
	 * @return the width of rectangles
	 */
	public double getH() {
		return m_h;
	}

	/**
	 * Clones this object and the array which stores the counters.
	 */
	public TallyHistogram clone() {
		TallyHistogram image = (TallyHistogram) super.clone();
		int[] coco = new int[2 + numBins];
		System.arraycopy(count, 0, coco, 0, 2 + numBins);
		image.count = coco;
		image.m_h = m_h;
		image.m_a = m_a;
		image.m_b = m_b;
		image.numBins = numBins;
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
		sb.append("Number of bins = " + numBins + " + 2" + PrintfFormat.NEWLINE);
		sb.append(PrintfFormat.NEWLINE + "Counters = {" + PrintfFormat.NEWLINE);
		sb.append("   (-inf, " + PrintfFormat.f(6, 3, m_a) + ")    " + count[0]
		        + PrintfFormat.NEWLINE);
		for (int i = 1; i <= numBins; i++) {
			double a = m_a + (i - 1) * m_h;
			double b = m_a + i * m_h;
			sb.append("   (" + PrintfFormat.f(6, 3, a) + ", " + PrintfFormat.f(6, 3, b) + ")    "
			        + count[i] + PrintfFormat.NEWLINE);
		}
		sb.append("   (" + PrintfFormat.f(6, 3, m_b) + ", inf)    " + count[numBins + 1]
		        + PrintfFormat.NEWLINE);
		sb.append("}" + PrintfFormat.NEWLINE);
		return sb.toString();
	}

}
