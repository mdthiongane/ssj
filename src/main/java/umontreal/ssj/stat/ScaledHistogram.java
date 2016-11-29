package umontreal.ssj.stat;
import java.util.logging.Logger;


/* This is a histogram for which the integer counts (frequencies)
*  are replaced real numbers (in double).
*  These numbers represent the heights of the rectangles.
*  They can be chosen for example so that the integral of the histogram over 
*  [a,b] is equal to a specific value. 
*  If this value is taken as equal to 1, i.e., the sum of heights equal to 1/h = numBins / (b-a),
*  then the histogram can be seen as a density estimator.  
*/
public class ScaledHistogram {
	protected int numBins; // number of bins
	protected double m_h; // width of each bin
	protected double m_a; // left boundary of first bin
	protected double m_b; // right boundary of last bin
	protected double[] height; // rescaled counters: height[j] is the height of bin j.
	protected double integral;  // Total area under the histogram, = (b-a) x sum of heights.
	private Logger log = Logger.getLogger("umontreal.ssj.stat");

	private ScaledHistogram() {}

	/**
	 * Constructs a `ScaledHistogram` over the interval @f$[a,b]@f$, which is divided into 
	 * @f$numBins@f$ bins of equal width. Initializes to 0 the size of each bin.
	 * @param a
	 *            left boundary of interval
	 * @param b
	 *            right boundary of interval
	 * @param numBins
	 *            number of bins
	 */
	public ScaledHistogram(double a, double b, int numBins) {
		init(a, b, numBins);
	}

	/*
	 * Constructs a `ScaledHistogram` from hist by normalizing the bin counts so
	 * that the integral of the histogram is equal to integral.
	 */
	public ScaledHistogram (TallyHistogram hist, double integral) {
		init (hist, integral);
	}

	/*
	 * Constructs a `ScaledHistogram` from hist by normalizing the bin counts so
	 * that the integral of the histogram is equal to integral.
	 */
	public ScaledHistogram (HistogramOnly hist, double integral) {
		init (hist, integral);
	}
	
	public void init (double a, double b, int numBins) {
		if (b <= a)
			throw new IllegalArgumentException("   b <= a");
		this.numBins = numBins;
		m_h = (b - a) / numBins;
		m_a = a;
		m_b = b;
		height = new double[numBins];
		for (int i = 0; i <= numBins + 1; i++)
			height[i] = 0;
		integral = 0.0;
	}

	public void init (TallyHistogram hist, double integral) {
		m_a = hist.getA();
		m_b = hist.getB();
		m_h = hist.getH();
		numBins = hist.numBins;
		this.integral = integral;
		height = new double[numBins];
		int count[] = hist.getCounters();
		double scaleFactor = integral / (hist.numberObs() * m_h);
		for (int i = 0; i < numBins; i++) 
			height[i] = count[i+1] * scaleFactor;  // With extra counters for values outside.
		    // height[i] = count[i] * scaleFactor;  // Without the extra counters.
	}

	public void init (HistogramOnly hist, double integral) {
		m_a = hist.getA();
		m_b = hist.getB();
		m_h = hist.m_h;
		numBins = hist.numBins;
		this.integral = integral;
		height = new double[numBins];
		int count[] = hist.getCounters();
		double scaleFactor = integral / (hist.numberObs() * m_h);
		for (int i = 0; i < numBins; i++)
			height[i] = count[i] * scaleFactor;
	}
	
	
	/*
	 * Rescales the histogram by renormalizing the heights so its integral has the specified value.
	 */
	public void rescale (double integral) {
		double scaleFactor = integral / this.integral;
		for (int i = 0; i < numBins; i++)
			height[i] *= scaleFactor;
		this.integral = integral;
	}

	/*
	 * 
	 * Returns an ASH-transformed version of this scaled histogram. The
	 * ASH-transformed histogram has the same bin size as the original. The new
	 * frequency (height) in any given bin is the weighted average of the frequencies in
	 * the neighboring bins, with weights $(r-d)/r^2$ given to bins that are at
	 * distance $d$ from the target bin, for all $d < r$. See \cite{tSCO85a}.
	 * 
	 */
	public ScaledHistogram averageShiftedHistogram(int r) {
		ScaledHistogram image = clone();
		double[] heightNew = image.getHeights();
		double rscale = 1.0 / (r * r);   // Rescaling to be made for each bin.
		double sum = 0.0;
		for (int k = 0; k < numBins; k++) {
			heightNew[k] = r * height[k]; //comment
			for (int ell = 1; ell < r; ell++) {
				if (k-ell >= 0) heightNew[k] += (r-ell) * height[k-ell];
				if (k+ell < numBins) heightNew[k] += (r-ell) * height[k+ell];	
			}
		
			heightNew[k] *= rscale;
			sum += heightNew[k];
   	    }
		image.height = heightNew;
		image.integral = sum*m_h;
		return image;
	}

/*	public ScaledHistogram averageShiftedHistogram1 (int r) {
		// if (numBins % r != 0)
			// throw new IllegalArgumentException("r must divide the number of bins.");
		ScaledHistogram image = clone();
		double[] heightNew = image.getHeights();
		double rscale = 1.0 / (r * r);   // Rescaling to be made for each bin.
		double sum = 0.0;
        double S1[] = new double[numBins];
        double S2[] = new double[numBins];
        S1[0] = height[0];
        S2[0] = 0.0;
        heightNew[0] = r * S1[0];
		for (int ell = 1; ell < Math.min(r, numBins); ell++) {
			S2[0] += height[ell];
			heightNew[0] += (r-ell) * height[ell];
		}
       	for (int k = 1; k < numBins; k++) {
      		for (int ell = 1; ell < r; ell++) {
      			S1[k] = S1[k-1] + height[k];
      			if (k >= r) S1[k] -= height[k-r];
      			S2[k] = S2[k-1] - height[k];
      			if (k+r-1 < numBins) S2[k] += height[k+r-1];	
      			heightNew[k] = heightNew[k-1] + S2[k-1] - S1[k-1];
			}
			heightNew[k] *= rscale;
			sum += heightNew[k];
   	    }
		image.height = heightNew;
		image.integral = sum*m_h;
		return image;
	}*/
	

	public ScaledHistogram averageShiftedHistogram1 (int r) {
		ScaledHistogram image = clone();
		double[] heightNew = image.getHeights();
		double rscale = 1.0 / (r * r);   // Rescaling to be made for each bin.
		double sum = 0.0;
        double S1[] = new double[numBins];
        double S2[] = new double[numBins];
        S1[0] = height[0];
        S2[0] = 0.0;
        heightNew[0] = r * S1[0];
		for (int ell = 1; ell <= Math.min(r, numBins); ell++){
			S2[0] += height[ell];	
			heightNew[0] += (r-ell) * height[ell];
		}
		
       	for (int k = 2; k <=numBins; k++) {
      			S1[k-1] = S1[k-2] + height[k-1];
      			if (k >= r) S1[k-1] -= height[k-r];
      			S2[k-1] = S2[k-2] - height[k-1];
      			if (k+r< numBins) S2[k-1] += height[k+r-1];	
      			heightNew[k-1] = heightNew[k-2] + S2[k-2] - S1[k-2];
		    }      	
    	for (int k = 0; k <numBins; k++) {
    		heightNew[k] *= rscale;
			sum += heightNew[k];	
    	 }
       	
		image.height = heightNew;
		image.integral = sum*m_h;
		return image;
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


	/*
	 * return the array count of the average shifted histogram
	 */
	public double[] getHeights() {
		return height;
	}

	public double getIntegral() {
		return integral;
	}

	/**
	 * Computes and returns the integrated square error (ISE) of a histogram 
	 * w.r.t. the U(0,1) distribution.  Assumes the histogram integrates to 1.
	 */
	public double ISEvsU01 () {
		double sum = 0.0;
		for (int j = 0; j < numBins; ++j)
		   sum += (height[j] - 1.0) * (height[j] - 1.0);
		return sum / numBins;
	}
	
	/**
	 * Computes and returns the ISE of a polygonal density w.r.t. the U(0,1)
	 * distribution.   Assumes the histogram integrates to 1.
	 */
	public double ISEvsU01polygonal () {
		double w0[] = new double[numBins];  // histogram shifted down to zero mean.
		for (int j = 0; j < numBins; ++j) {
			w0[j] = height[j] - 1.0;
		}
		double sum = 0.5 * (w0[1] * w0[1] + w0[numBins-1] * w0[numBins-1]);
		for (int j = 0; j < numBins-2; ++j)
		   sum += 0.333333333333333 * (w0[j+1] - w0[j]) * (w0[j+1] - w0[j]) 
		          + w0[j] * w0[j] + w0[j] * w0[j+1];
		return sum / (double)numBins;
	}
	
	/**
	 * Clones this object and the array which stores the counters.
	 */
	public ScaledHistogram clone() {
		ScaledHistogram image = new ScaledHistogram();
		image.numBins = numBins;
		image.m_h = m_h;
		image.m_a = m_a;
		image.m_b = m_b;
		image.height = new double[numBins]; 
		image.integral = integral;  
		for (int j = 1; j < numBins; ++j)
		   image.height[j] = height[j];
		return image;
	}
		

}
