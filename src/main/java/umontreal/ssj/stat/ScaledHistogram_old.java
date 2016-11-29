package umontreal.ssj.stat;


/*
 * 
 * This is a histogram for which the counts (frequencies)
*  are real numbers (in double) instead of integers.
*The frequencies can add up to one, for example, so the 
*histogram can be seen as a density estimator.
*/
public class ScaledHistogram_old extends TallyHistogram_old{
	  private double[] co2;       // counter: num of values in bin[i]
	  private double sum;
	  private TallyHistogram_old hist;
	
	     
	 /*
	  * 
	  * Constructs a `ScaledHistogram` from hist
	  * by normalizing the bin counts so they add up to sum.  
	  * 
	  */
	public ScaledHistogram_old (TallyHistogram_old hist, double sum)
	{ 
	    super(hist.getA(), hist.getB(),hist.getNumBins());
		this.sum=sum;
		this.hist=hist;
		init(hist, hist.getNumBins(),sum);
		
	}
	
	
	public ScaledHistogram_old (TallyHistogram_old hist)
	{ 
	    super(hist.getA(), hist.getB(),hist.getNumBins());
		
		this.hist=hist;
		this.sum=sumTotalCount(hist);
		init(hist, hist.getNumBins());
		
	}
	
	public void init(TallyHistogram_old hist, int s, double sum){
		 co2 = new double[s];	
		 int coCount[]=hist.getCounters();
		 double sumOb=sumTotalCount(hist);
		 for(int i = 0; i < s ; i++)
		       co2[i] = (coCount[i+1]*sum)/sumOb;

	}
	public void init(TallyHistogram_old hist, int s){
		 co2 = new double[s];	
		 int coCount[]=hist.getCounters();
		 for(int i = 0; i < s ; i++)
		     co2[i] = coCount[i+1];

	}
	
	public int sumTotalCount(TallyHistogram_old hist){
		int sum=0;
		int t[]=hist.getCounters();
		for(int i=1; i<hist.getCounters().length-1;i++)
			sum=sum+t[i];
		return sum;
	}
	
	/*
	 * 
	 * Returns an ASH-transformed version of this scaled histogram. The transformed histogram has 
     * the same bin size as the original.  The new frequency in any given bin is the weighted average
     * of the frequencies in the neighboring bins, with weights $(r-d)/r^2$ given to bins that are at distance 
     * $d$ from the target bin, for all $d < r$.  See \cite{tSCO85a}.
	 * 
	 */
	 public ScaledHistogram_old averageShiftedHistogram (int r){
		   ScaledHistogram_old image = (ScaledHistogram_old)super.clone();
		   image.sum=0;
		   double[] coco = new double[hist.getNumBins()];
		   System.arraycopy (co2, 0, coco, 0,hist.getNumBins());
		   for(int k=1;k<=coco.length-1;k++) 
			   coco[k]=getTK(k,r);
		   image.co2=coco;
		   return image;
		
	 }
	
	 /*
	  * return the array count of the average shifted histogram
	  */
	 
	 public double [] getCountersASH(){ 
		 return this.co2; 
	 }
	 
    public double getSum(){
    	return this.sum;
    }
   
    private double getTK(int k, int r)
    {  double d=sumTotalCount(this.hist)*( (hist.getB()-hist.getA())/hist.getNumBins())*r;
       for(int l=1;l<=(r-1);l++ )
        { double a=getCo2(k-l);
    	  double b=getCo2(k+l);
    	  sum=sum+(r-l)*(a+b) ;
    	}
       return (r*co2[k-1]+2*sum)/d;
    }
    
    public TallyHistogram_old getTallyHistogram(){
    	
    	return this.hist;}
    
    
    
    public ScaledHistogram_old averageShiftedHistogram2 (int r){
       ScaledHistogram_old image = (ScaledHistogram_old)super.clone();
   	   image.sum=this.getSum();
       image.co2=getTabTK(r);
   	   return image;	
    }  
    
    public double[] S1k( int r){  
    	double S1k[]= new double [co2.length];
    	S1k[0]=co2[0]; //is equal to s11
    	for(int j=2;j<=co2.length;j++)
    		S1k[j-1]=S1k[j-2]+ getCo2(j-1)-getCo2(j-1-r);
    	return S1k;
    		
    }
    
    
    public double[] S2k(int r){  
    	double S2k[]= new double [co2.length];
    	S2k[0]=getS21(r);
    	for(int j=2;j<=co2.length;j++)
    		S2k[j-1]=S2k[j-2]+ getCo2(j+r)-getCo2(j+1);
        return S2k; 		
    }
    
 public double getCo2(int j){
   double res=0;
  
   if(j>=1 && j<=co2.length)
      res=co2[j-1];
   return res; }

  public double getS21(int r){
	double res=0;
	for(int l=1;l<=(r-1);l++)
	   res=res+getCo2(l); //l+1
	return res;
   }

   public double getT1(int r)
   { double res=0;
      for(int l=1;l<=(r-1);l++)
	     res=res+(r-l)*(getCo2(1-l)+getCo2(1+l));
       return (r*co2[0])+(2*res);	
     }

 
 public double[] getTabTK(int r){
	 double d=sumTotalCount(this.hist)*( (hist.getB()-hist.getA())/hist.getNumBins())*r;
	 double s1k[]=S1k(r);
	 double s2k[]=S2k(r);
	 double tab[]=new double[co2.length];
	 tab[0]=getT1(r);
	 for(int j=2;j<=co2.length;j++)
		 tab[j-1]=(tab[j-2]+s2k[j-2]-s1k[j-2]); 
	 for(int j=0;j<co2.length;j++)
		 tab[j]=tab[j]/d;  
	 return tab; 
 }
 
 
 public ScaledHistogram_old averageShiftedHistogram3 (int r){
	   ScaledHistogram_old image = (ScaledHistogram_old)super.clone();
	   image.sum=0;
       image.co2=getTKbyScottAlgo(r,co2.length);
	   return image;	
 } 
 
 public double[] getTKbyScottAlgo(int r, int nbBin){
	 double d=sumTotalCount(this.hist)*( (hist.getB()-hist.getA())/hist.getNumBins())*r;
	 double fk[]=new double[co2.length];
	 for(int k=1;k<=co2.length;k++)
	     fk[k-1]=0;
	 for(int k=1;k<=co2.length;k++)
	 {	 if(co2[k-1]==0) continue;
	     for(int i=Math.max(1, k-r+1);i<=Math.min(nbBin, k+r-1);i++)
	    	 fk[i-1]=fk[i-1]+co2[k-1];	
	  }
	 for (int i=0;i<fk.length;i++)
		 fk[i]=fk[i]/(d);
	 return fk;
 }
 
}
