/* $Id: LocalParabolaT.h,v 1.2 2002-07-02 19:57:18 cjkimme Exp $ */
/* created: paklein (01/28/1997)                                          */
/* LocalParabolaT.h                                                       */
/* Interface for a successively fitted parabolic approximation.  Data     */
/* points can be added any number of times, but there should be at least  */
/* three data points before asking for properties of the parabola.  If    */
/* more than 3 data points are added, only the 3 most recent points are   */
/* used to make the parabolic approximation.                              */
/* Note: no protection has been added to test that the 3 most recent data */
/* points are non-colinear.                                               */

#ifndef _LOCALPARABOLA_T_H_
#define _LOCALPARABOLA_T_H_


namespace Tahoe {

class LocalParabolaT
{
public:

	/*
	 * Constructor
	 */
	LocalParabolaT(void);

	/*
	 * Reset - clear all data
	 */
	void Reset(void);

	/*
	 * Add point to fitting data, discarding the point in the current
	 * data farthest in y from the new point if at least 3 points
	 * have already been added.
	 */
	void AddPoint(double x, double y);

	/*
	 * Minima - return the critical point information through value.
	 * The functions return 0 if they are called before at least 3
	 * data points have been added, leaving value unmodified.
	 */
	int CriticalPointValue(double& value) const;
	int CriticalPointLocation(double& value) const;
	int CriticalPointConcavity(double& value) const;
	    	   	    	
	/*
	 * Returning values - returns 0.0 if less that 3 data points
	 * have been added.
	 */
	double operator()(double x) const;

private:

	/*
	 * Compute the parabolic fit
	 */
	void FitData(void) const;
	void FitDataWrapper(void);
		
private:

	/* status flags */
	int fPointCount;
	int fCurrentCoeffs;

	/* data points */
	double fx[3];
	double fy[3];
	
	/* coefficients y = fa0 + fa1 x + fa2 x^2 */
	double fa0;
	double fa1;
	double fa2;
	
	/* critical point values */
	double fxmin;  		    	   	    	
	double fymin;
};

} // namespace Tahoe 
#endif /* _LOCALPARABOLA_T_H_ */
