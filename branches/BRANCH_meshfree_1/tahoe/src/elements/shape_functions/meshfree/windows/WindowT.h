/* $Id: WindowT.h,v 1.5 2001-06-18 18:21:27 paklein Exp $ */

#ifndef _WINDOW_T_H_
#define _WINDOW_T_H_

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class dArrayT;
class dArray2DT;
template <class TYPE> class ArrayT;

/** base class for various support types, and hence different 
 * window functions */
class WindowT
{
  public:
  
	/** constants to identify derive types */
	enum TypeT { kGaussian = 0,
	          kCubicSpline = 1,
	                kBrick = 2};
	
	/** constants to identify the required neighbor search type */
	enum SearchTypeT { kSpherical = 0,
	                kConnectivity = 1};

	/** constructor */
	WindowT(void) { };  

	/** destructor */
	virtual ~WindowT(void) { };

	/** shared parameters.
	 * \return the number of nodal parameters associated
	 * with the window function */
	virtual int NumberOfNodalParameters(void) const = 0;

	/** window function name */
	virtual const char* Name(void) const = 0;

	/** write parameters to output stream */
	virtual void WriteParameters(ostream& out) const;

	/* single point evaluations */
	
	/** window function evaluation.
	 * \param x_n center of the window function
	 * \param param_n array of window function parameters associated with x_n
	 * \param x field point of evaluation
	 * \param order highest order derivative to be calculated
	 * \param w the value at x of the window function centered at x_n
	 * \param Dw window function derivatives: [nsd]
	 * \param DDw window function second derivatives: [nstr] */
	virtual void window(const dArrayT& x_n, const dArrayT& param_n, const dArrayT& x,
		int order, double& w, dArrayT& Dw, dArrayT& DDw) = 0;

	/* multi-point evaluations */
	
	/** window function evaluation. 
	 * \param x_n array of window function centers: [npts] x [nsd]
	 * \param param_n array of window function parameters: [npts] x [nparam]
	 * \param x field point of evaluation
	 * \param order highest order derivative to be calculated
	 * \param w values of the window function: [npts]
	 * \param w values of the window function derivaties: [npts] x [nsd]
	 * \param w values of the window function derivaties: [npts] x [nstr] 
	 * \return the number of points covered by the window function */
	virtual int window(const dArray2DT& x_n, const dArray2DT& param_n, const dArrayT& x,
		int order, dArrayT& w, dArray2DT& Dw, dArray2DT& DDw) = 0;	
	
	/** coverage test.
	 * \return true if the window function centered at x_n covers the
	 * point x. */
	virtual bool Covers(const dArrayT& x_n, const dArrayT& x, const dArrayT& param_n) = 0;

	/** coverage test for multiple points.
	 * \param x_n array of window function centers: [npts] x [nsd]
	 * \param x field point of evaluation
	 * \param covers array of coverage test results: [npts] */
	virtual void Covers(const dArray2DT& x_n, const dArrayT& x, 
			    const dArray2DT& param_n, ArrayT<bool>& covers) = 0;
};

#endif /* _WINDOW_T_H_ */
