/* $Id: WindowT.h,v 1.2 2001-06-14 22:08:15 paklein Exp $ */

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
	 * \param x_0 center of the window function
	 * \param param array of window function parameters associated with x_0
	 * \param x field point
	 * \param order highest order derivative to be calculated
	 * \param w the value at x of the window function centered at x_0
	 * \param Dw window function derivatives: [nsd]
	 * \param DDw window function second derivatives: [nstr] */
	virtual void window(const dArrayT& x_0, const dArrayT& param, const dArrayT& x,
		int order, double& w, dArrayT& Dw, dArrayT& DDw) = 0;

	/* multi-point evaluations */
	
	/** window function evaluation. 
	 * \param x_0 center of the window function
	 * \param param array of window function parameters: [npts] x [nparam]
	 * \param x array of field point coordinates: [npts] x [nsd]
	 * \param order highest order derivative to be calculated
	 * \param w values of the window function: [npts]
	 * \param w values of the window function derivaties: [npts] x [nsd]
	 * \param w values of the window function derivaties: [npts] x [nstr] 
	 * \return the number of points covered by the window function */
	virtual int window(const dArrayT& x_0, const dArray2DT& param, const dArray2DT& x,
		int order, dArrayT& w, dArray2DT& Dw, dArray2DT& DDw) = 0;	
	
	/** coverage test.
	 * \return true if the window function centered at x_0 covers the
	 * point x. */
	virtual bool Covers(const dArrayT& x_0, const dArrayT& x) = 0;

	/** coverage test for multiple points.
	 * \param x_0 center of the window function
	 * \param x array of field point coordinates: [npts] x [nsd]
	 * \param covers array of coverage test results: [npts] */
	virtual void Covers(const dArrayT& x_0, const dArray2DT& x, const ArrayT<bool>& covers) = 0;
};

#endif /* _WINDOW_T_H_ */
