/* $Id: WindowT.h,v 1.1 2001-06-13 20:53:57 paklein Exp $ */

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
	virtual int NumberOfParameters(void) const = 0;

	/** window function name */
	virtual const char* Name(void) const = 0;

	/** write parameters to output stream */
	virtual void WriteParameters(ostream& out) const;

	/* single point evaluations */
	
	/** window function evaluation.
	 * \param x_0 center of the window function
	 * \param param array of window function parameters associated with x_0
	 * \param x field point
	 * \param w the value at x of the window function centered at x_0 */
	virtual void window(const dArrayT& x_0, const dArrayT& param, const dArrayT& x,
		double& w) = 0;

	/** window function derivatives.
	 * \param x_0 center of the window function
	 * \param param array of window function parameters associated with x_0
	 * \param x field point
	 * \param Dw window function derivatives: [nsd] */
	virtual void Dwindow(const dArrayT& x_0, const dArrayT& param, const dArrayT& x,
		dArrayT& Dw) = 0;

	/** window function second derivatives.
	 * \param x_0 center of the window function
	 * \param param array of window function parameters associated with x_0
	 * \param x field point
	 * \param DDw window function second derivatives: [nstr] */
	virtual void DDwindow(const dArrayT& x_0, const dArrayT& param, const dArrayT& x,
		dArrayT& DDw) = 0;

	/* multi-point evaluations */
	
	/** window function evaluation. 
	 * \param x_0 center of the window function
	 * \param param array of window function parameters: [npts] x [nparam]
	 * \param x array of field point coordinates: [npts] x [nsd]
	 * \param w values of the window function: [npts] */
	virtual void window(const dArrayT& x_0, const dArray2DT& param, const dArray2DT& x,
		dArrayT& w) = 0;	

	/** window function derivatives. 
	 * \param x_0 center of the window function
	 * \param param array of window function parameters: [npts] x [nparam]
	 * \param x array of field point coordinates: [npts] x [nsd]
	 * \param w values of the window function derivaties: [npts] x [nsd] */
	virtual void Dwindow(const dArrayT& x_0, const dArray2DT& param, const dArray2DT& x,
		dArray2DT& Dw) = 0;	

	/** window function second derivatives. 
	 * \param x_0 center of the window function
	 * \param param array of window function parameters: [npts] x [nparam]
	 * \param x array of field point coordinates: [npts] x [nsd]
	 * \param w values of the window function derivaties: [npts] x [nstr] */
	virtual void DDwindow(const dArrayT& x_0, const dArray2DT& param, const dArray2DT& x,
		dArray2DT& DDw) = 0;	
	
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
