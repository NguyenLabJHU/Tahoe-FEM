/* $Id: FrontNodeT.h,v 1.3 2002-10-20 22:39:00 paklein Exp $ */
/* created: paklein (03/19/1999)                                          */

#ifndef _FRONT_NODE_T_H_
#define _FRONT_NODE_T_H_

/* direct members */
#include "dArray2DT.h"


namespace Tahoe {

class FrontNodeT
{
public:

	/* constructor */
	FrontNodeT(int dim, const double* x, const double* v_n, const double* v_t,
		double cone, double da, int num_pts);
	
	/* reset node parameters */
	void Reset(int dim, const double* x, const double* v_n, const double* v_t,
		double cone, double da, int num_pts);
	
	/* coordinates */
	const double* Coords(void) const;

	/* "upward" crack plane normal */
	const double* Normal(void) const;
	
	/* return the direction to the specified sampling point */
	const double* Direction(int i) const;
	
	/* sampling point info */
	const dArray2DT& SampleCoords(void) const;
	const dArray2DT& SampleTransForm(void) const;
	
	/* write point data to output */
	void Write(ostream& out) const;

private:

	/* reset functions */
	void Reset2D(const double* x, const double* v_n, const double* v_t,
		double cone, double da, int num_pts);
	void Reset3D(const double* x, const double* v_n, const double* v_t,
		double cone, double da, int num_pts);  	
	  	
private:

	/* dimension - 2 or 3 */
	int fdim;

	/* coordinates */
	double fx[3];
	
	/* local crack plane normal */
	double fN_nt[3];
	
	/* work space */
	dArray2DT ftip_pts;
	dArray2DT fQ; // rows are transformation tensors

};

/* inlines */
inline void FrontNodeT::Reset(int dim, const double* x, const double* v_n,
	const double* v_t, double cone, double da, int num_pts)
{
	/* check */
	fdim = dim;
	if (fdim == 2)
		Reset2D(x, v_n, v_t, cone, da, num_pts);
	else if (fdim == 3)
		Reset3D(x, v_n, v_t, cone, da, num_pts);
	else
		throw ExceptionT::kOutOfRange;
}

inline const double* FrontNodeT::Coords(void) const { return fx; }
inline const double* FrontNodeT::Normal(void) const { return fN_nt; }
inline const double* FrontNodeT::Direction(int i) const { return fQ(i); }

inline const dArray2DT& FrontNodeT::SampleCoords(void) const { return ftip_pts; }
inline const dArray2DT& FrontNodeT::SampleTransForm(void) const { return fQ; }

} // namespace Tahoe 
#endif /* _FRONT_NODE_T_H_ */
