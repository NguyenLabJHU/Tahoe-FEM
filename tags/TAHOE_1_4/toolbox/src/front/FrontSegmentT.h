/* $Id: FrontSegmentT.h,v 1.2 2002-07-02 19:57:10 cjkimme Exp $ */
/* created: paklein (03/19/1999)                                          */

#ifndef _FRONT_SEG_T_H_
#define _FRONT_SEG_T_H_


namespace Tahoe {

class FrontSegmentT
{
public:

	/* constructor */
	FrontSegmentT(const double* A, const double* B, const double* C);  	
	FrontSegmentT(const double* A, const double* B, const FrontSegmentT& source);
	
	/* reset coordinates */
	void Reset(const double* A, const double* B, const double* C);
	void Reset(const double* A, const double* B, const FrontSegmentT& source);
	
	/* endpoints */
	const double* x1(void) const;
	const double* x2(void) const;

	void x1(double* x) const; // x must be length 3
	void x2(double* x) const; // x must be length 3
	
	/* directions */
	const double* N_n(void) const;
	const double* N_t(void) const;

	/* segment length */
	double Length(void) const;

	/* compute midpoint */
	void MidPoint(double* x) const; // x must be length 3

private:

	/* endpoints */
	double fA[3];
	double fB[3];

	/* directions */
	double fN_n[3]; // normal
	double fN_t[3]; // tangent
};

/* inlines */
inline const double* FrontSegmentT::x1(void) const { return fA; }
inline const double* FrontSegmentT::x2(void) const { return fB; }
inline void FrontSegmentT::x1(double* x) const
{
	x[0] = fA[0];
	x[1] = fA[1];
	x[2] = fA[2];
}

inline void FrontSegmentT::x2(double* x) const
{
	x[0] = fB[0];
	x[1] = fB[1];
	x[2] = fB[2];
}

inline const double* FrontSegmentT::N_n(void) const { return fN_n; }
inline const double* FrontSegmentT::N_t(void) const { return fN_t; }

} // namespace Tahoe 
#endif /* _FRONT_SEG_T_H_ */
