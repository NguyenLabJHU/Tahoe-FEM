/* $Id: ThermalDilatationT.h,v 1.2 2001-09-15 01:17:36 paklein Exp $ */
/* created: paklein (08/25/1996) */

#ifndef _THERMALDILAT_H_
#define _THERMALDILAT_H_

#include "Environment.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class LoadTime;

class ThermalDilatationT
{
public:

	/* constructor */
	ThermalDilatationT(ifstreamT& in);
	
	/* to set LTf pointer */
	int LTfNumber(void) const;
	void SetLTfPtr(const LoadTime* LTf);

	/* returns true if active */
	bool IsActive(void) const;
	
	/* I/O functions */
	void Print(ostream& out) const;

	/* returns the current elongation factor */
	double PercentElongation(void) const;
							
private:
	
	double			fPercentElongation;
	int				LTfnum;
	const LoadTime*	LTfPtr;	
};

/* inline functions */

/* returns true if active */
inline bool ThermalDilatationT::IsActive(void) const { return fPercentElongation != 0.0; }

/* set LTf pointer */
inline int ThermalDilatationT::LTfNumber(void) const { return LTfnum; }
inline void ThermalDilatationT::SetLTfPtr(const LoadTime* LTf) { LTfPtr = LTf; }

#endif /* _THERMALDILAT_H_ */
