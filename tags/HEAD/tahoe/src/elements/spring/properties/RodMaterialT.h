/* $Id: RodMaterialT.h,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (11/20/1996)                                          */

#ifndef _RODMATERIALT_H_
#define _RODMATERIALT_H_

#include "Environment.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class ThermalDilatationT;
class LoadTime;

class RodMaterialT
{
public:

	/* constructor */
	RodMaterialT(ifstreamT& in);

	/* destructor */
	virtual ~RodMaterialT(void);

	/* I/O functions */
	virtual void Print(ostream& out) const = 0;
	virtual void PrintName(ostream& out) const = 0;

	/* print parameters */
	void PrintParameters(ostream& out) const;
	
	/* potential function and derivatives */
	virtual double Potential(double rmag, double Rmag) const = 0;
	virtual double DPotential(double rmag, double Rmag) const = 0;
	virtual double DDPotential(double rmag, double Rmag) const = 0;

	/* thermal accessors */
	int ThermalLTfNumber(void) const;
	void SetThermalLTfPtr(const LoadTime* LTfPtr);

	/* returns true if the material has internal forces in the unloaded
	 * configuration, ie thermal strains */
	virtual int HasInternalStrain(void) const;

protected:

	ThermalDilatationT* fThermal;

};

#endif /* _RODMATERIALT_H_ */
