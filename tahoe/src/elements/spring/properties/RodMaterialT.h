/* $Id: RodMaterialT.h,v 1.2.2.1 2002-06-27 18:03:54 cjkimme Exp $ */
/* created: paklein (11/20/1996)                                          */

#ifndef _RODMATERIALT_H_
#define _RODMATERIALT_H_

#include "Environment.h"

/* forward declarations */
#include "ios_fwd_decl.h"

namespace Tahoe {

class ifstreamT;
class ThermalDilatationT;
class ScheduleT;

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
	int ThermalScheduleNumber(void) const;
	void SetThermalSchedule(const ScheduleT* LTfPtr);

	/* returns true if the material has internal forces in the unloaded
	 * configuration, ie thermal strains */
	virtual int HasInternalStrain(void) const;

protected:

	ThermalDilatationT* fThermal;

};

} // namespace Tahoe 
#endif /* _RODMATERIALT_H_ */
