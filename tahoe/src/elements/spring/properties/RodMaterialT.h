/* $Id: RodMaterialT.h,v 1.5 2002-07-05 22:28:28 paklein Exp $ */
/* created: paklein (11/20/1996) */

#ifndef _RODMATERIALT_H_
#define _RODMATERIALT_H_

#include "Environment.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class ThermalDilatationT;
class ScheduleT;

/** pair interactions */
class RodMaterialT
{
public:

	/** constructor */
	RodMaterialT(ifstreamT& in);

	/** destructor */
	virtual ~RodMaterialT(void);

	/** I/O functions */
	virtual void Print(ostream& out) const = 0;
	virtual void PrintName(ostream& out) const = 0;

	/** print parameters */
	void PrintParameters(ostream& out) const;
	
	/** return the particle mass */
	double Mass(void) const { return fMass; };
	
	/** potential function and derivatives */
	virtual double Potential(double rmag, double Rmag) const = 0;
	virtual double DPotential(double rmag, double Rmag) const = 0;
	virtual double DDPotential(double rmag, double Rmag) const = 0;

	/** thermal accessors */
	int ThermalScheduleNumber(void) const;
	void SetThermalSchedule(const ScheduleT* LTfPtr);

	/** returns true if the material has internal forces in the unloaded
	 * configuration, ie thermal strains */
	virtual int HasInternalStrain(void) const;

protected:

	double fMass;
	ThermalDilatationT* fThermal;
};

} // namespace Tahoe 
#endif /* _RODMATERIALT_H_ */
