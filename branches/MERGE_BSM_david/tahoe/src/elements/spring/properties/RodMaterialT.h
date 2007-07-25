/* $Id: RodMaterialT.h,v 1.7 2005-11-06 00:37:58 paklein Exp $ */
/* created: paklein (11/20/1996) */
#ifndef _RODMATERIALT_H_
#define _RODMATERIALT_H_

#include "Environment.h"

namespace Tahoe {

/* forward declarations */
class ScheduleT;

/** pair interactions */
class RodMaterialT
{
public:

	/** constructor */
	RodMaterialT(double mass);

	/** destructor */
	virtual ~RodMaterialT(void);

	/** return the particle mass */
	double Mass(void) const { return fMass; };
	
	/** potential function and derivatives */
	virtual double Potential(double rmag, double Rmag) const = 0;
	virtual double DPotential(double rmag, double Rmag) const = 0;
	virtual double DDPotential(double rmag, double Rmag) const = 0;

#if 0
	/** thermal accessors */
	int ThermalScheduleNumber(void) const;
	void SetThermalSchedule(const ScheduleT* LTfPtr);

	/** returns true if the material has internal forces in the unloaded
	 * configuration, ie thermal strains */
	virtual int HasInternalStrain(void) const;
#endif

protected:

	double fMass;
};

} // namespace Tahoe 
#endif /* _RODMATERIALT_H_ */
