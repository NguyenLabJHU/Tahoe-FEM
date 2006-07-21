/* $Id: FieldMFAugLagMultT.h,v 1.1 2005-04-12 15:34:40 paklein Exp $ */
#ifndef _FIELD_MF_AUG_LAG_MULT_T_H_
#define _FIELD_MF_AUG_LAG_MULT_T_H_

#include "ElementsConfig.h"
#ifdef CONTINUUM_ELEMENT

/* base class */
#include "MFAugLagMultT.h"

namespace Tahoe {

/* forward declarations */
class ScheduleT;

/** prescribed field with augmented Lagrangian enforcement of
 * non-interpenetration */
class FieldMFAugLagMultT: public MFAugLagMultT
{
public:

	/** constructor */
	FieldMFAugLagMultT(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/* accumulate the constraint force vector fConstraintForce */
	virtual void ComputeConstraintValues(double kforce);

private:

	/** all constrained nodes */
	iArrayT fAllNodes;

	/** prescribed field */
	dArray2DT fPrescribedField;
	
	/** schedule to scale the prescribed field */
	const ScheduleT* fSchedule;
};

} // namespace Tahoe 

#endif /* CONTINUUM_ELEMENT */

#endif /* _FIELD_MF_AUG_LAG_MULT_T_H_ */
