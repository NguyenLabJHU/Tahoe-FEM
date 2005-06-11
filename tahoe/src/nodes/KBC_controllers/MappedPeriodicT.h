/* $Id: MappedPeriodicT.h,v 1.7 2004-07-15 08:31:21 paklein Exp $ */
/* created: paklein (04/07/1997) */

#ifndef _MAPPED_PERIODIC_T_H
#define _MAPPED_PERIODIC_T_H

/* base class */
#include "KBC_ControllerT.h"

/* direct members */
#include "dMatrixT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "ScheduleT.h"

namespace Tahoe {

/* forward declarations */
class BasicFieldT;

/** boundary condition class for finite deformation elasto-static with 2 
 * additional types of kinematic boundary conditions:
 * <ul>
 * <li> nodal position mapped forward using a prescribed deformation gradient.
 * <li> master-slave node pairs - applies to ALL the dof's of the nodes in each pair.
 * </ul>
 * The deformation gradient is specified by the perturbation from
 * an identity mapping:
 * \f[
 * \mathbf{F}(t) = \mathbf{1} + s(t) \mathbf{F}_{perturb}
 * \f]
 */
class MappedPeriodicT: public KBC_ControllerT
{
public:

	/** constructor */
	MappedPeriodicT(const BasicSupportT& support, BasicFieldT& field);

	/* initial condition */
	virtual void InitialCondition(void);

	/* set BC cards for current step */
	virtual void InitStep(void);

	/* output */
	virtual void WriteOutput(ostream& out) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
protected:

	/** the field */
	BasicFieldT& fField;

	/** schedule for fFperturb */
	const ScheduleT* fSchedule;   	
	
	/* specified deformation gradient */
	dMatrixT fFperturb;
	dMatrixT fF; /* F = 1 + LTf*Fperturb */
	  	
	/* list of mapped nodes */
	iArrayT fMappedNodeList;

	/* master-slave node/dof pairs */
	iArray2DT fSlaveMasterPairs;
	dArrayT   fD_sm; //used in SlaveNodes, (X_s - X_m)	
	
	/* dummy schedule for slave nodes */
	ScheduleT fDummySchedule;

	/* shallow copies to main list */
	ArrayT<KBC_CardT> fMappedCards;
	ArrayT<KBC_CardT> fSlaveCards;
};

} // namespace Tahoe 
#endif /* _MAPPED_PERIODIC_T_H */
