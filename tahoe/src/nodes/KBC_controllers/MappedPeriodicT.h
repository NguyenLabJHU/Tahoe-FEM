/* $Id: MappedPeriodicT.h,v 1.6 2002-07-05 22:28:31 paklein Exp $ */
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

	/* constructor */
	MappedPeriodicT(NodeManagerT& node_manager, BasicFieldT& field);

	/* initialize data - called immediately after construction */
	virtual void Initialize(ifstreamT& in);
	virtual void WriteParameters(ostream& out) const;

	/* initial condition */
	virtual void InitialCondition(void);

	/* set BC cards for current step */
	virtual void InitStep(void);

	/* output */
	virtual void WriteOutput(ostream& out) const;
	
protected:

	/** the field */
	BasicFieldT& fField;

	/* schedule for fFperturb */
	int fnumLTf;
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
