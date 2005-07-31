/* $Id: MappedPeriodicT.h,v 1.1.1.1 2001-01-29 08:20:40 paklein Exp $ */
/* created: paklein (04/07/1997)                                          */
/* Manager class for finite deformation elasto-static with 2 additional   */
/* types of kinematic boundary conditions:                                */
/* (1) nodal position mapped forward using a                              */
/* prescribed deformation gradient.                                       */
/* (2) master-slave node pairs - applies to ALL the dof's                 */
/* of the nodes in each pair.                                             */
/* The deformation gradient is specified by the perturbation from         */
/* an identity mapping:                                                   */
/* F = 1 + LTf*F_perturb                                                  */

#ifndef _MAPPED_PERIODIC_T_H
#define _MAPPED_PERIODIC_T_H

/* base class */
#include "KBC_ControllerT.h"

/* direct members */
#include "dMatrixT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "LoadTime.h"

class MappedPeriodicT: public KBC_ControllerT
{
public:

	/* constructor */
	MappedPeriodicT(NodeManagerT& node_manager);

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

	/* schedule for fFperturb */
	int fnumLTf;
	const LoadTime* fLTf;   	
	
	/* specified deformation gradient */
	dMatrixT fFperturb;
	dMatrixT fF; /* F = 1 + LTf*Fperturb */
	  	
/* list of mapped nodes */
iArrayT fMappedNodeList;

	/* master-slave node/dof pairs */
	iArray2DT fSlaveMasterPairs;
	dArrayT   fD_sm; //used in SlaveNodes, (X_s - X_m)	
	
	/* dummy schedule for slave nodes */
	LoadTime fDummySchedule;

	/* shallow copies to main list */
	ArrayT<KBC_CardT> fMappedCards;
	ArrayT<KBC_CardT> fSlaveCards;
};

#endif /* _MAPPED_PERIODIC_T_H */
