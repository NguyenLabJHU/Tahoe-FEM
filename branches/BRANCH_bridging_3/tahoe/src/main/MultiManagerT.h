/* $Id: MultiManagerT.h,v 1.1.2.1 2003-10-16 12:56:14 paklein Exp $ */
#ifndef _MULTI_MANAGER_H_
#define _MULTI_MANAGER_H_

/* base class  */
#include "FEManagerT.h"

namespace Tahoe {

class FEManagerT_bridging;

/** wrapper for running multiple FEManagerT with a single solver */
class MultiManagerT: public FEManagerT
{
public:

	/** constructor */
	MultiManagerT(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
		FEManagerT_bridging* fine, FEManagerT_bridging* coarse);

	/** destructor */
	virtual ~MultiManagerT(void);
	
	/** initialize members */
	virtual void Initialize(InitCodeT init = kFull);

	/** (re-)set the equation number for the given group */
	virtual void SetEquationSystem(int group);

	/** \name solution steps */
	/*@{*/
	/** initialize the current time increment for all groups */
	virtual ExceptionT::CodeT InitStep(void);

	/** close the current time increment for all groups */
	virtual ExceptionT::CodeT CloseStep(void);

	/** called for all groups if the solution procedure for any group fails */
	virtual ExceptionT::CodeT ResetStep(void);
	/*@}*/

	/** \name solution messaging */
	/*@{*/
	/** compute LHS-side matrix and assemble to solver */
	virtual void FormLHS(int group, GlobalT::SystemTypeT sys_type) const;

	/** compute RHS-side */
	virtual void FormRHS(int group) const;

	/** send update of the solution to the NodeManagerT */
	virtual void Update(int group, const dArrayT& update);

	/** system relaxation */
	virtual GlobalT::RelaxCodeT RelaxSystem(int group) const;
	/*@}*/

private:

	/** \name sub-managers */
	/*@{*/
	FEManagerT_bridging* fFine;
	FEManagerT_bridging* fCoarse;	
	/*@}*/
	
	/** \name equations for each sub-manager */
	/*@{*/
	iArrayT fEqnos1;
	iArrayT fEqnos2;
	/*@}*/

	/** work space */
	dArray2DT fFieldAtGhosts;
};

} /* namespace Tahoe */

#endif /* _MULTI_MANAGER_H_ */
