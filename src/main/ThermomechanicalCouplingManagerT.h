/* $Id  */
#ifndef _THERMOMECHANICAL_COUPLING_MANAGER_H_
#define _THERMOMECHANICAL_COUPLING_MANAGER_H_

/* element configuration header */
#include "ElementsConfig.h"
#ifdef BRIDGING_ELEMENT

/* base class  */
#include "MultiManagerT.h"

namespace Tahoe {

/* forward declarations */
class ParticleT;

/** manager for thermomechanical coupling  */
class ThermomechanicalCouplingManagerT: public MultiManagerT
{
public:

	/** constructor */
	ThermomechanicalCouplingManagerT(const StringT& input_file, ofstreamT& output, CommunicatorT& comm,
		const ArrayT<StringT>& argv, TaskT task);

	/** destructor */
	virtual ~ThermomechanicalCouplingManagerT(void);

	/** solve all the time sequences */
	virtual void Solve(void);

	/** (re-)set system to initial conditions */
	virtual ExceptionT::CodeT InitialCondition(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:	
	int fCoarseOutputID;
	ParticleT* fParticles;

};

} /* namespace Tahoe */

#endif /* BRIDGING_ELEMENT */
#endif /* _THERMOMECHANICAL_COUPLING_MANAGER_H_  */
