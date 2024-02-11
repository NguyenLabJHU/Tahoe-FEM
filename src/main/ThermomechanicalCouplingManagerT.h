/* $Id: ThermomechanicalCouplingManagerT.h,v 1.2 2005-04-06 15:42:00 paklein Exp $  */
#ifndef _THERMOMECHANICAL_COUPLING_MANAGER_H_
#define _THERMOMECHANICAL_COUPLING_MANAGER_H_

/* element configuration header */
#include "ElementsConfig.h"
#include "DevelopmentElementsConfig.h"
#if defined(BRIDGING_ELEMENT) && defined(BRIDGING_ELEMENT_DEV)

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

#endif /* BRIDGING_ELEMENT && BRIDGING_ELEMENT_DEV */
#endif /* _THERMOMECHANICAL_COUPLING_MANAGER_H_  */
