/* $Id: BridgingScaleManagerT.h,v 1.5 2005-04-06 15:42:00 paklein Exp $ */
#ifndef _BRIDGING_SCALE_MANAGER_H_
#define _BRIDGING_SCALE_MANAGER_H_

/* element configuration header */
#include "ElementsConfig.h"
#include "DevelopmentElementsConfig.h"
#if defined(BRIDGING_ELEMENT) && defined(BRIDGING_ELEMENT_DEV)

/* base class  */
#include "MultiManagerT.h"

namespace Tahoe {

/* forward declarations */
class FEManagerT_THK;

/** manager for dynamic bridging scale calculations */
class BridgingScaleManagerT: public MultiManagerT
{
public:

	/** constructor */
	BridgingScaleManagerT(const StringT& input_file, ofstreamT& output, CommunicatorT& comm,
		const ArrayT<StringT>& argv, TaskT task);

	/** destructor */
	virtual ~BridgingScaleManagerT(void);

	/** solve all the time sequences */
	virtual void Solve(void);

	/** (re-)set system to initial conditions */
	virtual ExceptionT::CodeT InitialCondition(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/** calculate internal+external force for the given displacement u */
	const dArray2DT& TotalForce(const StringT& field_name, const dArray2DT& field_values, 
		FEManagerT_bridging& bridging, dArray2DT& rhs_2D) const;

	/** calculate internal force for the given displacement u */
	const dArray2DT& InternalForce(const StringT& field_name, const dArray2DT& field_values, 
		FEManagerT_bridging& bridging) const;

private:

	/** cast of MultiManagerT::fFine to FEManagerT_THK */
	FEManagerT_THK* fFine_THK;	
};

} /* namespace Tahoe */

#endif /* BRIDGING_ELEMENT && BRIDGING_ELEMENT_DEV */
#endif /* _BRIDGING_SCALE_MANAGER_H_ */
