/* $Id: ConveyorAdhesiveT.h,v 1.1.2.1 2005-04-08 00:50:24 thao Exp $ */
#ifndef _CONVEYOR_ADHESIVE_T_H_
#define _CONVEYOR_ADHESIVE_T_H_

/* base class */
#include "ConveyorT.h"

/* direct members */
#include "AutoArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "iArrayT.h"
#include "ofstreamT.h"

namespace Tahoe {

/** forward declarations */
class FieldT;

/** conveyor belt */
class ConveyorAdhesiveT: public ConveyorT
{
public:

	/** constructor */
	ConveyorAdhesiveT(NodeManagerT& node_manager, FieldT& field);

	/** initialization */
	virtual void Initialize(ifstreamT& in);

	/** write parameters */
	virtual void WriteParameters(ostream& out) const;

	/** not implemented - there's no going back */
	virtual void Reset(void);

	/** set to initial conditions */
	virtual void InitialCondition(void);

	/** open time interva; */
	virtual void InitStep(void);

	/** computing residual force */
	virtual void FormRHS(void);

	/** apply the update to the solution. Does nothing by default. */
	virtual void Update(const dArrayT& update);

	/** signal that the solution has been found */
	virtual void CloseStep(void);

	/* returns true if the internal force has been changed since
	 * the last time step */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);
	virtual void ReadRestart(ifstreamT& in);
	virtual void WriteRestart(ofstreamT& out) const;

protected:

	/** reset system to new center
	 * \return true if the system focus has been changed */
	virtual bool SetSystemFocus(double focus);

	
	private:
		double fuY_interface;
};

} /* namespace Tahoe */

#endif /* _CONVEYOR_ADHESIVE_T_H_ */
