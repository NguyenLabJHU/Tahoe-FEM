/* $Id: ConveyorSymT.h,v 1.1 2004-12-29 16:00:58 thao Exp $ */
#ifndef _CONVEYOR_SYM_T_H_
#define _CONVEYOR_SYM_T_H_

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
class ConveyorSymT: public ConveyorT
{
public:

	/** constructor */
	ConveyorSymT(const BasicSupportT& support, FieldT& field);

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

	/** \name restart functions */
	/*@{*/
	virtual void ReadRestart(ifstreamT& in);
	virtual void WriteRestart(ofstreamT& out) const;
	/*@}*/

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
	/** reset system to new center
	 * \return true if the system focus has been changed */
	virtual bool SetSystemFocus(double focus);

private:
	double fuY_interface;
};

} /* namespace Tahoe */

#endif /* _CONVEYOR_T_H_ */
