/* $Id: ParameterInterfaceT.h,v 1.3 2003-04-27 07:30:12 paklein Exp $ */
#ifndef _PARAMETER_INTERFACE_T_H_
#define _PARAMETER_INTERFACE_T_H_

/* direct members */
#include "ParameterListT.h"

namespace Tahoe {

/* forward declarations */
class StringT;

/** abstract interface for classes which define and use parameters */
class ParameterInterfaceT
{
public:

	/** constructor */
	ParameterInterfaceT(void) {};

	/** identifier */
	virtual const StringT& Name(void) const = 0;

	/** \name parameters */
	/*@{*/
	/** accept completed parameter list */
	virtual void SetParameters(const ParameterListT& list);
	
	/** build complete parameter list description */
	virtual void DefineParameters(ParameterListT& list) const;
	/*@}*/

	/** \name sub-lists */
	/*@{*/
	/** return the list of sub-list names */
	virtual void SubListNames(ArrayT<StringT>& list, ArrayT<ParameterListT::OccurrenceT>& occur) const;
	
	/** a pointer to the ParameterInterfaceT of the given sublist
	 * or NULL if the name is invalid. The objects returned must remain
	 * valid as long as this. */
	virtual ParameterInterfaceT* SubList(const StringT& list_name);
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _PARAMETER_SOURCE_T_H_ */
