/* $Id: ParameterListT.h,v 1.2 2002-09-03 07:54:08 paklein Exp $ */
#ifndef _PARAMETER_LIST_T_H_
#define _PARAMETER_LIST_T_H_

/* direct members */
#include "ParameterT.h"
#include "AutoArrayT.h"

namespace Tahoe {

/** list of parameters */
class ParameterListT
{
public:

	/** parameter occurrence */
	enum OccurrenceT {
		Once,     /**< exactly once */
		ZeroPlus, /**< zero or more times */
		OnePlus   /**< one or more times */
	};

	/** constructor */
	ParameterListT(const StringT& name): fName(name) { };
	
	/** (re-)set namespace name */
	void SetNameSpace(const StringT& ns_name) { fNameSpace = ns_name; };
	
	/** add parameter. Returns true of there where no conflicts with
	 * existing parameters. The names of parameters cannot be repeated */
	bool AddParameter(const ParameterT& param, OccurrenceT occur); 
	
protected:

	/** parameters name space */
	StringT fNameSpace;

	/** list name */
	StringT fName;

	/** \name simple parameters */
	/*@{*/
	AutoArrayT<ParameterT>  fParameters;
	AutoArrayT<OccurrenceT> fParametersOccur;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _PARAMETER_LIST_T_H_ */
