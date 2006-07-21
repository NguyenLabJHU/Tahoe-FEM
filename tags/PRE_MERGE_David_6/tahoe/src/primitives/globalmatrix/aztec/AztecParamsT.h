/* $Id: AztecParamsT.h,v 1.1 2005-04-05 16:07:07 paklein Exp $ */
#ifndef _AZTEC_PARAMS_T_H_
#define _AZTEC_PARAMS_T_H_

/* library support options */
#ifdef __AZTEC__

/* base class */
#include "ParameterInterfaceT.h"

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"

namespace Tahoe {

class AztecParamsT: public ParameterInterfaceT
{
public:

	/** constructor */
	AztecParamsT(void);

	/** set solver options */
	void SetAztecOptions(int* options, double* params) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/** option name to index conversion */
	int OptionNameToIndex(const char* name) const;
	
	/** parameter name to index conversion */
	int ParamNameToIndex(const char* name) const;

private:

	/** \name Aztec option and parameter overrides */
	/*@{*/
	iArrayT AZ_options_dex;
	iArrayT AZ_options;
	iArrayT AZ_params_dex;
	dArrayT AZ_params;
	/*@}*/
};

} /* namespace Tahoe */
#endif /* __AZTEC__ */
#endif /* _AZTEC_PARAMS_T_H_ */
