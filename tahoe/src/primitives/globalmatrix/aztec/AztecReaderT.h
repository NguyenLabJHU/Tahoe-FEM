/* $Id: AztecReaderT.h,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (08/12/1998)                                          */
/* utility to read Aztec options and parameters.                          */
/* Option input is in the form:                                           */
/* [option] = [value]                                                     */
/* where [option] is either the name or integer constant of the           */
/* Aztec option, and [value] is either the name of integer                */
/* constant of the option.                                                */
/* ***** NOTE: only integers supported in rhs for now *******             */
/* Parameter input is in the form:                                        */
/* [parameter] = [value]                                                  */
/* where [parameter] is either the name or integer constant of the        */
/* Aztec parameter, and [value] is a floating point value.                */

#ifndef _AZTEC_READER_T_H_
#define _AZTEC_READER_T_H_

#include "Environment.h"

/* library support options */
#ifdef __AZTEC__

/* forward declarations */
class ifstreamT;

class AztecReaderT
{
public:

	/* constructor */
	AztecReaderT(void);
	
	/* read option and return integer value */
	void ReadOption(ifstreamT& in, int& dex, int& value);

	/* read parameter and return double value */
	void ReadParameter(ifstreamT& in, int& dex, double& value);

private:

	/* option name to index conversion */
	int OptionNameToIndex(const char* name) const;
	
	/* parameter name to index conversion */
	int ParamNameToIndex(const char* name) const;
};

/* library support options */
#endif /* __AZTEC__ */
#endif /* _AZTEC_READER_T_H_ */
