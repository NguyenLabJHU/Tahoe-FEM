/* $Id: ABAQUS_Ti.h,v 1.1 2004-08-03 00:40:04 paklein Exp $ */
#ifndef _ABAQUS_TI_H_
#define _ABAQUS_TI_H_

/* base class */
#include "ABAQUS_UMAT_BaseT.h"

namespace Tahoe {

class ABAQUS_Ti: public ABAQUS_UMAT_BaseT
{
public:

	/** constructor */
	ABAQUS_Ti(void);

private:

	/* UMAT function wrapper */
	virtual void UMAT(doublereal*, doublereal*, doublereal*, doublereal*,
		doublereal*, doublereal*, doublereal*, doublereal*,
		doublereal*, doublereal*, doublereal*, doublereal*,
		doublereal*, doublereal*, doublereal*, doublereal*,
		doublereal*, doublereal*, char*,
		integer*, integer*, integer*, integer*,
		doublereal*, integer*, doublereal*, doublereal*,
		doublereal*, doublereal*, doublereal*, doublereal*,
		integer*, integer*, integer*, integer*, integer*,
		integer*, ftnlen);

	/* set material output */
	virtual void SetOutputVariables(iArrayT& variable_index,
		ArrayT<StringT>& output_labels);
};

} /* namespace Tahoe */

#endif /* _ABAQUS_TI_H_ */
