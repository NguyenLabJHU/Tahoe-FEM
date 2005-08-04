/* $Id: ABAQUS_UMAT_IsoPlast.h,v 1.1 2005-08-04 07:09:17 paklein Exp $ */
#ifndef _ABAQUS_UMAT_ISOPLAST_H_
#define _ABAQUS_UMAT_ISOPLAST_H_

/* base class */
#include "ABAQUS_UMAT_SS_BaseT.h"

namespace Tahoe {

class ABAQUS_UMAT_IsoPlast: public ABAQUS_UMAT_SS_BaseT
{
public:

	/** constructor */
	ABAQUS_UMAT_IsoPlast(void);

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

#endif /* _ABAQUS_UMAT_ISOPLAST_H_ */
