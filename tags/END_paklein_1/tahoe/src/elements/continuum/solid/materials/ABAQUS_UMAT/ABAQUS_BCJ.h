/* $Id: ABAQUS_BCJ.h,v 1.3 2002-07-02 19:55:31 cjkimme Exp $ */
/* created: paklein (05/09/2000)                                          */

#ifndef _ABAQUS_BCJ_H_
#define _ABAQUS_BCJ_H_

/* base class */
#include "ABAQUS_UMAT_BaseT.h"

/* library support options */
#ifdef __F2C__


namespace Tahoe {

class ABAQUS_BCJ: public ABAQUS_UMAT_BaseT
{
public:

	/* constructor */
	ABAQUS_BCJ(ifstreamT& in, const FiniteStrainT& element);

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

} // namespace Tahoe 
#endif /* __F2C__ */
#endif /* _ABAQUS_BCJ_H_ */
