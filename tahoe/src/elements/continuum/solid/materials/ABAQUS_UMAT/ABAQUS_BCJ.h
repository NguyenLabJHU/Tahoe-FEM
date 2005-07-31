/* $Id: ABAQUS_BCJ.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (05/09/2000)                                          */

#ifndef _ABAQUS_BCJ_H_
#define _ABAQUS_BCJ_H_

/* base class */
#include "ABAQUS_UMAT_BaseT.h"

/* library support options */
#ifdef __F2C__

class ABAQUS_BCJ: public ABAQUS_UMAT_BaseT
{
public:

	/* constructor */
	ABAQUS_BCJ(ifstreamT& in, const ElasticT& element);

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

#endif /* __F2C__ */
#endif /* _ABAQUS_BCJ_H_ */
