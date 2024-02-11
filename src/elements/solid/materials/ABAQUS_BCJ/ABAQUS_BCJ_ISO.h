/* $Id: ABAQUS_BCJ_ISO.h,v 1.3 2004-08-01 20:42:35 paklein Exp $ */
#ifndef _ABAQUS_BCJ_ISO_H_
#define _ABAQUS_BCJ_ISO_H_

/* base class */
#include "ABAQUS_UMAT_BaseT.h"

/* library support */
#ifdef __F2C__

namespace Tahoe {

class ABAQUS_BCJ_ISO: public ABAQUS_UMAT_BaseT
{
public:

	/** constructor */
	ABAQUS_BCJ_ISO(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

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

#else /* __F2C__ */

#ifndef __MWERKS__
#error "ABAQUS_BCJ_ISO requires __F2C__"
#endif

#endif /* __F2C__ */

#endif /* _ABAQUS_BCJ_ISO_H_ */
