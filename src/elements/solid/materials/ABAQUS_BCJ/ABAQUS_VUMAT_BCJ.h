/* $Id: ABAQUS_VUMAT_BCJ.h,v 1.2 2003-09-06 08:42:54 paklein Exp $ */
#ifndef _ABAQUS_VUMAT_BCJ_H_
#define _ABAQUS_VUMAT_BCJ_H_

/* base class */
#include "ABAQUS_VUMAT_BaseT.h"

/* library support */
#ifdef __F2C__

namespace Tahoe {

class ABAQUS_VUMAT_BCJ: public ABAQUS_VUMAT_BaseT
{
public:

	/* constructor */
	ABAQUS_VUMAT_BCJ(ifstreamT& in, const FSMatSupportT& support);

private:

	/* VUMAT function wrapper */
	virtual void VUMAT(integer*, integer*, integer*, integer*, integer*, integer*, integer*, doublereal*,
		doublereal*, doublereal*, char*, doublereal*, doublereal*, doublereal*,
                doublereal*, doublereal*, doublereal*, doublereal*, doublereal*, doublereal*,
		doublereal*, doublereal*, doublereal*, doublereal*, doublereal*, doublereal*,
                doublereal*, doublereal*, doublereal*, doublereal*, doublereal*, doublereal*,
                doublereal*);

	/* set material output */
	virtual void SetOutputVariables(iArrayT& variable_index,
		ArrayT<StringT>& output_labels);
};

#else /* __F2C__ */

#ifndef __MWERKS__
#error "ABAQUS_VUMAT_BCJ requires __F2C__"
#endif

#endif /* __F2C__ */

} /* namespace Tahoe */

#endif /* _ABAQUS_VUMAT_BCJ_H_ */
