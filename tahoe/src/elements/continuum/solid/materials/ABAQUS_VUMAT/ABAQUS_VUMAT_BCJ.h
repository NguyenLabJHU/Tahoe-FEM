/* $Id: ABAQUS_VUMAT_BCJ.h,v 1.4 2002-07-02 19:55:33 cjkimme Exp $ */

#ifndef _ABAQUS_VUMAT_BCJ_H_
#define _ABAQUS_VUMAT_BCJ_H_

/* base class */
#include "ABAQUS_VUMAT_BaseT.h"

/* library support options */
#ifdef __F2C__


namespace Tahoe {

class ABAQUS_VUMAT_BCJ: public ABAQUS_VUMAT_BaseT
{
public:

	/* constructor */
	ABAQUS_VUMAT_BCJ(ifstreamT& in, const FiniteStrainT& element);

private:

	/* VUMAT function wrapper */

//this need to be changed to a VUMAT wrapper	
	
	virtual void VUMAT(integer*, integer*, integer*, integer*, integer*, integer*, integer*, doublereal*,
		doublereal*, doublereal*, char*, doublereal*, doublereal*, doublereal*,
                doublereal*, doublereal*, doublereal*, doublereal*, doublereal*, doublereal*,
		doublereal*, doublereal*, doublereal*, doublereal*, doublereal*, doublereal*,
                doublereal*, doublereal*, doublereal*, doublereal*, doublereal*, doublereal*,
                doublereal*);

	//	virtual void UMAT(doublereal*, doublereal*, doublereal*, doublereal*,
	//doublereal*, doublereal*, doublereal*, doublereal*,
	//doublereal*, doublereal*, doublereal*, doublereal*,
	//doublereal*, doublereal*, doublereal*, doublereal*,
	//doublereal*, doublereal*, char*,
	//integer*, integer*, integer*, integer*,
	//doublereal*, integer*, doublereal*, doublereal*,
	//doublereal*, doublereal*, doublereal*, doublereal*,
	//integer*, integer*, integer*, integer*, integer*,
	//integer*, ftnlen);

	/* set material output */
	virtual void SetOutputVariables(iArrayT& variable_index,
		ArrayT<StringT>& output_labels);
};

} // namespace Tahoe 
#endif /* __F2C__ */
#endif /* _ABAQUS_VUMAT_BCJ_H_ */
