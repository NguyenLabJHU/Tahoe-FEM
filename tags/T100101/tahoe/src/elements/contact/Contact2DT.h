/* $Id: Contact2DT.h,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (05/26/1999)                                          */

#ifndef _CONTACT2D_T_H_
#define _CONTACT2D_T_H_

/* base classes */
#include "ContactT.h"

/* direct members */
#include "AutoArrayT.h"

/* forward declarations */
class iGridManager2DT;

class Contact2DT: public ContactT
{
public:

	/* constructor */
	Contact2DT(FEManagerT& fe_manager);

	/* destructor */
	virtual ~Contact2DT(void);

	/* allocates space and reads connectivity data */
	virtual void Initialize(void);

protected:

	/* steps in setting contact configuration */
	virtual bool SetActiveInteractions(void); // "internal" data
	virtual void SetConnectivities(void); // "external" data - interface to FEManager

private:

	/* set working arrays */
	void SetShapeFunctionArrays(void);

	/* update by-body stored data */
	void SetSurfacesData(void);

	/* sets active striker data (based on current bodies data) */
	void SetActiveStrikers(void); // one contact per striker

protected:
	
	/* search grid */
	iGridManager2DT* fGrid2D; // not a general 2D/3D base class yet

	/* work space */
	dArrayT	fx1, fx2; // facet node coords (shallow)
	dArrayT fStriker; // striker node coords (shallow)
	dArrayT	fv1, fv2; // penetration vectors
	dArrayT	 fNEEvec;
	dMatrixT fNEEmat;
	
	/* derivative arrays */
	dMatrixT fdv1T;
	dMatrixT fdv2T;
	dArrayT  fColtemp1,fColtemp2;
	dMatrixT fdtanT;  	

private:

	/* by-body data */
	ArrayT<dArray2DT> fTanVecs;
	ArrayT<dArrayT>   fTanMags;
};

#endif /* _CONTACT2D_T_H_ */
