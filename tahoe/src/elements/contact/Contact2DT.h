/* $Id: Contact2DT.h,v 1.6 2004-07-15 08:26:08 paklein Exp $ */
/* created: paklein (05/26/1999) */
#ifndef _CONTACT2D_T_H_
#define _CONTACT2D_T_H_

/* base classes */
#include "ContactT.h"

/* direct members */
#include "AutoArrayT.h"

namespace Tahoe {

/* forward declarations */
class iGridManager2DT;

class Contact2DT: public ContactT
{
public:

	/** constructor */
	Contact2DT(const ElementSupportT& support);

	/** destructor */
	virtual ~Contact2DT(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/	
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** \name steps in setting contact configuration */
	/*@{*/
	/** set "internal" data */
	virtual bool SetActiveInteractions(void);

	/** set "external" data to send to FEManager */
	virtual void SetConnectivities(void);
	/*@}*/

	/** \name called by Contact2DT::SetActiveInteractions */
	/*@{*/
	/** update by-body stored data */
	void SetSurfacesData(void);

	/** sets active striker data (based on current bodies data). Produces
	 * one contact per striker */
	void SetActiveStrikers(void);
	/*@}*/

private:

	/* set working arrays */
	void SetShapeFunctionArrays(void);

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

} // namespace Tahoe 
#endif /* _CONTACT2D_T_H_ */
