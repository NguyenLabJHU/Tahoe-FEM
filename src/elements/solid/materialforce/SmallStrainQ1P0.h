/* $Id: SmallStrainQ1P0.h,v 1.1 2003-08-22 16:59:20 thao Exp $ */
#ifndef _SMALL_STRAIN_Q1P0_H_
#define _SMALL_STRAIN_Q1P0_H_

/* base class */
#include "SmallStrainT.h"

namespace Tahoe {

/* forward declarations */
class SSMatSupportT;

/** Interface for linear strain deformation and field gradients */
class SmallStrainQ1P0: public SmallStrainT
{
  public:
      
	/** constructor */
	SmallStrainQ1P0(const ElementSupportT& support, const FieldT& field);

	/** initialization. called immediately after constructor */
	virtual void Initialize(void);

	virtual void CloseStep(void);
	virtual void ResetStep(void);
	virtual void ReadRestart(istream& in);
	virtual void WriteRestart(ostream& out) const;

 protected:

	/** set local arrays */
	virtual void SetLocalArrays(void);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/** Set modified shape functions */
	void Set_Gradbar(void);
	
	/** Form Bbar */
	void Set_Bbar(const dArray2DT& DNa, const dArray2DT& DNabar, 
		      dMatrixT& Bbar, dMatrixT& B_dev, dMatrixT& Bbar_dil);

	/** calculate the internal force contribution ("-k*d") */
	void FormKd(double constK);

	/** form the element stiffness matrix */
	void FormStiffness(double constK);

  private:

	/** indicies of elements in the list of material needs */
	enum MaterialNeedsT {kstrain = 0,
	                kstrain_last = 1};
 protected:
 	dArray2DT fGradbar;
	dArrayT fTheta_List;      /*Dilation at each ip for each element*/ 
	dArrayT fTheta_List_last;
	dArrayT fPressure_List;   /*Dilation at each ip for each element*/ 

	const double fthird;

 private:	
	double fElemVol_inv;    /*inverse of element volume*/

	/*Bbar*/
	dMatrixT fBbar;
	dMatrixT fB_dev;
	dMatrixT fBbar_dil;
	dArrayT fGradTranspose;
	
  	/*@}*/
};

/* inlines */

} // namespace Tahoe 
#endif /* _SMALLSTRAIN_Q1P0_H_ */
