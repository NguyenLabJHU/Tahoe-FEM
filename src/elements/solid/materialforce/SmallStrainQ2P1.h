/* $Id: SmallStrainQ2P1.h,v 1.2 2003-08-08 22:57:28 thao Exp $ */
#ifndef _SMALL_STRAIN_Q2P1_H_
#define _SMALL_STRAIN_Q2P1_H_

/* base class */
#include "SmallStrainT.h"

namespace Tahoe {

/* forward declarations */
class SSMatSupportT;

/** Interface for linear strain deformation and field gradients */
class SmallStrainQ2P1: public SmallStrainT
{
  public:
      
	/** constructor */
	SmallStrainQ2P1(const ElementSupportT& support, const FieldT& field);

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

	void CalcPressure(void);

	/** calculate the internal force contribution ("-k*d") */
	void FormKd(double constK);

	/** form the element stiffness matrix */
	void FormStiffness(double constK);

  private:

	/** indicies of elements in the list of material needs */
	enum MaterialNeedsT {kstrain = 0,
	                kstrain_last = 1};

 private:
	dArray2DT fTheta_List;  /*Dilation at each ip for each element*/ 
	dArray2DT fTheta_List_last;
	dArray2DT fPressure_List;  /*Dilation at each ip for each element*/ 
	

	const int fpdof;            /*pressure dof*/
       	dArray2DT fMShapes; /*shape function for pressure dof*/
	dMatrixT fH_inv;    /*Mass matrix from fGamma*/

	/*workspace*/
	dArrayT fip_coords;
	dArrayT fGamma;
	dArrayT fAVec;
	dMatrixT fAMat;
	dMatrixT fBMat;

	/*Bbar*/
	dMatrixT fBbar;
	dMatrixT fB_dev;
	dMatrixT fBbar_dil;
	dArrayT fGradTranspose;
	ArrayT<dArray2DT> fGradbar;

	const double fthird;
  	/*@}*/
};

/* inlines */

} // namespace Tahoe 
#endif /* _SMALLSTRAIN_Q2P1_H_ */
