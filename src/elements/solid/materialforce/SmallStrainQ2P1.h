/* $Id: SmallStrainQ2P1.h,v 1.1 2003-06-28 00:40:29 thao Exp $ */
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

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

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
	

	const int fpdof;                /*pressure dof*/
	ArrayT<dArray2DT> fMShapes; /*shape function for pressure dof*/
	ArrayT<dSymMatrixT> fH_inv; /*Mass matrix from fGamma*/

	/*workspace*/
	dArrayT fip_coords;
	dArrayT fGamma;
	dSymMatrixT fH_i;
	dSymMatrixT fStress;

	dArray2DT fbabar;
	dMatrixT fMixMat;
	dArrayT fGradNa;

	dMatrixT fModulus;
	dMatrixT fPdev;
	dArrayT fOne;

	dArrayT fMEArray;
	dArrayT fNEArray;
	

	const double fthird;
  	/*@}*/
};

/* inlines */

} // namespace Tahoe 
#endif /* _SMALLSTRAIN_Q2P1_H_ */
