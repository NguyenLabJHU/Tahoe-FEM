/* $Id: BimaterialK_FieldT.h,v 1.6 2004-07-15 08:31:21 paklein Exp $ */
/* created: paklein (09/06/2000)*/

#ifndef _BIMATERIAL_K_FIELD_T_H_
#define _BIMATERIAL_K_FIELD_T_H_

/* base class */
#include "K_FieldT.h"

namespace Tahoe {

/** Displacements for a bimaterial K-field.
 * Displacement field taken from\\
 * P.P.L.Matos et al (1989), Int. J. of Fract. \b 40, 235-254. */
class BimaterialK_FieldT: public K_FieldT
{
public:

	/** constructor */
	BimaterialK_FieldT(const BasicSupportT& support);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* compute K-field displacement factors */
	virtual void ComputeDisplacementFactors(const dArrayT& tip_coords);

private:

	/* bimaterial displacement field factors */
	void SetFieldFactors(int side, double eps, double mu, double G,
		const dArrayT& tip_coords, const iArrayT& nodes, dArray2DT& K1_disp,
		dArray2DT& K2_disp);

	/* group in the "upper half plane" */
	int UpperHalfPlane(void) const;
	
protected:

	/* links to element groups */
	int fFarFieldGroupNum_2;
	int fFarFieldMaterialNum_2;

	/* BC nodes */
	ArrayT<StringT> fID_List_1;
	iArrayT fNodes_1;
	ArrayT<StringT> fID_List_2;
	iArrayT fNodes_2;

	dArray2DT fK1Disp_1;
	dArray2DT fK2Disp_1;
	dArray2DT fK1Disp_2;
	dArray2DT fK2Disp_2;

	/** \name elastic constants */
	/*@{*/
	double fmu_1; /**< shear modulus */
	double fnu_1; /**< Poisson's ratio */
	double fkappa_1; /**< function of nu */

	double fmu_2; /**< shear modulus */
	double fnu_2; /**< Poisson's ratio */
	double fkappa_2; /**< function of nu */
	
	int fGroupNumber_1;
	int fMaterialNumber_1;

	int fGroupNumber_2;
	int fMaterialNumber_2;
	/*@}*/

	/** group in "upper half plane" (t > 0) */
	int fUHP;
};

} // namespace Tahoe 
#endif /* _BIMATERIAL_K_FIELD_T_H_ */
