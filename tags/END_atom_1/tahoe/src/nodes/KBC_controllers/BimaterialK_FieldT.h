/* $Id: BimaterialK_FieldT.h,v 1.5 2002-07-05 22:28:31 paklein Exp $ */
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

	/* constructor */
	BimaterialK_FieldT(NodeManagerT& node_manager);

	/* initialize data - called immediately after construction */
	virtual void Initialize(ifstreamT& in);
	virtual void WriteParameters(ostream& out) const;

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
	
	/* external links */
	const IsotropicT*  fIsotropic_2;
	const Material2DT* fMaterial2D_2;

	dArray2DT fK1Disp_1;
	dArray2DT fK2Disp_1;
	dArray2DT fK1Disp_2;
	dArray2DT fK2Disp_2;
	
	/* group in "upper half plane" (t > 0) */
	int fUHP;
};

} // namespace Tahoe 
#endif /* _BIMATERIAL_K_FIELD_T_H_ */
