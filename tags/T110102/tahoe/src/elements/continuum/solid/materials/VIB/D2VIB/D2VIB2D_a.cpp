/* $Id: D2VIB2D_a.cpp,v 1.5 2002-10-20 22:48:58 paklein Exp $ */
/* created: paklein (10/23/1999)                                          */

#include "D2VIB2D_a.h"

#include <math.h>
#include <iostream.h>

#include "toolboxConstants.h"
#include "fstreamT.h"
#include "D2MeshFreeFDElasticT.h"
#include "D2MeshFreeShapeFunctionT.h"

/* constructors */

using namespace Tahoe;

D2VIB2D_a::D2VIB2D_a(ifstreamT& in, const D2MeshFreeFDElasticT& element):
	D2VIB2D(in, element),
	fLocDisp(element.Displacements()),
	/* work space */
	fPK2(NumSD()),
	fPK2mat(NumSD()),
	fPK1(NumSD()),
	fGradGradU(NumDOF(), dSymMatrixT::NumValues(NumSD()))
{
	/* gradient term coefficient */
	in >> fD2coeff;
	if (fD2coeff < 0.0) throw ExceptionT::kBadInputValue;
}

/* print parameters */
void D2VIB2D_a::Print(ostream& out) const
{
	/* inherited */
	D2VIB2D::Print(out);
	
	out << " Gradient coefficient. . . . . . . . . . . . . . = " << fD2coeff << '\n';
}

/* material internal stress terms */
void D2VIB2D_a::StressTerms(dMatrixT& DW, dMatrixT& DDW)
{
	/* strain */
	Compute_E(fE);

	/* compute the symetric 2nd Piola-Kirchhoff reduced index vector */
	ComputePK2(fE, fPK2);
	fPK2.ToMatrix(fPK2mat);

	/* PK1 */
	DW.MultAB(F(), fPK2mat);

	/* higher order gradient term */
	fD2MLSShape.GradGradU(fLocDisp, fGradGradU);
	DDW.SetToScaled(feps2*fD2coeff, fGradGradU);
}
