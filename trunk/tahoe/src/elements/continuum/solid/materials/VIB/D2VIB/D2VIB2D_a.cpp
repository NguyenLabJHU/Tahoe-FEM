/* $Id: D2VIB2D_a.cpp,v 1.9 2004-07-15 08:27:57 paklein Exp $ */
/* created: paklein (10/23/1999) */
#include "D2VIB2D_a.h"

#include <math.h>
#include <iostream.h>

#include "ifstreamT.h"
#include "D2FSMatSupportT.h"
#include "D2MeshFreeFSSolidT.h"
#include "D2MeshFreeShapeFunctionT.h"

using namespace Tahoe;

/* constructors */
D2VIB2D_a::D2VIB2D_a(ifstreamT& in, const D2FSMatSupportT& support):
	ParameterInterfaceT("gradient_VIB_a_2D"),
	D2VIB2D(in, support),
	fLocDisp(support.D2MeshFreeFDElastic()->Displacements()),
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
