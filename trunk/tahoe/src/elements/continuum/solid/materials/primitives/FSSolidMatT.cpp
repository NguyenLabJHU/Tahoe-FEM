/* $Id: FSSolidMatT.cpp,v 1.21 2006-10-24 00:24:26 tdnguye Exp $ */
/* created: paklein (06/09/1997) */
#include "FSSolidMatT.h"
#include "FSMatSupportT.h"
#include "ThermalDilatationT.h"
#include "iArray2DT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<FSSolidMatT>::fByteCopy = false;
DEFINE_TEMPLATE_STATIC const bool ArrayT<FSSolidMatT*>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
FSSolidMatT::FSSolidMatT(void):
	ParameterInterfaceT("large_strain_material"),
	fFSMatSupport(NULL),
	fTemperatureField(false)
{

}

/* set the material support or pass NULL to clear */
void FSSolidMatT::SetFSMatSupport(const FSMatSupportT* support)
{
	/* set inherited material support */
	SetMaterialSupport(support);

	fFSMatSupport = support;
	
	/* dimension */
	int nsd = NumSD();
	TensorTransformT::Dimension(nsd);
	fQ.Dimension(nsd);
	fF_therm_inv.Dimension(nsd);
	fF_therm_inv_last.Dimension(nsd);
	fF_mechanical.Dimension(nsd);

	/* initialize */
	fF_therm_inv.Identity();
	fF_therm_inv_last.Identity();	
	fF_mechanical.Identity();
}

/* spatial elastic modulus */
const dMatrixT& FSSolidMatT::ce_ijkl(void) {
	return c_ijkl();
}

const dMatrixT& FSSolidMatT::c_ijkl(void)
{
	/* basis vectors */
	int nsd = NumSD();
	double basis_1D[1*1] = {0.0};
	double basis_2D[2*2] = {1.0, 0.0, 0.0, 1.0};
	double basis_3D[3*3] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
	double* basis_list[4] = {NULL, basis_1D, basis_2D, basis_3D};
	dArray2DT basis(nsd, nsd, basis_list[nsd]);

	/* work space */
	dMatrixT& F = const_cast<dMatrixT&>(fFSMatSupport->DeformationGradient());
	F_0_ = F;
	double J_0 = F_0_.Det();
	dArrayT e_c, e_d;

	/* compute perturbed stress (in columns) */
	double eps = 1.0e-08;
	for (int i = 0; i < fModulus.Cols(); i++) {

		/* map column to symmetric indicies */
		int c, d;
		dSymMatrixT::ExpandIndex(nsd, i, c, d);
		basis.RowAlias(c, e_c);
		basis.RowAlias(d, e_d);

		/* perturbed deformation gradient (2.17) */
		F = F_0_;
		F_0_.MultTx(e_d, vec_);
		F.Outer(e_c, vec_, 0.5*eps, dMatrixT::kAccumulate); 
		F_0_.MultTx(e_c, vec_);
		F.Outer(e_d, vec_, 0.5*eps, dMatrixT::kAccumulate); 
		double J = F.Det();

		/* compute stress */
		stress_.SetToScaled(J, s_ij());
	
		/* write into modulus */
		fModulus.SetCol(i, stress_);
	}
	
	/* restore nominal state of deformation and stress */
	F = F_0_;
	stress_.SetToScaled(J_0, s_ij());
	
	/* compute finite difference and geometric contribution (2.18) */
	for (int i = 0; i < fModulus.Cols(); i++) {

		/* map column to symmetric indicies */
		int c, d;
		dSymMatrixT::ExpandIndex(nsd, i, c, d);
		basis.RowAlias(c, e_c);
		basis.RowAlias(d, e_d);

		/* geometric contribution */
		stress_.Multx(e_d, vec_);
		F_0_.Outer(e_c, vec_, 0.5, dMatrixT::kOverwrite); 
		F_0_.Outer(vec_, e_c, 0.5, dMatrixT::kAccumulate); 
		stress_.Multx(e_c, vec_);
		F_0_.Outer(e_d, vec_, 0.5, dMatrixT::kAccumulate); 
		F_0_.Outer(vec_, e_d, 0.5, dMatrixT::kAccumulate); 
		
		/* combine results */
		for (int j = 0; j < fModulus.Rows(); j++)
		{
			int a, b;
			dSymMatrixT::ExpandIndex(nsd, j, a, b);
			fModulus(j,i) = (fModulus(j,i) - stress_[j])/eps - F_0_(a,b);
		}
	}

	/* J factor */
	fModulus /= J_0;

	return fModulus;
}

/* material description */
const dMatrixT& FSSolidMatT::C_IJKL(void)
{
	/* spatial -> material */
	const dMatrixT& Fmat = F(); // NOTE: use F or F_mechanical?
	fModulus.SetToScaled(Fmat.Det(), PullBack(Fmat, c_ijkl()));	
	return fModulus;
}

const dSymMatrixT& FSSolidMatT::S_IJ(void)
{
	/* spatial -> material */
	const dMatrixT& Fmat = F(); // NOTE: use F or F_mechanical?
	fStress.SetToScaled(Fmat.Det(), PullBack(Fmat, s_ij()));	
	return fStress;
}

/* test for localization using "current" values for Cauchy
* stress and the spatial tangent moduli. Returns 1 if the
* determinant of the acoustic tensor is negative and returns
* the normal for which the determinant is minimum. Returns 0
* of the determinant is positive. */
bool FSSolidMatT::IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs, double &DetA)
{
#pragma unused(normals)
#pragma unused(slipdirs)
#pragma unused(DetA)

ExceptionT::GeneralFail("FSSolidMatT::IsLocalized", "broken");
return false;
//DEV
/*
#if 0
	if (FDContinuumT::IsLocalized(normal))
		return 1;
	else
		return 0;
#endif
*/
}

bool FSSolidMatT::IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs, 
                        AutoArrayT <double> &detAs, AutoArrayT <double>
                        &dissipations_fact)
{
#pragma unused(normals)
#pragma unused(slipdirs)
#pragma unused(detAs)
#pragma unused(dissipations_fact)

ExceptionT::GeneralFail("FSSolidMatT::IsLocalized", "broken");
return false;
}

bool FSSolidMatT::IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs)
{
  double dummyDetA = 0.0;
  return IsLocalized(normals, slipdirs, dummyDetA);
}


/* initialize current step. compute thermal dilatation */
void FSSolidMatT::InitStep(void)
{
	/* inherited */
	SolidMaterialT::InitStep();

	/* set multiplicative thermal transformation */
	SetInverseThermalTransformation(fF_therm_inv);
}

/* close current step. store thermal dilatation */
void FSSolidMatT::CloseStep(void)
{
	/* store the thermal dilatation used for the current time step */
	fF_therm_inv_last = fF_therm_inv;
}

/* deformation gradients */
const dMatrixT& FSSolidMatT::F(void) const
{
	return fFSMatSupport->DeformationGradient();
	//DEV - what about corrections for thermal strains?
}

const dMatrixT& FSSolidMatT::F(int ip) const
{
	return fFSMatSupport->DeformationGradient(ip);
	//DEV - what about corrections for thermal strains?
}

const dMatrixT& FSSolidMatT::F_total(void) const
{
	return fFSMatSupport->DeformationGradient();
}

const dMatrixT& FSSolidMatT::F_total(int ip) const
{
	return fFSMatSupport->DeformationGradient(ip);
}

const dMatrixT& FSSolidMatT::F_mechanical(void)
{
	/* has nodal temperature field */
	if (fTemperatureField)
	{
		/* integration point temperature */
		fFSMatSupport->Interpolate(*(fFSMatSupport->Temperatures()), fTemperature);
	
		/* expansion factor */
		double dilatation = 1.0 + fTemperature[0]*fCTE;

		/* remove thermal strain */
		fF_therm_inv.Identity(1.0/dilatation);
		fF_mechanical.MultAB(fFSMatSupport->DeformationGradient(), fF_therm_inv);
		return fF_mechanical;
	}
	else if (fThermal->IsActive())
	{	
		fF_mechanical.MultAB(fFSMatSupport->DeformationGradient(), fF_therm_inv);
		return fF_mechanical;
	}
	else /* no thermal strain */
		return fFSMatSupport->DeformationGradient();
}

const dMatrixT& FSSolidMatT::F_mechanical(int ip)
{
	/* has nodal temperature field */
	if (fTemperatureField)
	{
		/* integration point temperature */
		fFSMatSupport->Interpolate(*(fFSMatSupport->Temperatures()), fTemperature, ip);
	
		/* expansion factor */
		double dilatation = 1.0 + fTemperature[0]*fCTE;

		/* remove thermal strain */
		fF_therm_inv.Identity(1.0/dilatation);
		fF_mechanical.MultAB(fFSMatSupport->DeformationGradient(), fF_therm_inv);
		return fF_mechanical;
	}
	/* has thermal strain */
	else if (fThermal->IsActive())
	{
		fF_mechanical.MultAB(fFSMatSupport->DeformationGradient(ip), fF_therm_inv);
		return fF_mechanical;
	}
	else /* no thermal strain */
		return fFSMatSupport->DeformationGradient(ip);
}

/* deformation gradient from end of previous step */
const dMatrixT& FSSolidMatT::F_total_last(void) const
{
	return fFSMatSupport->DeformationGradient_last();
}

const dMatrixT& FSSolidMatT::F_total_last(int ip) const
{
	return fFSMatSupport->DeformationGradient_last(ip);
}

const dMatrixT& FSSolidMatT::F_mechanical_last(void)
{
	/* has nodal temperature field */
	if (fTemperatureField)
	{
		/* integration point temperature */
		fFSMatSupport->Interpolate(*(fFSMatSupport->LastTemperatures()), fTemperature);
	
		/* expansion factor */
		double dilatation = 1.0 + fTemperature[0]*fCTE;

		/* remove thermal strain */
		fF_therm_inv.Identity(1.0/dilatation);
		fF_mechanical.MultAB(fFSMatSupport->DeformationGradient_last(), fF_therm_inv);
		return fF_mechanical;
	}
	/* has thermal strain */
	else if (fThermal->IsActive())
	{
		fF_mechanical.MultAB(fFSMatSupport->DeformationGradient_last(), 
			fF_therm_inv_last);
		return fF_mechanical;
	}
	else /* no thermal strain */
		return fFSMatSupport->DeformationGradient_last();
}

const dMatrixT& FSSolidMatT::F_mechanical_last(int ip)
{
	/* has nodal temperature field */
	if (fTemperatureField)
	{
		/* integration point temperature */
		fFSMatSupport->Interpolate(*(fFSMatSupport->LastTemperatures()), fTemperature, ip);
	
		/* expansion factor */
		double dilatation = 1.0 + fTemperature[0]*fCTE;

		/* remove thermal strain */
		fF_therm_inv.Identity(1.0/dilatation);
		fF_mechanical.MultAB(fFSMatSupport->DeformationGradient(), fF_therm_inv);
		return fF_mechanical;
	}
	/* has thermal strain */
	else if (fThermal->IsActive())
	{
		fF_mechanical.MultAB(fFSMatSupport->DeformationGradient_last(ip), 
			fF_therm_inv_last);
		return fF_mechanical;
	}
	else /* no thermal strain */
		return fFSMatSupport->DeformationGradient_last(ip);
}

/* accept parameter list */
void FSSolidMatT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SolidMaterialT::TakeParameterList(list);

	/* dimension return values */
	int nsd = NumSD();
	fStress.Dimension(nsd);
	fModulus.Dimension(dSymMatrixT::NumValues(nsd));

	/* FSSolidMatT::c_ijkl work space */
	F_0_.Dimension(nsd);
	vec_.Dimension(nsd);
	stress_.Dimension(nsd);

	/* set multiplicative thermal transformation */
	SetInverseThermalTransformation(fF_therm_inv);
	fF_therm_inv_last = fF_therm_inv;
	
	/* check for temperature field */
	if (fFSMatSupport->Temperatures() && fFSMatSupport->LastTemperatures()) {
		/* set flag */
		fTemperatureField = true;
		
		/* work space */
		fTemperature.Dimension(fFSMatSupport->Temperatures()->MinorDim());
		
		/* disable prescribed dilatation */
		fThermal->SetSchedule(NULL);
	}	
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* left stretch tensor */
void FSSolidMatT::Compute_b(dSymMatrixT& b) const
{
	Compute_b(fFSMatSupport->DeformationGradient(), b);
}

/* right stretch tensor */
void FSSolidMatT::Compute_C(dSymMatrixT& C) const
{
	Compute_C(fFSMatSupport->DeformationGradient(), C);
}

/* Green-Lagrangian strain */
void FSSolidMatT::Compute_E(dSymMatrixT& E) const
{
	Compute_E(fFSMatSupport->DeformationGradient(), E);
}

/* left stretch tensor */
void FSSolidMatT::Compute_b(const dMatrixT& F, dSymMatrixT& b) const
{
	const char caller[] = "FSSolidMatT::Compute_b";
	int nsd = F.Rows();
#if __option(extended_errorcheck)
	if (F.Cols() != nsd || b.Rows() != nsd) ExceptionT::SizeMismatch(caller);	
#endif
	if (nsd == 2)
	{
		const double* f = F.Pointer();
		double* a = b.Pointer();
		
		/* unrolled */
		a[0] = f[0]*f[0] + f[2]*f[2];
		a[1] = f[1]*f[1] + f[3]*f[3];
		a[2] = f[0]*f[1] + f[2]*f[3];
	}
	else if (nsd == 3)
	{
		const double* f = F.Pointer();
		double* a = b.Pointer();
	
		/* unrolled */
		a[0] = f[0]*f[0] + f[3]*f[3] + f[6]*f[6];
		a[1] = f[1]*f[1] + f[4]*f[4] + f[7]*f[7];
		a[2] = f[2]*f[2] + f[5]*f[5] + f[8]*f[8];
		a[3] = f[1]*f[2] + f[4]*f[5] + f[7]*f[8];
		a[4] = f[0]*f[2] + f[3]*f[5] + f[6]*f[8];
		a[5] = f[0]*f[1] + f[3]*f[4] + f[6]*f[7];
	}
	else
		ExceptionT::GeneralFail(caller, "unsupported dimension %d", nsd);
}

/* right stretch tensor */
void FSSolidMatT::Compute_C(const dMatrixT& F, dSymMatrixT& C) const
{
	const char caller[] = "FSSolidMatT::Compute_C";
	int nsd = F.Rows();
#if __option(extended_errorcheck)
	if (F.Cols() != nsd || C.Rows() != nsd) ExceptionT::SizeMismatch(caller);	
#endif
	if (nsd == 2)
	{
		const double* f = F.Pointer();
		double* c = C.Pointer();
		
		/* unrolled */
		c[0] = f[0]*f[0] + f[1]*f[1];
		c[1] = f[2]*f[2] + f[3]*f[3];
		c[2] = f[0]*f[2] + f[1]*f[3];
	}
	else if (nsd == 3)
	{
		const double* f = F.Pointer();
		double* c = C.Pointer();
	
		/* unrolled */
		c[0] = f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
		c[1] = f[3]*f[3] + f[4]*f[4] + f[5]*f[5];
		c[2] = f[6]*f[6] + f[7]*f[7] + f[8]*f[8];
		c[3] = f[3]*f[6] + f[4]*f[7] + f[5]*f[8];
		c[4] = f[0]*f[6] + f[1]*f[7] + f[2]*f[8];
		c[5] = f[0]*f[3] + f[1]*f[4] + f[2]*f[5];
	}
	else
		ExceptionT::GeneralFail(caller, "unsupported dimension %d", nsd);
}

/* Green-Lagrangian strain */
void FSSolidMatT::Compute_E(const dMatrixT& F, dSymMatrixT& E) const
{
	const char caller[] = "FSSolidMatT::Compute_E";
	int nsd = F.Rows();
#if __option(extended_errorcheck)
	if (F.Cols() != nsd || E.Rows() != nsd) ExceptionT::SizeMismatch(caller);
#endif
	if (nsd == 2)
	{
		const double* f = F.Pointer();
		double* e = E.Pointer();
		
		/* unrolled */
		e[0] = (f[0]*f[0] + f[1]*f[1] - 1.0)*0.5;
		e[1] = (f[2]*f[2] + f[3]*f[3] - 1.0)*0.5;
		e[2] = (f[0]*f[2] + f[1]*f[3])*0.5;
	}
	else if (nsd == 3)
	{
		const double* f = F.Pointer();
		double* e = E.Pointer();
	
		/* unrolled */
		e[0] = (f[0]*f[0] + f[1]*f[1] + f[2]*f[2] - 1.0)*0.5;
		e[1] = (f[3]*f[3] + f[4]*f[4] + f[5]*f[5] - 1.0)*0.5;
		e[2] = (f[6]*f[6] + f[7]*f[7] + f[8]*f[8] - 1.0)*0.5;
		e[3] = (f[3]*f[6] + f[4]*f[7] + f[5]*f[8])*0.5;
		e[4] = (f[0]*f[6] + f[1]*f[7] + f[2]*f[8])*0.5;
		e[5] = (f[0]*f[3] + f[1]*f[4] + f[2]*f[5])*0.5;
	}
	else if (nsd == 1)
		E[0] = 0.5*(F[0]*F[0] - 1.0);
	else
		ExceptionT::GeneralFail(caller, "unsupported dimension %d", nsd);
}

/* return the acoustical tensor and wave speeds */
const dSymMatrixT& FSSolidMatT::AcousticalTensor(const dArrayT& normal)
{
	const char caller[] = "FSSolidMatT::AcousticalTensor";
#if __option(extended_errorcheck)
	if (fQ.Rows() != normal.Length()) ExceptionT::SizeMismatch(caller);
#endif

	/* collect matrices */
	/* material description */
	const dMatrixT&    C_ = C_IJKL(); // material tangent moduli
	const dSymMatrixT& S_ = S_IJ();   // 2nd Piola-Kirchhoff stress
	const dMatrixT&    F_ = F();      // (current) deformation gradient
	
	/* dispatch */
	if (normal.Length() == 2)
		ComputeQ_2D(C_, S_, F_, normal, fQ);
	else if (normal.Length() == 3)
		ComputeQ_3D(C_, S_, F_, normal, fQ);
	else
		ExceptionT::GeneralFail(caller);

	return fQ;
}

/*************************************************************************
 * Private
 *************************************************************************/

/* set inverse of thermal transformation - return true if active */
bool FSSolidMatT::SetInverseThermalTransformation(dMatrixT& F_trans_inv)
{
	if (fThermal->IsActive())
	{
		/* assuming isotropic expansion */
		double Fii_inv = 1.0/(1.0 + fThermal->PercentElongation());
		F_trans_inv.Identity(Fii_inv);
		return true;
	}
	else
	{
		F_trans_inv.Identity(1.0);
		return false;
	}
}

/* acoustical tensor routines */
void FSSolidMatT::ComputeQ_2D(const dMatrixT& CIJKL, const dSymMatrixT& SIJ,
	const dMatrixT& FkK, const dArrayT& N, dSymMatrixT& Q) const
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24;
	double z25, z26, z27, z28, z29, z30, z31, z32, z33, z34, z35, z36;
	double z37, z38, z39;

	/* assign temp names */
	double C11 = CIJKL[0];
	double C22 = CIJKL[4];
	double C33 = CIJKL[8];
	double C23 = CIJKL[7];
	double C13 = CIJKL[6];
	double C12 = CIJKL[3];

	double S11 = SIJ[0];
	double S22 = SIJ[1];
	double S12 = SIJ[2];

	double F11 = FkK[0];
	double F21 = FkK[1];
	double F12 = FkK[2];
	double F22 = FkK[3];

	double N1 = N[0];
	double N2 = N[1];

	z1 = F11*F11;
	z2 = F12*F12;
	z3 = F21*F21;
	z4 = F22*F22;
	z5 = N1*N1;
	z6 = N1*N2;
	z7 = N2*N2;
	z8 = 2.*z6;
	z9 = C12*z6;
	z10 = C33*z6;
	z11 = F11*z10;
	z10 = F12*z10;
	z11 = F22*z11;
	z10 = F21*z10;
	z12 = F11*z6;
	z6 = F12*z6;
	z6 = C12*z8;
	z12 = F11*F12*z6;
	z6 = F21*F22*z6;
	z13 = C33*z8;
	z14 = F11*F12*z13;
	z13 = F21*F22*z13;
	z15 = F11*z8;
	z15 = C13*F21*z15;
	z16 = F12*z8;
	z16 = C23*F22*z16;
	z17 = F21*z8;
	z17 = S12*z8;
	z18 = F11*z9;
	z18 = F22*z18;
	z9 = F12*z9;
	z9 = F21*z9;
	z19 = 2.*F11*F12;
	z20 = F11*F21;
	z21 = F12*F21;
	z22 = F11*F22;
	z23 = F12*F22;
	z24 = 2.*F21*F22;
	z25 = C11*z5;
	z26 = C13*z5;
	z27 = C33*z5;
	z5 = S11*z5;
	z28 = z1*z25;
	z29 = z25*z3;
	z25 = z20*z25;
	z30 = z19*z26;
	z31 = z21*z26;
	z32 = z22*z26;
	z26 = z24*z26;
	z33 = z2*z27;
	z34 = z27*z4;
	z27 = z23*z27;
	z35 = C22*z7;
	z36 = C23*z7;
	z37 = C33*z7;
	z7 = S22*z7;
	z38 = z2*z35;
	z39 = z35*z4;
	z23 = z23*z35;
	z19 = z19*z36;
	z21 = z21*z36;
	z22 = z22*z36;
	z24 = z24*z36;
	z35 = z1*z37;
	z36 = z3*z37;
	z20 = z20*z37;
	z37 = C13*z8;
	z8 = C23*z8;
	z1 = z1*z37;
	z3 = z3*z37;
	z2 = z2*z8;
	z4 = z4*z8;
	z5 = z17 + z5 + z7;
	z7 = z10 + z11 + z15 + z16 + z18 + z25 + z27 + z31 + z32 + z9;
	z7 = z20 + z21 + z22 + z23 + z7;
	z2 = z12 + z14 + z19 + z2 + z28 + z30 + z33 + z35 + z38 + z5;
	z1 = z1 + z2;
	z2 = z13 + z24 + z26 + z29 + z34 + z36 + z39 + z4 + z5 + z6;
	z2 = z2 + z3;

	// {{z1, z7},
	//  {z7, z2}}

	Q[0] = z1;
	Q[1] = z2;
	Q[2] = z7;
}

void FSSolidMatT::ComputeQ_3D(const dMatrixT& CIJKL, const dSymMatrixT& SIJ,
	const dMatrixT& FkK, const dArrayT& N, dSymMatrixT& Q) const
{	
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24;
	double z25, z26, z27, z28, z29, z30, z31, z32, z33, z34, z35, z36;
	double z37, z38, z39, z40, z41, z42, z43, z44, z45, z46, z47, z48;
	double z49, z50, z51, z52, z53, z54, z55, z56, z57, z58, z59, z60;
	double z61, z62, z63, z64, z65, z66, z67, z68, z69, z70, z71, z72;
	double z73, z74, z75, z76, z77, z78, z79, z80, z81, z82, z83, z84;
	double z85, z86, z87, z88, z89, z90, z91, z92, z93, z94, z95, z96;
	double z97, z98, z99, z100, z101, z102, z103, z104, z105, z106, z107, z108;
	double z109, z110, z111, z112, z113, z114, z115, z116, z117, z118, z119, z120;
	double z121, z122, z123, z124, z125, z126, z127, z128, z129, z130, z131, z132;
	double z133, z134, z135, z136, z137, z138, z139, z140, z141, z142, z143, z144;
	double z145, z146, z147, z148, z149, z150, z151, z152, z153, z154, z155, z156;
	double z157, z158, z159, z160, z161, z162, z163, z164, z165, z166, z167, z168;
	double z169, z170, z171, z172, z173, z174, z175, z176, z177, z178, z179, z180;
	double z181, z182, z183, z184, z185, z186, z187, z188, z189, z190, z191, z192;
	double z193, z194, z195, z196, z197, z198, z199, z200, z201, z202, z203, z204;
	double z205, z206, z207, z208, z209, z210, z211, z212, z213, z214, z215, z216;
	double z217, z218, z219, z220, z221, z222, z223, z224, z225, z226, z227, z228;
	double z229, z230, z231, z232, z233, z234, z235, z236, z237, z238, z239, z240;
	double z241, z242, z243, z244, z245, z246, z247, z248, z249, z250, z251, z252;
	double z253, z254, z255, z256, z257, z258, z259, z260, z261, z262, z263, z264;
	double z265, z266, z267, z268, z269, z270, z271, z272, z273, z274, z275, z276;
	double z277, z278, z279, z280, z281, z282, z283, z284, z285, z286, z287, z288;
	double z289, z290, z291, z292, z293, z294, z295, z296, z297, z298, z299, z300;
	double z301, z302, z303, z304, z305, z306, z307, z308, z309, z310, z311, z312;
	double z313, z314, z315, z316, z317, z318, z319, z320, z321, z322, z323, z324;
	double z325, z326, z327, z328, z329, z330, z331, z332, z333, z334, z335, z336;
	double z337, z338, z339, z340, z341, z342, z343, z344, z345, z346, z347, z348;
	double z349, z350, z351, z352, z353, z354, z355, z356, z357, z358, z359, z360, z361;

	/* assign temp names */
	double C11 = CIJKL[0];
	double C12 = CIJKL[1];
	double C13 = CIJKL[2];
	double C14 = CIJKL[3];
	double C15 = CIJKL[4];
	double C16 = CIJKL[5];
	double C22 = CIJKL[7];
	double C23 = CIJKL[8];
	double C24 = CIJKL[9];
	double C25 = CIJKL[10];
	double C26 = CIJKL[11];
	double C33 = CIJKL[14];
	double C34 = CIJKL[15];
	double C35 = CIJKL[16];
	double C36 = CIJKL[17];
	double C44 = CIJKL[21];
	double C45 = CIJKL[22];
	double C46 = CIJKL[23];
	double C55 = CIJKL[28];
	double C56 = CIJKL[29];
	double C66 = CIJKL[35];

	double S11 = SIJ[0];
	double S22 = SIJ[1];
	double S33 = SIJ[2];
	double S23 = SIJ[3];
	double S13 = SIJ[4];
	double S12 = SIJ[5];

	double F11 = FkK[0];
	double F21 = FkK[1];
	double F31 = FkK[2];
	double F12 = FkK[3];
	double F22 = FkK[4];
	double F32 = FkK[5];
	double F13 = FkK[6];
	double F23 = FkK[7];
	double F33 = FkK[8];

	double N1 = N[0];
	double N2 = N[1];
	double N3 = N[2];


	z1 = F11*F11;
	z2 = F12*F12;
	z3 = F13*F13;
	z4 = F21*F21;
	z5 = F22*F22;
	z6 = F23*F23;
	z7 = F31*F31;
	z8 = F32*F32;
	z9 = F33*F33;
	z10 = 2.*N1;
	z11 = C14*N1;
	z12 = C45*N1;
	z13 = C46*N1;
	z14 = C56*N1;
	z15 = F11*N1;
	z16 = F12*N1;
	z17 = F13*N1;
	z18 = F21*N1;
	z19 = F22*N1;
	z20 = F23*N1;
	z21 = F31*N1;
	z22 = F32*N1;
	z23 = F33*N1;
	z24 = N1*N1;
	z25 = 2.*N2;
	z26 = C25*N2;
	z27 = C45*N2;
	z28 = C46*N2;
	z29 = C56*N2;
	z30 = F11*N2;
	z31 = F12*N2;
	z32 = F13*N2;
	z33 = F21*N2;
	z34 = F22*N2;
	z35 = F23*N2;
	z36 = F31*N2;
	z37 = F32*N2;
	z38 = F33*N2;
	z39 = N1*N2;
	z39 = N2*N2;
	z40 = 2.*N3;
	z41 = C36*N3;
	z42 = C45*N3;
	z43 = C46*N3;
	z44 = C56*N3;
	z45 = F11*N3;
	z46 = F12*N3;
	z47 = F13*N3;
	z48 = F21*N3;
	z49 = F22*N3;
	z50 = F23*N3;
	z51 = F31*N3;
	z52 = F32*N3;
	z53 = F33*N3;
	z54 = N1*N3;
	z54 = N2*N3;
	z54 = N3*N3;
	z55 = C14*z10;
	z56 = C45*z10;
	z56 = C46*z10;
	z56 = C56*z10;
	z56 = F11*z10;
	z57 = F12*z10;
	z58 = F13*z10;
	z59 = F21*z10;
	z60 = F22*z10;
	z61 = F23*z10;
	z62 = F31*z10;
	z63 = F32*z10;
	z64 = F33*z10;
	z65 = N2*z10;
	z66 = N3*z10;
	z67 = F21*z11;
	z68 = F31*z11;
	z69 = F23*z12;
	z70 = F33*z12;
	z71 = F22*z13;
	z72 = F32*z13;
	z73 = F21*z14;
	z74 = F31*z14;
	z75 = N2*z15;
	z75 = N3*z15;
	z75 = F23*z16;
	z76 = F33*z16;
	z77 = N2*z16;
	z77 = F22*z17;
	z78 = F32*z17;
	z79 = N3*z17;
	z79 = N2*z18;
	z79 = N3*z18;
	z79 = F33*z19;
	z80 = N2*z19;
	z80 = F32*z20;
	z81 = N3*z20;
	z81 = N2*z21;
	z21 = N3*z21;
	z21 = N2*z22;
	z21 = N3*z23;
	z21 = C25*z25;
	z22 = F12*z21;
	z23 = F22*z21;
	z21 = F32*z21;
	z81 = C45*z25;
	z81 = C46*z25;
	z81 = C56*z25;
	z81 = F11*z25;
	z82 = F12*z81;
	z83 = F13*z81;
	z84 = F21*z81;
	z81 = F31*z81;
	z85 = F12*z25;
	z85 = F13*z25;
	z86 = F21*z25;
	z87 = F22*z86;
	z88 = F23*z86;
	z86 = F31*z86;
	z89 = F22*z25;
	z90 = F23*z25;
	z91 = F31*z25;
	z92 = F32*z91;
	z91 = F33*z91;
	z93 = F32*z25;
	z94 = F33*z25;
	z95 = N3*z25;
	z96 = S23*z95;
	z97 = F22*z26;
	z98 = F32*z26;
	z99 = z26*z76;
	z100 = z26*z77;
	z101 = z26*z78;
	z102 = z26*z79;
	z103 = z26*z80;
	z104 = F23*z27;
	z105 = F33*z27;
	z106 = F22*z28;
	z107 = F32*z28;
	z108 = F23*z30;
	z109 = z108*z11;
	z110 = z108*z14;
	z111 = F33*z30;
	z112 = z11*z111;
	z113 = z111*z14;
	z114 = z10*z30;
	z114 = F23*z31;
	z114 = z114*z13;
	z115 = F33*z31;
	z115 = z115*z13;
	z116 = N3*z31;
	z116 = z10*z31;
	z116 = F21*z32;
	z117 = F31*z32;
	z118 = N3*z32;
	z118 = F33*z33;
	z119 = z11*z118;
	z120 = z118*z14;
	z121 = z10*z33;
	z121 = z16*z33;
	z122 = C12*z121;
	z121 = C66*z121;
	z123 = F33*z34;
	z13 = z123*z13;
	z123 = N3*z34;
	z123 = z10*z34;
	z123 = z15*z34;
	z124 = C12*z123;
	z123 = C66*z123;
	z125 = F31*z35;
	z126 = N3*z35;
	z126 = z10*z36;
	z16 = z16*z36;
	z126 = C12*z16;
	z16 = C66*z16;
	z19 = z19*z36;
	z127 = C12*z19;
	z19 = C66*z19;
	z128 = N3*z37;
	z128 = z10*z37;
	z128 = z15*z37;
	z129 = C12*z128;
	z128 = C66*z128;
	z37 = z18*z37;
	z130 = C12*z37;
	z37 = C66*z37;
	z38 = N3*z38;
	z38 = C36*z40;
	z38 = C45*z40;
	z38 = C46*z40;
	z38 = C56*z40;
	z38 = F11*z40;
	z38 = F12*z40;
	z38 = F13*z40;
	z38 = F21*z40;
	z38 = F22*z40;
	z38 = F23*z40;
	z38 = F31*z40;
	z38 = F32*z40;
	z38 = F13*z41;
	z40 = F23*z41;
	z131 = F33*z41;
	z76 = z41*z76;
	z77 = z41*z77;
	z78 = z41*z78;
	z79 = z41*z79;
	z80 = z41*z80;
	z132 = z41*z83;
	z133 = z41*z88;
	z134 = z41*z91;
	z108 = z108*z41;
	z111 = z111*z41;
	z116 = z116*z41;
	z117 = z117*z41;
	z118 = z118*z41;
	z125 = z125*z41;
	z135 = F13*z42;
	z136 = F23*z42;
	z137 = F33*z42;
	z83 = z42*z83;
	z88 = z42*z88;
	z42 = z42*z91;
	z91 = F22*z43;
	z138 = F32*z43;
	z82 = z43*z82;
	z87 = z43*z87;
	z92 = z43*z92;
	z84 = z44*z84;
	z81 = z44*z81;
	z86 = z44*z86;
	z139 = F22*z45;
	z140 = z11*z139;
	z139 = z139*z14;
	z141 = F32*z45;
	z142 = z11*z141;
	z141 = z14*z141;
	z143 = z10*z45;
	z22 = z22*z45;
	z97 = z45*z97;
	z143 = z45*z98;
	z144 = z104*z45;
	z145 = z105*z45;
	z106 = z106*z45;
	z146 = z107*z45;
	z147 = F21*z46;
	z148 = z147*z26;
	z147 = z147*z28;
	z149 = F31*z46;
	z150 = z149*z26;
	z149 = z149*z28;
	z85 = z46*z85;
	z89 = C24*z46*z89;
	z151 = C24*z46*z93;
	z152 = z25*z46;
	z152 = C23*z85;
	z85 = C44*z85;
	z153 = F21*z47;
	z154 = F22*z47;
	z155 = F31*z47;
	z156 = F32*z47;
	z157 = z10*z47;
	z157 = C34*z47*z90;
	z158 = C34*z47*z94;
	z159 = z25*z47;
	z153 = z153*z27;
	z154 = z12*z154;
	z155 = z155*z27;
	z156 = z12*z156;
	z159 = F32*z48;
	z160 = z10*z48;
	z23 = z23*z48;
	z160 = z17*z48;
	z98 = z48*z98;
	z161 = z105*z48;
	z107 = z107*z48;
	z162 = F31*z49;
	z90 = z49*z90;
	z93 = C24*z49*z93;
	z163 = z25*z49;
	z163 = z32*z49;
	z11 = z11*z159;
	z14 = z14*z159;
	z159 = C13*z160;
	z160 = C55*z160;
	z164 = F31*z50;
	z165 = F32*z50;
	z166 = z10*z50;
	z166 = z15*z50;
	z167 = C34*z50*z94;
	z168 = z25*z50;
	z168 = z31*z50;
	z169 = z10*z51;
	z21 = z21*z51;
	z17 = z17*z51;
	z20 = z20*z51;
	z169 = z162*z26;
	z162 = z162*z28;
	z170 = C23*z90;
	z90 = C44*z90;
	z94 = z52*z94;
	z171 = z25*z52;
	z171 = z32*z52;
	z52 = z35*z52;
	z172 = C23*z163;
	z163 = C44*z163;
	z173 = z10*z53;
	z15 = z15*z53;
	z18 = z18*z53;
	z173 = z25*z53;
	z173 = z31*z53;
	z53 = z34*z53;
	z164 = z164*z27;
	z12 = z12*z165;
	z165 = C13*z166;
	z166 = C55*z166;
	z174 = C23*z168;
	z168 = C44*z168;
	z175 = C13*z17;
	z17 = C55*z17;
	z176 = C13*z20;
	z20 = C55*z20;
	z177 = F13*z30*z55;
	z178 = F23*z33*z55;
	z179 = F33*z36*z55;
	z180 = F12*z45*z55;
	z181 = F22*z48*z55;
	z55 = F32*z51*z55;
	z182 = C23*z94;
	z94 = C44*z94;
	z183 = C23*z171;
	z171 = C44*z171;
	z184 = C23*z52;
	z52 = C44*z52;
	z185 = C13*z15;
	z15 = C55*z15;
	z186 = C13*z18;
	z18 = C55*z18;
	z187 = C23*z173;
	z173 = C44*z173;
	z188 = C23*z53;
	z53 = C44*z53;
	z189 = F12*z56;
	z56 = F13*z56;
	z190 = F13*z57;
	z191 = z30*z57;
	z38 = z38*z57;
	z135 = z135*z57;
	z91 = z57*z91;
	z57 = z138*z57;
	z104 = z104*z58;
	z192 = z105*z58;
	z58 = z45*z58;
	z193 = F22*z59;
	z194 = F23*z59;
	z195 = C16*z30*z59;
	z59 = C15*z45*z59;
	z189 = z189*z44;
	z56 = z29*z56;
	z196 = F23*z60;
	z197 = C26*z31*z60;
	z198 = z33*z60;
	z40 = z40*z60;
	z136 = z136*z60;
	z60 = z138*z60;
	z138 = z190*z26;
	z190 = z190*z28;
	z199 = C12*z191;
	z191 = C66*z191;
	z105 = z105*z61;
	z200 = C35*z47*z61;
	z61 = z48*z61;
	z201 = C13*z58;
	z58 = C55*z58;
	z202 = F32*z62;
	z203 = F33*z62;
	z30 = C16*z30*z62;
	z33 = C16*z33*z62;
	z45 = C15*z45*z62;
	z48 = C15*z48*z62;
	z62 = z193*z44;
	z193 = z194*z29;
	z194 = z196*z26;
	z196 = z196*z28;
	z204 = F33*z63;
	z31 = C26*z31*z63;
	z34 = C26*z34*z63;
	z36 = z36*z63;
	z131 = z131*z63;
	z63 = z137*z63;
	z137 = C12*z198;
	z198 = C66*z198;
	z47 = C35*z47*z64;
	z50 = C35*z50*z64;
	z51 = z51*z64;
	z64 = C13*z61;
	z61 = C55*z61;
	z205 = S12*z65;
	z202 = z202*z44;
	z29 = z203*z29;
	z203 = S13*z66;
	z206 = z204*z26;
	z28 = z204*z28;
	z204 = C12*z36;
	z36 = C66*z36;
	z207 = C13*z51;
	z51 = C55*z51;
	z208 = z32*z67;
	z67 = z46*z67;
	z209 = z32*z68;
	z210 = z35*z68;
	z211 = z46*z68;
	z68 = z49*z68;
	z69 = z46*z69;
	z212 = z46*z70;
	z70 = z49*z70;
	z71 = z32*z71;
	z213 = z32*z72;
	z72 = z35*z72;
	z214 = z32*z73;
	z73 = z46*z73;
	z32 = z32*z74;
	z35 = z35*z74;
	z46 = z46*z74;
	z49 = z49*z74;
	z26 = z26*z75;
	z41 = z41*z75;
	z74 = 2.*F11;
	z75 = 2.*F12;
	z215 = 2.*F13;
	z215 = 2.*F21;
	z216 = F11*F21;
	z217 = F12*F21;
	z218 = F13*F21;
	z219 = 2.*F22;
	z220 = F11*F22;
	z221 = F12*F22;
	z222 = F13*F22;
	z223 = 2.*F23;
	z223 = F11*F23;
	z224 = F12*F23;
	z225 = F13*F23;
	z226 = 2.*F31;
	z227 = F11*F31;
	z228 = F12*F31;
	z229 = F13*F31;
	z230 = F21*F31;
	z231 = F22*F31;
	z232 = F23*F31;
	z233 = 2.*F32;
	z234 = F11*F32;
	z235 = F12*F32;
	z236 = F13*F32;
	z237 = F21*F32;
	z238 = F22*F32;
	z239 = F23*F32;
	z240 = 2.*F33;
	z240 = F11*F33;
	z241 = F12*F33;
	z242 = F13*F33;
	z243 = F21*F33;
	z244 = F22*F33;
	z245 = F23*F33;
	z246 = C24*z95;
	z95 = C34*z95;
	z247 = 2.*z24;
	z247 = C11*z24;
	z248 = C15*z24;
	z249 = C16*z24;
	z250 = C55*z24;
	z251 = C56*z24;
	z252 = C66*z24;
	z24 = S11*z24;
	z27 = z10*z27;
	z253 = 2.*z39;
	z253 = C22*z39;
	z254 = C24*z39;
	z255 = C26*z39;
	z256 = C44*z39;
	z257 = C46*z39;
	z258 = C66*z39;
	z39 = S22*z39;
	z10 = z10*z43;
	z25 = z25*z44;
	z43 = 2.*z54;
	z43 = C33*z54;
	z44 = C34*z54;
	z259 = C35*z54;
	z260 = C44*z54;
	z261 = C45*z54;
	z262 = C55*z54;
	z54 = S33*z54;
	z263 = C16*z65;
	z65 = C26*z65;
	z264 = C15*z66;
	z66 = C35*z66;
	z265 = F12*z74;
	z74 = F13*z74;
	z75 = F13*z75;
	z266 = F22*z215;
	z215 = F23*z215;
	z219 = F23*z219;
	z267 = F32*z226;
	z226 = F33*z226;
	z233 = F33*z233;
	z268 = z2*z246;
	z269 = z246*z5;
	z246 = z246*z8;
	z270 = z3*z95;
	z271 = z6*z95;
	z272 = z1*z247;
	z273 = z247*z4;
	z274 = z247*z7;
	z275 = z216*z247;
	z276 = z227*z247;
	z247 = z230*z247;
	z277 = z218*z248;
	z278 = z223*z248;
	z279 = z229*z248;
	z280 = z232*z248;
	z281 = z240*z248;
	z282 = z243*z248;
	z283 = z217*z249;
	z284 = z220*z249;
	z285 = z228*z249;
	z286 = z231*z249;
	z287 = z234*z249;
	z288 = z237*z249;
	z289 = z250*z3;
	z290 = z250*z6;
	z291 = z225*z250;
	z292 = z242*z250;
	z293 = z245*z250;
	z294 = z222*z251;
	z295 = z224*z251;
	z296 = z236*z251;
	z297 = z239*z251;
	z298 = z241*z251;
	z299 = z244*z251;
	z300 = z2*z252;
	z301 = z252*z5;
	z302 = z221*z252;
	z303 = z252*z8;
	z304 = z235*z252;
	z252 = z238*z252;
	z305 = z27*z3;
	z306 = z27*z6;
	z307 = z2*z253;
	z308 = z253*z5;
	z309 = z221*z253;
	z310 = z253*z8;
	z311 = z235*z253;
	z253 = z238*z253;
	z312 = z222*z254;
	z313 = z224*z254;
	z314 = z236*z254;
	z315 = z239*z254;
	z316 = z241*z254;
	z317 = z244*z254;
	z318 = z217*z255;
	z319 = z220*z255;
	z320 = z228*z255;
	z321 = z231*z255;
	z322 = z234*z255;
	z323 = z237*z255;
	z324 = z256*z3;
	z325 = z256*z6;
	z326 = z225*z256;
	z327 = z242*z256;
	z328 = z245*z256;
	z329 = z218*z257;
	z330 = z223*z257;
	z331 = z229*z257;
	z332 = z232*z257;
	z333 = z240*z257;
	z334 = z243*z257;
	z335 = z1*z258;
	z336 = z258*z4;
	z337 = z258*z7;
	z338 = z216*z258;
	z339 = z227*z258;
	z258 = z230*z258;
	z340 = z10*z2;
	z341 = z10*z5;
	z10 = z10*z8;
	z342 = z1*z25;
	z343 = z25*z4;
	z25 = z25*z7;
	z344 = z3*z43;
	z345 = z43*z6;
	z225 = z225*z43;
	z242 = z242*z43;
	z245 = z245*z43;
	z222 = z222*z44;
	z224 = z224*z44;
	z236 = z236*z44;
	z239 = z239*z44;
	z241 = z241*z44;
	z244 = z244*z44;
	z218 = z218*z259;
	z223 = z223*z259;
	z229 = z229*z259;
	z232 = z232*z259;
	z240 = z240*z259;
	z243 = z243*z259;
	z346 = z2*z260;
	z347 = z260*z5;
	z221 = z221*z260;
	z348 = z260*z8;
	z235 = z235*z260;
	z238 = z238*z260;
	z217 = z217*z261;
	z220 = z220*z261;
	z228 = z228*z261;
	z231 = z231*z261;
	z234 = z234*z261;
	z237 = z237*z261;
	z260 = z1*z262;
	z349 = z262*z4;
	z350 = z262*z7;
	z216 = z216*z262;
	z227 = z227*z262;
	z230 = z230*z262;
	z262 = z1*z263;
	z351 = z263*z4;
	z263 = z263*z7;
	z2 = z2*z65;
	z5 = z5*z65;
	z8 = z65*z8;
	z1 = z1*z264;
	z4 = z264*z4;
	z7 = z264*z7;
	z3 = z3*z66;
	z6 = z6*z66;
	z65 = z249*z265;
	z264 = z255*z265;
	z265 = z261*z265;
	z352 = z248*z74;
	z353 = z257*z74;
	z74 = z259*z74;
	z354 = z251*z75;
	z355 = z254*z75;
	z75 = z44*z75;
	z356 = z249*z266;
	z357 = z255*z266;
	z266 = z261*z266;
	z358 = z215*z248;
	z359 = z215*z257;
	z215 = z215*z259;
	z360 = z219*z251;
	z361 = z219*z254;
	z219 = z219*z44;
	z249 = z249*z267;
	z255 = z255*z267;
	z261 = z261*z267;
	z248 = z226*z248;
	z257 = z226*z257;
	z226 = z226*z259;
	z95 = z9*z95;
	z250 = z250*z9;
	z27 = z27*z9;
	z256 = z256*z9;
	z43 = z43*z9;
	z9 = z66*z9;
	z66 = z233*z251;
	z251 = z233*z254;
	z44 = z233*z44;
	z24 = z203 + z205 + z24 + z39 + z54 + z96;
	z1 = z1 + z2 + z260 + z262 + z324 + z335 + z340 + z342 + z344 + z346;
	z2 = z264 + z265 + z3 + z352 + z353 + z354 + z355 + z65 + z74 + z75;
	z1 = z1 + z132 + z152 + z177 + z2 + z22 + z24 + z82 + z83 + z85;
	z1 = z1 + z135 + z138 + z180 + z189 + z190 + z191 + z199 + z38 + z56;
	z1 = z1 + z201 + z268 + z270 + z272 + z289 + z300 + z305 + z307 + z58;
	z2 = z325 + z336 + z341 + z343 + z345 + z347 + z349 + z351 + z4 + z5;
	z3 = z215 + z219 + z266 + z356 + z357 + z358 + z359 + z360 + z361 + z6;
	z2 = z133 + z170 + z178 + z2 + z23 + z24 + z3 + z87 + z88 + z90;
	z2 = z136 + z137 + z181 + z193 + z194 + z196 + z198 + z2 + z40 + z62;
	z2 = z2 + z269 + z271 + z273 + z290 + z301 + z306 + z308 + z61 + z64;
	z3 = z222 + z225 + z312 + z313 + z318 + z319 + z326 + z329 + z330 + z338;
	z3 = z100 + z109 + z216 + z217 + z218 + z220 + z221 + z223 + z224 + z3;
	z3 = z108 + z110 + z114 + z116 + z121 + z122 + z123 + z124 + z3 + z77;
	z3 = z106 + z139 + z140 + z144 + z147 + z148 + z3 + z84 + z89 + z97;
	z3 = z153 + z154 + z157 + z159 + z160 + z163 + z165 + z166 + z172 + z3;
	z3 = z104 + z168 + z174 + z195 + z197 + z200 + z208 + z3 + z59 + z91;
	z3 = z214 + z26 + z275 + z277 + z3 + z41 + z67 + z69 + z71 + z73;
	z3 = z278 + z283 + z284 + z291 + z294 + z295 + z3 + z302 + z309;
	z4 = z10 + z249 + z25 + z255 + z263 + z337 + z348 + z350 + z7 + z8;
	z5 = z226 + z248 + z250 + z256 + z257 + z261 + z27 + z43 + z9 + z95;
	z4 = z134 + z21 + z24 + z251 + z4 + z42 + z44 + z5 + z66 + z92;
	z4 = z131 + z179 + z182 + z202 + z206 + z29 + z4 + z55 + z63 + z94;
	z4 = z204 + z207 + z246 + z274 + z28 + z303 + z310 + z36 + z4 + z51;
	z5 = z236 + z242 + z314 + z316 + z320 + z322 + z327 + z331 + z333 + z339;
	z5 = z101 + z227 + z228 + z229 + z234 + z235 + z240 + z241 + z5 + z99;
	z5 = z112 + z113 + z115 + z126 + z128 + z129 + z16 + z5 + z76 + z78;
	z5 = z111 + z117 + z141 + z142 + z143 + z145 + z146 + z150 + z5 + z81;
	z5 = z149 + z151 + z155 + z156 + z158 + z17 + z171 + z175 + z183 + z5;
	z5 = z15 + z173 + z185 + z187 + z192 + z30 + z31 + z45 + z5 + z57;
	z5 = z209 + z211 + z212 + z213 + z276 + z279 + z32 + z46 + z47 + z5;
	z5 = z281 + z285 + z287 + z292 + z296 + z298 + z304 + z311 + z5;
	z6 = z239 + z245 + z258 + z315 + z317 + z321 + z323 + z328 + z332 + z334;
	z6 = z102 + z103 + z230 + z231 + z232 + z237 + z238 + z243 + z244 + z6;
	z6 = z119 + z120 + z127 + z13 + z130 + z19 + z37 + z6 + z79 + z80;
	z6 = z107 + z11 + z118 + z125 + z14 + z161 + z6 + z86 + z93 + z98;
	z6 = z12 + z162 + z164 + z167 + z169 + z176 + z184 + z20 + z52 + z6;
	z6 = z105 + z18 + z186 + z188 + z33 + z34 + z48 + z53 + z6 + z60;
	z6 = z210 + z247 + z280 + z35 + z49 + z50 + z6 + z68 + z70 + z72;
	z6 = z252 + z253 + z282 + z286 + z288 + z293 + z297 + z299 + z6;

	// {{z1, z3, z5},
	//  {z3, z2, z6},
	//  {z5, z6, z4}}

	Q[0] = z1;
	Q[1] = z2;
	Q[2] = z4;

	Q[3] = z6;
	Q[4] = z5;
	Q[5] = z3;
}
