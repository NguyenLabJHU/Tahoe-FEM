/* $Id: SSSolidMatT.cpp,v 1.8.28.1 2004-09-22 06:13:33 paklein Exp $ */
/* created: paklein (06/09/1997) */
#include "SSSolidMatT.h"
#include <iostream.h>
#include "SSMatSupportT.h"
#include "dSymMatrixT.h"
#include "ThermalDilatationT.h"

using namespace Tahoe;

/* perturbation used to compute c_ijkl from finite difference */
const double strain_perturbation = 1.0e-08;

/* constructor */
SSSolidMatT::SSSolidMatT(ifstreamT& in, const SSMatSupportT& support):
	SolidMaterialT(in, support),
	fSSMatSupport(&support),
	fModulus(dSymMatrixT::NumValues(NumSD())),
	fStrainTemp(NumSD()),
	fQ(NumSD()),
	fThermalStrain(NumSD()),
	fHasThermalStrain(false)
{
	SetName("small_strain_solid_material");
}

SSSolidMatT::SSSolidMatT(void):
	fSSMatSupport(NULL),
	fHasThermalStrain(false)
{
	SetName("small_strain_solid_material");
}

/* I/O */
void SSSolidMatT::PrintName(ostream& out) const
{
	/* inherited */
	SolidMaterialT::PrintName(out);
	
	out << "    Small strain\n";
}

/* strain - returns the elastic strain, ie. thermal removed */
const dSymMatrixT& SSSolidMatT::e(void)
{
	/* remove thermal strain */
	if (fHasThermalStrain)
	{
		/* thermal strain is purely dilatational */
		fStrainTemp  = fSSMatSupport->LinearStrain();
		fStrainTemp -= fThermalStrain;
		return fStrainTemp;
	}
	else
		return fSSMatSupport->LinearStrain();
}

/* elastic strain at the given integration point */
const dSymMatrixT& SSSolidMatT::e(int ip)
{
	/* remove thermal strain */
	if (fHasThermalStrain)
	{
		/* thermal strain is purely dilatational */
		fStrainTemp  = fSSMatSupport->LinearStrain(ip);
		fStrainTemp -= fThermalStrain;
		return fStrainTemp;
	}
	else
		return fSSMatSupport->LinearStrain(ip);
}

/* strain - returns the elastic strain, ie. thermal removed */
const dSymMatrixT& SSSolidMatT::e_last(void)
{
	/* cannot have thermal strain */
	if (fHasThermalStrain)
		ExceptionT::GeneralFail("SSSolidMatT::e_last", "not available with thermal strains");
	return fSSMatSupport->LinearStrain_last();
}

/* elastic strain at the given integration point */
const dSymMatrixT& SSSolidMatT::e_last(int ip)
{
	/* cannot have thermal strain */
	if (fHasThermalStrain)
		ExceptionT::GeneralFail("SSSolidMatT::e_last", "not available with thermal strains");
	return fSSMatSupport->LinearStrain_last(ip);
}

/* material description */
const dMatrixT& SSSolidMatT::C_IJKL(void)  { return c_ijkl(); }
const dSymMatrixT& SSSolidMatT::S_IJ(void) { return s_ij();   }

/* return modulus */
const dMatrixT& SSSolidMatT::c_ijkl(void)
{
	/* get the strain tensor for the current ip - use the strain
	 * from the material support since the return values from ensure e()
	 * is recomputed when there are thermal strains */
	dSymMatrixT& strain = const_cast<dSymMatrixT&>(fSSMatSupport->LinearStrain());

	/* compute columns of modulus */
	for (int i = 0; i < fModulus.Cols(); i++) {

		/* perturb strain */
		strain[i] += strain_perturbation;
	
		/* compute stress */
		const dSymMatrixT& stress = s_ij();
	
		/* write into modulus */
		fModulus.SetCol(i, stress);
		
		/* undo perturbation */
		strain[i] -= strain_perturbation;
	}
	
	/* restore stress to unperturbed state */
	const dSymMatrixT& stress = s_ij();
	
	/* compute modulus from finite difference */
	int nsd = NumSD();
	double den = strain_perturbation;
	for (int i = 0; i < fModulus.Cols(); i++) {

		/* shear strains */
		if (i == nsd) den *= 2.0;

		for (int j = 0; j < fModulus.Rows(); j++)
			fModulus(j,i) = (fModulus(j,i) - stress[j])/den;
	}

	return fModulus;
}

/* apply pre-conditions at the current time step */
void SSSolidMatT::InitStep(void)
{
	/* inherited */
	SolidMaterialT::InitStep();

	/* thermal strain */
	fHasThermalStrain = SetThermalStrain(fThermalStrain);
}

/*************************************************************************
* Protected
*************************************************************************/

/* set the internal thermal strain */
bool SSSolidMatT::SetThermalStrain(dSymMatrixT& thermal_strain)
{
	thermal_strain = 0;
	if (fThermal->IsActive())
	{
		thermal_strain.PlusIdentity(fThermal->PercentElongation());
		return true;
	}
	else
		return false;
}

/* return the acoustical tensor and wave speeds */
const dSymMatrixT& SSSolidMatT::AcousticalTensor(const dArrayT& normal)
{
#if __option(extended_errorcheck)
	if (fQ.Rows() != normal.Length()) throw ExceptionT::kSizeMismatch;
#endif

	/* fetch modulus */
	const dMatrixT& c_ = c_ijkl();

	if (normal.Length() == 2)
		Q_2D(c_, normal, fQ);
	else if (normal.Length() == 3)
		Q_3D(c_, normal, fQ);
	else
		throw ExceptionT::kGeneralFail;

	return fQ;
}

/*************************************************************************
* Private
*************************************************************************/

/* acoustical tensor routines */
void SSSolidMatT::Q_2D(const dMatrixT& c_ijkl, const dArrayT& n, dSymMatrixT& Q) const
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10;

	double n1 = n[0];
	double n2 = n[1];
	
	double c11 = c_ijkl[0];
	double c22 = c_ijkl[4];
	double c33 = c_ijkl[8];
	double c23 = c_ijkl[7];
	double c13 = c_ijkl[6];
	double c12 = c_ijkl[3];

	z1 = n1*n1;
	z2 = n1*n2;
	z3 = n2*n2;
	z4 = c12*z2;
	z5 = 2.*c13*z2;
	z6 = 2.*c23*z2;
	z2 = c33*z2;
	z7 = c11*z1;
	z8 = c13*z1;
	z1 = c33*z1;
	z9 = c22*z3;
	z10 = c23*z3;
	z3 = c33*z3;
	z1 = z1 + z6 + z9;
	z3 = z3 + z5 + z7;
	z2 = z10 + z2 + z4 + z8;

	//{{z3, z2},
	// {z2, z1}}
	
	Q[0] = z3;
	Q[1] = z1;
	Q[2] = z2;
}

void SSSolidMatT::Q_3D(const dMatrixT& c_ijkl, const dArrayT& n, dSymMatrixT& Q) const
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24;
	double z25, z26, z27, z28, z29, z30, z31, z32, z33, z34, z35, z36;
	double z37, z38, z39, z40, z41, z42, z43, z44, z45;

	double n1 = n[0];
	double n2 = n[1];
	double n3 = n[2];
	
	double c11 = c_ijkl[0];
	double c12 = c_ijkl[1];
	double c13 = c_ijkl[2];
	double c14 = c_ijkl[3];
	double c15 = c_ijkl[4];
	double c16 = c_ijkl[5];
	double c22 = c_ijkl[7];
	double c23 = c_ijkl[8];
	double c24 = c_ijkl[9];
	double c25 = c_ijkl[10];
	double c26 = c_ijkl[11];
	double c33 = c_ijkl[14];
	double c34 = c_ijkl[15];
	double c35 = c_ijkl[16];
	double c36 = c_ijkl[17];
	double c44 = c_ijkl[21];
	double c45 = c_ijkl[22];
	double c46 = c_ijkl[23];
	double c55 = c_ijkl[28];
	double c56 = c_ijkl[29];
	double c66 = c_ijkl[35];

	z1 = 2.*n1;
	z2 = c14*n1;
	z3 = c45*n1;
	z4 = c46*n1;
	z5 = c56*n1;
	z6 = n1*n1;
	z7 = 2.*n2;
	z8 = c25*n2;
	z9 = c45*n2;
	z10 = c46*n2;
	z11 = n1*n2;
	z12 = n2*n2;
	z13 = c36*n3;
	z14 = c46*n3;
	z15 = c56*n3;
	z16 = n1*n3;
	z17 = n2*n3;
	z18 = n3*n3;
	z19 = n2*z1;
	z20 = n3*z1;
	z10 = n3*z10;
	z21 = c12*z11;
	z11 = c66*z11;
	z22 = n1*z13;
	z13 = n2*z13;
	z14 = z1*z14;
	z23 = c13*z16;
	z16 = c55*z16;
	z24 = n2*z2;
	z2 = n3*z2;
	z25 = c23*z17;
	z17 = c44*z17;
	z26 = c16*z19;
	z19 = c26*z19;
	z27 = c15*z20;
	z20 = c35*z20;
	z3 = n3*z3;
	z4 = n2*z4;
	z28 = n2*z5;
	z5 = n3*z5;
	z29 = n3*z7;
	z30 = c24*z29;
	z29 = c34*z29;
	z7 = z15*z7;
	z15 = n1*z8;
	z8 = n3*z8;
	z31 = n3*z9;
	z1 = z1*z9;
	z9 = c22*z12;
	z32 = c24*z12;
	z33 = c26*z12;
	z34 = c44*z12;
	z35 = c46*z12;
	z12 = c66*z12;
	z36 = c33*z18;
	z37 = c34*z18;
	z38 = c35*z18;
	z39 = c44*z18;
	z40 = c45*z18;
	z18 = c55*z18;
	z41 = c11*z6;
	z42 = c15*z6;
	z43 = c16*z6;
	z44 = c55*z6;
	z45 = c56*z6;
	z6 = c66*z6;
	z7 = z12 + z18 + z26 + z27 + z41 + z7;
	z12 = z13 + z16 + z23 + z24 + z28 + z31 + z35 + z38 + z42;
	z2 = z10 + z11 + z2 + z21 + z33 + z40 + z43 + z5 + z8;
	z1 = z1 + z20 + z29 + z34 + z36 + z44;
	z3 = z15 + z17 + z22 + z25 + z3 + z32 + z37 + z4 + z45;
	z4 = z14 + z19 + z30 + z39 + z6 + z9;

	//{{z7, z2, z12},
	// {z2, z4, z3},
	// {z12, z3, z1}}

	Q[0] = z7;	
	Q[1] = z4;	
	Q[2] = z1;	
	Q[3] = z3;	
	Q[4] = z12;	
	Q[5] = z2;	
}
