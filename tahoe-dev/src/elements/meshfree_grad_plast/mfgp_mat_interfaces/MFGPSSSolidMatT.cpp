/* $Id: MFGPSSSolidMatT.cpp,v 1.2 2005-05-11 23:10:05 kyonten Exp $  */
#include "MFGPSSSolidMatT.h"
#include "MFGPMatSupportT.h"
#include "dSymMatrixT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<MFGPSSSolidMatT>::fByteCopy = false;
DEFINE_TEMPLATE_STATIC const bool ArrayT<MFGPSSSolidMatT*>::fByteCopy = true;
} /* namespace Tahoe */

/* perturbation used to compute c_ijkl from finite difference */
const double strain_perturbation = 1.0e-08;

/* constructor */
MFGPSSSolidMatT::MFGPSSSolidMatT(void):
	ParameterInterfaceT("mfgp_ss_material"),
	fMFGPMatSupport(NULL)
{

}

/* set the material support or pass NULL to clear */
void MFGPSSSolidMatT::SetMFGPMatSupport(const MFGPMatSupportT* support)
{
	/* set inherited material support */
	MFGPMaterialT::SetMFGPMatSupport(support);

	fMFGPMatSupport = support;

	/* dimension */
	int nsd = NumSD();
	fModulus.Dimension(dSymMatrixT::NumValues(nsd));
}

/* material description */
const dMatrixT& MFGPSSSolidMatT::C_IJKL(void)  { return c_ijkl(); }
const dSymMatrixT& MFGPSSSolidMatT::S_IJ(void) { return s_ij();   }

/* return modulus */
const dMatrixT& MFGPSSSolidMatT::c_ijkl(void)
{
	/* get the strain tensor for the current ip - use the strain
	 * from the material support since the return values from ensure e()
	 * is recomputed when there are thermal strains */
	dSymMatrixT& strain = const_cast<dSymMatrixT&>(fMFGPMatSupport->LinearStrain());

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

/* spatial elastic modulus */
const dMatrixT& MFGPSSSolidMatT::ce_ijkl(void) {
	return c_ijkl();
}

/* apply pre-conditions at the current time step */
void MFGPSSSolidMatT::InitStep(void)
{
	/* inherited */
	MFGPMaterialT::InitStep();
}

/* assumes small strains. returns true if the strain localization condition is satisfied,
* .ie if the acoustic tensor has zero (or negative eigenvalues),
* for the current conditions (current integration point and strain
* state). If localization is detected, the normals (current config)
* to the surface and slip directions are returned */
bool MFGPSSSolidMatT::IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs, 
							AutoArrayT <double> &detAs, AutoArrayT <double> &dissipations_fact)
{
	/* elastic modulus */
	/* this uses same space as c_ijkl(), so save separatley first */
	const dMatrixT modulus_e = ce_ijkl();

	/* localization condition checker */
	DetCheckT checker(s_ij(), c_ijkl(), modulus_e);
	normals.Dimension(NumSD());
	slipdirs.Dimension(NumSD());
	return checker.IsLocalized_SS(normals, slipdirs, detAs);
}

bool MFGPSSSolidMatT::IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs, double &detA)
{
	/* elastic modulus */
	/* this uses same space as c_ijkl(), so save separatley first */
	const dMatrixT modulus_e = ce_ijkl();

	/* localization condition checker */
	DetCheckT checker(s_ij(), c_ijkl(), modulus_e);
	normals.Dimension(NumSD());
	slipdirs.Dimension(NumSD());
	return checker.IsLocalized_SS(normals, slipdirs, detA);
}

bool MFGPSSSolidMatT::IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs)
{
	double dummyDetA = 0.0;
	return IsLocalized(normals, slipdirs, dummyDetA);
}

