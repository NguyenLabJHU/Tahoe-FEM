/* $Id: EAMFCC3D.cpp,v 1.6.2.1 2004-07-06 06:53:26 paklein Exp $ */
/* created: paklein (12/02/1996) */
#include "EAMFCC3D.h"

#include "ParameterContainerT.h"
#include "dMatrixT.h"
#include "StringT.h"

/* EAM glue functions */
#include "ErcolessiAdamsAl.h"
#include "VoterChenAl.h"
#include "VoterChenCu.h"
#include "FBD_EAMGlue.h"
#include "EAM_particle.h"

using namespace Tahoe;

/* bond parameters */
const int kEAMFCC3DNumBonds			= 54;
const int kEAMFCC3DNumLatticeDim 	=  3;
const int kEAMFCC3DNumAtomsPerCell	=  4;

//TEMP
#pragma message("rename me to indicate this is a Cauchy-Born solver")

/* constructor */
EAMFCC3D::EAMFCC3D(void):
	FCCLatticeT(0), /* number of shells is not used by this class */
	fEAM(NULL),
	fLatticeParameter(0.0),
	fCellVolume(0.0)	
{

}

#if 0
EAMFCC3D::EAMFCC3D(ifstreamT& in, int EAMcode, int nsd):
	FCCLatticeT(0), /* number of shells is not used by this class */
	fEAM(NULL),
	fLatticeParameter(0.0),
	fCellVolume(0.0)
{
	/* set EAM solver functions */
	SetGlueFunctions(in, EAMcode, nsd);

	/* lattice parameter and cell volume */
	fLatticeParameter = fEAM->LatticeParameter();
	fCellVolume = fLatticeParameter*fLatticeParameter*fLatticeParameter;
}
#endif

#if 0
EAMFCC3D::EAMFCC3D(ifstreamT& in, const dMatrixT& Q, int EAMcode, int nsd):
	CBLatticeT(Q, numspatialdim, numbonds), fEAMcode(EAMcode)
{
	/* dimension check */
	if( Q.Rows() != kEAMFCC3DNumLatticeDim ||
	    Q.Cols() != kEAMFCC3DNumLatticeDim) throw ExceptionT::kGeneralFail;

	/* set EAM solver functions */
	SetGlueFunctions(in);
	
	/* lattice parameter and cell volume */
	if (!fEAM)
		fLatticeParameter = fEAM_particle->LatticeParameter();
	else
		fLatticeParameter = fEAM->LatticeParameter();
		
	fCellVolume       = fLatticeParameter*fLatticeParameter*fLatticeParameter;
}
#endif

/* destructor */
EAMFCC3D::~EAMFCC3D(void) { delete fEAM; }

/* strain energy density */
double EAMFCC3D::EnergyDensity(const dSymMatrixT& strain)
{
	/* compute deformed lattice geometry */
	ComputeDeformedLengths(strain);

	return (kEAMFCC3DNumAtomsPerCell/fCellVolume)*fEAM->ComputeUnitEnergy();
}

/* return the material tangent moduli in Cij */
void EAMFCC3D::Moduli(dMatrixT& Cij, const dSymMatrixT& strain)
{
	/* compute deformed lattice geometry */
	ComputeDeformedLengths(strain);

	/* unit moduli */
	if (!fEAM)
		fEAM_particle->ComputeUnitModuli(Cij);
	else
		fEAM->ComputeUnitModuli(Cij);
	
	/* scale by atoms per cell/volume per cell */
	Cij	*= kEAMFCC3DNumAtomsPerCell/fCellVolume;
}

/* return the symmetric 2nd PK stress tensor */
void EAMFCC3D::SetStress(const dSymMatrixT& strain, dSymMatrixT& stress)
{
	/* compute deformed lattice geometry */
	ComputeDeformedLengths(strain);

	/* unit stress */
	if (!fEAM)
		fEAM_particle->ComputeUnitStress(stress);
	else
		fEAM->ComputeUnitStress(stress);
	
	/* scale by atoms per cell/volume per cell */
	stress *= kEAMFCC3DNumAtomsPerCell/fCellVolume;
}

/* compute electron density at ghost atom */
void EAMFCC3D::ElectronDensity(const dSymMatrixT& strain, double& edensity, double& embforce)
{
	/* compute deformed lattice geometry */
	ComputeDeformedLengths(strain);
	
	/* get electron density */
	if (!fEAM)
	{
		edensity = fEAM_particle->TotalElectronDensity();
		embforce = fEAM_particle->ReturnEmbeddingForce(edensity);
	}
	else
		edensity = fEAM->TotalElectronDensity();	
}

/**********************************************************************
 * Protected
 **********************************************************************/
	
void EAMFCC3D::LoadBondTable(void)
{
	/* dimension work space */
	fBondCounts.Dimension(kEAMFCC3DNumBonds);
	fDefLength.Dimension(kEAMFCC3DNumBonds);
	fBonds.Dimension(kEAMFCC3DNumBonds, kEAMFCC3DNumLatticeDim);

	/* all bonds appear once */
	fBondCounts = 1;
	
	/* clear deformed lengths for now */
	fDefLength = 0.0;

	/* undeformed bond data for unit cube to 4th nearest neighbors */
	double bonddata[kEAMFCC3DNumBonds][kEAMFCC3DNumLatticeDim] = {
	
		{0, 0, -1.},
		{0, 0, 1.},
		{0, -1., 0},
		{0, 1., 0},
		{-1., 0, 0},
		{1., 0, 0},
		{-0.5, 0, -0.5},
		{-0.5, 0, 0.5},
		{0.5, 0, -0.5},
		{0.5, 0, 0.5},
		{0, -0.5, -0.5},
		{0, -0.5, 0.5},
		{0, 0.5, -0.5},
		{0, 0.5, 0.5},
		{-0.5, -0.5, 0},
		{-0.5, 0.5, 0},
		{0.5, -0.5, 0},
		{0.5, 0.5, 0},
		{-1., 0, -1.},
		{-1., 0, 1.},
		{1., 0, -1.},
		{1., 0, 1.},
		{0, -1., -1.},
		{0, -1., 1.},
		{0, 1., -1.},
		{0, 1., 1.},
		{-1., -1., 0},
		{-1., 1., 0},
		{1., -1., 0},
		{1., 1., 0},
		{-0.5, -0.5, -1.},
		{-0.5, -0.5, 1.},
		{-0.5, 0.5, -1.},
		{-0.5, 0.5, 1.},
		{0.5, -0.5, -1.},
		{0.5, -0.5, 1.},
		{0.5, 0.5, -1.},
		{0.5, 0.5, 1.},
		{-0.5, -1., -0.5},
		{-0.5, -1., 0.5},
		{-0.5, 1., -0.5},
		{-0.5, 1., 0.5},
		{0.5, -1., -0.5},
		{0.5, -1., 0.5},
		{0.5, 1., -0.5},
		{0.5, 1., 0.5},
		{-1., -0.5, -0.5},
		{-1., -0.5, 0.5},
		{-1., 0.5, -0.5},
		{-1., 0.5, 0.5},
		{1., -0.5, -0.5},
		{1., -0.5, 0.5},
		{1., 0.5, -0.5},
		{1., 0.5, 0.5}
	};
				     		
	/* copy into reference table */
	for (int i = 0; i < kEAMFCC3DNumBonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			fBonds(i,j) = bonddata[i][j];
			
	/* scale to correct lattice parameter */				     		
	fBonds *= fLatticeParameter;
}

/* describe the parameters needed by the interface */
void EAMFCC3D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FCCLatticeT::DefineParameters(list);

	/* number of spatial dimensions */
	ParameterT nsd(ParameterT::Integer, "dimensions");
	nsd.AddLimit(2, LimitT::Only);
	nsd.AddLimit(3, LimitT::Only);
	nsd.SetDefault(3);
	list.AddParameter(nsd);
}

/* information about subordinate parameter lists */
void EAMFCC3D::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FCCLatticeT::DefineSubs(sub_list);
	
	/* choice of EAM Cauchy-Born glue functions */
	sub_list.AddSub("EAM_FCC_glue_choice", ParameterListT::Once, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* EAMFCC3D::NewSub(const StringT& list_name) const
{
	if (list_name == "EAM_FCC_glue_choice")
	{
		ParameterContainerT* choice = new ParameterContainerT(list_name);
		choice->SetListOrder(ParameterListT::Choice);

		/* choices */
		ParameterContainerT EA_Al("Ercolessi-Adams_Al");
		ParameterContainerT VC_Al("Voter-Chen_Al");
		ParameterContainerT VC_Cu("Voter-Chen_Cu");
		ParameterContainerT FBD("Paradyn_EAM");
		FBD.AddParameter(ParameterT::Word, "parameter_file");

		choice->AddSub(EA_Al);
		choice->AddSub(VC_Al);
		choice->AddSub(VC_Cu);
		choice->AddSub(FBD);
		
		return choice;
	}
	else /* inherited */
		return FCCLatticeT::NewSub(list_name);
}

/* accept parameter list */
void EAMFCC3D::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "EAMFCC3D::TakeParameterList";

	/* NOTE: initialization is a bit convoluted because fEAM has information needed
	 * in the call to FCCLatticeT::TakeParameterList, when the lattice vectors are
	 * computed, while fEAM can't dimension its internal work space until the vectors
	 * are set. Therefore, the procedure is:
	 * (1) construct fEAM and set fLatticeParameter
	 * (2) call FCCLatticeT::TakeParameterList which sets the bond vectors
	 * (3) dimension fEAM with the number of bond vectors
	 */

	/* construct glue */
	int nsd = list.GetParameter("dimensions");
	const char glue_name[] = "EAM_FCC_glue_choice";
	const ParameterListT& glue = list.GetListChoice(*this, glue_name);
	if (glue.Name() == "Ercolessi-Adams_Al")
		fEAM = new ErcolessiAdamsAl(*this);
	else if (glue.Name() == "Voter-Chen_Al")
		fEAM = new VoterChenAl(*this);
	else if (glue.Name() == "Voter-Chen_Cu")
		fEAM = new VoterChenCu(*this);
	else if (glue.Name() == "Paradyn_EAM")
	{
		/* data file */
		StringT data_file = glue.GetParameter("parameter_file");
		data_file.ToNativePathName();

#if 0
		/* path to source file */
		StringT path;
		path.FilePath(in.filename());
		data_file.Prepend(path);
#endif

		ifstreamT data(data_file);
		if (!data.is_open())
			ExceptionT::GeneralFail(caller, "could not open file \"%s\"", data_file.Pointer());

		fEAM = new FBD_EAMGlue(*this, data);
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized glue function \"%s\"",
			glue.Name().Pointer());

	/* lattice parameter and cell volume */
	fLatticeParameter = fEAM->LatticeParameter();
	fCellVolume = fLatticeParameter*fLatticeParameter*fLatticeParameter;

	/* inherited - lattice parameter needs to be set first */
	FCCLatticeT::TakeParameterList(list);

	/* initialize glue functions */
	fEAM->Initialize(nsd, NumberOfBonds());
}
