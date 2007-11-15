/* $Id: EAMFCC3D_surf.cpp,v 1.14 2007-11-15 15:26:08 hspark Exp $ */
/* created: paklein (12/02/1996) */
#include "EAMFCC3D_surf.h"

#include "ParameterContainerT.h"
#include "dMatrixT.h"
#include "StringT.h"
#include "ifstreamT.h"

/* EAM glue functions */
#include "ErcolessiAdamsAl.h"
#include "VoterChenAl.h"
#include "VoterChenCu.h"
#include "FBD_EAMGlue.h"
#include "EAM_particle.h"

using namespace Tahoe;

/* bond parameters */
const int kEAMFCC3DSurfBonds        = 62;	// 78
const int kEAMFCC3DNumBonds			= 42;	// 54
const int kEAMFCC3DSurf1Bonds       = 25;	// 33
const int kEAMFCC3DSurf2Bonds       = 37;	// 45
const int kEAMFCC3DNumLatticeDim 	=  3;
const int kEAMFCC3DNumAtomsPerCell	=  4;
const int kEAMFCC3DNumAtomsPerArea  =  2;
const double piby2 = 4.0 * atan(1.0) / 2.0;

// WHAT ABOUT FACTOR OF 0.5 IN PAIR POTENTIAL TERMS IN EAM.CPP?

/* constructor */
EAMFCC3D_surf::EAMFCC3D_surf(int nshells, int normal):
	FCCLatticeT_Surf(nshells, normal),
	fEAM(NULL),
	fEAM_particle(NULL),
	fLatticeParameter(0.0),
	fSurfaceThickness(-1),
	fCellArea(0),
	fNormalCode(normal),
	fCellVolume(0.0)	
{

}

/* destructor */
EAMFCC3D_surf::~EAMFCC3D_surf(void) { 
	delete fEAM; 
	delete fEAM_particle;	
}

/* strain energy density */
double EAMFCC3D_surf::EnergyDensity(const dSymMatrixT& strain)
{
	/* compute deformed lattice geometry */
	ComputeDeformedLengths(strain);

	/* scale by atoms per cell/AREA per cell, split energy by one-half */
	return (kEAMFCC3DNumAtomsPerArea/fCellArea)*fEAM->ComputeUnitSurfaceEnergy();
}

/* return the material tangent moduli in Cij */
void EAMFCC3D_surf::Moduli(dMatrixT& Cij, const dSymMatrixT& strain)
{
	/* compute deformed lattice geometry */
	ComputeDeformedLengths(strain);

	/* unit moduli */
	if (!fEAM)
		fEAM_particle->ComputeUnitModuli(Cij);
	else
		fEAM->ComputeUnitModuli(Cij);
	
	/* scale by atoms per cell/AREA per cell, split by one half for counting all bonds */
	Cij	*= kEAMFCC3DNumAtomsPerArea/fCellArea;
}

/* return the symmetric 2nd PK surface stress tensor */
void EAMFCC3D_surf::SetStress(const dSymMatrixT& strain, dSymMatrixT& stress)
{
	/* compute deformed lattice geometry */
	ComputeDeformedBulkBonds(strain);

	/* Compute deformed lattice geometry for surface/bulk unit cells */
	ComputeDeformedLengths(strain);
	ComputeDeformedSurf1Bonds(strain);
	ComputeDeformedSurf2Bonds(strain);

	/* unit stress */
	if (!fEAM)
		fEAM_particle->ComputeUnitStress(stress);	// IMPLEMENT THIS?
	else
		fEAM->ComputeUnitSurfaceStress(stress);
	
	/* scale by atoms per cell/AREA per cell, split by one half for counting all bonds */
	/* MAY NOT NEED TO SPLIT BONDS BY ONE HALF */
	stress *= kEAMFCC3DNumAtomsPerArea/fCellArea;
}

/* compute electron density at ghost atom */
//void EAMFCC3D_surf::ElectronDensity(const dSymMatrixT& strain, double& edensity, double& embforce)
//{
	/* compute deformed lattice geometry */
//	ComputeDeformedLengths(strain);
	
	/* get electron density */
//	if (!fEAM)
//	{
//		edensity = fEAM_particle->TotalElectronDensity();
//		embforce = fEAM_particle->ReturnEmbeddingForce(edensity);
//	}
//	else
//		edensity = fEAM->TotalElectronDensity();	
//}

/* density computed from the potential parameters */
double EAMFCC3D_surf::Density(void) const {
	if (fEAM_particle) {
		return kEAMFCC3DNumAtomsPerCell*fEAM_particle->Mass()/fCellVolume;
	} else {
		return kEAMFCC3DNumAtomsPerCell*fEAM->Mass()/fCellVolume;
	}
}

/**********************************************************************
 * Protected
 **********************************************************************/
	
void EAMFCC3D_surf::LoadBondTable(void)
{
	/* dimension work space - ARE THESE DIMENSIONS CORRECT? */
	fBondCounts.Dimension(kEAMFCC3DSurfBonds);
	fBulkCounts.Dimension(kEAMFCC3DNumBonds);
	fSurf1Counts.Dimension(kEAMFCC3DSurf1Bonds);
	fSurf2Counts.Dimension(kEAMFCC3DSurf2Bonds);
	fDefLength.Dimension(kEAMFCC3DSurfBonds);
	fDefBulk.Dimension(kEAMFCC3DNumBonds);
	fDefSurf1.Dimension(kEAMFCC3DSurf1Bonds);
	fDefSurf2.Dimension(kEAMFCC3DSurf2Bonds);
	fBonds.Dimension(kEAMFCC3DSurfBonds, kEAMFCC3DNumLatticeDim);
	fBulkBonds.Dimension(kEAMFCC3DNumBonds,3);
	fSurf1Bonds.Dimension(kEAMFCC3DSurf1Bonds,3);
	fSurf2Bonds.Dimension(kEAMFCC3DSurf2Bonds,3);

	dArray2DT temp_bonds, temp_bonds2;
	temp_bonds.Dimension(kEAMFCC3DSurfBonds, 3);	// temporary bond table before rotation
	temp_bonds2.Dimension(kEAMFCC3DSurfBonds, 3);	// Currently have # of bonds for {100} surfaces

	/* all bonds appear once */
	fBondCounts = 1;
	
	/* clear deformed lengths for now */
	fDefLength = 0.0;
	fDefBulk = 0.0;
	fDefSurf1 = 0.0;
	fDefSurf2 = 0.0;

	/* undeformed bond data for bulk atom with 4th neighbor interactions */
	double bulkbond[kEAMFCC3DNumBonds][kEAMFCC3DNumLatticeDim] = {
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

	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DNumBonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			fBulkBonds(i,j) = bulkbond[i][j];

	/* Bond table for an atom on the surface */
	double surf1bond[kEAMFCC3DSurf1Bonds][kEAMFCC3DNumLatticeDim] = {
		{0.5, 0.5, 0.0}, // Surface cluster (8 nearest neighbors)
		{0.5, -0.5, 0.0},
		{0.5, 0.0, 0.5},
		{0.5, 0.0, -0.5},
		{0.0, 0.5, 0.5},
		{0.0, -0.5, 0.5},
		{0.0, 0.5, -0.5},
		{0.0, -0.5, -0.5},
		{1.0, 0.0, 0.0}, // Surface cluster (5 2nd shell neighbors)
		{0.0, 1.0, 0.0},
		{0.0, 0.0, 1.0},
		{0.0, -1.0, 0.0},
		{0.0, 0.0, -1.0},
		{1.0, 0.5, 0.5}, // Surface cluster (12 3rd shell neighbors)
		{0.5, 1.0, 0.5},
		{0.5, 0.5, 1.0},
		{1.0, 0.5, -0.5},
		{0.5, 1.0, -0.5},
		{0.5, 0.5, -1.0},
		{1.0, -0.5, 0.5},
		{0.5, -1.0, 0.5},
		{0.5, -0.5, 1.0},
		{1.0, -0.5, -0.5},
		{0.5, -1.0, -0.5},
		{0.5, -0.5, -1.0}
//		{0.0, 1.0, -1.0}, // Surface cluster (8 4th shell neighbors)
//		{0.0, 1.0, 1.0},
//		{1.0, 1.0, 0.0},
//		{1.0, -1.0, 0.0},
//		{1.0, 0.0, 1.0},
//		{1.0, 0.0, -1.0},
//		{0.0, -1.0, -1.0},
//		{0.0, -1.0, 1.0}
	};

	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DSurf1Bonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			fSurf1Bonds(i,j) = surf1bond[i][j];

	/* Bond table for an atom 1 layer into the bulk */
	double surf2bond[kEAMFCC3DSurf2Bonds][kEAMFCC3DNumLatticeDim] = {
		{0.5, 0.5, 0.0}, // Surface cluster (8 nearest neighbors)
		{0.5, -0.5, 0.0},
		{0.5, 0.0, 0.5},
		{0.5, 0.0, -0.5},
		{0.0, 0.5, 0.5},
		{0.0, -0.5, 0.5},
		{0.0, 0.5, -0.5},
		{0.0, -0.5, -0.5},
		{-0.5, -0.5, 0.0}, // New bonds for second surface cluster begin here
		{-0.5, 0.5, 0.0},
		{-0.5, 0.0, 0.5},
		{-0.5, 0.0, -0.5},
		{1.0, 0.0, 0.0}, // Surface cluster (5 2nd shell neighbors)
		{0.0, 1.0, 0.0},
		{0.0, 0.0, 1.0},
		{0.0, -1.0, 0.0},
		{0.0, 0.0, -1.0},
		{1.0, 0.5, 0.5}, // Surface cluster (12 3rd shell neighbors)
		{0.5, 1.0, 0.5},
		{0.5, 0.5, 1.0},
		{1.0, 0.5, -0.5},
		{0.5, 1.0, -0.5},
		{0.5, 0.5, -1.0},
		{1.0, -0.5, 0.5},
		{0.5, -1.0, 0.5},
		{0.5, -0.5, 1.0},
		{1.0, -0.5, -0.5},
		{0.5, -1.0, -0.5},
		{0.5, -0.5, -1.0},
		{-0.5, 1.0, 0.5}, // New bonds for second surface cluster begin here
		{-0.5, 1.0, -0.5},
		{-0.5, 0.5, 1.0},
		{-0.5, 0.5, -1.0},
		{-0.5, -0.5, 1.0},
		{-0.5, -0.5, -1.0},
		{-0.5, -1.0, 0.5},
		{-0.5, -1.0, -0.5}
//		{0.0, 1.0, -1.0}, // Surface cluster (8 4th shell neighbors)
//		{0.0, 1.0, 1.0},
//		{1.0, 1.0, 0.0},
//		{1.0, -1.0, 0.0},
//		{1.0, 0.0, 1.0},
//		{1.0, 0.0, -1.0},
//		{0.0, -1.0, -1.0},
//		{0.0, -1.0, 1.0}
	};

	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DSurf2Bonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			fSurf2Bonds(i,j) = surf2bond[i][j];

	/* New bond table for surface clusters - change dimensions! */
	/* AVOID HARD CODING NUMBER OF BONDS SPECIFIC TO {100} */
	double bonddata[kEAMFCC3DSurfBonds][kEAMFCC3DNumLatticeDim] = {
		{0.5, 0.5, 0.0}, // Surface cluster (8 nearest neighbors)
		{0.5, -0.5, 0.0},
		{0.5, 0.0, 0.5},
		{0.5, 0.0, -0.5},
		{0.0, 0.5, 0.5},
		{0.0, -0.5, 0.5},
		{0.0, 0.5, -0.5},
		{0.0, -0.5, -0.5},
		{0.5, 0.5, 0.0}, // Repeat cluster here - one atomic thickness into bulk
		{0.5, -0.5, 0.0}, // Total of 12 nearest neighbors
		{0.5, 0.0, 0.5},
		{0.5, 0.0, -0.5},
		{0.0, 0.5, 0.5},
		{0.0, -0.5, 0.5},
		{0.0, 0.5, -0.5},
		{0.0, -0.5, -0.5},
		{-0.5, -0.5, 0.0}, // New bonds for second surface cluster begin here
		{-0.5, 0.5, 0.0},
		{-0.5, 0.0, 0.5},
		{-0.5, 0.0, -0.5},
		{1.0, 0.0, 0.0}, // Surface cluster (5 2nd shell neighbors)
		{0.0, 1.0, 0.0},
		{0.0, 0.0, 1.0},
		{0.0, -1.0, 0.0},
		{0.0, 0.0, -1.0},
		{1.0, 0.0, 0.0}, // Repeat cluster here - one atomic thickness into bulk
		{0.0, 1.0, 0.0}, // Total of 5 2nd shell neighbors
		{0.0, 0.0, 1.0},
		{0.0, -1.0, 0.0},
		{0.0, 0.0, -1.0},
		{1.0, 0.5, 0.5}, // Surface cluster (12 3rd shell neighbors)
		{0.5, 1.0, 0.5},
		{0.5, 0.5, 1.0},
		{1.0, 0.5, -0.5},
		{0.5, 1.0, -0.5},
		{0.5, 0.5, -1.0},
		{1.0, -0.5, 0.5},
		{0.5, -1.0, 0.5},
		{0.5, -0.5, 1.0},
		{1.0, -0.5, -0.5},
		{0.5, -1.0, -0.5},
		{0.5, -0.5, -1.0},
		{1.0, 0.5, 0.5}, // Repeat cluster here - one atomic thickness into bulk
		{0.5, 1.0, 0.5}, // Total of 20 3rd shell neighbors
		{0.5, 0.5, 1.0},
		{1.0, 0.5, -0.5},
		{0.5, 1.0, -0.5},
		{0.5, 0.5, -1.0},
		{1.0, -0.5, 0.5},
		{0.5, -1.0, 0.5},
		{0.5, -0.5, 1.0},
		{1.0, -0.5, -0.5},
		{0.5, -1.0, -0.5},
		{0.5, -0.5, -1.0},
		{-0.5, 1.0, 0.5}, // New bonds for second surface cluster begin here
		{-0.5, 1.0, -0.5},
		{-0.5, 0.5, 1.0},
		{-0.5, 0.5, -1.0},
		{-0.5, -0.5, 1.0},
		{-0.5, -0.5, -1.0},
		{-0.5, -1.0, 0.5},
		{-0.5, -1.0, -0.5}
//		{0.0, 1.0, -1.0}, // Surface cluster (8 4th shell neighbors)
//		{0.0, 1.0, 1.0},
//		{1.0, 1.0, 0.0},
//		{1.0, -1.0, 0.0},
//		{1.0, 0.0, 1.0},
//		{1.0, 0.0, -1.0},
//		{0.0, -1.0, -1.0},
//		{0.0, -1.0, 1.0},
//		{0.0, 1.0, -1.0}, // Repeat cluster here - one atomic thickness into the bulk
//		{0.0, 1.0, 1.0},
//		{1.0, 1.0, 0.0},
//		{1.0, -1.0, 0.0},
//		{1.0, 0.0, 1.0},
//		{1.0, 0.0, -1.0},
//		{0.0, -1.0, -1.0},
//		{0.0, -1.0, 1.0}
	};
	
	/* Rotate Bond Tables based on fNormalCode and rotation matrices */
	/* Create temporary bond table temp_bonds that combines bonddata */
	for (int i = 0; i < kEAMFCC3DSurfBonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_bonds(i,j) = bonddata[i][j];
	
	/* Now manipulate temp_bonds */
	dMatrixT blah1(3);
	dArrayT asdf(3), prod(3);
	if (fNormalCode == 0)	// normal is [1,0,0]
	{
		temp_bonds2 = temp_bonds;
		fBonds = temp_bonds2;
		fBonds *= -1.0;
	}
	else if (fNormalCode == 1)
		fBonds = temp_bonds;	// this table is the default, i.e. [-1,0,0]
	else if (fNormalCode == 2)	// rotate [-1,0,0] to [0,1,0]
	{
		temp_bonds2 = temp_bonds;
		blah1 = RotationMatrixA(piby2);
		for (int i = 0; i < kEAMFCC3DSurfBonds; i++)
		{
			temp_bonds2.RowCopy(i,asdf);	// take bond
			blah1.Multx(asdf,prod);		// rotate bond via rotation matrix
			temp_bonds2.SetRow(i,prod);	// place new bond back into temp_bonds
		}
		fBonds = temp_bonds2;
	}
	else if (fNormalCode == 3)	// rotate [-1,0,0] to [0,-1,0]
	{
		temp_bonds2 = temp_bonds;
		blah1 = RotationMatrixA(-piby2);
		for (int i = 0; i < kEAMFCC3DSurfBonds; i++)
		{
			temp_bonds2.RowCopy(i,asdf);	// take bond
			blah1.Multx(asdf,prod);		// rotate bond via rotation matrix
			temp_bonds2.SetRow(i,prod);	// place new bond back into temp_bonds
		}	
		fBonds = temp_bonds2;
	}
	else if (fNormalCode == 4)	// rotate [-1,0,0] to [0,0,1]
	{
		temp_bonds2 = temp_bonds;
		blah1 = RotationMatrixB(-piby2);
		for (int i = 0; i < kEAMFCC3DSurfBonds; i++)
		{
			temp_bonds2.RowCopy(i,asdf);	// take bond
			blah1.Multx(asdf,prod);		// rotate bond via rotation matrix
			temp_bonds2.SetRow(i,prod);	// place new bond back into temp_bonds
		}	
		fBonds = temp_bonds2;
	}	
	else if (fNormalCode == 5)	// rotate [-1,0,0] to [0,0,-1]
	{
		temp_bonds2 = temp_bonds;
		blah1 = RotationMatrixB(piby2);
		for (int i = 0; i < kEAMFCC3DSurfBonds; i++)
		{
			temp_bonds2.RowCopy(i,asdf);	// take bond
			blah1.Multx(asdf,prod);		// rotate bond via rotation matrix
			temp_bonds2.SetRow(i,prod);	// place new bond back into temp_bonds
		}	
		fBonds = temp_bonds2;
	}	
	
	/* copy into reference table */
	/* AVOID HARD CODING NUMBER OF BONDS DUE TO {100} SURFACES */
	/* OBVIATED BY fBonds ABOVE
	for (int i = 0; i < 78; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			fBonds(i,j) = bonddata[i][j];
	
	*/
	/* scale to correct lattice parameter */				     		
	fBonds *= fLatticeParameter;
	fBulkBonds *= fLatticeParameter;
	fSurf1Bonds *= fLatticeParameter;
	fSurf2Bonds *= fLatticeParameter;
}

/* describe the parameters needed by the interface */
void EAMFCC3D_surf::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FCCLatticeT_Surf::DefineParameters(list);

	/* number of spatial dimensions */
	ParameterT nsd(ParameterT::Integer, "dimensions");
	nsd.AddLimit(2, LimitT::Only);
	nsd.AddLimit(3, LimitT::Only);
	nsd.SetDefault(3);
	list.AddParameter(nsd);
}

/* information about subordinate parameter lists */
void EAMFCC3D_surf::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FCCLatticeT_Surf::DefineSubs(sub_list);
	
	/* choice of EAM Cauchy-Born glue functions */
	sub_list.AddSub("EAM_FCC_glue_choice", ParameterListT::Once, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* EAMFCC3D_surf::NewSub(const StringT& name) const
{
	if (name == "EAM_FCC_glue_choice")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		/* choices */
		ParameterContainerT EA_Al("Ercolessi-Adams_Al");
		ParameterContainerT VC_Al("Voter-Chen_Al");
		ParameterContainerT VC_Cu("Voter-Chen_Cu");
		ParameterContainerT FBD("Paradyn_EAM_glue");
		FBD.AddParameter(ParameterT::Word, "parameter_file");
		ParameterContainerT particle_Paradyn_EAM("particle_Paradyn_EAM_glue");
		particle_Paradyn_EAM.AddParameter(ParameterT::Word, "parameter_file");

		choice->AddSub(EA_Al);
		choice->AddSub(VC_Al);
		choice->AddSub(VC_Cu);
		choice->AddSub(FBD);
		choice->AddSub(particle_Paradyn_EAM);
		
		return choice;
	}
	else /* inherited */
		return FCCLatticeT_Surf::NewSub(name);
}

/* accept parameter list */
void EAMFCC3D_surf::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "EAMFCC3D_surf::TakeParameterList";

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
	else if (glue.Name() == "Paradyn_EAM_glue")
	{
		/* data file */
		StringT data_file = glue.GetParameter("parameter_file");
		data_file.ToNativePathName();

		/* path to source file */
		data_file.Prepend(fstreamT::Root());

		ifstreamT data(data_file);
		if (!data.is_open())
			ExceptionT::GeneralFail(caller, "could not open file \"%s\"", data_file.Pointer());

		fEAM = new FBD_EAMGlue(*this, data);
	}
	else if (glue.Name() == "particle_Paradyn_EAM_glue")
	{
		/* data file */
		StringT data_file = glue.GetParameter("parameter_file");
		data_file.ToNativePathName();

		/* path to source file */
		data_file.Prepend(fstreamT::Root());

		fEAM_particle = new EAM_particle(*this, data_file);
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized glue function \"%s\"",
			glue.Name().Pointer());

	/* lattice parameter and cell volume */
	fLatticeParameter = (fEAM) ? fEAM->LatticeParameter() : fEAM_particle->LatticeParameter();
	fCellVolume = fLatticeParameter*fLatticeParameter*fLatticeParameter;
	fCellArea = fLatticeParameter*fLatticeParameter; // OK because dividing by 2 in Stress/Modulus calls
	//fSurfaceThickness = 0.75*fLatticeParameter; // CHECK THIS
	fSurfaceThickness = 1.75*fLatticeParameter;	// Accounting for S3 and S4 layers

	/* inherited - lattice parameter needs to be set first */
	FCCLatticeT_Surf::TakeParameterList(list);

	/* initialize glue functions */
	if (fEAM)
		fEAM->Initialize(nsd, NumberOfBonds());
	else
		fEAM_particle->Initialize(nsd, NumberOfBonds());
}

/*************************************************************************
 * Private
 *************************************************************************/
 
 /* Rotate bonds with [-1,0,0] normal to bonds with [0,1,0]-type normals */
dMatrixT EAMFCC3D_surf::RotationMatrixA(const double angle)
 {
	dMatrixT rmatrix(3);
	rmatrix = 0.0;
    rmatrix(0,0) = cos(angle);
	rmatrix(0,1) = sin(angle);
	rmatrix(1,0) = -sin(angle);
	rmatrix(1,1) = cos(angle);
	rmatrix(2,2) = 1.0;
	
	return rmatrix;
 }
 
/* Rotate bonds with [-1,0,0] normal to bonds with [0,0,1]-type normals */
dMatrixT EAMFCC3D_surf::RotationMatrixB(const double angle)
{
	dMatrixT rmatrix(3);
	rmatrix = 0.0;
    rmatrix(0,0) = cos(angle);
	rmatrix(0,2) = -sin(angle);
	rmatrix(1,1) = 1.0;
	rmatrix(2,0) = sin(angle);
	rmatrix(2,2) = cos(angle);
	
	return rmatrix;
}