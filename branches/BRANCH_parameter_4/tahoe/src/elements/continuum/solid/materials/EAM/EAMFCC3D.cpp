/* $Id: EAMFCC3D.cpp,v 1.6 2004-06-17 07:40:19 paklein Exp $ */
/* created: paklein (12/02/1996)                                          */
/* EAMFCC3D.cpp                                                           */

#include "EAMFCC3D.h"

#include <iostream.h>

#include "ifstreamT.h"
#include "dMatrixT.h"
#include "StringT.h"

/* EAM glue functions */
#include "ErcolessiAdamsAl.h"
#include "VoterChenAl.h"
#include "VoterChenCu.h"
#include "FBD_EAMGlue.h"
#include "EAM_particle.h"

/* constructor */

using namespace Tahoe;

EAMFCC3D::EAMFCC3D(ifstreamT& in, int EAMcode, int numspatialdim, int numbonds):
	CBLatticeT(kEAMFCC3DNumLatticeDim, numspatialdim, numbonds),
	fEAMcode(EAMcode)
{
	/* set EAM solver functions */
	SetGlueFunctions(in);

	/* lattice parameter and cell volume */
	if (!fEAM)
		fLatticeParameter = fEAM_particle->LatticeParameter();
	else
		fLatticeParameter = fEAM->LatticeParameter();
		
	fCellVolume       = fLatticeParameter*fLatticeParameter*fLatticeParameter;
}

EAMFCC3D::EAMFCC3D(ifstreamT& in, const dMatrixT& Q, int EAMcode, int numspatialdim,
	int numbonds):
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

/* destructor */
EAMFCC3D::~EAMFCC3D(void) { delete fEAM; }

/* strain energy density */
double EAMFCC3D::EnergyDensity(const dSymMatrixT& strain)
{
	/* compute deformed lattice geometry */
	ComputeDeformedLengths(strain);

	if (!fEAM)
		return (kEAMFCC3DNumAtomsPerCell/fCellVolume)*fEAM_particle->ComputeUnitEnergy();
	else
		return  (kEAMFCC3DNumAtomsPerCell/fCellVolume)*fEAM->ComputeUnitEnergy();
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

/* I/O functions */
void EAMFCC3D::Print(ostream& out) const
{
	out << " EAM glue function code. . . . . . . . . . . . . = " << fEAMcode << '\n';
	out << "    eq. " << kErcolessiAdamsAl << ", Ercolessi & Adams Al\n";
	out << "    eq. " << kVoterChenAl      << ", Voter & Chen Al\n";   	
	out << "    eq. " << kVoterChenCu      << ", Voter & Chen Cu\n";   	
	out << "    eq. " << kFoilesBaskesDaw  << ", Foiles, Baskes, Daw\n";   	
	out << "    eq. " << kEAMParticle      << ", Foiles, Baskes, Daw(Particle)\n";
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

/* initialize bond tables */
void EAMFCC3D::InitBondTables(void)
{
	BondLatticeT::Initialize();
}

/**********************************************************************
* Protected
**********************************************************************/
	
void EAMFCC3D::LoadBondTable(void)
{
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

	/* dimension check */				     		
	if (fBonds.MajorDim() != fNumBonds ||
	    fBonds.MinorDim() != kEAMFCC3DNumLatticeDim) throw ExceptionT::kGeneralFail;
				     		
	/* copy into reference table */
	for (int i = 0; i < fNumBonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			fBonds(i,j) = bonddata[i][j];
			
	/* scale to correct lattice parameter */				     		
	fBonds *= fLatticeParameter;
}

/**********************************************************************
* Private
**********************************************************************/

/* Set glue functions */
void EAMFCC3D::SetGlueFunctions(ifstreamT& in)
{
	switch (fEAMcode)
	{
		case kErcolessiAdamsAl:
			fEAM = new ErcolessiAdamsAl(*this);
			break;
			
		case kVoterChenAl:
			fEAM = new VoterChenAl(*this);
			break;
	
		case kVoterChenCu:
			fEAM = new VoterChenCu(*this);
			break;
	
		case kFoilesBaskesDaw:
		{
			/* data file */
			StringT data_file;
			in >> data_file;
			data_file.ToNativePathName();

			/* path to source file */
			StringT path;
			path.FilePath(in.filename());
			data_file.Prepend(path);

			ifstreamT data(data_file);
			if (!data.is_open())
			{
				cout << "\n EAMFCC3D::SetGlueFunctions: could not open file: "
				     << data_file << endl;
				throw ExceptionT::kBadInputValue;
			}
			fEAM = new FBD_EAMGlue(*this, data);
			break;
		}
		case kEAMParticle:
		{
			fEAM_particle = new EAM_particle(*this);
			fEAM = NULL;	// initialize fEAM to NULL if using kEAMParticle
			break;
		}			
		
		default:
		
			cout << "\nEAMFCC3D::EAMFCC3D: unknown EAM code: " << fEAMcode << endl;
			throw ExceptionT::kBadInputValue;
	}
	
	//if (!fEAM) throw ExceptionT::kOutOfMemory;
	if (!fEAM)
	{
		StringT path;
		path.FilePath(in.filename());
		StringT data_file;
		in >> data_file;
		data_file.ToNativePathName();
		data_file.Prepend(path);
		fEAM_particle->SetGlueFunctions(data_file);
	}
	else
		fEAM->SetGlueFunctions();		
}
