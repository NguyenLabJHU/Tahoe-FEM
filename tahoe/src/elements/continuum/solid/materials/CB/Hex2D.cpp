/* $Id: Hex2D.cpp,v 1.1.2.1 2003-02-19 01:12:59 paklein Exp $ */
/* created: paklein (07/01/1996) */
#include "Hex2D.h"
#include "ElementsConfig.h"
#include "HexLattice2DT.h"

#include <math.h>
#include <iostream.h>

#include "fstreamT.h"

/* pair properties */
#ifdef PARTICLE_ELEMENT
#include "HarmonicPairT.h"
#else
#pragma message("Hex2D requires PARTICLE_ELEMENT")
#endif

const double sqrt3 = sqrt(3.0);

using namespace Tahoe;

/* constructor */
Hex2D::Hex2D(ifstreamT& in, const FSMatSupportT& support):
	NL_E_Mat2DT(in, support, kPlaneStress),
	fNearestNeighbor(-1),
	fQ(2),
	fHexLattice2D(NULL),
	fPairProperty(NULL),
	fBondTensor4(dSymMatrixT::NumValues(2)),
	fBondTensor2(dSymMatrixT::NumValues(2))	
{
	const char caller[] = "Hex2D::Hex2D";

	/* read the number of shells */
	int nshells;
	in >> nshells;

	/* construct pair property */
	ParticlePropertyT::TypeT property;
	in >> property;
	switch (property)
	{
		case ParticlePropertyT::kHarmonicPair:
		{
			double mass, K;
			in >> mass >> fNearestNeighbor >> K;
			fPairProperty = new HarmonicPairT(mass, fNearestNeighbor, K);
			break;
		}
		default:
			ExceptionT::BadInputValue(caller, "unrecognized property type %d", property);
	}
	
	/* construct the bond tables */
	fQ.Identity();
	fHexLattice2D = new HexLattice2DT(fQ, nshells);
	fHexLattice2D->Initialize();
	
	/* check */
	if (fNearestNeighbor < kSmall)
		ExceptionT::BadInputValue(caller, "nearest bond ! (%g > 0)", fNearestNeighbor);
		
	/* compute the cell volume */
	fCellVolume = fNearestNeighbor*fNearestNeighbor*sqrt3/2.0;
}

/* destructor */
Hex2D::~Hex2D(void)
{
	delete fHexLattice2D;
	delete fPairProperty;
}

/* I/O functions */
void Hex2D::PrintName(ostream& out) const
{
	NL_E_Mat2DT::PrintName(out);
	out << "    2D Hexagonal lattice\n";
}

void Hex2D::Print(ostream& out) const
{
	/* inherited */
	NL_E_Mat2DT::Print(out);

	/* lattice parameters */
	out << " Number of neighbor shells . . . . . . . . . . . = " << fHexLattice2D->NumShells() << '\n';
	out << " Number of neighbors . . . . . . . . . . . . . . = " << fHexLattice2D->NumberOfBonds() << '\n';

	/* write pair properties to output */
	fPairProperty->Write(out);
}

/*************************************************************************
 * Protected
 *************************************************************************/

void Hex2D::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{	
	fHexLattice2D->ComputeDeformedLengths(E);
	const dArrayT& bond_length = fHexLattice2D->DeformedLengths();

	/* fetch function pointers */
	PairPropertyT::ForceFunction force = fPairProperty->getForceFunction();
	PairPropertyT::StiffnessFunction stiffness = fPairProperty->getStiffnessFunction();
	
	moduli = 0.0; 
	int nb = fHexLattice2D->NumberOfBonds();
	for (int i = 0; i < nb; i++)
	{
		double ri = bond_length[i]*fNearestNeighbor;
		double coeff = (stiffness(ri, NULL, NULL) - force(ri, NULL, NULL)/ri)/ri/ri;
		fHexLattice2D->BondComponentTensor4(i, fBondTensor4);
		moduli.AddScaled(coeff/fCellVolume, fBondTensor4);
	}
	
	/* symmetric */
	moduli.CopySymmetric();
}

/* 2nd Piola-Kirchhoff stress vector */
void Hex2D::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	fHexLattice2D->ComputeDeformedLengths(E);
	const dArrayT& bond_length = fHexLattice2D->DeformedLengths();

	/* fetch function pointer */
	PairPropertyT::ForceFunction force = fPairProperty->getForceFunction();
	
	PK2 = 0.0;
	int nb = fHexLattice2D->NumberOfBonds();
	for (int i = 0; i < nb; i++)
	{
		double ri = bond_length[i]*fNearestNeighbor;
		double coeff = force(ri, NULL, NULL)/ri;
		fHexLattice2D->BondComponentTensor2(i, fBondTensor2);
		PK2.AddScaled(coeff/fCellVolume, fBondTensor2);
	}
}

/* strain energy density */
double Hex2D::ComputeEnergyDensity(const dSymMatrixT& E)
{
	fHexLattice2D->ComputeDeformedLengths(E);
	const dArrayT& bond_length = fHexLattice2D->DeformedLengths();

	/* fetch function pointer */
	PairPropertyT::EnergyFunction energy = fPairProperty->getEnergyFunction();

	double tmpSum  = 0.;	
	int nb = fHexLattice2D->NumberOfBonds();
	for (int i = 0; i < nb; i++)
	{
		double r = bond_length[i]*fNearestNeighbor;
		tmpSum += energy(r, NULL, NULL);
	}
	tmpSum /= fCellVolume;
	
	return tmpSum;
}
