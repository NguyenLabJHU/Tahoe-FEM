/* $Id: Hex2D.cpp,v 1.4 2004-06-26 06:02:24 paklein Exp $ */
/* created: paklein (07/01/1996) */
#include "Hex2D.h"
#include "ElementsConfig.h"
#include "HexLattice2DT.h"

#include <math.h>
#include <iostream.h>

#include "ifstreamT.h"

/* pair properties */
#ifdef PARTICLE_ELEMENT
#include "HarmonicPairT.h"
#include "LennardJonesPairT.h"
#else
#pragma message("Hex2D requires PARTICLE_ELEMENT")
#error "Hex2D requires PARTICLE_ELEMENT"
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
	fBondTensor2(dSymMatrixT::NumValues(2)),
	fFullDensityForStressOutput(true)
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
		case ParticlePropertyT::kLennardJonesPair:
		{
			double mass, eps, sigma, alpha;
			in >> mass >> eps >> sigma >> alpha;
			fPairProperty = new LennardJonesPairT(mass, eps, sigma, alpha);

			/* equilibrium length of a single unmodified LJ bond */
			fNearestNeighbor = pow(2.0,1.0/6.0)*sigma;
			break;
		}
		default:
			ExceptionT::BadInputValue(caller, "unrecognized property type %d", property);
	}
	
	/* construct the bond tables */
	fQ.Identity();
	fHexLattice2D = new HexLattice2DT(fQ, nshells);
	fHexLattice2D->Initialize();

	/* construct default bond density array */
	fFullDensity.Dimension(fHexLattice2D->NumberOfBonds());
	fFullDensity = 1.0;
	
	/* check */
	if (fNearestNeighbor < kSmall)
		ExceptionT::BadInputValue(caller, "nearest bond ! (%g > 0)", fNearestNeighbor);
		
	/* compute the cell volume */
	fCellVolume = fNearestNeighbor*fNearestNeighbor*sqrt3/2.0;

	/* compute stress-free dilatation */
	double stretch = ZeroStressStretch();
	fNearestNeighbor *= stretch;
	fCellVolume = fNearestNeighbor*fNearestNeighbor*sqrt3/2.0;
	
	/* reset the continuum density (2 atoms per unit cell) */
	fDensity = 2*fPairProperty->Mass()/fCellVolume;
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

	/* higher precision */
	int prec = out.precision();
	out.precision(12);

	/* lattice parameters */
	out << " Number of neighbor shells . . . . . . . . . . . = " << fHexLattice2D->NumShells() << '\n';
	out << " Number of neighbors . . . . . . . . . . . . . . = " << fHexLattice2D->NumberOfBonds() << '\n';
	out << " Nearest neighbor distance . . . . . . . . . . . = " << fNearestNeighbor << '\n';

	/* write pair properties to output */
	out << " Interaction potential parameters:\n";
	fPairProperty->Write(out);

	/* restore precision */
	out.precision(prec);
}

/* return a reference to the bond lattice */
const BondLatticeT& Hex2D::BondLattice(void) const {
	if (!fHexLattice2D) ExceptionT::GeneralFail("Hex2D::BondLattice", "pointer not set");
	return *fHexLattice2D;
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

	/* bond density */
	const double* density = fFullDensity.Pointer();
	int nb = fHexLattice2D->NumberOfBonds();
	const ElementCardT* element = MaterialSupport().CurrentElement();
	if (element && element->IsAllocated()) {
		const dArrayT& d_array = element->DoubleData();
		density = d_array.Pointer(CurrIP()*nb);
	}
	
	/* sum over bonds */
	moduli = 0.0; 
	double R4byV = fNearestNeighbor*fNearestNeighbor*fNearestNeighbor*fNearestNeighbor/fCellVolume;
	for (int i = 0; i < nb; i++) 
	{
		double ri = bond_length[i]*fNearestNeighbor;
		double coeff = (*density++)*(stiffness(ri, NULL, NULL) - force(ri, NULL, NULL)/ri)/ri/ri;
		fHexLattice2D->BondComponentTensor4(i, fBondTensor4);
		moduli.AddScaled(R4byV*coeff, fBondTensor4);
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

	/* bond density */
	const double* density = fFullDensity.Pointer();
	int nb = fHexLattice2D->NumberOfBonds();
	const ElementCardT* element = MaterialSupport().CurrentElement();
	bool keep_full_density = MaterialSupport().RunState() == GlobalT::kWriteOutput && fFullDensityForStressOutput;
	if (!keep_full_density && element && element->IsAllocated()) {
		const dArrayT& d_array = element->DoubleData();
		density = d_array.Pointer(CurrIP()*nb);
	}
	
	/* sum over bonds */
	PK2 = 0.0;
	double R2byV = fNearestNeighbor*fNearestNeighbor/fCellVolume;
	for (int i = 0; i < nb; i++)
	{
		double ri = bond_length[i]*fNearestNeighbor;
		double coeff = (*density++)*force(ri, NULL, NULL)/ri;
		fHexLattice2D->BondComponentTensor2(i, fBondTensor2);
		PK2.AddScaled(R2byV*coeff, fBondTensor2);
	}
}

/* strain energy density */
double Hex2D::ComputeEnergyDensity(const dSymMatrixT& E)
{
	fHexLattice2D->ComputeDeformedLengths(E);
	const dArrayT& bond_length = fHexLattice2D->DeformedLengths();

	/* fetch function pointer */
	PairPropertyT::EnergyFunction energy = fPairProperty->getEnergyFunction();

	/* bond density */
	const double* density = fFullDensity.Pointer();
	int nb = fHexLattice2D->NumberOfBonds();
	const ElementCardT* element = MaterialSupport().CurrentElement();
	if (element && element->IsAllocated()) {
		const dArrayT& d_array = element->DoubleData();
		density = d_array.Pointer(CurrIP()*nb);
	}

	/* sum over bonds */
	double tmpSum  = 0.;	
	for (int i = 0; i < nb; i++) 
	{
		double r = bond_length[i]*fNearestNeighbor;
		tmpSum += (*density++)*energy(r, NULL, NULL);
	}
	tmpSum /= fCellVolume;
	
	return tmpSum;
}

/* return the equitriaxial stretch at which the stress is zero */
double Hex2D::ZeroStressStretch(void)
{
	const char caller[] = "Hex2D::ZeroStress";

	int nsd = 2;
	dSymMatrixT E(nsd), PK2(nsd);
	dMatrixT C(dSymMatrixT::NumValues(nsd));

	E = 0.0;
	ComputePK2(E, PK2);
	
	/* Newton iteration */
	int count = 0;
	double error, error0;
	error = error0 = fabs(PK2(0,0));
	while (count++ < 10 && error0 > kSmall && error/error0 > kSmall)
	{
		ComputeModuli(E, C);
		double dE = -PK2(0,0)/(C(0,0) + C(0,1));
		E.PlusIdentity(dE);
		
		ComputePK2(E, PK2);
		error = fabs(PK2(0,0));
	}

	/* check convergence */
	if (error0 > kSmall && error/error0 > kSmall) {
		cout << "\n " << caller << ":\n";
		cout << " E =\n" << E << '\n';
		cout << " PK2 =\n" << PK2 << endl;
		ExceptionT::GeneralFail(caller, "failed to find stress-free state");
	}
	
	return sqrt(2.0*E(0,0) + 1.0);
}
