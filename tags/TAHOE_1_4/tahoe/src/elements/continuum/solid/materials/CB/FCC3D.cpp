/* $Id: FCC3D.cpp,v 1.7 2004-07-07 18:03:28 paklein Exp $ */
/* created: paklein (07/01/1996) */
#include "FCC3D.h"
#include "ElementsConfig.h"
#include "FCCLatticeT.h"

#include <math.h>
#include <iostream.h>

#include "ifstreamT.h"

/* pair properties */
#ifdef PARTICLE_ELEMENT
#include "HarmonicPairT.h"
#include "LennardJonesPairT.h"
#else
#pragma message("FCC3D requires PARTICLE_ELEMENT")
#error "FCC3D requires PARTICLE_ELEMENT"
#endif

using namespace Tahoe;

/* constructor */
FCC3D::FCC3D(ifstreamT& in, const FSMatSupportT& support):
	NL_E_MatT(in, support),
	fNearestNeighbor(-1),
	fQ(3),
	fFCCLattice(NULL),
	fPairProperty(NULL),
	fAtomicVolume(0),
	fBondTensor4(dSymMatrixT::NumValues(3)),
	fBondTensor2(dSymMatrixT::NumValues(3)),
	fFullDensityForStressOutput(true)	
{
	const char caller[] = "FCC3D::FCC3D";

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
	fFCCLattice = new FCCLatticeT(fQ, nshells);
	fFCCLattice->Initialize();

	/* construct default bond density array */
	fFullDensity.Dimension(fFCCLattice->NumberOfBonds());
	fFullDensity = 1.0;

	/* check */
	if (fNearestNeighbor < kSmall)
		ExceptionT::BadInputValue(caller, "nearest bond ! (%g > 0)", fNearestNeighbor);
		
	/* compute the (approx) cell volume */
	double cube_edge = fNearestNeighbor*sqrt(2.0);
	fAtomicVolume = cube_edge*cube_edge*cube_edge/4.0;

	/* compute stress-free dilatation */
	double stretch = ZeroStressStretch();
	fNearestNeighbor *= stretch;
	cube_edge = fNearestNeighbor*sqrt(2.0);
	fAtomicVolume = cube_edge*cube_edge*cube_edge/4.0;

	/* reset the continuum density (4 atoms per unit cell) */
	fDensity = fPairProperty->Mass()/fAtomicVolume;
}

/* destructor */
FCC3D::~FCC3D(void)
{
	delete fFCCLattice;
	delete fPairProperty;
}

/* I/O functions */
void FCC3D::PrintName(ostream& out) const
{
	NL_E_MatT::PrintName(out);
	out << "    FCC lattice\n";
}

void FCC3D::Print(ostream& out) const
{
	/* inherited */
	NL_E_MatT::Print(out);

	/* higher precision */
	int prec = out.precision();
	out.precision(12);

	/* lattice parameters */
	out << " Number of neighbor shells . . . . . . . . . . . = " << fFCCLattice->NumShells() << '\n';
	out << " Number of neighbors . . . . . . . . . . . . . . = " << fFCCLattice->NumberOfBonds() << '\n';
	out << " Nearest neighbor distance . . . . . . . . . . . = " << fNearestNeighbor << '\n';

	/* write pair properties to output */
	out << " Interaction potential parameters:\n";
	fPairProperty->Write(out);

	/* restore precision */
	out.precision(prec);
}

/* return a reference to the bond lattice */
const BondLatticeT& FCC3D::BondLattice(void) const {
	if (!fFCCLattice) ExceptionT::GeneralFail("FCC3D::BondLattice", "pointer not set");
	return *fFCCLattice;
}

/*************************************************************************
 * Protected
 *************************************************************************/

void FCC3D::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{	
	fFCCLattice->ComputeDeformedLengths(E);
	const dArrayT& bond_length = fFCCLattice->DeformedLengths();

	/* fetch function pointers */
	PairPropertyT::ForceFunction force = fPairProperty->getForceFunction();
	PairPropertyT::StiffnessFunction stiffness = fPairProperty->getStiffnessFunction();

	/* bond density */
	const double* density = fFullDensity.Pointer();
	int nb = fFCCLattice->NumberOfBonds();
	const ElementCardT* element = MaterialSupport().CurrentElement();
	if (element && element->IsAllocated()) {
		const dArrayT& d_array = element->DoubleData();
		density = d_array.Pointer(CurrIP()*nb);
	}
	
	/* sum over bonds */
	moduli = 0.0; 
	double R4byV = fNearestNeighbor*fNearestNeighbor*fNearestNeighbor*fNearestNeighbor/fAtomicVolume;
	for (int i = 0; i < nb; i++)
	{
		double ri = bond_length[i]*fNearestNeighbor;
		double coeff = (*density++)*(stiffness(ri, NULL, NULL) - force(ri, NULL, NULL)/ri)/ri/ri;
		fFCCLattice->BondComponentTensor4(i, fBondTensor4);
		moduli.AddScaled(R4byV*coeff, fBondTensor4);
	}
	
	/* symmetric */
	moduli.CopySymmetric();
}

/* 2nd Piola-Kirchhoff stress vector */
void FCC3D::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	fFCCLattice->ComputeDeformedLengths(E);
	const dArrayT& bond_length = fFCCLattice->DeformedLengths();

	/* fetch function pointer */
	PairPropertyT::ForceFunction force = fPairProperty->getForceFunction();
	
	/* bond density */
	const double* density = fFullDensity.Pointer();
	int nb = fFCCLattice->NumberOfBonds();
	const ElementCardT* element = MaterialSupport().CurrentElement();
	bool keep_full_density = MaterialSupport().RunState() == GlobalT::kWriteOutput && fFullDensityForStressOutput;
	if (!keep_full_density && element && element->IsAllocated()) {
		const dArrayT& d_array = element->DoubleData();
		density = d_array.Pointer(CurrIP()*nb);
	}

	/* sum over bonds */
	PK2 = 0.0;
	double R2byV = fNearestNeighbor*fNearestNeighbor/fAtomicVolume;
	for (int i = 0; i < nb; i++)
	{
		double ri = bond_length[i]*fNearestNeighbor;
		double coeff = (*density++)*force(ri, NULL, NULL)/ri;
		fFCCLattice->BondComponentTensor2(i, fBondTensor2);
		PK2.AddScaled(R2byV*coeff, fBondTensor2);
	}
}

/* strain energy density */
double FCC3D::ComputeEnergyDensity(const dSymMatrixT& E)
{
	fFCCLattice->ComputeDeformedLengths(E);
	const dArrayT& bond_length = fFCCLattice->DeformedLengths();

	/* fetch function pointer */
	PairPropertyT::EnergyFunction energy = fPairProperty->getEnergyFunction();

	/* bond density */
	const double* density = fFullDensity.Pointer();
	int nb = fFCCLattice->NumberOfBonds();
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
	tmpSum /= fAtomicVolume;
	
	return tmpSum;
}

/* return the equitriaxial stretch at which the stress is zero */
double FCC3D::ZeroStressStretch(void)
{
	const char caller[] = "FCC3D::ZeroStress";

	int nsd = 3;
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
		double dE = -PK2(0,0)/(C(0,0) + 2.0*C(0,1));
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
