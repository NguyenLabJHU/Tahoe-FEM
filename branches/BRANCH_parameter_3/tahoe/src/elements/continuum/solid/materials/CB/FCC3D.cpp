/* $Id: FCC3D.cpp,v 1.3.12.4 2004-06-19 23:27:58 paklein Exp $ */
/* created: paklein (07/01/1996) */
#include "FCC3D.h"
#include "ElementsConfig.h"
#include "FCCLatticeT.h"

#include <math.h>
#include <iostream.h>

#include "fstreamT.h"

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
	ParameterInterfaceT("FCC_3D"),
	NL_E_MatT(in, support),
	fNearestNeighbor(-1),
	fFCCLattice(NULL),
	fPairProperty(NULL),
	fAtomicVolume(0),
	fBondTensor4(dSymMatrixT::NumValues(3)),
	fBondTensor2(dSymMatrixT::NumValues(3))	
{
#if 0
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
#endif
}

FCC3D::FCC3D(void):
	ParameterInterfaceT("FCC_3D"),
	fNearestNeighbor(-1),
	fFCCLattice(NULL),
	fPairProperty(NULL),
	fAtomicVolume(0),
	fBondTensor4(dSymMatrixT::NumValues(3)),
	fBondTensor2(dSymMatrixT::NumValues(3))	
{

}

/* destructor */
FCC3D::~FCC3D(void)
{
	delete fFCCLattice;
	delete fPairProperty;
}

/* describe the parameters needed by the interface */
void FCC3D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	NL_E_MatT::DefineParameters(list);

	/* number of neighbor shells */
	ParameterT n_shells(ParameterT::Integer, "shells");
	n_shells.AddLimit(1, LimitT::LowerInclusive);
	list.AddParameter(n_shells);
}

/* information about subordinate parameter lists */
void FCC3D::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	NL_E_MatT::DefineSubs(sub_list);

	/* FCC lattice */
	sub_list.AddSub("CB_lattice_FCC");

	/* pair potential choice */
	sub_list.AddSub("FCC_3D_potential_choice", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void FCC3D::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
	SubListT& sub_sub_list) const
{
	if (sub == "FCC_3D_potential_choice")
	{
		order = ParameterListT::Choice;

		/* choice of potentials */
		sub_sub_list.AddSub("harmonic");
		sub_sub_list.AddSub("Lennard_Jones");
		sub_sub_list.AddSub("Paradyn_pair");
		sub_sub_list.AddSub("Matsui");
	}
	else /* inherited */
		NL_E_MatT::DefineInlineSub(sub, order, sub_sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FCC3D::NewSub(const StringT& list_name) const
{
	/* try to construct pair property */
	PairPropertyT* pair_prop = PairPropertyT::New(list_name, fMaterialSupport);
	if (pair_prop)
		return pair_prop;
	else if (list_name == "CB_lattice_FCC")	
		return new FCCLatticeT(0);
	else /* inherited */
		return NL_E_MatT::NewSub(list_name);
}

/* accept parameter list */
void FCC3D::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "FCC3D::TakeParameterList";

	/* inherited */
	NL_E_MatT::TakeParameterList(list);

	/* number of shells */
	int nshells = list.GetParameter("shells");

	/* construct pair property */
	const ParameterListT& pair_prop = list.GetListChoice(*this, "FCC_3D_potential_choice");
	fPairProperty = PairPropertyT::New(pair_prop.Name(), &(MaterialSupport()));
	if (!fPairProperty) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", pair_prop.Name().Pointer());
	fPairProperty->TakeParameterList(pair_prop);
	fNearestNeighbor = fPairProperty->NearestNeighbor();

	/* construct lattice */
	fFCCLattice = new FCCLatticeT(nshells);
	fFCCLattice->TakeParameterList(list.GetList("CB_lattice_FCC"));

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
	
	moduli = 0.0; 
	int nb = fFCCLattice->NumberOfBonds();
	double R4byV = fNearestNeighbor*fNearestNeighbor*fNearestNeighbor*fNearestNeighbor/fAtomicVolume;
	for (int i = 0; i < nb; i++)
	{
		double ri = bond_length[i]*fNearestNeighbor;
		double coeff = (stiffness(ri, NULL, NULL) - force(ri, NULL, NULL)/ri)/ri/ri;
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
	
	PK2 = 0.0;
	int nb = fFCCLattice->NumberOfBonds();
	double R2byV = fNearestNeighbor*fNearestNeighbor/fAtomicVolume;
	for (int i = 0; i < nb; i++)
	{
		double ri = bond_length[i]*fNearestNeighbor;
		double coeff = force(ri, NULL, NULL)/ri;
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

	double tmpSum  = 0.;	
	int nb = fFCCLattice->NumberOfBonds();
	for (int i = 0; i < nb; i++)
	{
		double r = bond_length[i]*fNearestNeighbor;
		tmpSum += energy(r, NULL, NULL);
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
