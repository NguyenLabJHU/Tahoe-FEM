/* $Id: Hex2D.cpp,v 1.2.42.6 2004-06-19 23:27:58 paklein Exp $ */
/* created: paklein (07/01/1996) */
#include "Hex2D.h"
#include "ElementsConfig.h"
#include "HexLattice2DT.h"
#include "MaterialSupportT.h"

#include <math.h>
#include <iostream.h>

#include "fstreamT.h"

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
	ParameterInterfaceT("hex_2D"),
	NL_E_MatT(in, support),
	fNearestNeighbor(-1),
	fQ(2),
	fHexLattice2D(NULL),
	fPairProperty(NULL),
	fBondTensor4(dSymMatrixT::NumValues(2)),
	fBondTensor2(dSymMatrixT::NumValues(2))	
{
#if 0
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
#endif
}

Hex2D::Hex2D(void):
	ParameterInterfaceT("hex_2D"),
	fNearestNeighbor(-1),
	fQ(2),
	fHexLattice2D(NULL),
	fPairProperty(NULL),
	fBondTensor4(dSymMatrixT::NumValues(2)),
	fBondTensor2(dSymMatrixT::NumValues(2)),
	fCellVolume(0.0)
{

}
	
/* destructor */
Hex2D::~Hex2D(void)
{
	delete fHexLattice2D;
	delete fPairProperty;
}

/* describe the parameters needed by the interface */
void Hex2D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	NL_E_MatT::DefineParameters(list);
	
	/* 2D option must be plain stress */
	ParameterT& constraint = list.GetParameter("constraint_2D");
	constraint.SetDefault(kPlaneStress);

	/* number of neighbor shells */
	ParameterT n_shells(ParameterT::Integer, "shells");
	n_shells.AddLimit(1, LimitT::LowerInclusive);
	list.AddParameter(n_shells);
}

/* information about subordinate parameter lists */
void Hex2D::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	NL_E_MatT::DefineSubs(sub_list);

	/* pair potential choice */
	sub_list.AddSub("hex_2D_potential_choice", ParameterListT::Once, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* Hex2D::NewSub(const StringT& list_name) const
{
	/* try to construct pair property */
	PairPropertyT* pair_prop = PairPropertyT::New(list_name, fMaterialSupport);
	if (pair_prop)
		return pair_prop;
	else /* inherited */
		return NL_E_MatT::NewSub(list_name);
}

/* return the description of the given inline subordinate parameter list */
void Hex2D::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
	SubListT& sub_sub_list) const
{
	if (sub == "hex_2D_potential_choice")
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

/* accept parameter list */
void Hex2D::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "Hex2D::TakeParameterList";

	/* inherited */
	NL_E_MatT::TakeParameterList(list);

	/* number of shells */
	int nshells = list.GetParameter("shells");

	/* construct pair property */
	const ParameterListT& pair_prop = list.GetListChoice(*this, "hex_2D_potential_choice");
	fPairProperty = PairPropertyT::New(pair_prop.Name(), &(MaterialSupport()));
	if (!fPairProperty) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", pair_prop.Name().Pointer());
	fPairProperty->TakeParameterList(pair_prop);

	/* construct the bond tables */
	fHexLattice2D = new HexLattice2DT(nshells);
	fHexLattice2D->Initialize();
	fNearestNeighbor = fPairProperty->NearestNeighbor();
	
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
	double R4byV = fNearestNeighbor*fNearestNeighbor*fNearestNeighbor*fNearestNeighbor/fCellVolume;
	for (int i = 0; i < nb; i++)
	{
		double ri = bond_length[i]*fNearestNeighbor;
		double coeff = (stiffness(ri, NULL, NULL) - force(ri, NULL, NULL)/ri)/ri/ri;
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
	
	PK2 = 0.0;
	int nb = fHexLattice2D->NumberOfBonds();
	double R2byV = fNearestNeighbor*fNearestNeighbor/fCellVolume;
	for (int i = 0; i < nb; i++)
	{
		double ri = bond_length[i]*fNearestNeighbor;
		double coeff = force(ri, NULL, NULL)/ri;
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
