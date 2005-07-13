/* $Id: FCC3D_Surf.cpp,v 1.8 2005-07-13 05:26:42 hspark Exp $ */
/* created: paklein (07/01/1996) */
#include "FCC3D_Surf.h"

#include "ElementsConfig.h"
#include "FCCLatticeT_Surf.h"
#include "ContinuumElementT.h"
#include "dSymMatrixT.h"
#include "ModelManagerT.h"
#include <math.h>

/* pair properties */
#ifdef PARTICLE_ELEMENT
#include "HarmonicPairT.h"
#include "LennardJonesPairT.h"
#else
#pragma message("FCC3D_Surf requires PARTICLE_ELEMENT")
#error "FCC3D_Surf requires PARTICLE_ELEMENT"
#endif

using namespace Tahoe;

/* constructor */
FCC3D_Surf::FCC3D_Surf(void):
	ParameterInterfaceT("FCC_3D_Surf"),
	fNearestNeighbor(-1),
	fFCCLattice_Surf(NULL),
	fPairProperty(NULL),
	fAtomicVolume(0),
	fAtomicArea(0),
	fBondTensor4(dSymMatrixT::NumValues(3)),
	fBondTensor2(dSymMatrixT::NumValues(3)),
	fFullDensityForStressOutput(true)	
{

}

/* destructor */
FCC3D_Surf::~FCC3D_Surf(void)
{
	delete fFCCLattice_Surf;
	delete fPairProperty;
}

/* describe the parameters needed by the interface */
void FCC3D_Surf::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	NL_E_MatT::DefineParameters(list);

	/* number of neighbor shells */
	ParameterT n_shells(ParameterT::Integer, "shells");
	n_shells.AddLimit(1, LimitT::LowerInclusive);
	list.AddParameter(n_shells);
	
	/* surface normal */
	ParameterT normal(ParameterT::Integer, "normal_code");
	normal.AddLimit(0, LimitT::LowerInclusive);
	normal.AddLimit(5, LimitT::UpperInclusive);
	list.AddParameter(normal);
	
	/* bulk nearest neighbor distance */
	ParameterT nearest_neighbor(ParameterT::Double, "bulk_nearest_neighbor");
	nearest_neighbor.AddLimit(0, LimitT::Lower);
	list.AddParameter(nearest_neighbor);
}

/* information about subordinate parameter lists */
void FCC3D_Surf::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	NL_E_MatT::DefineSubs(sub_list);

	/* FCC lattice */
	sub_list.AddSub("CB_lattice_FCC");

	/* pair potential choice */
	sub_list.AddSub("FCC_3D_potential_choice", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void FCC3D_Surf::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "FCC_3D_potential_choice")
	{
		order = ParameterListT::Choice;

		/* choice of potentials */
		sub_lists.AddSub("harmonic");
		sub_lists.AddSub("Lennard_Jones");
		sub_lists.AddSub("Paradyn_pair");
		sub_lists.AddSub("Matsui");
	}
	else /* inherited */
		NL_E_MatT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FCC3D_Surf::NewSub(const StringT& name) const
{
	/* try to construct pair property */
	PairPropertyT* pair_prop = PairPropertyT::New(name, fMaterialSupport);
	if (pair_prop)
		return pair_prop;
	else if (name == "CB_lattice_FCC")	
		return new FCCLatticeT_Surf(0,0);
	else /* inherited */
		return NL_E_MatT::NewSub(name);
}

/* accept parameter list */
void FCC3D_Surf::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "FCC3D_Surf::TakeParameterList";

	/* inherited */
	NL_E_MatT::TakeParameterList(list);

	/* number of shells */
	int nshells = list.GetParameter("shells");

	/* construct pair property */
	const ParameterListT& pair_prop = list.GetListChoice(*this, "FCC_3D_potential_choice");
	fPairProperty = PairPropertyT::New(pair_prop.Name(), &(MaterialSupport()));
	if (!fPairProperty) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", pair_prop.Name().Pointer());
	fPairProperty->TakeParameterList(pair_prop);

	/* Obtain surface normal, use it to rotate Bond Table to correct orientation */
	/* Pass normal in as input to new FCCLatticeT_Surf(nshells,normal) */
	int normal = list.GetParameter("normal_code");

	/* construct lattice */
	fFCCLattice_Surf = new FCCLatticeT_Surf(nshells,normal);
	fFCCLattice_Surf->TakeParameterList(list.GetList("CB_lattice_FCC"));
	
	/* construct default bond density array */
	/* THIS IS AREA/VOLUME NORMALIZATION FACTOR, i.e. fFullDensity */
	fFullDensity.Dimension(fFCCLattice_Surf->NumberOfBonds());
	fFullDensity = 1.0;
		
	/* compute stress-free dilatation */
	fNearestNeighbor = list.GetParameter("bulk_nearest_neighbor");
	double cube_edge = fNearestNeighbor*sqrt(2.0);
	fAtomicVolume = cube_edge*cube_edge*cube_edge/4.0;
	fAtomicArea = .5*cube_edge*cube_edge;	// area normalization same for all surface cluster atoms

	/* reset the continuum density (4 atoms per unit cell) */
	/* DOES THIS NEED TO BE CHANGED? */
	fDensity = fPairProperty->Mass()/fAtomicVolume;
	
	/* CALL STRESS/MODULUS TO TEST VALUES - GET VALUES FOR EACH OF 6 SURFACES */
	//dSymMatrixT E(3), PK2(3);
	//dMatrixT C(6);
	//double rho;
	//E = 0.0;
	//E(0,1)=E(1,0)=0.001;
	//ComputePK2(E,PK2);
	//ComputeModuli(E,C);
	//rho = ComputeEnergyDensity(E);
}

/* return a reference to the bond lattice */
const BondLatticeT_Surf& FCC3D_Surf::BondLattice_Surf(void) const {
	if (!fFCCLattice_Surf) ExceptionT::GeneralFail("FCC3D_Surf::BondLattice_Surf", "pointer not set");
	return *fFCCLattice_Surf;
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* Need to modify AtomicVolume, force/energy to account for full bonds/split energy */
void FCC3D_Surf::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{	
	fFCCLattice_Surf->ComputeDeformedLengths(E);
	const dArrayT& bond_length = fFCCLattice_Surf->DeformedLengths();

	/* fetch function pointers */
	PairPropertyT::ForceFunction force = fPairProperty->getForceFunction();
	PairPropertyT::StiffnessFunction stiffness = fPairProperty->getStiffnessFunction();

	/* bond density */
	const double* density = fFullDensity.Pointer();
	int nb = fFCCLattice_Surf->NumberOfBonds();
	const ElementCardT* element = MaterialSupport().CurrentElement();
	if (element && element->IsAllocated()) {
		const dArrayT& d_array = element->DoubleData();
		density = d_array.Pointer(CurrIP()*nb);
	}
	
	/* sum over bonds */
	moduli = 0.0; 
	/* HAVE NOT NORMALIZED MODULUS BY AREA YET */
	double R4byV = fNearestNeighbor*fNearestNeighbor*fNearestNeighbor*fNearestNeighbor/fAtomicVolume;
	for (int i = 0; i < nb; i++)
	{
		double ri = bond_length[i]*fNearestNeighbor;
		double coeff = (*density++)*(stiffness(ri, NULL, NULL) - force(ri, NULL, NULL)/ri)/ri/ri;
		fFCCLattice_Surf->BondComponentTensor4(i, fBondTensor4);
		moduli.AddScaled(R4byV*coeff, fBondTensor4);
	}
	
	/* Multiply modulus by half due to splitting bond energies */
	moduli*=.5;
	
	/* symmetric */
	moduli.CopySymmetric();
}

/* 2nd Piola-Kirchhoff stress vector */
void FCC3D_Surf::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
#if __option(extended_errorcheck)
	if (E.Rows() != PK2.Rows() ||
	   (E.Rows() != 2 && E.Rows() != 3))
	   ExceptionT::GeneralFail("FCC3D_Surf::ComputePK2");
#endif

	/* bond density */
	const double* density = fFullDensity.Pointer();
	int nb = fFCCLattice_Surf->NumberOfBonds();
	const ElementCardT* element = MaterialSupport().CurrentElement();
	bool keep_full_density = MaterialSupport().RunState() == GlobalT::kWriteOutput && fFullDensityForStressOutput;
	if (!keep_full_density && element && element->IsAllocated()) {
		const dArrayT& d_array = element->DoubleData();
		density = d_array.Pointer(CurrIP()*nb);
	}

	/* fetch function pointer */
	PairPropertyT::ForceFunction force = fPairProperty->getForceFunction();

	/* lattice properties */
	dArrayT& bond_length = fFCCLattice_Surf->DeformedLengths();
	const dArray2DT& bonds = fFCCLattice_Surf->Bonds();
	
	/* NORMALIZE BY AREA INSTEAD OF VOLUME FOR SURFACE CB */
	double R2byV = fNearestNeighbor*fNearestNeighbor/fAtomicArea;

	double* pC = fFCCLattice_Surf->Stretch().Pointer();
	const double* pE = E.Pointer();
	double* pPK2 = PK2.Pointer();
	
	/* plane strain */
	if (E.Rows() == 2)
	{
		/* compute the stretch */
		pC[0] = 2.0*pE[0] + 1.0;
		pC[1] = 2.0*pE[1] + 1.0;
		pC[2] = 1.0;
		pC[3] = 0.0; /* 23 */
		pC[4] = 0.0; /* 13 */
		pC[5] = 2.0*pE[2]; /* 12 */	

		/* initialize */
		pPK2[0] = 0.0;
		pPK2[1] = 0.0;
		pPK2[2] = 0.0;

		/* sum over bonds */
		for (int i = 0; i < nb; i++)
		{
			/* bond vector */
			const double* R = bonds(i);
	
			/* deformed bond length */
			double& l = bond_length[i];
			l = sqrt(
				(pC[0]*R[0] + pC[5]*R[1] + pC[4]*R[2])*R[0] +
				(pC[5]*R[0] + pC[1]*R[1] + pC[3]*R[2])*R[1] +
				(pC[4]*R[0] + pC[3]*R[1] + pC[2]*R[2])*R[2]
			);

			/* accumulate */
			double ri = l*fNearestNeighbor;
			double coeff = R2byV*(*density++)*force(ri, NULL, NULL)/ri;
			pPK2[0] += coeff*R[0]*R[0];
			pPK2[1] += coeff*R[1]*R[1];
			pPK2[2] += coeff*R[0]*R[1];
		}
	}
	else /* 3D */ 
	{
		/* compute the stretch */
		pC[0] = 2.0*pE[0] + 1.0;
		pC[1] = 2.0*pE[1] + 1.0;
		pC[2] = 2.0*pE[2] + 1.0;
		pC[3] = 2.0*pE[3]; /* 23 */
		pC[4] = 2.0*pE[4]; /* 13 */
		pC[5] = 2.0*pE[5]; /* 12 */	

		/* initialize */
		pPK2[0] = 0.0;
		pPK2[1] = 0.0;
		pPK2[2] = 0.0;
		pPK2[3] = 0.0;
		pPK2[4] = 0.0;
		pPK2[5] = 0.0;

		/* sum over bonds */
		for (int i = 0; i < nb; i++)
		{
			/* bond vector */
			const double* R = bonds(i);
	
			/* deformed bond length */
			double& l = bond_length[i];
			l = sqrt(
				(pC[0]*R[0] + pC[5]*R[1] + pC[4]*R[2])*R[0] +
				(pC[5]*R[0] + pC[1]*R[1] + pC[3]*R[2])*R[1] +
				(pC[4]*R[0] + pC[3]*R[1] + pC[2]*R[2])*R[2]
			);

			/* accumulate */
			double ri = l*fNearestNeighbor;
			double coeff = R2byV*(*density++)*force(ri, NULL, NULL)/ri;
			
			/* multiply PK2 by half because of splitting bond energies */
			pPK2[0] += coeff*R[0]*R[0]*0.5;
			pPK2[1] += coeff*R[1]*R[1]*0.5;
			pPK2[2] += coeff*R[2]*R[2]*0.5;
			pPK2[3] += coeff*R[1]*R[2]*0.5;
			pPK2[4] += coeff*R[0]*R[2]*0.5;
			pPK2[5] += coeff*R[0]*R[1]*0.5;
		}
	}
}

/* strain energy density */
double FCC3D_Surf::ComputeEnergyDensity(const dSymMatrixT& E)
{
	fFCCLattice_Surf->ComputeDeformedLengths(E);
	const dArrayT& bond_length = fFCCLattice_Surf->DeformedLengths();

	/* fetch function pointer */
	PairPropertyT::EnergyFunction energy = fPairProperty->getEnergyFunction();

	/* bond density */
	const double* density = fFullDensity.Pointer();
	int nb = fFCCLattice_Surf->NumberOfBonds();
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
		/* Split energy by half since counting all bonds unlike bulk */
		tmpSum += (*density++)*energy(r, NULL, NULL)*0.5;
	}
	/* MODIFIED FOR SURFACE CB CALCULATIONS, I.E. NORMALIZE BY AREA */
	tmpSum /= fAtomicArea;
	
	return tmpSum;
}

/* return the equitriaxial stretch at which the stress is zero */
double FCC3D_Surf::ZeroStressStretch(void)
{
	const char caller[] = "FCC3D_Surf::ZeroStress";

	int nsd = 3;
	dSymMatrixT E(nsd), PK2(nsd);
	dMatrixT C(dSymMatrixT::NumValues(nsd));

	E = 0.0;
	ComputePK2(E, PK2);
	
	/* Newton iteration */
	int count = 0;
	double error, error0;
	error = error0 = fabs(PK2(0,0));
	while (count++ < 10 && error > kSmall && error/error0 > kSmall)
	{
		ComputeModuli(E, C);
		double dE = -PK2(0,0)/(C(0,0) + C(0,1) + C(0,2));
		E.PlusIdentity(dE);
		
		ComputePK2(E, PK2);
		error = fabs(PK2(0,0));
	}

	/* check convergence */
	if (error > kSmall && error/error0 > kSmall) {
		cout << "\n " << caller << ":\n";
		cout << " E =\n" << E << '\n';
		cout << " PK2 =\n" << PK2 << endl;
		ExceptionT::GeneralFail(caller, "failed to find stress-free state");
	}
	
	return sqrt(2.0*E(0,0) + 1.0);
}
