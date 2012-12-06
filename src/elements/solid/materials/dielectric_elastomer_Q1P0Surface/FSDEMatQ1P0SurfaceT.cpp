#include "ExceptionT.h"
#include "FSDEMatQ1P0SurfaceT.h"
#include "FSDE_incQ1P0Surface.h"

/* FCC3D_Surf */
#include "ElementsConfig.h"
#include "FCCLatticeT_Q1P0Surf.h"
#include "ContinuumElementT.h"
#include "dSymMatrixT.h"
#include "ModelManagerT.h"
#include <cmath>

/* pair properties */
#ifdef PARTICLE_ELEMENT
#include "HarmonicPairT.h"
#include "LennardJonesPairT.h"
#else
#pragma message("FSDEMatQ1P0SurfaceT requires PARTICLE_ELEMENT")
#error "FSDEMatQ1P0SurfaceT requires PARTICLE_ELEMENT"
#endif

namespace Tahoe {

  FSDEMatQ1P0SurfaceT::FSDEMatQ1P0SurfaceT() :
    ParameterInterfaceT("Dielectric_Elastomer_Q1P0Surface"),
        fFSDEMatSupportQ1P0Surface(0),
	/* FCC3D_Surf */
	fNearestNeighbor(-1),
	fSurfaceThickness(-1),
	fFCCLattice_Q1P0Surf(NULL),
	fPairProperty(NULL),
	fAtomicVolume(0),
	fAtomicArea(0),
	fBondTensor4(dSymMatrixT::NumValues(3)),
	fBondTensor2(dSymMatrixT::NumValues(3)),
	fFullDensityForStressOutput(true)		        
  {
    SetName(FSDEMatQ1P0SurfaceT::Name);
    Initialize();
  }

FSDEMatQ1P0SurfaceT::~FSDEMatQ1P0SurfaceT(void)
{
	delete fFCCLattice_Q1P0Surf;
	delete fPairProperty;	
}

  //
  static const char DE[] = "Dielectric_Elastomer_Q1P0Surface";
  const char* FSDEMatQ1P0SurfaceT::Name = DE;
  const int kNSD       = 3;
  const int kNumDOF    = 3;
  const int kStressDim =dSymMatrixT::NumValues(kNSD);

  void FSDEMatQ1P0SurfaceT::Initialize()
  {
    fMu = 0.0;
    fElectricPermittivity = 0.0;
    fNrig = 0.0;
    fLambda = 0.0;
  }

  //
  void FSDEMatQ1P0SurfaceT::DefineParameters(ParameterListT& list) const
  {
	FSSolidMatT::DefineParameters(list);
	
	list.AddParameter(fMu, "mu");
	list.AddParameter(fElectricPermittivity, "epsilon");
 	list.AddParameter(fNrig, "Nrig");
 	list.AddParameter(fLambda, "lambda");
 	
 	/* FCC3D_Surf */
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

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSDEMatQ1P0SurfaceT::NewSub(const StringT& name) const
{
	/* try to construct pair property */
	PairPropertyT* pair_prop = PairPropertyT::New(name, fMaterialSupport);
	if (pair_prop)
		return pair_prop;
	else if (name == "CB_lattice_FCC")	
		return new FCCLatticeT_Q1P0Surf(0,0);
	else /* inherited */
		return FSSolidMatT::NewSub(name);
}

  //
  void FSDEMatQ1P0SurfaceT::TakeParameterList(const ParameterListT& list)
  {
	FSSolidMatT::TakeParameterList(list);
	fMu = list.GetParameter("mu");
	fElectricPermittivity = list.GetParameter("epsilon");
 	fNrig = list.GetParameter("Nrig");
 	fLambda = list.GetParameter("lambda");

	/* write into vector to pass to C code for stress/modulus calculations */
	fParams.Dimension(4);
	fParams[0] = fMu;
	fParams[1] = fLambda;
 	fParams[2] = fElectricPermittivity;
 	fParams[3] = fNrig;
	
	/* dimension work space */
	fTangentMechanical.Dimension(kStressDim);
	fTangentMechanicalElec.Dimension(kStressDim);
	fStress.Dimension(kNumDOF);
	fTangentElectrical.Dimension(kNumDOF);
	fTangentElectromechanical.Dimension(kStressDim, kNumDOF);
	fTangentElectromechanicalSpatial.Dimension(kStressDim, kNumDOF);	
	fElectricDisplacement.Dimension(kNumDOF);
	fElectricField.Dimension(kNumDOF);
	
	/* FCC3D_Surf */
	/* number of shells */
	int nshells = list.GetParameter("shells");

	/* construct pair property */
	const ParameterListT& pair_prop = list.GetListChoice(*this, "FCC_3D_potential_choice");
	fPairProperty = PairPropertyT::New(pair_prop.Name(), &(MaterialSupport()));
//	if (!fPairProperty) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", pair_prop.Name().Pointer());
	fPairProperty->TakeParameterList(pair_prop);

	/* Obtain surface normal, use it to rotate Bond Table to correct orientation */
	/* Pass normal in as input to new FCCLatticeT_Surf(nshells,normal) */
	int normal = list.GetParameter("normal_code");

	/* construct lattice */
	fFCCLattice_Q1P0Surf = new FCCLatticeT_Q1P0Surf(nshells,normal);
	fFCCLattice_Q1P0Surf->TakeParameterList(list.GetList("CB_lattice_FCC"));
	
	/* construct default bond density array */
	/* THIS IS AREA/VOLUME NORMALIZATION FACTOR, i.e. fFullDensity */
	fFullDensity.Dimension(fFCCLattice_Q1P0Surf->NumberOfBonds());
	fFullDensity = 1.0;
		
	/* compute stress-free dilatation */
	fNearestNeighbor = list.GetParameter("bulk_nearest_neighbor");
	double cube_edge = fNearestNeighbor*sqrt(2.0);
	fAtomicVolume = cube_edge*cube_edge*cube_edge/4.0;
	fAtomicArea = .5*cube_edge*cube_edge;	// area normalization same for all surface cluster atoms

	/* set surface thickness - should be right */
	fSurfaceThickness = 0.0;

	/* reset the continuum density (4 atoms per unit cell) */
	/* DOES THIS NEED TO BE CHANGED? */
	fDensity = fPairProperty->Mass()/fAtomicVolume;		
  }

/* return the description of the given inline subordinate parameter list */
void FSDEMatQ1P0SurfaceT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
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
		FSSolidMatT::DefineInlineSub(name, order, sub_lists);
}

  // information about subordinate parameter lists
  void FSDEMatQ1P0SurfaceT::DefineSubs(SubListT& sub_list) const
  {
	FSSolidMatT::DefineSubs(sub_list);
	
	/* FCC3D_Surf */
	/* FCC lattice */
	sub_list.AddSub("CB_lattice_FCC");

	/* pair potential choice */
	sub_list.AddSub("FCC_3D_potential_choice", ParameterListT::Once, true);
			
    return;
  }

  // Set electrical permittivity
   void FSDEMatQ1P0SurfaceT::SetElectricPermittivity(double epsilon)
  {	
    fElectricPermittivity = epsilon;
  }

  // Get electrical permittivity
   double FSDEMatQ1P0SurfaceT::GetElectricPermittivity() const
  {
    return fElectricPermittivity;
  }

  //
   void FSDEMatQ1P0SurfaceT::SetFSDEMatSupportQ1P0Surface(
      const FSDEMatSupportQ1P0SurfaceT* support)
  {
    fFSDEMatSupportQ1P0Surface = support;
  }

  //
   const dArrayT FSDEMatQ1P0SurfaceT::ElectricField()
  {
    fElectricField = fFSDEMatSupportQ1P0Surface->ElectricField();
    return fElectricField;
  }

  //
   const dArrayT FSDEMatQ1P0SurfaceT::ElectricField(int ip)
  {
    fElectricField = fFSDEMatSupportQ1P0Surface->ElectricField(ip);
    return fElectricField;
  }

  //
   const dMatrixT FSDEMatQ1P0SurfaceT::RightCauchyGreenDeformation()
  {
    const dMatrixT& F = F_mechanical();
    dMatrixT FTF(3);
    FTF.MultATB(F, F);

    return FTF;
  }

  // material energy density
   double FSDEMatQ1P0SurfaceT::StrainEnergyDensity()
  {
	  return 0.0;
  }

  // material mechanical tangent modulus
   const dMatrixT&
  FSDEMatQ1P0SurfaceT::C_IJKL()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
    const dArrayT& E = ElectricField();
    const dMatrixT& F = F_mechanical();

	double I1 = C(0,0)+C(1,1)+C(2,2);
	double J = C.Det();
	J = sqrt(J);
	
	/* call C function for mechanical part of tangent modulus */
// 	mech_tanmod_Q1P0Surface(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, I1, fTangentMechanical.Pointer()); 
// 	me_tanmod_Q1P0Surface(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, fTangentMechanicalElec.Pointer());
 	fTangentMechanical+=fTangentMechanicalElec;
    return fTangentMechanical;
  }

  // Second Piola-Kirchhoff stress (mechanical)
   const dSymMatrixT&
  FSDEMatQ1P0SurfaceT::S_IJ()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();
    const dMatrixT& F = F_mechanical();
	
	dMatrixT stress_temp(3);
	dMatrixT stress_temp2(3);
	double I1 = C(0,0)+C(1,1)+C(2,2);
	double J = C.Det();
	J = sqrt(J);
	
	/* call C function for mechanical part of PK2 stress */
// 	mech_pk2_Q1P0Surface(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, I1, stress_temp.Pointer()); 
//	me_pk2_Q1P0Surface(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, stress_temp2.Pointer());
	stress_temp+=stress_temp2;
	
	fStress.FromMatrix(stress_temp);
    return fStress;
  }

  // material electromechanical tangent modulus
   const dMatrixT&
  FSDEMatQ1P0SurfaceT::E_IJK()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();
    const dMatrixT& F = F_mechanical();	
	double J = C.Det();
	J = sqrt(J);

	/* call C function for electromechanical tangent modulus */
// 	me_mixedmodulus_Q1P0Surface(fParams.Pointer(), E.Pointer(),  
// 		C.Pointer(), F.Pointer(), J, fTangentElectromechanical.Pointer()); 
 
    return fTangentElectromechanical;
  }

  // spatial electromechanical tangent modulus
   const dMatrixT&
  FSDEMatQ1P0SurfaceT::e_ijk()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();
	double J = C.Det();
	J = sqrt(J);  
    const dMatrixT& F = F_mechanical();
	
	/* call C function for (spatial) electromechanical tangent modulus */
// 	me_mixedmodulus_Q1P0spatialSurface(fParams.Pointer(), E.Pointer(),  
// 		C.Pointer(), F.Pointer(), J, fTangentElectromechanicalSpatial.Pointer()); 
 	
 	fTangentElectromechanicalSpatial /= J;
    return fTangentElectromechanicalSpatial;
  }

  // material electric tangent modulus
   const dMatrixT&
  FSDEMatQ1P0SurfaceT::B_IJ()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	double J = C.Det();
	J = sqrt(J);

	dMatrixT Cinv(3);
	Cinv.Inverse(C);
	fTangentElectrical = Cinv;
	fTangentElectrical *= fElectricPermittivity;
	fTangentElectrical *= J;
    return fTangentElectrical;
  }

  // spatial electric tangent modulus
   const dMatrixT&
  FSDEMatQ1P0SurfaceT::b_ij()
  {
    const dMatrixT& F = F_mechanical();
    const dMatrixT& C = RightCauchyGreenDeformation();    
    const double J = F.Det();
	
	// repeat B_IJ first, then push forward
	dMatrixT Cinv(3), bij(3);
	Cinv.Inverse(C);
	bij = Cinv;
	bij *= fElectricPermittivity;
	bij *= J;
	
    // prevent aliasing
//    const dMatrixT b = B_IJ();
    fTangentElectrical.MultABCT(F, bij, F);
    fTangentElectrical /= J;
    return fTangentElectrical;
  }

  // Electric displacement 
   const dArrayT&
  FSDEMatQ1P0SurfaceT::D_I()
  {
  	const dMatrixT& C = RightCauchyGreenDeformation();
  	const dArrayT& E = ElectricField();
    const dMatrixT& F = F_mechanical();  	
  	
  	double J = C.Det();
  	J = sqrt(J);
  
	/* call C function for electric stress (i.e. electric displacement D_{I}) */
// 	elec_pk2_Q1P0Surface(fParams.Pointer(), E.Pointer(),  
// 		C.Pointer(), F.Pointer(), J, fElectricDisplacement.Pointer()); 
 		
  	return fElectricDisplacement;
  }

  // spatial electric tangent modulus
   const dArrayT&
  FSDEMatQ1P0SurfaceT::d_i()
  {
    const dMatrixT& F = F_mechanical();
    const double J = F.Det();
	
    // prevent aliasing
    const dArrayT D = D_I();
	F.Multx(D, fElectricDisplacement);
	fElectricDisplacement /= J;
    return fElectricDisplacement;		// need to divide by J
  }

  // Electric field
   const dArrayT&
  FSDEMatQ1P0SurfaceT::E_I()
  {
    return fElectricField;
  }

  // spatial tangent modulus
   const dMatrixT&
  FSDEMatQ1P0SurfaceT::c_ijkl()
  {
//     const dMatrixT& F = F_mechanical();
//     const double J = F.Det();
// 
//     // prevent aliasing
//     const dMatrixT CIJKL = C_IJKL();
//     fTangentMechanical.SetToScaled(1.0 / J, PushForward(F, CIJKL));

	fTangentMechanical = FSSolidMatT::c_ijkl();

    return fTangentMechanical;

  }

  // Cauchy stress
   const dSymMatrixT&
  FSDEMatQ1P0SurfaceT::s_ij()
  {
    const dMatrixT& F = F_mechanical();
    const double J = F.Det();
	
    // prevent aliasing
    const dSymMatrixT S = S_IJ();

    fStress.SetToScaled(1.0 / J, PushForward(F, S));
    return fStress;
  }

  // pressure associated with the last computed stress
   double FSDEMatQ1P0SurfaceT::Pressure() const
  {

    return 0.0;

  }

/* FCC3D_Surf */
/* Need to modify AtomicVolume, force/energy to account for full bonds/split energy */
void FSDEMatQ1P0SurfaceT::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{	
	fFCCLattice_Q1P0Surf->ComputeDeformedLengths(E);
	const dArrayT& bond_length = fFCCLattice_Q1P0Surf->DeformedLengths();
	const dArray2DT& bonds = fFCCLattice_Q1P0Surf->Bonds();
	
	/* fetch function pointers */
	PairPropertyT::ForceFunction force = fPairProperty->getForceFunction();
	PairPropertyT::StiffnessFunction stiffness = fPairProperty->getStiffnessFunction();

	/* bond density */
	const double* density = fFullDensity.Pointer();
	int nb = fFCCLattice_Q1P0Surf->NumberOfBonds();
	const ElementCardT* element = MaterialSupport().CurrentElement();
	if (element && element->IsAllocated()) {
		const dArrayT& d_array = element->DoubleData();
		density = d_array.Pointer(CurrIP()*nb);
	}
	
	/* sum over bonds */
	moduli = 0.0; 
	
	/* Normalize modulus by area instead of volume */
	double R4byV = fNearestNeighbor*fNearestNeighbor*fNearestNeighbor*fNearestNeighbor/fAtomicArea;
	for (int i = 0; i < nb; i++)
	{
		double ri = bond_length[i]*fNearestNeighbor;
		double coeff = (*density++)*(stiffness(ri, NULL, NULL) - force(ri, NULL, NULL)/ri)/ri/ri;
		fFCCLattice_Q1P0Surf->BondComponentTensor4(i, fBondTensor4);
		moduli.AddScaled(R4byV*coeff, fBondTensor4);
	}
	
	/* Multiply modulus by half due to splitting bond energies */
	moduli*=.5;
	
	/* symmetric */
	moduli.CopySymmetric();
}

/* 2nd Piola-Kirchhoff stress vector */
void FSDEMatQ1P0SurfaceT::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
#if __option(extended_errorcheck)
	if (E.Rows() != PK2.Rows() ||
	   (E.Rows() != 2 && E.Rows() != 3))
	   ExceptionT::GeneralFail("FSDEMatQ1P0SurfaceT::ComputePK2");
#endif

	/* bond density */
	const double* density = fFullDensity.Pointer();
	int nb = fFCCLattice_Q1P0Surf->NumberOfBonds();
	const ElementCardT* element = MaterialSupport().CurrentElement();
	bool keep_full_density = MaterialSupport().RunState() == GlobalT::kWriteOutput && fFullDensityForStressOutput;
	if (!keep_full_density && element && element->IsAllocated()) {
		const dArrayT& d_array = element->DoubleData();
		density = d_array.Pointer(CurrIP()*nb);
	}

	/* fetch function pointer */
	PairPropertyT::ForceFunction force = fPairProperty->getForceFunction();

	/* lattice properties */
	dArrayT& bond_length = fFCCLattice_Q1P0Surf->DeformedLengths();
	const dArray2DT& bonds = fFCCLattice_Q1P0Surf->Bonds();
	
	/* NORMALIZE BY AREA INSTEAD OF VOLUME FOR SURFACE CB */
	double R2byV = fNearestNeighbor*fNearestNeighbor/fAtomicArea;

	double* pC = fFCCLattice_Q1P0Surf->Stretch().Pointer();
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
			/* Split by 1/2 for surface */
			double ri = l*fNearestNeighbor;
			double coeff = R2byV*(*density++)*force(ri, NULL, NULL)/ri;
			pPK2[0] += coeff*R[0]*R[0]*0.5;
			pPK2[1] += coeff*R[1]*R[1]*0.5;
			pPK2[2] += coeff*R[0]*R[1]*0.5;
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

/* return the equitriaxial stretch at which the stress is zero */
double FSDEMatQ1P0SurfaceT::ZeroStressStretch(void)
{
	const char caller[] = "FSDEMatQ1P0SurfaceT::ZeroStress";

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

} //namespace Tahoe
