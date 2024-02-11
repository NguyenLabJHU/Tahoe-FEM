// ISSUES: (1) Where call Initialize() -> in constructor!
// (2) 

#include "ExceptionT.h"
#include "FSDEMatViscoT.h"
#include "FSDEVisco_inc.h"
// BELOW HEADERS FROM RGSplitT2.cpp
#include "ParameterContainerT.h"
#include "ifstreamT.h"
#include "ExceptionT.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

#include "MooneyRivlin.h"
#include "NeoHookean.h"
#include "VWPotentialT.h"
#include "FungPotentialT.h"
#include "ArrudaBoyce.h"

#include "LinearExponentialT.h"
#ifdef __DEVELOPMENT__
#include "ScaledCsch.h"
#endif

namespace Tahoe {

  FSDEMatViscoT::FSDEMatViscoT():
    ParameterInterfaceT("Finite Strain Dielectric Elastomer Viscoelastic"),
	fFSDEMatSupportVisco(0), 
    fSpectralDecompRef(3),
    fSpectralDecompSpat(3)
  {
    SetName(FSDEMatViscoT::Name);
    Initialize();	// MAY CALL LATER IN TAKEPARAMETERLIST?
	/*set default value*/ 
	/*overide in derived element classes before calling *
	 *RGViscoelasticity::TakeParameterLis               */  
  }

  //
  static const char DE[] = "Dielectric_Elastomer_Viscoelastic";
  const char* FSDEMatViscoT::Name = DE;
  const int kNSD       = 3;
  const int kNumDOF    = 3;
  const int kStressDim =dSymMatrixT::NumValues(kNSD);
  const double third = 1.0/3.0;   

  //
  void FSDEMatViscoT::Initialize()
  {
    fMu = 0.0;
    fElectricPermittivity = 0.0;
    fNrig = 0.0;
    fLambda = 0.0;
    
 	/* dimension work space - from RGSplitT2.cpp */
//  	fF_M.Dimension(NumSD());
//   	fF_T_inv.Dimension(NumSD());

//   	fF3D.Dimension(3);
//   	fInverse.Dimension(3);
// 
//   	fb.Dimension(3);
//   	fbe.Dimension(3);
//   	fb_tr.Dimension(3);
// 
//   	fEigs_dev.Dimension(3);
//   	fEigs.Dimension(3);
//   	fEigs_e.Dimension(3);
//   	fEigs_tr.Dimension(3);
// 
//   	ftau_EQ.Dimension(3);
//   	ftau_NEQ.Dimension(3);
// 
//   	fStressNEQ.Dimension(NumSD());
//   	fStress3D.Dimension(3);
// 
//   	fDtauDe_EQ.Dimension(3);
//   	fDtauDe_NEQ.Dimension(3);
// 
//   	fiKAB.Dimension(3);
//   	fCalg.Dimension(3);
// 
//   	fModulus3D.Dimension(6);
//   	fModMat.Dimension(6);
//   	fModulusNEQ.Dimension(dSymMatrixT::NumValues(NumSD()));    
  }

  //
  void FSDEMatViscoT::DefineParameters(ParameterListT& list) const
  {
	FSSolidMatT::DefineParameters(list);
	
	list.AddParameter(fMu, "mu");
	list.AddParameter(fElectricPermittivity, "epsilon");
 	list.AddParameter(fNrig, "Nrig");
 	list.AddParameter(fLambda, "lambda");
  }

  //
  void FSDEMatViscoT::TakeParameterList(const ParameterListT& list)
  {
  	/* inherited */
	FSSolidMatT::TakeParameterList(list);
	fMu = list.GetParameter("mu");
	fElectricPermittivity = list.GetParameter("epsilon");
 	fNrig = list.GetParameter("Nrig");
 	fLambda = list.GetParameter("lambda");

	/* write into vector to pass to C code for stress/modulus calculations */
	fParams.Dimension(3);
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
	fElectricDisplacement.Dimension(kNumDOF);
	fElectricField.Dimension(kNumDOF);
		
	/* Code from RGSplitT2.cpp */
	StringT caller = "FSDEMatViscoT::TakeParameterList";
	int num_neq =  list.NumLists("rg_neq_potential");
	int num_shear_visc = list.NumLists("rg_shear_viscosity");
	int num_bulk_visc = list.NumLists("rg_bulk_viscosity");
	if (num_neq != num_shear_visc || num_neq != num_bulk_visc)
		ExceptionT::GeneralFail("FSDEMatViscoT::TakeParameterList", 
			"number of matrix viscosity functions does not match number of matrix nonequilibrium potentials");
	fNumProcess = list.NumLists("rg_shear_viscosity");

	/* inherited from RGViscoelasticityT.cpp */
	/* Default Dimension state variable arrays - from RGViscoelasticityT.cpp */
	if (fNumProcess > 0) SetStateVariables(fNumProcess);	

	fPot.Dimension(fNumProcess+1);
	fVisc_s.Dimension(fNumProcess);
	fVisc_b.Dimension(fNumProcess);

	/* For some reason, dimensioning these in TakeParameterList works better */
 	fF_M.Dimension(NumSD());
  	fF_T_inv.Dimension(NumSD());
  	
  	fF3D.Dimension(3);
  	fInverse.Dimension(3);

  	fb.Dimension(3);
  	fbe.Dimension(3);
  	fb_tr.Dimension(3);

  	fEigs_dev.Dimension(3);
  	fEigs.Dimension(3);
  	fEigs_e.Dimension(3);
  	fEigs_tr.Dimension(3);

  	ftau_EQ.Dimension(3);
  	ftau_NEQ.Dimension(3);

  	fStressNEQ.Dimension(NumSD());
  	fStress3D.Dimension(3);

  	fDtauDe_EQ.Dimension(3);
  	fDtauDe_NEQ.Dimension(3);

  	fiKAB.Dimension(3);
  	fCalg.Dimension(3);

  	fModulus3D.Dimension(6);
  	fModMat.Dimension(6);
  	fModulusNEQ.Dimension(dSymMatrixT::NumValues(NumSD()));  
	
	const ParameterListT& eq_pot = list.GetListChoice(*this, "rg_eq_potential");
	if(eq_pot.Name() == "neo-hookean")
		fPot[0] = new NeoHookean;
	else if(eq_pot.Name() == "mooney-rivlin")
		fPot[0] = new MooneyRivlin;
	else if(eq_pot.Name() == "veronda-westmann")
		fPot[0] = new VWPotentialT;
	else if(eq_pot.Name() == "fung-potential")
		fPot[0] = new FungPotentialT;
	else if(eq_pot.Name() == "arruda-boyce")
		fPot[0] = new ArrudaBoyce;
	else 
		ExceptionT::GeneralFail(caller, "no such potential");
	if (!fPot[0]) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", eq_pot.Name().Pointer());			
	fPot[0]->TakeParameterList(eq_pot);
	

	for (int i = 0; i < fNumProcess; i++)
	{
		const ParameterListT& pot_neq = list.GetListChoice(*this, "rg_neq_potential",i);
		if(pot_neq.Name() == "mooney-rivlin")
			fPot[i+1] = new MooneyRivlin;
		else if(pot_neq.Name() == "neo-hookean")
			fPot[i+1] = new NeoHookean;
		else if(pot_neq.Name() == "veronda-westmann")
			fPot[i+1] = new VWPotentialT;
		else if(pot_neq.Name() == "fung-potential")
			fPot[i+1] = new FungPotentialT;
		else if(pot_neq.Name() == "arruda-boyce")
			fPot[i+1] = new ArrudaBoyce;
		else 
			ExceptionT::GeneralFail(caller, "no such potential");
		if (!fPot[i+1]) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", pot_neq.Name().Pointer());			
		fPot[i+1]->TakeParameterList(pot_neq);

		const ParameterListT& shear_visc = list.GetListChoice(*this, "rg_shear_viscosity", i);
		if (shear_visc.Name() == "linear_exponential")
			fVisc_s[i] = new LinearExponentialT;

#ifdef __DEVELOPMENT__
		else if (shear_visc.Name() == "scaled-csch")
			fVisc_s[i] = new ScaledCsch;
#endif

		else 
			ExceptionT::GeneralFail(caller, "no such potential");
		if (!fVisc_s[i]) throw ExceptionT::kOutOfMemory;
		fVisc_s[i]->TakeParameterList(shear_visc);

		const ParameterListT& bulk_visc = list.GetListChoice(*this, "rg_bulk_viscosity", i);
		if (bulk_visc.Name() == "linear_exponential")
			fVisc_b[i] = new LinearExponentialT;

#ifdef __DEVELOPMENT__
		else if (bulk_visc.Name() == "scaled-csch")
			fVisc_b[i] = new ScaledCsch;
#endif

		else 
			ExceptionT::GeneralFail(caller, "no such potential");
		if (!fVisc_b[i]) throw ExceptionT::kOutOfMemory;
		fVisc_b[i]->TakeParameterList(bulk_visc);
	}
	
	/* set dimension of workspaces */
//	Initialize();	
  }

  // information about subordinate parameter lists
  void FSDEMatViscoT::DefineSubs(SubListT& sub_list) const
  {
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);

	/*material parameters for matrix*/
	sub_list.AddSub("rg_eq_potential", ParameterListT::Once);
	sub_list.AddSub("rg_neq_potential", ParameterListT::Any);

	/* choice of viscosity */
	sub_list.AddSub("rg_shear_viscosity", ParameterListT::Any);
	sub_list.AddSub("rg_bulk_viscosity", ParameterListT::Any);  
  
    return;
  }

  // Set electrical permittivity
  void FSDEMatViscoT::SetElectricPermittivity(double epsilon)
  {	
    fElectricPermittivity = epsilon;
  }

  // Get electrical permittivity
  double FSDEMatViscoT::GetElectricPermittivity() const
  {
    return fElectricPermittivity;
  }

  //
  void FSDEMatViscoT::SetFSDEMatSupportVisco(
      const FSDEMatSupportViscoT* support)
  {
    fFSDEMatSupportVisco = support;
  }

  //
  const dArrayT FSDEMatViscoT::ElectricField()
  {
    fElectricField = fFSDEMatSupportVisco->ElectricField();
    return fElectricField;
  }

  //
  const dArrayT FSDEMatViscoT::ElectricField(int ip)
  {
    fElectricField = fFSDEMatSupportVisco->ElectricField(ip);
    return fElectricField;
  }

  //
  const dMatrixT FSDEMatViscoT::RightCauchyGreenDeformation()
  {
    const dMatrixT& F = F_mechanical();
    dMatrixT FTF(3);
    FTF.MultATB(F, F);

    return FTF;
  }

  // material energy density
  double FSDEMatViscoT::StrainEnergyDensity()
  {

//    fEnergyDensity = EnergyDensity(C, D);
	  return 0.0;
//    return fEnergyDensity;
  }

  // material mechanical tangent modulus
  const dMatrixT& FSDEMatViscoT::C_IJKL()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
    const dArrayT& E = ElectricField();

	double I1 = C(0,0)+C(1,1)+C(2,2);
	double J = C.Det();
	J = sqrt(J);
	
	/* call C function for mechanical part of tangent modulus */
 	mech_tanmod_ab_visco(fParams.Pointer(), E.Pointer(), C.Pointer(), J, I1, fTangentMechanical.Pointer()); 
 	me_tanmod_ab_visco(fParams.Pointer(), E.Pointer(), C.Pointer(), J, fTangentMechanicalElec.Pointer());
 	fTangentMechanical+=fTangentMechanicalElec;
    return fTangentMechanical;
  }

  // Second Piola-Kirchhoff stress (mechanical)
  const dSymMatrixT& FSDEMatViscoT::S_IJ()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();
	
	dMatrixT stress_temp(3);
	dMatrixT stress_temp2(3);
	double I1 = C(0,0)+C(1,1)+C(2,2);
	double J = C.Det();
	J = sqrt(J);
	
	/* call C function for mechanical part of PK2 stress */
 	mech_pk2_ab_visco(fParams.Pointer(), E.Pointer(), C.Pointer(), J, I1, stress_temp.Pointer()); 
	me_pk2_ab_visco(fParams.Pointer(), E.Pointer(), C.Pointer(), J, stress_temp2.Pointer());
	stress_temp+=stress_temp2;
	
	fStress.FromMatrix(stress_temp);
    return fStress;
  }

  // material electromechanical tangent modulus
  const dMatrixT& FSDEMatViscoT::E_IJK()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();
	double J = C.Det();
	J = sqrt(J);

	/* call C function for electromechanical tangent modulus */
 	me_mixedmodulus_ab_visco(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), J, fTangentElectromechanical.Pointer()); 
 
    return fTangentElectromechanical;
  }

  // material electric tangent modulus
  const dMatrixT& FSDEMatViscoT::B_IJ()
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

  // Electric displacement 
  const dArrayT& FSDEMatViscoT::D_I()
  {
  	const dMatrixT& C = RightCauchyGreenDeformation();
  	const dArrayT& E = ElectricField();
  	
  	double J = C.Det();
  	J = sqrt(J);
  
	/* call C function for electric stress (i.e. electric displacement D_{I}) */
 	elec_pk2_visco(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), J, fElectricDisplacement.Pointer()); 
 		
  	return fElectricDisplacement;
  }

  // Electric field
  const dArrayT& FSDEMatViscoT::E_I()
  {
    return fElectricField;
  }

  // spatial tangent modulus
  const dMatrixT& FSDEMatViscoT::c_ijkl()
  {
    const dMatrixT F = F_mechanical();
    const double J = F.Det();

    // prevent aliasing
    const dMatrixT CIJKL = C_IJKL();
    fTangentMechanical.SetToScaled(1.0 / J, PushForward(F, CIJKL));

    return fTangentMechanical;
  }

  // Cauchy stress
  const dSymMatrixT& FSDEMatViscoT::s_ij()
  {
    const dMatrixT F = F_mechanical();
    const double J = F.Det();
	
    // prevent aliasing
    const dSymMatrixT S = S_IJ();
    fStress.SetToScaled(1.0 / J, PushForward(F, S));
    return fStress;
  }

  // pressure associated with the last computed stress
  double FSDEMatViscoT::Pressure() const
  {
    return 0.0;
  }

	/* BELOW ARE FUNCTIONS COPIED FROM RGViscoelasticityT.cpp */
	/*initializes history variable */
void  FSDEMatViscoT::PointInitialize(void)
{
	/* allocate element storage */
	ElementCardT& element = CurrentElement();	
	if (CurrIP() == 0 && fNumProcess > 0)
	{
		ElementCardT& element = CurrentElement();
		element.Dimension(0, fnstatev*NumIP());
	
		/* initialize internal variables to identity*/
		for (int ip = 0; ip < NumIP(); ip++)
		{
		      /* load state variables */
		      Load(element, ip);
		      
			  for (int i = 0; i < fNumProcess; i++)
			  {
				fC_vn[i].Identity();
				fC_v[i].Identity();
			  }

		      /* write to storage */
		      Store(element, ip);
		}
	}
}
 
void FSDEMatViscoT::UpdateHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP() && fNumProcess > 0; ip++)
	{
		/* load state variables */
		Load(element, ip);
	
		/* assign "current" to "last" */	
		for (int i = 0; i < fNumProcess; i++)
			fC_vn[i] = fC_v[i];

		/* write to storage */
		Store(element, ip);
	}
}

void FSDEMatViscoT::ResetHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP() && fNumProcess > 0; ip++)
	{
		/* load state variables*/
		Load(element, ip);
	
		/* assign "last" to "current" */
		for (int i = 0; i < fNumProcess; i++)
			fC_v[i] = fC_vn[i];
		
		/* write to storage */
		Store(element, ip);
	}
}

/* form of tangent matrix */
GlobalT::SystemTypeT FSDEMatViscoT::TangentType(void) const
{
	/* symmetric by default */
	return GlobalT::kNonSymmetric;
}

const dArrayT& FSDEMatViscoT::Compute_Eigs_v(const int process_id)
{
	fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v[process_id], false);
	return(fSpectralDecompRef.Eigenvalues());
} 

const dArrayT& FSDEMatViscoT::Compute_Eigs_vn(const int process_id)
{
	fSpectralDecompRef.SpectralDecomp_Jacobi(fC_vn[process_id], false);
	return(fSpectralDecompRef.Eigenvalues());
} 

void FSDEMatViscoT::Load(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pd = d_array.Pointer(fnstatev*ip);
	double* pdr = fstatev.Pointer();
	for (int i = 0; i < fnstatev; i++)
		*pdr++ = *pd++;
}

void FSDEMatViscoT::Store(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pdr = fstatev.Pointer();
	double* pd = d_array.Pointer(fnstatev*ip);
	for (int i = 0; i < fnstatev; i++)
		*pd++ = *pdr++;
}

/* BELOW ARE FUNCTIONS COPIED FROM RGSplitT2.cpp */
const dMatrixT& FSDEMatViscoT::MechanicalDeformation(void)
{
	const dMatrixT& F_T_inv = ThermalDeformation_Inverse();
	const dMatrixT& F = F_total();

	fF_M.MultAB(F, F_T_inv);
	return(fF_M);
}

const dMatrixT& FSDEMatViscoT::ThermalDeformation_Inverse(void)
{
	/*calculates mechanical and thermal strains in FSSolidMat*/
	const dMatrixT& F_mech = F_mechanical();

	/*retrieves thermal strains*/
	fF_T_inv = F_thermal_inverse();
	return(fF_T_inv);
}

/* Calculate the NEQ stress and stiffness in both spatial and material configurations */
const dMatrixT& FSDEMatViscoT::c_ijkl_neq()
{
	const dMatrixT& F = MechanicalDeformation();
	if (NumSD() == 2)
	{
	    fF3D[0] = F[0];
	    fF3D[1] = F[1];
	    fF3D[2] = 0.0;
	    
	    fF3D[3] = F[2];
	    fF3D[4] = F[3];
	    fF3D[5] = 0.0;
	    
	    fF3D[6] = 0.0;
	    fF3D[7] = 0.0;
	    fF3D[8] = 1.0;
	}
	else fF3D = F;

	/*calcualte total stretch*/
    fb.MultAAT(fF3D);
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
    fEigs = fSpectralDecompSpat.Eigenvalues();
    const ArrayT<dArrayT>& eigenvectors=fSpectralDecompSpat.Eigenvectors();

	/*calc EQ component of stress and moduli*/
    double J = sqrt(fEigs.Product());
    fEigs_dev = fEigs;
    fEigs_dev *= pow(J, -2.0*third);
	
    fPot[0]->DevStress(fEigs_dev, ftau_EQ);
    ftau_EQ += fPot[0]->MeanStress(J);    
    fPot[0]->DevMod(fEigs_dev,fDtauDe_EQ);
    fDtauDe_EQ += fPot[0]->MeanMod(J);

    dSymMatrixT& Gamma = fDtauDe_EQ;
    Gamma(0,0) -= 2.0*ftau_EQ[0];
    Gamma(1,1) -= 2.0*ftau_EQ[1];
    Gamma(2,2) -= 2.0*ftau_EQ[2];
   
//	cout << "\nGamma: "<<Gamma;
	
//	fModulus3D = fSpectralDecompSpat.EigsToRank4(Gamma);	
	double dl, coeff;
//	cout << "\nfModulusNEQ3D: "<<fModulusNEQ3D;

    double& l0 = fEigs[0];
    double& l1 = fEigs[1];
    double& l2 = fEigs[2];
	
	dl = l0 - l1;
    if (fabs(dl) > kSmall)
		coeff = (ftau_EQ[0]*l1 - ftau_EQ[1]*l0)/dl;
    else 
		coeff = 0.5*(Gamma(0,0)-Gamma(0,1))-ftau_EQ[0];
    MixedRank4_3D(eigenvectors[0], eigenvectors[1], fModMat);
//    fModulus3D.AddScaled(2.0*coeff, fModMat);
    
    dl = l0 - l2;
    if (fabs(dl) > kSmall)
      coeff = (ftau_EQ[0]*l2 - ftau_EQ[2]*l0)/dl;
    else 
      coeff = 0.5*(Gamma(0,0)-Gamma(0,2))-ftau_EQ[2];	
    MixedRank4_3D(eigenvectors[0], eigenvectors[2], fModMat);
//    fModulus3D.AddScaled(2.0*coeff, fModMat);
    
    dl = l1 - l2;
   if (fabs(dl) > kSmall)
		coeff  = (ftau_EQ[1]*l2 - ftau_EQ[2]*l1)/dl;
    else
      coeff = 0.5*(Gamma(1,1)-Gamma(1,2))-ftau_EQ[1];	
    MixedRank4_3D(eigenvectors[1], eigenvectors[2], fModMat);
//    fModulus3D.AddScaled(2.0*coeff, fModMat);
	fModulus3D = 0.0;	// Don't need EQ part of modulus

//   cout << "\nc_eq: "<<fModulusNEQ3D;
	/*calc NEQ component of stress and moduli*/
	/*calcualte principal values of elastic stretch*/

if (fNumProcess > 0)
{
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    
	for (int i = 0; i < fNumProcess; i++)
	{
		fInverse.Inverse(fC_vn[i]);
		fb_tr.MultQBQT(fF3D, fInverse);

		fSpectralDecompSpat.SpectralDecomp_Jacobi(fb_tr, false);	
		fEigs_tr = fSpectralDecompSpat.Eigenvalues(); 		

		fInverse.Inverse(fC_v[i]);
		fbe.MultQBQT(fF3D, fInverse);
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
		fEigs_e = fSpectralDecompSpat.Eigenvalues(); 
		const ArrayT<dArrayT>& eigenvectors_e=fSpectralDecompSpat.Eigenvectors();
		
		double Je = sqrt(fEigs_e.Product());
		fEigs_dev = fEigs_e;
		fEigs_dev *= pow(Je,-2.0*third);
    
		/*stresses*/
		fPot[i+1]->DevStress(fEigs_dev, ftau_NEQ);
		double sm =  fPot[i+1]->MeanStress(Je);    

		fPot[i+1]->DevMod(fEigs_dev, fDtauDe_NEQ);
		double cm = fPot[i+1]->MeanMod(Je);
		
		/*Calculate Calg_AB*/
		Compute_Calg(ftau_NEQ, fDtauDe_NEQ, sm, cm, fCalg, i);

		ftau_NEQ += sm;
		fDtauDe_NEQ += cm;
			   
//		cout << "\nCalg: "<<fCalg;
		fModulus3D += fSpectralDecompSpat.NonSymEigsToRank4(fCalg);
    
		double dl_tr;

		double& l0_tr = fEigs_tr[0];
		double& l1_tr = fEigs_tr[1];
		double& l2_tr = fEigs_tr[2];
	
	
		dl_tr = l0_tr - l1_tr;
		if (fabs(dl_tr) > kSmall)
			coeff = (ftau_NEQ[0]*l1_tr - ftau_NEQ[1]*l0_tr)/dl_tr;
		else 
			coeff = 0.5*(fCalg(0,0)-fCalg(0,1))-ftau_NEQ[0];
		MixedRank4_3D(eigenvectors_e[0], eigenvectors_e[1], fModMat);
		fModulus3D.AddScaled(2.0*coeff, fModMat);
    
		dl_tr = l0_tr - l2_tr;
		if (fabs(dl_tr) > kSmall)
			coeff =(ftau_NEQ[0]*l2_tr - ftau_NEQ[2]*l0_tr)/dl_tr;
		else 
			coeff = 0.5*(fCalg(0,0)-fCalg(0,2))-ftau_NEQ[2];	
		MixedRank4_3D(eigenvectors_e[0], eigenvectors_e[2], fModMat);
		fModulus3D.AddScaled(2.0*coeff, fModMat);
    
		dl_tr = l1_tr - l2_tr;
		if (fabs(dl_tr) > kSmall)
			coeff  = (ftau_NEQ[1]*l2_tr - ftau_NEQ[2]*l1_tr)/dl_tr;
		else
			coeff = 0.5*(fCalg(1,1)-fCalg(1,2))-ftau_NEQ[1];	
		MixedRank4_3D(eigenvectors_e[1], eigenvectors_e[2], fModMat);
		fModulus3D.AddScaled(2.0*coeff, fModMat);
    }
}
//   cout << "\nc_tot: "<<fModulusNEQ3D;
	if (NumSD() == 2)
	{
		fModulusNEQ[0] = fModulus3D[0];
		fModulusNEQ[1] = fModulus3D[1];
		fModulusNEQ[2] = fModulus3D[5];

		fModulusNEQ[3] = fModulus3D[6];
		fModulusNEQ[4] = fModulus3D[7];
		fModulusNEQ[5] = fModulus3D[11];
		fModulusNEQ[6] = fModulus3D[30];
		fModulusNEQ[7] = fModulus3D[31];
		fModulusNEQ[8] = fModulus3D[35];
	}
	else fModulusNEQ = fModulus3D;

	const dMatrixT& Ftotal = F_total();	
	fModulusNEQ *= 1.0/Ftotal.Det();
//	cout << "\nfModulusNEQ: "<<fModulusNEQ;

    return fModulusNEQ;
}

const dMatrixT& FSDEMatViscoT::C_IJKL_NEQ()
{
    /* deformation gradient */
    const dMatrixT& Fmat = F_total();
  
    /* transform */
    fModulusNEQ.SetToScaled(Fmat.Det(), PullBack(Fmat, c_ijkl_neq()));
    return fModulusNEQ;
}

const dSymMatrixT& FSDEMatViscoT::s_ij_neq()
{
	const dMatrixT& F = MechanicalDeformation();
	
/*	cout << "\nfF_T_inv: "<<fF_T_inv;
	cout << "\nFm: "<<F;
*/
//	cout << "\nnsd: "<<NumSD();
//	cout << "\nF: "<<F;
	if (NumSD() == 2)
	{
		fF3D[0] = F[0];
		fF3D[1] = F[1];
		fF3D[2] = 0.0;
	    
		fF3D[3] = F[2];
		fF3D[4] = F[3];
		fF3D[5] = 0.0;
	    
		fF3D[6] = 0.0;
		fF3D[7] = 0.0;
		fF3D[8] = 1.0;
	}
	else fF3D = F;
	/*calculate EQ part of the stress*/
	fb.MultAAT(fF3D);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();

	/*jacobian determinant*/
	double J = sqrt(fEigs.Product());
	
	fEigs_dev = fEigs;
	fEigs_dev *= pow(J,-2.0*third);
	
	fPot[0]->DevStress(fEigs_dev, ftau_EQ);
	ftau_EQ += fPot[0]->MeanStress(J);
	
/*		const double mu_eq = fPot[0]->GetMu();
		const double kappa_eq = fPot[0]->GetKappa();
		cout << "\neq mu: "<< mu_eq;
		cout << "\neq kappa: "<< kappa_eq;
*/	
//	fStress3D = fSpectralDecompSpat.EigsToRank2(ftau_EQ);
	fStress3D = 0.0;	// Don't need EQ part of stress

    /*load the viscoelastic principal stretches from state variable arrays*/
if (fNumProcess > 0 )
{
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    if (fFSMatSupport->RunState() == GlobalT::kFormRHS)
    {		
		/*calc NEQ component of stress and moduli*/
		for (int i = 0; i < fNumProcess; i++)
		{
			/*calc trial state*/
			fInverse.Inverse(fC_vn[i]);
			fb_tr.MultQBQT(fF3D, fInverse);
			
			fSpectralDecompSpat.SpectralDecomp_Jacobi(fb_tr, false);	
			fEigs_tr = fSpectralDecompSpat.Eigenvalues(); // OK to here

			/*calc elastic stretch*/
			fEigs_e = fEigs_tr; /*initial condition*/	// OK to here
			ComputeEigs_e(fEigs, fEigs_e, ftau_NEQ, fDtauDe_NEQ, i);
			// After ComputeEigs_e:  fEigs OK, fEigs_e=NAN, ftau_NEQ OK, fDtauDe_NEQ OK

			double Je = sqrt(fEigs_e.Product());
			fEigs_dev = fEigs_e;
			fEigs_dev *= pow(Je,-2.0*third);
	
			fPot[i+1]->DevStress(fEigs_dev, ftau_NEQ);
			ftau_NEQ += fPot[i+1]->MeanStress(Je);
			fStress3D += fSpectralDecompSpat.EigsToRank2(ftau_NEQ);
	
			/*Calculate Cv*/
			fInverse = fSpectralDecompSpat.EigsToRank2(fEigs_e); /*be which is colinear with btr*/
			fInverse.Inverse();
			fC_v[i].MultQTBQ(fF3D, fInverse); 
		}
		Store(element, CurrIP());
	}	
    else 
    {
		/*calc NEQ component of stress and moduli*/
		for (int i = 0; i < fNumProcess; i++)
		{
			/*calc elastic stretch*/
			fInverse.Inverse(fC_v[i]);
			fbe.MultQBQT(fF3D, fInverse);
			fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
			fEigs_e = fSpectralDecompSpat.Eigenvalues(); 

	//		fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	

			double Je = sqrt(fEigs_e.Product());
			fEigs_dev = fEigs_e;
			fEigs_dev *= pow(Je,-2.0*third);
		
			fPot[i+1]->DevStress(fEigs_dev, ftau_NEQ);
			ftau_NEQ += fPot[i+1]->MeanStress(Je);
			fStress3D += fSpectralDecompSpat.EigsToRank2(ftau_NEQ);
		}
    }
}

	if (NumSD() == 2)
    {
        fStressNEQ[0] = fStress3D[0];
        fStressNEQ[1] = fStress3D[1];
        fStressNEQ[2] = fStress3D[5];
    }
    else fStressNEQ = fStress3D;
	
	const dMatrixT& Ftotal = F_total();	
//	cout << "\nFtot: "<<Ftotal;
    fStressNEQ *= 1.0/Ftotal.Det();
//	cout << "\nstress: "<<fStress;
	return fStressNEQ;

}

const dSymMatrixT& FSDEMatViscoT::S_IJ_NEQ()
{
    /* deformation gradient */
    const dMatrixT& Fmat = F_total();
  
    /* transform */
    fStressNEQ.SetToScaled(Fmat.Det(), PullBack(Fmat, s_ij_neq()));
    return fStressNEQ;
}

/*************************************************************************
 * Private
 *************************************************************************/
/* construct symmetric rank-4 mixed-direction tensor (6.1.44) */
void FSDEMatViscoT::MixedRank4_2D(const dArrayT& a, const dArrayT& b, 
	dMatrixT& rank4_ab) const
{
#if __option(extended_errorcheck)
	if (a.Length() != 2 ||
	    b.Length() != 2 ||
	    rank4_ab.Rows() != 3 ||
	    rank4_ab.Cols() != 3) throw ExceptionT::kSizeMismatch;
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11;

	z1 = a[0];
	z2 = a[1];
	z3 = b[0];
	z4 = b[1];

	z5 = z1*z1;
	z6 = z2*z2;
	z7 = z3*z3;
	z8 = 2.*z1*z2*z3*z4;
	z9 = z4*z4;
	z3 = 2.*z3*z4;
	z4 = z3*z5;
	z3 = z3*z6;
	z10 = 2.*z1*z2*z7;
	z11 = 2.*z5*z7;
	z7 = z6*z7;
	z1 = 2.*z1*z2*z9;
	z2 = z5*z9;
	z5 = 2.*z6*z9;
	z4 = z10 + z4;
	z1 = z1 + z3;
	z2 = z2 + z7 + z8;
	z3 = 0.5*z4;
	z1 = 0.5*z1;
	z2 = 0.5*z2;

	//{{z11, z8, z3}, 
	// {z8, z5, z1}, 
	// {z3, z1, z2}}

	double* p = rank4_ab.Pointer();
	*p++ = z11;
    *p++ = z8;
    *p++ = z3;
    *p++ = z8;
    *p++ = z5;
    *p++ = z1;
    *p++ = z3;
    *p++ = z1;
    *p   = z2;
}

void FSDEMatViscoT::MixedRank4_3D(const dArrayT& a, const dArrayT& b, 
	dMatrixT& rank4_ab) const
{
#if __option(extended_errorcheck)
	if (a.Length() != 3 ||
	    b.Length() != 3 ||
	    rank4_ab.Rows() != 6 ||
	    rank4_ab.Cols() != 6) throw ExceptionT::kSizeMismatch;
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24;
	double z25, z26, z27, z28, z29, z30, z31, z32, z33, z34, z35, z36;
	double z37, z38, z39, z40, z41;

	z1 = a[0];	
	z2 = a[1];
	z3 = a[2];
	z4 = b[0];
	z5 = b[1];
	z6 = b[2];
	z7 = z1*z1;
	z8 = z2*z2;
	z9 = z3*z3;
	z10 = 2.*z1*z4;
	z11 = z10*z2;
	z12 = z2*z4;
	z13 = z1*z12;
	z14 = z4*z4;
	z15 = z10*z5;
	z15 = z15*z3;
	z16 = z11*z5;
	z17 = z12*z5;
	z18 = 2.*z17;
	z17 = z17*z3;
	z18 = z18*z3;
	z19 = z1*z4*z5;
	z19 = z19*z3;
	z20 = z5*z5;
	z11 = z11*z6;
	z13 = z13*z6;
	z10 = z10*z3*z6;
	z12 = z12*z3*z6;
	z21 = 2.*z12;
	z22 = z1*z2*z5*z6;
	z23 = 2.*z22;
	z24 = z1*z3*z5*z6;
	z25 = 2.*z24;
	z26 = 2.*z2*z3*z5*z6;
	z27 = z6*z6;
	z28 = 2.*z1*z14;
	z29 = 2.*z1*z2;
	z30 = z14*z2;
	z11 = z11 + z15;
	z15 = z1*z20*z3;
	z31 = 2.*z2*z20*z3;
	z18 = z18 + z23;
	z21 = z21 + z25;
	z23 = z1*z2*z27;
	z1 = 2.*z1*z27*z3;
	z25 = 2.*z2*z27*z3;
	z2 = z2*z28;
	z28 = z28*z3;
	z29 = z20*z29;
	z3 = z3*z30;
	z30 = 2.*z14*z7;
	z32 = z20*z7;
	z33 = z27*z7;
	z34 = 2.*z4*z5*z7;
	z35 = 2.*z4*z6*z7;
	z7 = z5*z6*z7;
	z36 = z14*z8;
	z37 = 2.*z20*z8;
	z38 = z27*z8;
	z39 = 2.*z4*z5*z8;
	z40 = z4*z6*z8;
	z8 = 2.*z5*z6*z8;
	z14 = z14*z9;
	z20 = z20*z9;
	z27 = 2.*z27*z9;
	z41 = z4*z5*z9;
	z4 = 2.*z4*z6*z9;
	z5 = 2.*z5*z6*z9;
	z6 = 0.5*z11;
	z9 = 0.5*z18;
	z11 = 0.5*z21;
	z2 = z2 + z34;
	z18 = z28 + z35;
	z3 = z13 + z19 + z3 + z7;
	z7 = z16 + z32 + z36;
	z13 = z29 + z39;
	z15 = z15 + z17 + z22 + z40;
	z8 = z31 + z8;
	z14 = z10 + z14 + z33;
	z17 = z20 + z26 + z38;
	z12 = z12 + z23 + z24 + z41;
	z1 = z1 + z4;
	z4 = z25 + z5;
	z2 = 0.5*z2;
	z5 = 0.5*z18;
	z3 = 0.5*z3;
	z7 = 0.5*z7;
	z13 = 0.5*z13;
	z15 = 0.5*z15;
	z8 = 0.5*z8;
	z14 = 0.5*z14;
	z17 = 0.5*z17;
	z12 = 0.5*z12;
	z1 = 0.5*z1;
	z4 = 0.5*z4;
	
	//{{z30, z16, z10,  z6,  z5,  z2}, 
	// {z16, z37, z26,  z8,  z9, z13}, 
	// {z10, z26, z27,  z4,  z1, z11}, 
	// { z6,  z8,  z4, z17, z12, z15}, 
	// { z5,  z9,  z1, z12, z14,  z3},
	// { z2, z13, z11, z15,  z3,  z7}}
	
	double* p = rank4_ab.Pointer();
    *p++ = z30;
    *p++ = z16;
    *p++ = z10;
    *p++ = z6;
    *p++ = z5;
    *p++ = z2;
    *p++ = z16;
    *p++ = z37;
    *p++ = z26;
    *p++ = z8;
    *p++ = z9;
    *p++ = z13;
    *p++ = z10;
    *p++ = z26;
    *p++ = z27;
    *p++ = z4;
    *p++ = z1;
    *p++ = z11;
    *p++ = z6;
    *p++ = z8;
    *p++ = z4;
    *p++ = z17;
    *p++ = z12;
    *p++ = z15;
    *p++ = z5;
    *p++ = z9;
    *p++ = z1;
    *p++ = z12;
    *p++ = z14;
    *p++ = z3;
    *p++ = z2;
    *p++ = z13;
    *p++ = z11;
    *p++ = z15;
    *p++ = z3;
    *p  = z7;
}

/* accept parameter list */
void FSDEMatViscoT::SetStateVariables(const int numprocess)
{
	fC_v.Dimension(numprocess);
	fC_vn.Dimension(numprocess);

	int ndof = 3;
	int numstress = dSymMatrixT::NumValues(ndof);

	fnstatev = 0;
	fnstatev += numstress;   /*current C_v*/
	fnstatev += numstress;   /*last C_vn*/

	fnstatev *= numprocess;
	
	fstatev.Dimension(fnstatev);
	double* pstatev = fstatev.Pointer();
		
	/* assign pointers to current and last blocks of state variable array */
	for (int i = 0; i < numprocess; i++)
	{
		fC_v[i].Set(ndof, pstatev);
		pstatev += numstress;
		fC_vn[i].Set(ndof, pstatev);
		pstatev += numstress;
	}
}

/* BELOW ARE FUNCTIONS FROM RGSplitT2.cpp */

 void FSDEMatViscoT::Compute_Calg(const dArrayT& tau_dev, const dSymMatrixT& dtau_dev, const double& tau_m, 
	const double& dtau_m, dMatrixT& Calg, const int type)
 {
		double s0 = tau_dev[0];
		double s1 = tau_dev[1];
		double s2 = tau_dev[2];
		
		double sm = tau_m;
		
		double c0 = dtau_dev(0,0);
		double c1 = dtau_dev(1,1);
		double c2 = dtau_dev(2,2);

		double c12 = dtau_dev(1,2);
		double c02 = dtau_dev(0,2);
		double c01 = dtau_dev(0,1);
	    
		double cm = dtau_m;
		
		/*calculates  KAB = 1+dt*D(dWdE_Idev/nD+isostress/nV)/Dep_e*/
		double dt = fFSMatSupport->TimeStep();
		
		double stress_mag = sqrt(s0*s0 + s1*s1 + s2*s2);
		fietaS = 1.0/fVisc_s[type]->Function(stress_mag);
		fietaB = 1.0/fVisc_b[type]->Function(tau_m);

		fiKAB(0,0) = 1+0.5*fietaS*dt*c0+third*fietaB*dt*cm;
		fiKAB(1,1) = 1+0.5*fietaS*dt*c1+third*fietaB*dt*cm;
		fiKAB(2,2) = 1+0.5*fietaS*dt*c2+third*fietaB*dt*cm;

		fiKAB(1,2) = 0.5*fietaS*dt*c12+third*fietaB*dt*cm;
		fiKAB(0,2) = 0.5*fietaS*dt*c02+third*fietaB*dt*cm;
		fiKAB(0,1) = 0.5*fietaS*dt*c01+third*fietaB*dt*cm;
       
		fiKAB(2,1) = fiKAB(1,2);
		fiKAB(2,0) = fiKAB(0,2);
		fiKAB(1,0) = fiKAB(0,1);
		
		if (stress_mag > kSmall)
		{
			double coeffs = fietaS*fVisc_s[type]->DFunction(stress_mag);
			double coeffb = fietaB*fVisc_b[type]->DFunction(sm);
			
			fiKAB(0,0) -= 0.5*fietaS*dt*coeffs*(s0*c0+s1*c01+s2*c02)/stress_mag*s0 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(1,1) -= 0.5*fietaS*dt*coeffs*(s0*c01+s1*c1+s2*c12)/stress_mag*s1 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(2,2) -= 0.5*fietaS*dt*coeffs*(s0*c02+s1*c12+s2*c2)/stress_mag*s2 
						- third*fietaB*dt*coeffb*(cm*sm);

			fiKAB(1,2) -= 0.5*fietaS*dt*coeffs*(s0*c02+s1*c12+s2*c2)/stress_mag*s1 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(0,2) -= 0.5*fietaS*dt*coeffs*(s0*c02+s1*c12+s2*c2)/stress_mag*s0 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(0,1) -= 0.5*fietaS*dt*coeffs*(s0*c01+s1*c1+s2*c12)/stress_mag*s0 
						- third*fietaB*dt*coeffb*(cm*sm);
						
			fiKAB(2,1) -= 0.5*fietaS*dt*coeffs*(s0*c01+s1*c1+s2*c12)/stress_mag*s2 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(2,0) -= 0.5*fietaS*dt*coeffs*(s0*c0+s1*c01+s2*c02)/stress_mag*s2 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(1,0) -= 0.5*fietaS*dt*coeffs*(s0*c0+s1*c01+s2*c02)/stress_mag*s1 
						- third*fietaB*dt*coeffb*(cm*sm);
		}
	
		/*inverts KAB*/
		fiKAB.Inverse();

//		dSymMatrixT& DAB = fDtauDe_NEQ;
//		DAB += cm; 
	
		Calg(0,0) = (c0+cm)*fiKAB(0,0) + (c01+cm)*fiKAB(1,0) + (c02+cm)*fiKAB(2,0) - 2.0*(tau_dev[0]+tau_m);
		Calg(1,0) = (c01+cm)*fiKAB(0,0) + (c1+cm)*fiKAB(1,0) + (c12+cm)*fiKAB(2,0);
		Calg(2,0) = (c02+cm)*fiKAB(0,0) + (c12+cm)*fiKAB(1,0) + (c2+cm)*fiKAB(2,0);
		Calg(0,1) = (c0+cm)*fiKAB(0,1) + (c01+cm)*fiKAB(1,1) + (c02+cm)*fiKAB(2,1);
		Calg(1,1) = (c01+cm)*fiKAB(0,1) + (c1+cm)*fiKAB(1,1) + (c12+cm)*fiKAB(2,1) - 2.0*(tau_dev[1]+tau_m);
		Calg(2,1) = (c02+cm)*fiKAB(0,1) + (c12+cm)*fiKAB(1,1) + (c2+cm)*fiKAB(2,1);
		Calg(0,2) = (c0+cm)*fiKAB(0,2) + (c01+cm)*fiKAB(1,2) + (c02+cm)*fiKAB(2,2);
		Calg(1,2) = (c01+cm)*fiKAB(0,2) + (c1+cm)*fiKAB(1,2) + (c12+cm)*fiKAB(2,2);
		Calg(2,2) = (c02+cm)*fiKAB(0,2) + (c12+cm)*fiKAB(1,2) + (c2+cm)*fiKAB(2,2) - 2.0*(tau_dev[2]+tau_m);
}

void FSDEMatViscoT::ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
			     dArrayT& eigenstress, dSymMatrixT& eigenmodulus, const int type) 
{		
	const double ctol = 1.00e-14;
		
	/*set references to principle stretches*/
     
	double& le0 = eigenstretch_e[0];
	double& le1 = eigenstretch_e[1];
	double& le2 = eigenstretch_e[2];
  
	double tol;

	/*initialize principle elastic and trial elastic log strains */
	double ep_tr0 = 0.5*log(le0);
	double ep_tr1 = 0.5*log(le1);
	double ep_tr2 = 0.5*log(le2);

	double ep_e0 = ep_tr0;		
	double ep_e1 = ep_tr1;	
	double ep_e2 = ep_tr2;
	
	/*initializes principle viscous stretch*/
	int iteration  = 0;	
	do 
	{
		iteration ++;
	    double Je=sqrt(le0*le1*le2);
	    fEigs_dev = eigenstretch_e;
	    fEigs_dev *= pow(Je,-2.0*third);
		
	    /*calculate stresses and moduli*/
	    fPot[type+1]->DevStress(fEigs_dev, eigenstress);
	    
	    double& s0 = eigenstress[0];
	    double& s1 = eigenstress[1];
	    double& s2 = eigenstress[2];
		
		double stress_mag = sqrt(s0*s0 + s1*s1 + s2*s2);
		fietaS = 1.0/fVisc_s[type]->Function(stress_mag);
		
	    fPot[type+1]->DevMod(fEigs_dev,eigenmodulus);
		/*deviatoric values*/
		double& c0 = eigenmodulus(0,0);
		double& c1 = eigenmodulus(1,1);
		double& c2 = eigenmodulus(2,2);

		double& c12 = eigenmodulus(1,2);
		double& c02 = eigenmodulus(0,2);
		double& c01 = eigenmodulus(0,1);
	    
	    /*caculate means*/
	    double sm = fPot[type+1]->MeanStress(Je);
		fietaB = 1.0/fVisc_b[type]->Function(sm);

	    double cm = fPot[type+1]->MeanMod(Je);
	    
		double dt = fFSMatSupport->TimeStep();
		fiKAB(0,0) = 1+0.5*fietaS*dt*c0+third*fietaB*dt*cm;
		fiKAB(1,1) = 1+0.5*fietaS*dt*c1+third*fietaB*dt*cm;
		fiKAB(2,2) = 1+0.5*fietaS*dt*c2+third*fietaB*dt*cm;

		fiKAB(1,2) = 0.5*fietaS*dt*c12+third*fietaB*dt*cm;
		fiKAB(0,2) = 0.5*fietaS*dt*c02+third*fietaB*dt*cm;
		fiKAB(0,1) = 0.5*fietaS*dt*c01+third*fietaB*dt*cm;
       
		fiKAB(2,1) = fiKAB(1,2);
		fiKAB(2,0) = fiKAB(0,2);
		fiKAB(1,0) = fiKAB(0,1);
	
		if (stress_mag > kSmall)
		{
			double coeffs = fietaS*fVisc_s[type]->DFunction(stress_mag);
			double coeffb = fietaB*fVisc_b[type]->DFunction(sm);
			
			fiKAB(0,0) -= 0.5*fietaS*dt*coeffs*(s0*c0+s1*c01+s2*c02)/stress_mag*s0 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(1,1) -= 0.5*fietaS*dt*coeffs*(s0*c01+s1*c1+s2*c12)/stress_mag*s1 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(2,2) -= 0.5*fietaS*dt*coeffs*(s0*c02+s1*c12+s2*c2)/stress_mag*s2 
						- third*fietaB*dt*coeffb*(cm*sm);

			fiKAB(1,2) -= 0.5*fietaS*dt*coeffs*(s0*c02+s1*c12+s2*c2)/stress_mag*s1 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(0,2) -= 0.5*fietaS*dt*coeffs*(s0*c02+s1*c12+s2*c2)/stress_mag*s0 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(0,1) -= 0.5*fietaS*dt*coeffs*(s0*c01+s1*c1+s2*c12)/stress_mag*s0 
						- third*fietaB*dt*coeffb*(cm*sm);
						
			fiKAB(2,1) -= 0.5*fietaS*dt*coeffs*(s0*c01+s1*c1+s2*c12)/stress_mag*s2 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(2,0) -= 0.5*fietaS*dt*coeffs*(s0*c0+s1*c01+s2*c02)/stress_mag*s2 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(1,0) -= 0.5*fietaS*dt*coeffs*(s0*c0+s1*c01+s2*c02)/stress_mag*s1 
						- third*fietaB*dt*coeffb*(cm*sm);
		}
	
		/*inverts KAB*/
		fiKAB.Inverse();
	    
	    /*calculate the residual*/
	    double res0 = ep_e0 + dt*(0.5*fietaS*s0 +
			  third*fietaB*sm) - ep_tr0;
	    double res1 = ep_e1 + dt*(0.5*fietaS*s1 +
			  third*fietaB*sm) - ep_tr1;
	    double res2 = ep_e2 + dt*(0.5*fietaS*s2 +
			  third*fietaB*sm) - ep_tr2;
		
	    /*solve for the principal strain increments*/
	    double dep_e0=-fiKAB(0,0)*res0-fiKAB(0,1)*res1-fiKAB(0,2)*res2;
	    double dep_e1=-fiKAB(1,0)*res0-fiKAB(1,1)*res1-fiKAB(1,2)*res2;
	    double dep_e2=-fiKAB(2,0)*res0-fiKAB(2,1)*res1-fiKAB(2,2)*res2;
	    
	    /*updates principal elastic stretches*/ 
	    ep_e0 += dep_e0;
	    ep_e1 += dep_e1;
	    ep_e2 += dep_e2;
	    
	    le0 = exp(2.0*ep_e0);
	    le1 = exp(2.0*ep_e1);
	    le2 = exp(2.0*ep_e2);
	    
	    /*Check that the L2 norm of the residual is less than tolerance*/
	    tol = sqrt(res0*res0 + res1*res1+res2*res2);
	}while (tol>ctol && iteration < 10); 
	if (iteration >= 10) 
		ExceptionT::GeneralFail("RGSplitT2::ComputeEigs_e", 
			"number of iteration exceeds maximum of 10");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSDEMatViscoT::NewSub(const StringT& name) const
{
	PotentialT* pot = NULL;
	if (name == "neo-hookean")
		pot = new NeoHookean;
	else if (name == "mooney-rivlin")
		pot = new MooneyRivlin;
	else if (name == "veronda-westmann")
		pot = new  VWPotentialT;
	else if (name == "fung-potential")
		pot = new  FungPotentialT;
	else if (name == "arruda-boyce")
		pot = new ArrudaBoyce;
	if (pot)
		return pot;

	C1FunctionT* func = NULL;
	if (name == "linear_exponential")
		func = new LinearExponentialT;
#ifdef __DEVELOPMENT__
	else if (name == "scaled-csch")
		func = new ScaledCsch;
#endif

	if (func)
		return func;

	if (name == "rg_eq_potential" || name == "rg_neq_potential")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
	
		/* choice of parameters */
		choice->AddSub("neo-hookean");
		choice->AddSub("mooney-rivlin");
		choice->AddSub("veronda-westmann");
		choice->AddSub("fung-potential");
		choice->AddSub("arruda-boyce");
		return(choice);
	}
	else if (name == "rg_shear_viscosity")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetDescription("eta_S(|sig_dev|)");	
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
		
#ifdef __DEVELOPMENT__
		choice->AddSub("scaled-csch");
#endif
		choice->AddSub("linear_exponential");
		return(choice);
	}
	else if (name == "rg_bulk_viscosity")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetDescription("eta_B(sig_m)");	
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
		
#ifdef __DEVELOPMENT__
		choice->AddSub("scaled-csch");
#endif
		choice->AddSub("linear_exponential");
		return(choice);
	}
	else return(FSDEMatViscoT::NewSub(name));
}

} //namespace Tahoe
