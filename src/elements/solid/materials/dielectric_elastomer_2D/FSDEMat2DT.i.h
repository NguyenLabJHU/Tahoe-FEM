#include "FSDE_inc_2D.h"

namespace Tahoe {

  inline FSDEMat2DT::FSDEMat2DT() :
    ParameterInterfaceT("Finite Strain Dielectric Elastomer 2D"),
        fFSDEMatSupport2D(0)
  {
    SetName(FSDEMat2DT::Name);
    Initialize();
  }

  //
  // Set electrical permittivity
  //
  inline void FSDEMat2DT::SetElectricPermittivity(double epsilon)
  {	
    fElectricPermittivity = epsilon;
  }

  //
  // Get electrical permittivity
  //
  inline double FSDEMat2DT::GetElectricPermittivity() const
  {
    return fElectricPermittivity;
  }

  //
  //
  //
  inline void FSDEMat2DT::SetFSDEMatSupport2D(
      const FSDEMatSupport2DT* support)
  {
    fFSDEMatSupport2D = support;
  }


  //
  //
  //
  inline const dArrayT FSDEMat2DT::ElectricField()
  {
    fElectricField = fFSDEMatSupport2D->ElectricField();
    return fElectricField;
  }

  //
  //
  //
  inline const dArrayT FSDEMat2DT::ElectricField(int ip)
  {
    fElectricField = fFSDEMatSupport2D->ElectricField(ip);
    return fElectricField;
  }

  //
  //
  //
  inline const dMatrixT FSDEMat2DT::RightCauchyGreenDeformation()
  {
    const dMatrixT F = F_mechanical();
    dMatrixT FTF(2);
    FTF.MultATB(F, F);

    return FTF;

  }

  //
  // material energy density
  //
  inline double FSDEMat2DT::StrainEnergyDensity()
  {

//    fEnergyDensity = EnergyDensity(C, D);
	  return 0.0;
//    return fEnergyDensity;

  }


  /* Mechanical and electromechanical tangent modulus */
  inline void FSDEMat2DT::C_Mech_Elec(dMatrixT& mech, dMatrixT& elec)
  {
//  	cout << "FSDEMat2T::C_Mech_Elec" << endl;
  	const dMatrixT& C = RightCauchyGreenDeformation();
  	const dArrayT& E = ElectricField();
  
  	/* call C function for both tangent moduli */
  	get_ddCmech_elec_2D(fParams.Pointer(), E.Pointer(), C.Pointer(),
  		mech.Pointer(), elec.Pointer());
  
  	mech *= 4.0;
  	elec *= -2.0;
//  	cout << "mech = " << mech << endl;
//  	cout << "elec = " << elec << endl;
  }

  /* Electrical displacement and tangent modulus */
  inline void FSDEMat2DT::S_C_Elec(dArrayT& D, dMatrixT& CE)
  {
//  	cout << "FSDEMat2T::S_C_Elec" << endl;
  	const dMatrixT& C = RightCauchyGreenDeformation();
  	const dArrayT& E = ElectricField();
  
  	/* call C function for both tangent moduli */
  	get_ddC_sc_elec_2D(fParams.Pointer(), E.Pointer(), C.Pointer(),
  		D.Pointer(), CE.Pointer());
  
  	D *= -1.0;
  	CE *= -1.0;
//  	cout << "D = " << D << endl;
//  	cout << "CE = " << CE << endl;
  }


  inline void 
  FSDEMat2DT::ComputeAllLHS(dMatrixT& Cmech, dMatrixT& Cemech, dMatrixT& elec)
  {
//  	cout << "FSDEMat2DT::C_IJKL" << endl;
//    const dMatrixT& C = RightCauchyGreenDeformation();
    const dArrayT& E = ElectricField();
    const dMatrixT F = F_mechanical();
    dMatrixT C(2);
    C.MultATB(F, F);  
  
  	double det_C = C.Det();
  	double Invar_1 = C(0,0) + C(1,1) + 1.0/det_C;
 	double I1 = Invar_1 / fNrig;	
  	double fp1 = 0.5+I1/10.0 + 11.0*I1*I1*I1/350.0 + 19.0*I1*I1*I1/1750.0 + 519.0*I1*I1*I1*I1/134750.0;	
  	double fp2 = 1/(10*fNrig) + 11*I1/(175*fNrig) + 57*I1*I1/(1750*fNrig) + 1038*I1*I1*I1/(67375*fNrig);
  	double coef = fElectricPermittivity;
  	dMatrixT CInv = C.Inverse();
  	double CI11 = CInv(0,0);
  	double CI12 = CInv(0,1);
  	double CI22 = CInv(1,1);
  	double EE1 = E[0];
  	double EE2 = E[1];
  	double lamda_new = 4.0*fp1*fMu/det_C;
  	double mu_new = 2.0*fp1*fMu/det_C;
  	double mufp2 = 4.0*fMu*fp2;

    double D1111=mufp2*(1-CI11/det_C)*(1-CI11/det_C);
    double D1122=mufp2*(1-CI11/det_C)*(1-CI22/det_C);
    double D1112=mufp2*(1-CI11/det_C)*(0-CI12/det_C);
    double D2222=mufp2*(1-CI22/det_C)*(1-CI22/det_C);
    double D2212=mufp2*(1-CI22/det_C)*(0-CI12/det_C);
    double D1212=mufp2*(0-CI12/det_C)*(0-CI12/det_C);

	D1111=D1111+lamda_new*CI11*CI11+mu_new*2.0*CI11*CI11;
  	D1122=D1122+lamda_new*CI11*CI22+mu_new*2.0*CI12*CI12;
  	D1112=D1112+lamda_new*CI11*CI12+mu_new*2.0*CI11*CI12;
  	D2222=D2222+lamda_new*CI22*CI22+mu_new*2.0*CI22*CI22;
  	D2212=D2212+lamda_new*CI22*CI12+mu_new*2.0*CI12*CI22;
  	D1212=D1212+lamda_new*CI12*CI12+mu_new*(CI11*CI22+CI12*CI12);

	double G111111=8.0*CI11*CI11*CI11; 
	double G111122=8.0*CI11*CI12*CI12; 
	double G111112=8.0*CI11*CI11*CI12; 
  	double G112211=G111122; 
  	double G112222=8.0*CI12*CI12*CI22; 
  	double G112212=4.0*(CI11*CI22+CI12*CI12)*CI12; 
  	double G111211=G111112; 
  	double G111222=G112212; 
  	double G111212=2.0*CI11*(CI11*CI22+CI12*CI12)+4.0*CI11*CI12*CI12;
  	double G222211=G112222; 
  	double G222222=8.0*CI22*CI22*CI22; 
  	double G222212=8.0*CI12*CI22*CI22;
  	double G221211=G112212; 
  	double G221222=G222212; 
  	double G221212=2.0*CI22*(CI11*CI22+CI12*CI12)+4.0*CI12*CI12*CI22;
  	double G121211=G111212; 
  	double G121222=G221212; 
  	double G121212=2.0*CI12*(CI11*CI22+CI12*CI12)+4.0*CI11*CI12*CI22;

	Cmech(0,0) =D1111-0.5*coef*(G111111*EE1*EE1+G111122*EE2*EE2+2.0*G111112*EE1*EE2);
  	Cmech(0,1) =D1122-0.5*coef*(G112211*EE1*EE1+G112222*EE2*EE2+2.0*G112212*EE1*EE2);
  	Cmech(0,2) =D1112-0.5*coef*(G111211*EE1*EE1+G111222*EE2*EE2+2.0*G111212*EE1*EE2);
  	Cmech(1,0) = Cmech(0,1);
  	Cmech(1,1) = D2222-0.5*coef*(G222211*EE1*EE1+G222222*EE2*EE2+2.0*G222212*EE1*EE2);
  	Cmech(1,2) =D2212-0.5*coef*(G221211*EE1*EE1+G221222*EE2*EE2+2.0*G221212*EE1*EE2);
  	Cmech(2,0) = Cmech(0,2);
  	Cmech(2,1) = Cmech(1,2);
  	Cmech(2,2) =D1212-0.5*coef*(G121211*EE1*EE1+G121222*EE2*EE2+2.0*G121212*EE1*EE2);
//	cout << "Cmech = " << Cmech << endl;
	
    Cemech(0,0) =-2.0*coef*(CI11*CI11*EE1+CI11*CI12*EE2);
 	Cemech(0,1) =-2.0*coef*(CI12*CI11*EE1+CI12*CI12*EE2);
  	Cemech(1,0) =-2.0*coef*(CI12*CI12*EE1+CI12*CI22*EE2);
  	Cemech(1,1) =-2.0*coef*(CI22*CI12*EE1+CI22*CI22*EE2);
  	Cemech(2,0) =-1.0*coef*(2.0*CI11*CI12*EE1+(CI11*CI22+CI12*CI12)*EE2);
  	Cemech(2,1) =-1.0*coef*((CI11*CI22+CI12*CI12)*EE1+2.0*CI12*CI22*EE2);
// 	cout << "Cemech = " << Cemech << endl;
  
  	elec = CInv;
  	elec *= coef;
//  	cout << "Elec = " << elec << endl;
  }


  //
  // material mechanical tangent modulus
  //
  inline const dMatrixT&
  FSDEMat2DT::C_IJKL()
  {
	cout << "FSDEMat2DT::C_IJKL" << endl;
    const dMatrixT& C = RightCauchyGreenDeformation();
    const dArrayT& E = ElectricField();
 	
	/* call C function for mechanical tangent modulus */
 	get_ddCmech_2D(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), fTangentMechanical.Pointer()); 
	fTangentMechanical*=4.0;
	cout << "C_IJKL = " << fTangentMechanical << endl;
    return fTangentMechanical;
  }

  //
  // material electromechanical tangent modulus
  //
  inline const dMatrixT&
  FSDEMat2DT::E_IJK()
  {
	cout << "FSDEMat2DT::E_IJK" << endl;
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();

	/* call C function for electromechanical tangent modulus */
 	get_ddCE_2D(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), fTangentElectromechanical.Pointer()); 
 
 	fTangentElectromechanical*=-2.0;
 	cout << "E_IJK = " << fTangentElectromechanical << endl;
    return fTangentElectromechanical;

  }

  //
  // material electric tangent modulus
  //
  inline const dMatrixT&
  FSDEMat2DT::B_IJ()
  {
	cout << "FSDEMat2DT::B_IJ" << endl;
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();

	/* call C function for electrical tangent modulus */
 	get_ddE_2D(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), fTangentElectrical.Pointer()); 
 
 	fTangentElectrical *= -1.0;
 	cout << "B_IJ = " << fTangentElectrical << endl;
    return fTangentElectrical;

  }

  //
  // Second Piola-Kirchhoff stress (mechanical)
  //
  inline const dSymMatrixT&
  FSDEMat2DT::S_IJ()
  {
//	cout << "FSDEMat2DT::S_IJ" << endl;
//     const dMatrixT& C = RightCauchyGreenDeformation();
//     const dArrayT& E = ElectricField();
//      
//  	/* call C function for mechanical stress */
//  	get_dUdCmech_2D(fParams.Pointer(), E.Pointer(),  
//  		C.Pointer(), stress_temp.Pointer()); 
//  
//     fStress.FromMatrix(stress_temp);
//  	fStress*=2.0;

    const dMatrixT F = F_mechanical();
    dMatrixT C(2);
    C.MultATB(F, F); 
	const dArrayT& E = ElectricField();

	double det_C = C.Det();
	double Invar_1 = C(0,0) + C(1,1) + 1.0/det_C;
 	double I1 = Invar_1 / fNrig;	// fill in
  	double fp1 = 0.5+I1/10.0 + 11.0*I1*I1*I1/350.0 + 19.0*I1*I1*I1/1750.0 + 519.0*I1*I1*I1*I1/134750.0;	
	double mu = fMu*2.0*fp1;
	dMatrixT CInv = C.Inverse();
	double CI11 = CInv(0,0);
	double CI22 = CInv(1,1);
	double CI12 = CInv(0,1);
	double coef = fElectricPermittivity;
	double EE1 = E[0];
	double EE2 = E[1];
	
	fStress[0]=(-mu/det_C)*CI11+mu+coef*(CI11*CI11*EE1*EE1+CI12*CI12*EE2*EE2+2.0*CI11*CI12*EE1*EE2);
    fStress[1]=(-mu/det_C)*CI22+mu+coef*(CI12*CI12*EE1*EE1+CI22*CI22*EE2*EE2+2.0*CI12*CI22*EE1*EE2);
    fStress[2]=(-mu/det_C)*CI12+coef*(CI11*CI12*EE1*EE1+CI12*CI22*EE2*EE2+(CI11*CI22+CI12*CI12)*EE1*EE2); 
	
	/* ZHou also defines Efield here */
	
//	cout << "S_IJ = " << fStress << endl;
    return fStress;

  }

  //
  // Electric displacement - is it necessary to pass Efield?
  //
  inline const dArrayT&
  FSDEMat2DT::D_I()
  {
//  	cout << "FSDEMat2DT::D_I" << endl;
//     const dMatrixT& C = RightCauchyGreenDeformation();
 	const dArrayT& E = ElectricField();
 
//     dMatrixT CE(2);
//   
//   	/* call C function for both tangent moduli */
//   	get_ddC_sc_elec_2D(fParams.Pointer(), E.Pointer(), C.Pointer(),
//   		fElectricDisplacement.Pointer(), CE.Pointer());
//   
//   	fElectricDisplacement *= -1.0;
//   	CE *= -1.0;
	
// 	/* call C function for electric displacement */
//  	get_dUdE_2D(fParams.Pointer(), E.Pointer(),  
//  		C.Pointer(), fElectricDisplacement.Pointer());     
//     
//     fElectricDisplacement *= -1.0;
//     cout << "fElectricDisplacement = " << fElectricDisplacement << endl;\\

	const dMatrixT F = F_mechanical();
    dMatrixT C(2);
    C.MultATB(F, F); 
  	dMatrixT CInv = C.Inverse();
	CInv *= fElectricPermittivity;
	CInv.Multx(E, fElectricDisplacement);
//	cout << "fElectricDisplacement = " << fElectricDisplacement << endl;

    return fElectricDisplacement;
  }

  //
  // Electric field
  //
  inline const dArrayT&
  FSDEMat2DT::E_I()
  {
    return fElectricField;
  }

  //
  // spatial tangent modulus
  //
  inline const dMatrixT&
  FSDEMat2DT::c_ijkl()
  {
//	cout << "FSDEMat2DT::c_ijkl" << endl;
    const dMatrixT F = F_mechanical();
    const double J = F.Det();

    // prevent aliasing
    const dMatrixT CIJKL = C_IJKL();
    fTangentMechanical.SetToScaled(1.0 / J, PushForward(F, CIJKL));

    return fTangentMechanical;

  }

  //
  // Cauchy stress
  //
  inline const dSymMatrixT&
  FSDEMat2DT::s_ij()
  {
//	cout << "FSDEMat2DT::s_ij" << endl;
    const dMatrixT F = F_mechanical();
    const double J = F.Det();
	
    // prevent aliasing
    const dSymMatrixT S = S_IJ();
    fStress.SetToScaled(1.0 / J, PushForward(F, S));
    return fStress;

  }

  //
  // pressure associated with the last computed stress
  //
  inline double FSDEMat2DT::Pressure() const
  {

    return 0.0;

  }

  //
  // compute symmetric Cij reduced index matrix */
  //
  inline void FSDEMat2DT::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli) 
  {


  }

  //
  // compute symmetric 2nd PK2 reduced index vector */
  //
  inline void FSDEMat2DT::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2) 
  {


  }

  //
  // compute strain energy density for the specified strain */
  //
  inline double FSDEMat2DT::ComputeEnergyDensity(const dSymMatrixT& E) 
  {
	 return 0.0;

  }

} //namespace Tahoe
