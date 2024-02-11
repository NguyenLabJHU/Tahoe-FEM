
namespace Tahoe {

  inline FSDEMat2DViscoT::FSDEMat2DViscoT() :
    ParameterInterfaceT("Finite Strain Dielectric Elastomer 2D Visco"),
        fFSDEMatSupport2DVisco(0)
  {
    SetName(FSDEMat2DViscoT::Name);
    Initialize();
  }

  //
  // Set electrical permittivity
  //
  inline void FSDEMat2DViscoT::SetElectricPermittivity(double epsilon)
  {	
    fElectricPermittivity = epsilon;
  }

  //
  // Get electrical permittivity
  //
  inline double FSDEMat2DViscoT::GetElectricPermittivity() const
  {
    return fElectricPermittivity;
  }

  //
  //
  //
  inline void FSDEMat2DViscoT::SetFSDEMatSupport2DVisco(
      const FSDEMatSupport2DViscoT* support)
  {
    fFSDEMatSupport2DVisco = support;
  }


  //
  //
  //
  inline const dArrayT FSDEMat2DViscoT::ElectricField()
  {
    fElectricField = fFSDEMatSupport2DVisco->ElectricField();
    return fElectricField;
  }

  //
  //
  //
  inline const dArrayT FSDEMat2DViscoT::ElectricField(int ip)
  {
    fElectricField = fFSDEMatSupport2DVisco->ElectricField(ip);
    return fElectricField;
  }

  //
  //
  //
  inline const dMatrixT FSDEMat2DViscoT::RightCauchyGreenDeformation()
  {
    const dMatrixT F = F_mechanical();
    dMatrixT FTF(2);
    FTF.MultATB(F, F);

    return FTF;

  }

  //
  // material energy density
  //
  inline double FSDEMat2DViscoT::StrainEnergyDensity()
  {

//    fEnergyDensity = EnergyDensity(C, D);
	  return 0.0;
//    return fEnergyDensity;

  }

  //
  // material mechanical tangent modulus
  //
  inline const dMatrixT&
  FSDEMat2DViscoT::C_IJKL()
  {
    const dMatrixT F = F_mechanical();
    dMatrixT C(2);
    C.MultATB(F, F);  
    const dArrayT& E = ElectricField();
  
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

	fTangentMechanical(0,0) =D1111-0.5*coef*(G111111*EE1*EE1+G111122*EE2*EE2+2.0*G111112*EE1*EE2);
  	fTangentMechanical(0,1) =D1122-0.5*coef*(G112211*EE1*EE1+G112222*EE2*EE2+2.0*G112212*EE1*EE2);
  	fTangentMechanical(0,2) =D1112-0.5*coef*(G111211*EE1*EE1+G111222*EE2*EE2+2.0*G111212*EE1*EE2);
  	fTangentMechanical(1,0) = fTangentMechanical(0,1);
  	fTangentMechanical(1,1) = D2222-0.5*coef*(G222211*EE1*EE1+G222222*EE2*EE2+2.0*G222212*EE1*EE2);
  	fTangentMechanical(1,2) =D2212-0.5*coef*(G221211*EE1*EE1+G221222*EE2*EE2+2.0*G221212*EE1*EE2);
  	fTangentMechanical(2,0) = fTangentMechanical(0,2);
  	fTangentMechanical(2,1) = fTangentMechanical(1,2);
  	fTangentMechanical(2,2) = D1212-0.5*coef*(G121211*EE1*EE1+G121222*EE2*EE2+2.0*G121212*EE1*EE2); 	

    return fTangentMechanical;
  }

  //
  // Second Piola-Kirchhoff stress (mechanical)
  //
  inline const dSymMatrixT&
  FSDEMat2DViscoT::S_IJ()
  {
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

    return fStress;
  }

  //
  // material electromechanical tangent modulus
  //
  inline const dMatrixT&
  FSDEMat2DViscoT::E_IJK()
  {
    const dMatrixT F = F_mechanical();
    dMatrixT C(2);
    C.MultATB(F, F);  
	const dArrayT& E = ElectricField();
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
	
    fTangentElectromechanical(0,0) =-2.0*coef*(CI11*CI11*EE1+CI11*CI12*EE2);
 	fTangentElectromechanical(0,1) =-2.0*coef*(CI12*CI11*EE1+CI12*CI12*EE2);
  	fTangentElectromechanical(1,0) =-2.0*coef*(CI12*CI12*EE1+CI12*CI22*EE2);
  	fTangentElectromechanical(1,1) =-2.0*coef*(CI22*CI12*EE1+CI22*CI22*EE2);
  	fTangentElectromechanical(2,0) =-1.0*coef*(2.0*CI11*CI12*EE1+(CI11*CI22+CI12*CI12)*EE2);
  	fTangentElectromechanical(2,1) =-1.0*coef*((CI11*CI22+CI12*CI12)*EE1+2.0*CI12*CI22*EE2);
 
    return fTangentElectromechanical;
  }

  //
  // material electric tangent modulus
  //
  inline const dMatrixT&
  FSDEMat2DViscoT::B_IJ()
  {
    const dMatrixT F = F_mechanical();
    dMatrixT C(2);
    C.MultATB(F, F);  

	dMatrixT Cinv(2);
	Cinv.Inverse(C);
	fTangentElectrical = Cinv;
	fTangentElectrical *= fElectricPermittivity;

    return fTangentElectrical;
  }

  //
  // Electric displacement - is it necessary to pass Efield?
  //
  inline const dArrayT&
  FSDEMat2DViscoT::D_I()
  {
    const dMatrixT F = F_mechanical();
    dMatrixT C(2);
    C.MultATB(F, F);  
  	const dArrayT& E = ElectricField();
	dMatrixT CInv = C.Inverse();
	CInv *= fElectricPermittivity;
	CInv.Multx(E, fElectricDisplacement);

    return fElectricDisplacement;
  }

  //
  // Electric field
  //
  inline const dArrayT&
  FSDEMat2DViscoT::E_I()
  {
    return fElectricField;
  }

  //
  // spatial tangent modulus
  //
  inline const dMatrixT&
  FSDEMat2DViscoT::c_ijkl()
  {
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
  FSDEMat2DViscoT::s_ij()
  {
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
  inline double FSDEMat2DViscoT::Pressure() const
  {

    return 0.0;

  }


} //namespace Tahoe
