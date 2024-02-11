#include "FSDE_incQ1P02D.h"

#include <iostream>
#include <fstream>
using namespace std;

namespace Tahoe {

  inline FSDEMatQ1P02DT::FSDEMatQ1P02DT() :
    ParameterInterfaceT("Finite Strain Dielectric Elastomer Q1P02D"),
        fFSDEMatSupportQ1P02D(0)
  {
    SetName(FSDEMatQ1P02DT::Name);
    Initialize();
  }

  // Set electrical permittivity
  inline void FSDEMatQ1P02DT::SetElectricPermittivity(double epsilon)
  {	
    fElectricPermittivity = epsilon;
  }

  // Get electrical permittivity
  inline double FSDEMatQ1P02DT::GetElectricPermittivity() const
  {
    return fElectricPermittivity;
  }

  //
  inline void FSDEMatQ1P02DT::SetFSDEMatSupportQ1P02D(
      const FSDEMatSupportQ1P02DT* support)
  {
    fFSDEMatSupportQ1P02D = support;
  }

  //
  inline const dArrayT FSDEMatQ1P02DT::ElectricField()
  {
    fElectricField = fFSDEMatSupportQ1P02D->ElectricField();
	/*---------- Prescribing Electric Field ---------- */
	//fElectricField[0] = 0.0;
    //fElectricField[1] = -0.1;
    //fElectricField[2] = 0.0;
     /*-----------------------------------------------*/
    return fElectricField;
  }

  //
  inline const dArrayT FSDEMatQ1P02DT::ElectricField(int ip)
  {
    //fElectricField = fFSDEMatSupportQ1P02D->ElectricField(ip);
	/*---------- Prescribing Electric Field ----------*/
 	//fElectricField[0] = 0.0;
    //fElectricField[1] = -0.1;
    //fElectricField[2] = 0.0;
    /*-----------------------------------------------*/
    return fElectricField;
  }

  /* ---------- Prescribed Deformation Tensor ---------- */
  inline const dMatrixT FSDEMatQ1P02DT::DeformationMatrix()
  {
	  dMatrixT F = F_mechanical();
	  //dMatrixT F(3);
	  //F(0,0) = 1.05; F(1,0) = 0.0;  F(2,0) = 0.0;
	  //F(0,1) = 0.0;  F(1,1) = 0.95; F(2,1) = 0.0;
	  //F(0,2) = 0.0;  F(1,2) = 0.0;  F(2,2) = 1.0;
	  return F;
  }
  /* -----------------------------------------------------*/
  //
  inline const dMatrixT FSDEMatQ1P02DT::RightCauchyGreenDeformation()
  {
    const dMatrixT& F = F_mechanical();
    dMatrixT FTF(2);
    FTF.MultATB(F, F);

    return FTF;
  }

  // material energy density
  inline double FSDEMatQ1P02DT::StrainEnergyDensity()
  {
	  return 0.0;
  }

  // material mechanical tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P02DT::C_IJKL()
  {
    const dMatrixT& F2D = F_mechanical();
    const dMatrixT& C2D = RightCauchyGreenDeformation();
    double J = F2D.Det();


    dMatrixT C3D(3), F3D(3), mechtan3D(6), metan3D(6);
    dArrayT E3D(3); 

    metan3D = 0.0;
    mechtan3D = 0.0;

    C3D[0] = C2D[0];
    C3D[1] = C2D[1];
    C3D[2] = 0.0;
    
    C3D[3] = C2D[2];
    C3D[4] = C2D[3];
    C3D[5] = 0.0;
    
    C3D[6] = 0.0;
    C3D[7] = 0.0;
    C3D[8] = 1.0;
    
    F3D[0] = F2D[0];
    F3D[1] = F2D[1];
    F3D[2] = 0.0;
    
    F3D[3] = F2D[2];
    F3D[4] = F2D[3];
    F3D[5] = 0.0;
    
    F3D[6] = 0.0;
    F3D[7] = 0.0;
    F3D[8] = 1.0;
    
    const dArrayT& E = ElectricField();
    E3D[0] = E[0];
    E3D[1] = E[1];
    E3D[2] = 0.0;
    
  	double I1 = C2D(0,0) + C2D(1,1) + 1.0;	// plane strain constraint
	
	/* call C function for mechanical part of tangent modulus */
 	mech_tanmod_q1p02D(fParams.Pointer(), E3D.Pointer(), C3D.Pointer(), F3D.Pointer(), J, I1, mechtan3D.Pointer());
 	me_tanmod_q1p02D(fParams.Pointer(), E3D.Pointer(), C3D.Pointer(), F3D.Pointer(), J, metan3D.Pointer());
 	mechtan3D+=metan3D;
 
 	fTangentMechanical(0,0) = mechtan3D(0,0);
	fTangentMechanical(0,1) = mechtan3D(0,1);
	fTangentMechanical(0,2) = mechtan3D(0,5);
	fTangentMechanical(1,0) = mechtan3D(1,0);
	fTangentMechanical(1,1) = mechtan3D(1,1);
	fTangentMechanical(1,2) = mechtan3D(1,5);
	fTangentMechanical(2,0) = mechtan3D(5,0);
	fTangentMechanical(2,1) = mechtan3D(5,1);
	fTangentMechanical(2,2) = mechtan3D(5,5);

 	/* -------------- Writing into a file -------------
 	ofstream myCIJKL;
 	myCIJKL.open("C_IJKL.txt");
 	for (int i = 0; i < fTangentMechanical.Rows(); i++) {
 		for (int j = 0; j < fTangentMechanical.Cols(); j++) {
 			myCIJKL << fTangentMechanical(i,j) << " ";
 		}
 		myCIJKL << endl;
 	}
 	myCIJKL.close();
 	 ------------------------------------------------- */

    return fTangentMechanical;
  }

  // Second Piola-Kirchhoff stress (mechanical)
  inline const dSymMatrixT&
  FSDEMatQ1P02DT::S_IJ()
  {
    const dMatrixT& F2D = F_mechanical();
    const dMatrixT& C2D = RightCauchyGreenDeformation();
    double J = F2D.Det();


    dMatrixT C3D(3), F3D(3), stress_temp(3), stress_temp2(3);
    dArrayT E3D(3); 

    stress_temp = 0.0;
    stress_temp2 = 0.0;

    C3D[0] = C2D[0];
    C3D[1] = C2D[1];
    C3D[2] = 0.0;
    
    C3D[3] = C2D[2];
    C3D[4] = C2D[3];
    C3D[5] = 0.0;
    
    C3D[6] = 0.0;
    C3D[7] = 0.0;
    C3D[8] = 1.0;
    
    F3D[0] = F2D[0];
    F3D[1] = F2D[1];
    F3D[2] = 0.0;
    
    F3D[3] = F2D[2];
    F3D[4] = F2D[3];
    F3D[5] = 0.0;
    
    F3D[6] = 0.0;
    F3D[7] = 0.0;
    F3D[8] = 1.0;
    
    const dArrayT& E = ElectricField();
    E3D[0] = E[0];
    E3D[1] = E[1];
    E3D[2] = 0.0;
    
  	double I1 = C2D(0,0) + C2D(1,1) + 1.0;	// plane strain constraint
	
	/* call C function for mechanical part of PK2 stress */
 	mech_pk2_q1p02D(fParams.Pointer(), E3D.Pointer(), C3D.Pointer(), F3D.Pointer(), J, I1, stress_temp.Pointer()); 
	me_pk2_q1p02D(fParams.Pointer(), E3D.Pointer(), C3D.Pointer(), F3D.Pointer(), J, stress_temp2.Pointer());
	stress_temp+=stress_temp2;
	
	fStress(0,0) = stress_temp(0,0);
    fStress(0,1) = stress_temp(0,1);
    fStress(1,0) = stress_temp(1,0);
    fStress(1,1) = stress_temp(1,1);
    
 	/* -------------- Writing into a file -------------
 	ofstream mySIJ;
 	mySIJ.open("S_IJ.txt");
 	for (int i = 0; i < fStress.Rows(); i++) {
 		for (int j = 0; j < fStress.Cols(); j++) {
 			mySIJ << fStress(i,j) << " ";
 		}
 		mySIJ << endl;
 	}
 	mySIJ.close();
 	 ------------------------------------------------- */

    return fStress;
  }

  // material electromechanical tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P02DT::E_IJK()
  {
    const dMatrixT& F2D = F_mechanical();
    const dMatrixT& C2D = RightCauchyGreenDeformation();
    double J = F2D.Det();


    dMatrixT C3D(3), F3D(3), tangentEM3D(6,3);
    dArrayT E3D(3); 
    
    tangentEM3D = 0.0;

    C3D[0] = C2D[0];
    C3D[1] = C2D[1];
    C3D[2] = 0.0;
    
    C3D[3] = C2D[2];
    C3D[4] = C2D[3];
    C3D[5] = 0.0;
    
    C3D[6] = 0.0;
    C3D[7] = 0.0;
    C3D[8] = 1.0;
    
    F3D[0] = F2D[0];
    F3D[1] = F2D[1];
    F3D[2] = 0.0;
    
    F3D[3] = F2D[2];
    F3D[4] = F2D[3];
    F3D[5] = 0.0;
    
    F3D[6] = 0.0;
    F3D[7] = 0.0;
    F3D[8] = 1.0;
    
    const dArrayT& E = ElectricField();
    E3D[0] = E[0];
    E3D[1] = E[1];
    E3D[2] = 0.0;

	/* call C function for electromechanical tangent modulus */
 	me_mixedmodulus_q1p02D(fParams.Pointer(), E3D.Pointer(), C3D.Pointer(), F3D.Pointer(), J, tangentEM3D.Pointer());

	fTangentElectromechanical(0,0) = tangentEM3D(0,0);
	fTangentElectromechanical(1,0) = tangentEM3D(1,0);
	fTangentElectromechanical(2,0) = tangentEM3D(5,0);
	fTangentElectromechanical(0,1) = tangentEM3D(0,1);
	fTangentElectromechanical(1,1) = tangentEM3D(1,1);
	fTangentElectromechanical(2,1) = tangentEM3D(5,1);

 	/* -------------- Writing into a file -------------
 	ofstream myEIJK;
 	myEIJK.open("E_IJK.txt");
 	for (int i = 0; i < fTangentElectromechanical.Rows(); i++) {
 		for (int j = 0; j < fTangentElectromechanical.Cols(); j++) {
 			myEIJK << fTangentElectromechanical(i,j) << " ";
 		}
 		myEIJK << endl;
 	}
 	myEIJK.close();
 	------------------------------------------------- */

    return fTangentElectromechanical;
  }

  // spatial electromechanical tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P02DT::e_ijk()
  {
    const dMatrixT& F2D = F_mechanical();
    const dMatrixT& C2D = RightCauchyGreenDeformation();
    double J = F2D.Det();


    dMatrixT C3D(3), F3D(3), tangentEM3Ds(6,3);
    dArrayT E3D(3);

    tangentEM3Ds = 0.0;

    C3D[0] = C2D[0];
    C3D[1] = C2D[1];
    C3D[2] = 0.0;
    
    C3D[3] = C2D[2];
    C3D[4] = C2D[3];
    C3D[5] = 0.0;
    
    C3D[6] = 0.0;
    C3D[7] = 0.0;
    C3D[8] = 1.0;
    
    F3D[0] = F2D[0];
    F3D[1] = F2D[1];
    F3D[2] = 0.0;
    
    F3D[3] = F2D[2];
    F3D[4] = F2D[3];
    F3D[5] = 0.0;
    
    F3D[6] = 0.0;
    F3D[7] = 0.0;
    F3D[8] = 1.0;
    
    const dArrayT& E = ElectricField();
    E3D[0] = E[0];
    E3D[1] = E[1];
    E3D[2] = 0.0;
	
	/* call C function for (spatial) electromechanical tangent modulus */
 	me_mixedmodulus_q1p02Dspatial(fParams.Pointer(), E3D.Pointer(), C3D.Pointer(), F3D.Pointer(), J, tangentEM3Ds.Pointer());
 
 	fTangentElectromechanicalSpatial(0,0) = tangentEM3Ds(0,0);
	fTangentElectromechanicalSpatial(1,0) = tangentEM3Ds(1,0);
	fTangentElectromechanicalSpatial(2,0) = tangentEM3Ds(5,0);
	fTangentElectromechanicalSpatial(0,1) = tangentEM3Ds(0,1);
	fTangentElectromechanicalSpatial(1,1) = tangentEM3Ds(1,1);
	fTangentElectromechanicalSpatial(2,1) = tangentEM3Ds(5,1);
 
 	fTangentElectromechanicalSpatial /= J;

 	/* -------------- Writing into a file -------------
 	ofstream myeijk;
 	myeijk.open("e_ijk_2.txt");
 	for (int i = 0; i < fTangentElectromechanicalSpatial.Rows(); i++) {
 		for (int j = 0; j < fTangentElectromechanicalSpatial.Cols(); j++) {
 			myeijk << fTangentElectromechanicalSpatial(i,j) << " ";
 		}
 		myeijk << endl;
 	}
 	myeijk.close();
 	------------------------------------------------- */

    return fTangentElectromechanicalSpatial;
  }

  // material electric tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P02DT::B_IJ()
  {
    const dMatrixT& F2D = F_mechanical();
    const dMatrixT& C2D = RightCauchyGreenDeformation();
    double J = F2D.Det();


    dMatrixT C3D(3), ElecTan3D(3);
    ElecTan3D = 0.0;	// initialize
    
    C3D[0] = C2D[0];
    C3D[1] = C2D[1];
    C3D[2] = 0.0;
    
    C3D[3] = C2D[2];
    C3D[4] = C2D[3];
    C3D[5] = 0.0;
    
    C3D[6] = 0.0;
    C3D[7] = 0.0;
    C3D[8] = 1.0;

	dMatrixT Cinv(3);
	Cinv.Inverse(C3D);
	ElecTan3D = Cinv;
	ElecTan3D *= fElectricPermittivity;
	ElecTan3D *= J;
	
	fTangentElectrical(0,0) = ElecTan3D(0,0);
	fTangentElectrical(1,0) = ElecTan3D(1,0);
	fTangentElectrical(0,1) = ElecTan3D(0,1);
	fTangentElectrical(1,1) = ElecTan3D(1,1);
	
 	/* -------------- Writing into a file -------------
 	ofstream myBIJ;
 	myBIJ.open("B_IJ.txt");
 	for (int i = 0; i < fTangentElectrical.Rows(); i++) {
 		for (int j = 0; j < fTangentElectrical.Cols(); j++) {
 			myBIJ << fTangentElectrical(i,j) << " ";
 		}
 		myBIJ << endl;
 	}
 	myBIJ.close();
 	------------------------------------------------- */

    return fTangentElectrical;
  }

  // spatial electric tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P02DT::b_ij()
  {
    const dMatrixT& F = F_mechanical();
    double J = F.Det();

    // prevent aliasing
    const dMatrixT b = B_IJ();
    fTangentElectrical.MultABCT(F, b, F);
    fTangentElectrical /= J;

 	/* -------------- Writing into a file -------------
 	ofstream mybij;
 	mybij.open("b_ij_2.txt");
 	for (int i = 0; i < fTangentElectrical.Rows(); i++) {
 		for (int j = 0; j < fTangentElectrical.Cols(); j++) {
 			mybij << fTangentElectrical(i,j) << " ";
 		}
 		mybij << endl;
 	}
 	mybij.close();
 	------------------------------------------------- */

    return fTangentElectrical;
  }

  // Electric displacement 
  inline const dArrayT&
  FSDEMatQ1P02DT::D_I()
  {
    const dMatrixT& F2D = F_mechanical();
    const dMatrixT& C2D = RightCauchyGreenDeformation();
    double J = F2D.Det();


    dMatrixT C3D(3), F3D(3);
    dArrayT E3D(3), ED(3); 

    ED = 0.0;

    C3D[0] = C2D[0];
    C3D[1] = C2D[1];
    C3D[2] = 0.0;
    
    C3D[3] = C2D[2];
    C3D[4] = C2D[3];
    C3D[5] = 0.0;
    
    C3D[6] = 0.0;
    C3D[7] = 0.0;
    C3D[8] = 1.0;
    
    F3D[0] = F2D[0];
    F3D[1] = F2D[1];
    F3D[2] = 0.0;
    
    F3D[3] = F2D[2];
    F3D[4] = F2D[3];
    F3D[5] = 0.0;
    
    F3D[6] = 0.0;
    F3D[7] = 0.0;
    F3D[8] = 1.0;
    
    const dArrayT& E = ElectricField();
    E3D[0] = E[0];
    E3D[1] = E[1];
    E3D[2] = 0.0;
  
	/* call C function for electric stress (i.e. electric displacement D_{I}) */
 	elec_pk2_q1p02D(fParams.Pointer(), E3D.Pointer(),  
 		C3D.Pointer(), F3D.Pointer(), J, ED.Pointer()); 
 		
 	fElectricDisplacement[0] = ED[0];
 	fElectricDisplacement[1] = ED[1];

 	/* -------------- Writing into a file -------------
 	ofstream myDI;
 	myDI.open("D_I.txt");
 	for (int i = 0; i < fElectricDisplacement.Length(); i++) {
 		myDI << fElectricDisplacement[i] << " ";
 	}
 	myDI.close();
 	------------------------------------------------- */

  	return fElectricDisplacement;
  }

  // spatial electric tangent modulus
  inline const dArrayT&
  FSDEMatQ1P02DT::d_i()
  {
    const dMatrixT& F = F_mechanical();
    const double J = F.Det();

    // prevent aliasing
    const dArrayT D = D_I();
	F.Multx(D, fElectricDisplacement);
	fElectricDisplacement /= J;

 	/* -------------- Writing into a file -------------
 	ofstream mydi;
 	mydi.open("D_I_2.txt");
 	for (int i = 0; i < fElectricDisplacement.Length(); i++) {
 		mydi << fElectricDisplacement[i] << " ";
 	}
 	mydi.close();
 	------------------------------------------------- */

    return fElectricDisplacement;		// need to divide by J
  }

  // Electric field
  inline const dArrayT&
  FSDEMatQ1P02DT::E_I()
  {
    return fElectricField;
  }

  // spatial tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P02DT::c_ijkl()
  {
//     const dMatrixT& F = F_mechanical();
//     const double J = F.Det();
// 
//     // prevent aliasing
//     const dMatrixT CIJKL = C_IJKL();
//     fTangentMechanical.SetToScaled(1.0 / J, PushForward(F, CIJKL));
    
	fTangentMechanical = FSSolidMatT::c_ijkl();

 	/* -------------- Writing into a file -------------
 	ofstream mycijkl;
 	mycijkl.open("c_ijkl_2.txt");
 	for (int i = 0; i < fTangentMechanical.Rows(); i++) {
 		for (int j = 0; j < fTangentMechanical.Cols(); j++) {
 			mycijkl << fTangentMechanical(i,j) << " ";
 		}
 		mycijkl << endl;
 	}
 	mycijkl.close();
 	------------------------------------------------- */

    return fTangentMechanical;

  }

  inline const dMatrixT&
  FSDEMatQ1P02DT::a_ijkl()
  {
	  const dSymMatrixT& sigma = s_ij();
	  const dMatrixT& cijkl = C_IJKL();

	  cout << sigma << endl;

	  fa(0, 0) = cijkl(0, 0) + 4*sigma(0, 0);
	  fa(0, 1) = cijkl(0, 1);
	  fa(0, 2) = cijkl(0, 2);
	  fa(1, 0) = cijkl(1, 0);
	  fa(1, 1) = cijkl(1, 1) + 4*sigma(1, 1);
	  fa(1, 2) = cijkl(1, 2) + 4*sigma(0, 1);
	  fa(2, 0) = cijkl(2, 0);
	  fa(2, 1) = cijkl(2, 1) + 4*sigma(1, 0);
	  fa(2, 2) = cijkl(2, 2) + 4*sigma(0, 0);
	  cout << "hello" <<endl;
	  return fa;
  }

  // Cauchy stress
  inline const dSymMatrixT&
  FSDEMatQ1P02DT::s_ij()
  {
    const dMatrixT& F = F_mechanical();
    const double J = F.Det();

    // prevent aliasing
    const dSymMatrixT S = S_IJ();
    fStress.SetToScaled(1.0 / J, PushForward(F, S));

 	/* -------------- Writing into a file -------------
 	ofstream mysij;
 	mysij.open("s_ij_2.txt");
 	for (int i = 0; i < fStress.Rows(); i++) {
 		for (int j = 0; j < fStress.Cols(); j++) {
 			mysij << fStress(i,j) << " ";
 		}
 		mysij << endl;
 	}
 	mysij.close();
 	------------------------------------------------- */

    return fStress;
  }

  // pressure associated with the last computed stress
  inline double FSDEMatQ1P02DT::Pressure() const
  {

    return 0.0;

  }

}	// namespace Tahoe
