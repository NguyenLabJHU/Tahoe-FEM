/* $Id: FDCrystalElast.cpp,v 1.10 2004-07-15 08:28:16 paklein Exp $ */
#include "FDCrystalElast.h"

#include <stdlib.h>
#include "CrystalElastLat.h"
#include "CrystalElastMat.h"
#include "ElementCardT.h"

#include "Utils.h"
#include "UpLagr_ExternalFieldT.h"
#include "FSMatSupportT.h"
#include "SpectralDecompT.h"

/* spatial dimensions of the problem */
const int kNSD = 3;
/* element output data */
const int kNumOutput = 4;
static const char* Labels[kNumOutput] = {"VM_stress","s1","s2","sn"};
/* Numerical constants */
const double sqrt23 = sqrt(2.0/3.0);

using namespace Tahoe;

FDCrystalElast::FDCrystalElast(ifstreamT& in, const FSMatSupportT& support):
	ParameterInterfaceT("large_strain_crystal_elasticity"),
	CrystalElast(in, support),
  // deformation gradient 
	fF    (kNSD,kNSD),
  // elastic and thermal (inverse) deformation gradients
        fFe   (kNSD,kNSD),
        fFthi (kNSD,kNSD),
  // symmetric tensors in intermediate configuration
  	fCe   (kNSD),
	fBe   (kNSD),
  	fEe   (kNSD),
  	fS    (kNSD),
  // tensors in polar decomposition of fe
        fUe   (kNSD),
        fRe   (kNSD,kNSD),
  // crystal orientation matrix
  	fRotMat (kNSD,kNSD),
  // crystal Cauchy stress
	fs_ij (kNSD),
  // crystal elasticity matrix
  	fc_ijkl  (dSymMatrixT::NumValues(kNSD)),
  // anisotropic -cubic- part of elastic matrix
  	fCanisoLat  (dSymMatrixT::NumValues(kNSD)),  // lattice axes
  	fCanisoSpl  (dSymMatrixT::NumValues(kNSD)),  // sample axes (B0)
  // 2nd order identity tensor
  	fISym (kNSD),
  // principal values of Cauchy stress
  	fsEigs (kNSD),
  // tensor of thermal expansion
        falpha (kNSD),
  // work spaces
  	fmatx1    (kNSD,kNSD),
  	fmatx2    (kNSD,kNSD),
  	fmatx3    (kNSD,kNSD),
  	fRank4    (dSymMatrixT::NumValues(kNSD)),
  	fsymmatx1 (kNSD),
	fvector1  (kNSD)
{
  // check number of grains
  if (fNumGrain != 1)
    throwRunTimeError("FDCrystalElast::FDCrystalElast: NumGrain != 1");

  // set anisotropic -cubic- part of elasticity matrix (lattice axes)
  fCanisoLat = 0.;
  fCanisoLat(0,0) = fCanisoLat(1,1) = fCanisoLat(2,2) = 1.;

  // set 2nd order unit tensor (sym matrix)
  fISym.Identity();

  fExFieldElement = TB_DYNAMIC_CAST(const UpLagr_ExternalFieldT*, FSMatSupport().FiniteStrain());
  if (!fExFieldElement)
    {
      cout << "\n FDCrystalElast::FDCrystalElast: could not cast element group to \n"
           << "UpLagr_ExternalFieldT" << endl;
      throw ExceptionT::kGeneralFail;
    }
}

FDCrystalElast::~FDCrystalElast() {}

int FDCrystalElast::NumVariablesPerElement()
{
  int d_size = 0;
  int dim = dSymMatrixT::NumValues(kNSD);

  // variables per crystal per ip
  d_size += kNSD*kNSD;                  // fRotMat (const)
  d_size += kNSD*kNSD;                  // fFe    
  d_size += dim ;                       // fs_ij  
  d_size += dim * dim;                  // fc_ijkl

  // total # variables per element
  d_size *= NumIP();

  return d_size;
}

const dSymMatrixT& FDCrystalElast::s_ij()
{
  // fetch current element and int point
  ElementCardT& element = CurrentElement();
  int intpt = CurrIP();
  int elem  = CurrElementNumber();

  // recover local data
  LoadCrystalData(element, intpt);

  array1.Dimension(1);

  // compute elastic moduli and stress due to thermal strain
  if (MaterialSupport().RunState() == GlobalT::kFormRHS)
    {
      // deformation gradient
      fF = F();

      // temperature at current IP
      fExFieldElement->IP_Interpolate(fExFieldElement->ExternalField(), array1);
      fTemp_DegC = array1[0];

      // temperature dependent properties
      fCrystalElastMat->ElasticityProps(fMatProp, fTemp_DegC, elem, intpt);
      fCrystalElastMat->ThermalProps(falpha, fTemp_DegC);

      // fCanisoLat -> fCanisoSpl (for anisotropic elasticity)
      if (!fCrystalElastMat->IsIsotropic()) 
                  FFFFC_3D(fCanisoSpl, fCanisoLat, fRotMat);

      // compute Cauchy stress
      CrystalS_ij();

      // compute consistent tangent
      CrystalC_ijkl();
    }

  // return Cauchy stress
  return fs_ij;
}

const dMatrixT& FDCrystalElast::c_ijkl()
{
  // fetch current element and int point
  ElementCardT& element = CurrentElement();
  int intpt = CurrIP();

  // recover local data
  LoadCrystalData(element, intpt);

  // return crystal moduli
  return fc_ijkl;
}

int FDCrystalElast::NumOutputVariables() const {return kNumOutput;}

void FDCrystalElast::OutputLabels(ArrayT<StringT>& labels) const
{
  // allocate space for labels
  labels.Dimension(kNumOutput);

  // copy labels
  for (int i = 0; i < kNumOutput; i++)
    labels[i] = Labels[i];
}

void FDCrystalElast::ComputeOutput(dArrayT& output)
{
  // gather element/integ point information
  ElementCardT& element = CurrentElement();
  int elem  = CurrElementNumber();
  int intpt = CurrIP();

  // fetch crystal data
  LoadCrystalData(element, intpt);

  // polar decomposition of Fe
  SpectralDecompT fSpecD(kNSD);
  fSpecD.PolarDecomp(fFe, fRe, fUe, false);

  // Von Mises stress 
  output[0] = sqrt(fsymmatx1.Deviatoric(fs_ij).ScalarProduct())/sqrt23;

  // principal values of stress
  fs_ij.PrincipalValues(fsEigs);
  output[1] = fsEigs[0];
  output[2] = fsEigs[1];

  // normal to slip plane crystal orientation
  fVecNormC.Dimension(3);
  fVecNormC[0] = 1.; fVecNormC[1] = 1.; fVecNormC[2] = 0.; 
  fVecNormC.UnitVector();

  // normal to slip plane specimen orientation
  fVecNorm.Dimension(3);
  fRotMat.Multx(fVecNormC, fvector1);
  fRe.Multx(fvector1, fVecNorm);

  // compute stress normal to slip plane
  fs_ij.Multx(fVecNorm, fvector1);
  a_i_b_i(fsnorm,fvector1,fVecNorm);
  output[3] = fsnorm;

  // write texture at IP/ELE 
  const int& step = ContinuumElement().ElementSupport().StepNumber();
  fmatx1.MultAB(fRe, fRotMat);
  dArrayT& angles = fangles[0];
  fCrystalElastLat->RotMatrixToAngles(fmatx1, angles);
  fCrystalElastLat->WriteTexture(elem, intpt, fNumGrain, step, fangles);
}

/* PROTECTED MEMBER FUNCTIONS */

void FDCrystalElast::InitializeCrystalVariables()
{
  // initialize state at each element and ...
  for (int elem = 0; elem < NumElements(); elem++)
    {
      // get pointer to element elem
      ElementCardT& element = ElementCard(elem);

      // ... at each integration point
      for (int intpt = 0; intpt < NumIP(); intpt++)
	{
          // fetch crystal data
          LoadCrystalData(element, intpt);

	  int igrn = 0;

          // fetch euler angles
          dArrayT& angles = fEuler[elem](intpt, igrn);
 
	  // storage rotation matrix from Euler angles
	  fCrystalElastLat->AnglesToRotMatrix(angles, fRotMat);

	  // elastic deformation gradient
	  fFe = 0.;

	  // crystal Cauchy stress
	  fs_ij = 0.;

          // anisotropic crystal elasticity matrix
          fc_ijkl = 0.;
	}
    }
}

void FDCrystalElast::LoadCrystalData(ElementCardT& element, int intpt)
{
  // fetch internal variable array
  dArrayT& d_array = element.DoubleData();

  // decode
  int dim   = dSymMatrixT::NumValues(kNSD);
  int block = 2*kNSD*kNSD + dim + dim * dim;
  int dex   = intpt*block;

  fRotMat.Set (kNSD,kNSD, &d_array[dex             ]);
  fFe.Set     (kNSD,kNSD, &d_array[dex += kNSD*kNSD]);
  fs_ij.Set   (kNSD,       &d_array[dex += kNSD*kNSD]);
  fc_ijkl.Set (dim, dim,  &d_array[dex += dim      ]);
}

void FDCrystalElast::CrystalS_ij()
{
  // inverse of thermal deformation gradient
  InverseFthermal();

  // elastic deformation gradient
  fFe.MultAB(fF, fFthi); 

  // elastic right Cauchy-Green tensor (Ce)
  fCe.MultATA(fFe);

  // elastic Green strain
  fEe.SetToCombination(0.5, fCe, -0.5, fISym);

  // intermediate 2nd P-K stress
  // ... isotropic contribution: 2*mu*Ee + lambda*tr(Ee)*I
  fS.SetToCombination(2.*fMatProp[0], fEe, fMatProp[1]*fEe.Trace(), fISym);

  // ... anisotropic contribution: -2*beta*Caniso[Ee]
  if (!fCrystalElastMat->IsIsotropic()) 
    {
      fsymmatx1.A_ijkl_B_kl(fCanisoSpl, fEe);
      fS.AddScaled(-2.*fMatProp[2], fsymmatx1 );
    }

  // Cauchy Stress
  fs_ij.MultQBQT(fFe, fS);
  fs_ij /= fFe.Det();
}

void FDCrystalElast::InverseFthermal()
{
  // Fthi = (I+0.5*dtemp*A)^(-1) * (I-0.5*dtheta*A)
  
  double dtheta = fTemp_DegC - fInit_Temp_DegC;

  fmatx1.SetToScaled(dtheta/2.0, falpha);

  fmatx2.Identity();
  fmatx2.AddScaled(1.0, fmatx1);
  fmatx2.Inverse();

  fmatx3.Identity();
  fmatx3.AddScaled(-1.0, fmatx1);
  fFthi.MultAB(fmatx2, fmatx3);
}

void FDCrystalElast::CrystalC_ijkl()
{
  // elastic left Cauchy-Green tensor (Be)
  fBe.MultAAT(fFe);

  // I_b tensor
  Set_I_b_Tensor(fBe, fc_ijkl);

  // Be (x) Be
  fRank4.Outer(fBe, fBe);

  // isotropic elastic contribution to elasticity tensor
  fc_ijkl *= 2.*fMatProp[0];
  fc_ijkl.AddScaled(fMatProp[1], fRank4);

  // anisotropic (cubic) contribution to elasticity tensor
  if (!fCrystalElastMat->IsIsotropic())
    {
      FFFFC_3D(fRank4, fCanisoSpl, fFe);
      fc_ijkl.AddScaled(-2.*fMatProp[2], fRank4);
    }
}

void FDCrystalElast::Set_I_b_Tensor(const dSymMatrixT& b, dMatrixT& c)
{
  double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
  double z13, z14, z15, z16, z17, z18, z19, z20, z21;
  
  z1 = b[0];
  z2 = b[1];
  z3 = b[2];
  z4 = b[3];
  z5 = b[4];
  z6 = b[5];
  z7 = z1*z1;
  z8 = z1*z2;
  z9 = z2*z2;
  z10 = z1*z3;
  z11 = z2*z3;
  z12 = z3*z3;
  z13 = z1*z4;
  z14 = z2*z4;
  z15 = z3*z4;
  z16 = z4*z4;
  z17 = z1*z5;
  z18 = z2*z5;
  z19 = z3*z5;
  z20 = z4*z5;
  z21 = z5*z5;
  z1 = z1*z6;
  z2 = z2*z6;
  z3 = z3*z6;
  z4 = z4*z6;
  z5 = z5*z6;
  z6 = z6*z6;
  z11 = z11 + z16;
  z10 = z10 + z21;
  z3 = z20 + z3;
  z18 = z18 + z4;
  z13 = z13 + z5;
  z8 = z6 + z8;
  z11 = 0.5*z11;
  z10 = 0.5*z10;
  z3 = 0.5*z3;
  z18 = 0.5*z18;
  z13 = 0.5*z13;
  z8 = 0.5*z8;
  
  //{{z7, z6, z21, z5, z17, z1}, 
  // {z6, z9, z16, z14, z4, z2}, 
  // {z21, z16, z12, z15, z19, z20}, 
  // {z5, z14, z15, z11, z3, z18}, 
  // {z17, z4, z19, z3, z10, z13}, 
  // {z1, z2, z20, z18, z13, z8}}
  
  c(0,0) = z7;
  c(0,1) = c(1,0) = z6;
  c(0,2) = c(2,0) = z21;
  c(0,3) = c(3,0) = z5;
  c(0,4) = c(4,0) = z17;
  c(0,5) = c(5,0) = z1;
  c(1,1) = z9;
  c(1,2) = c(2,1) = z16;
  c(1,3) = c(3,1) = z14;
  c(1,4) = c(4,1) = z4;
  c(1,5) = c(5,1) = z2;
  c(2,2) = z12;
  c(2,3) = c(3,2) = z15;
  c(2,4) = c(4,2) = z19;
  c(2,5) = c(5,2) = z20;
  c(3,3) = z11;
  c(3,4) = c(4,3) = z3;
  c(3,5) = c(5,3) = z18;
  c(4,4) = z10;
  c(4,5) = c(5,4) = z13;
  c(5,5) = z8;
}

void FDCrystalElast::FFFFC_3D(dMatrixT& Co, dMatrixT& Ci, const dMatrixT& F)
{
  int nsd = F.Rows();
  dSymMatrixT coltemp;
  dMatrixT outer(nsd);
  dMatrixT transform(dSymMatrixT::NumValues(nsd));

  // compute tranformation matrix
  dArrayT Fi1(nsd,F(0));
  dArrayT Fi2(nsd,F(1));
  dArrayT Fi3(nsd,F(2));

  coltemp.Set(nsd,transform(0));
  outer.Outer(Fi1,Fi1);
  coltemp.FromMatrix(outer);

  coltemp.Set(nsd,transform(1));
  outer.Outer(Fi2,Fi2);
  coltemp.FromMatrix(outer);

  coltemp.Set(nsd,transform(2));
  outer.Outer(Fi3,Fi3);
  coltemp.FromMatrix(outer);

  coltemp.Set(nsd,transform(3));
  outer.Outer(Fi2,Fi3);
  outer.Symmetrize();
  coltemp.FromMatrix(outer);
  coltemp *= 2.0;

  coltemp.Set(nsd,transform(4));
  outer.Outer(Fi1,Fi3);
  outer.Symmetrize();
  coltemp.FromMatrix(outer);
  coltemp *= 2.0;

  coltemp.Set(nsd,transform(5));
  outer.Outer(Fi1,Fi2);
  outer.Symmetrize();
  coltemp.FromMatrix(outer);
  coltemp *= 2.0;
	
  // compute transformed tensor
  Co.MultQBQT(transform, Ci);
}

void FDCrystalElast::a_i_b_i(double& C, dArrayT a, dArrayT b)
{
  C = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
