/*
  File: LocalCrystalPlast_C.cpp
*/

#include "LocalCrystalPlast_C.h"
#include "LatticeOrient.h"
#include "VoceHardening.h"
#include "ElementCardT.h"
#include "ifstreamT.h"
#include "Utils.h"

#include "ElasticT.h"
#include "FEManagerT.h"

/* spatial dimensions of the problem */
const int kNSD = 3; 

const double sqrt23 = sqrt(2.0/3.0);

/* element output data */
const int kNumOutput = 2;
static const char* Labels[kNumOutput] = {"VM_stress", "Hardness"};

LocalCrystalPlast_C::LocalCrystalPlast_C(ifstreamT& in, const ElasticT& element) :
  LocalCrystalPlast(in, element),
  fLocInitX (element.InitialCoordinates()),
  //fNNodes   (element.NumElemNodes()),
  fNNodes   (fLocInitX.NumberOfNodes()),
  fLNa      (1, fNNodes),
  fLDNa     (NumSD(), fNNodes),
  fGDNa     (NumSD(), fNNodes),
  fGradU    (NumSD())
{

}

LocalCrystalPlast_C::~LocalCrystalPlast_C() {} 

void LocalCrystalPlast_C::Initialize()
{
  // inherited
  LocalCrystalPlast::Initialize();
  
  // shape functions and derivatives at center
  SetLocalShape_C(fLNa, fLDNa);
}

int LocalCrystalPlast_C::NumVariablesPerElement()
{
  int d_size = 0;
  int dim = dSymMatrixT::NumValues(kNSD);

  // variables per crystal 
  d_size += kNSD*kNSD;                       // fRotMat   (const)
  d_size += kNSD*kNSD;                       // fFpi_n    (t_n)
  d_size += fNumSlip;                        // fDGamma_n (t_n)
  d_size += kNSD*kNSD;                       // fFpi      (t)
  d_size += fNumSlip;                        // fDGamma   (t)
  d_size += kNSD*kNSD;                       // fFe       (t)
  d_size += dim;                             // fs_ij     (t)
  d_size += fHardening->NumberOfVariables(); // hard vars at t and t_n

  // total # crystal variables per element (at the center)
  d_size *= 1 * fNumGrain;

  // averaged (aggregate) stress and moduli per element
  d_size += 1 * dim;                         // fsavg_ij   (t)
  d_size += 1 * dim * dim;                   // fcavg_ijkl (t)

  return d_size;
}

const dSymMatrixT& LocalCrystalPlast_C::s_ij()
{
  // fetch current element
  ElementCardT& element = CurrentElement();
  int intpt = 0;

  // recover aggregate data
  LoadAggregateData(element, intpt);

  // compute state, stress and moduli at center of element
  if (fStatus == GlobalT::kFormRHS && CurrIP() == 0)
    {
      // reset iteration counter to check NLCSolver
      fIterCount = 0;

      // initialize average stress and moduli at center
      fsavg_ij = 0.0;
      fcavg_ijkl = 0.0;

      // deformation gradient at center of element
      fFtot = DeformationGradient(fLocDisp);
      // fFtot = fContinuumElement.FEManager().DeformationGradient();

      // crystal state/stress
      for (int igrn = 0; igrn < fNumGrain; igrn++)
	{
	  // recover local data
	  LoadCrystalData(element, intpt, igrn);
	  
	  // Schmidt tensors in sample coordinates
	  for (int i = 0; i < fNumSlip; i++)
	    {
	      fZ[i].MultQBQT(fRotMat, fZc[i]);
	      fP[i].Symmetrize(fZ[i]);
	    }
	      
	  // compute crystal state
	  SolveCrystalState();
	  
	  // compute crystal Cauchy stress
	  CrystalS_ij();
	  
	  // compute plastic deformation gradient
	  FPInverse();

	  // compute crystal moduli
	  fc_ijkl = 0.;
          CrystalC_ijkl_Elastic();    // elastic part
	  CrystalC_ijkl_Plastic();    // plastic part	  

	  // add stress and moduli to corresponding averaged quantities
	  fsavg_ij.AddScaled(1./fNumGrain, fs_ij);
	  fcavg_ijkl.AddScaled(1./fNumGrain, fc_ijkl);
	}
    }

  // return averaged stress
  return fsavg_ij;
}

const dMatrixT& LocalCrystalPlast_C::c_ijkl()
{
  // fetch current element
  ElementCardT& element = CurrentElement();
  int intpt = 0;

  // recover aggregate data
  LoadAggregateData(element, intpt);

  // return averaged moduli
  return fcavg_ijkl;
}

void LocalCrystalPlast_C::UpdateHistory()
{
  // get element
  ElementCardT& element = CurrentElement();

  // update state at each crystal
  int intpt = 0;
  for (int igrn = 0; igrn < fNumGrain; igrn++)
    {
      // recover local data
      LoadCrystalData(element, intpt, igrn);
      
      // update state
      fFpi_n      = fFpi;
      fDGamma_n   = fDGamma;
      fHardening->UpdateHistory();
    }
}

void LocalCrystalPlast_C::ResetHistory()
{
  // get element
  ElementCardT& element = CurrentElement();

  // reset history at each crystal
  int intpt = 0;
  for (int igrn = 0; igrn < fNumGrain; igrn++)
    {
      // recover local data
      LoadCrystalData(element, intpt, igrn);
      
      // reset state
      fFpi      = fFpi_n;
      fDGamma   = fDGamma_n;
      fHardening->ResetHistory();
    }
}

int LocalCrystalPlast_C::NumOutputVariables() const {return kNumOutput;}

void LocalCrystalPlast_C::OutputLabels(ArrayT<StringT>& labels) const
{
  // allocate space for labels
  labels.Allocate(kNumOutput);

  // copy labels
  for (int i = 0; i < kNumOutput; i++)
    labels[i] = Labels[i];
}

void LocalCrystalPlast_C::ComputeOutput(dArrayT& output)
{
  if (CurrIP() == 0)
    {
      // gather element information
      ElementCardT& element = CurrentElement();
      int elem  = CurrElementNumber();
      int intpt = 0;

      // aggregate Cauchy stress
      LoadAggregateData(element, intpt);

      // aggregate Von Mises stress
      output[0] = sqrt(fsymmatx1.Deviatoric(fsavg_ij).ScalarProduct())/sqrt23;
      cerr << " S_eq = " << output[0] << endl;
      output[0] /= (fHardening->MaterialProperties())[1];  // Voce Hardening model

      // for hardness (?)
      output[1] = fIterCount;

      // compute texture of aggregate, if requested
      const int& step = fContinuumElement.FEManager().StepNumber();
      const int& nsteps = fContinuumElement.FEManager().NumberOfSteps();

      if (fmod(step, fODFOutInc) == 0 || step == nsteps)
	{
	  for (int igrn = 0; igrn < fNumGrain; igrn++)
	    {
	      // fetch crystal data 
	      LoadCrystalData(element, intpt, igrn);
	  
	      // texture: rotation tensor from fFe
	      PolarDecomp();
	  
	      // texture: compute new crystal orientation matrix
	      fmatx1.MultAB(fRe, fRotMat);

	      // texture: compute Euler angles (in radians)
	      dArrayT& angles = fangles[igrn];
	      fLatticeOrient->RotMatrixToAngles(fmatx1, angles);
	    }
	  
	  // write texture at center of element
	  fLatticeOrient->WriteTexture(elem, intpt, fNumGrain, step, fangles);
	}
    }
}

void LocalCrystalPlast_C::Print(ostream& out) const
{
  // inherited
  LocalCrystalPlast::Print(out);
}

void LocalCrystalPlast_C::PrintName(ostream& out) const
{
  // inherited
  LocalCrystalPlast::PrintName(out);

  // output model feature
  out << "    (Using one aggregate per element - center)\n";
}

/* PROTECTED MEMBER FUNCTIONS */

void LocalCrystalPlast_C::InitializeCrystalVariables()
{
  // initialize state at each element and ...
  for (int elem = 0; elem < NumElements(); elem++)
    {
      // get pointer to element elem
      ElementCardT& element = ElementCard(elem);
      int intpt = 0;

      // fetch aggregate quantities
      LoadAggregateData(element, intpt);
      
      // initialize averaged stress and moduli
      fsavg_ij = 0.;
      fcavg_ijkl = 0.;

      // ... at each crystal (center)
      for (int igrn = 0; igrn < fNumGrain; igrn++)
	{
	  // fetch crystal data 
	  LoadCrystalData(element, intpt, igrn);
	      
	  // fetch euler angles
	  dArrayT& angles = fEuler[elem](intpt, igrn);
	      
	  // storage rotation matrix from Euler angles
	  fLatticeOrient->AnglesToRotMatrix(angles, fRotMat);
	      
	  // plastic deformation gradients 
	  fFpi_n.Identity();
	  fFpi.Identity();
	  
	  // shear rates on slip systems
	  fDGamma_n = 0.;
	  fDGamma   = 0.;
	  
	  // elastic deformation gradient
	  fFe.Identity();

	  // crystal Cauchy stress
	  fs_ij = 0.;

	  // hardening variables
	  fHardening->InitializeHardVariables();
	}
    }
}

const dMatrixT& LocalCrystalPlast_C::DeformationGradient(const LocalArrayT& disp)
{
  return  DefGradientAtCenter(disp);
}

const dMatrixT& LocalCrystalPlast_C::DefGradientAtCenter(const LocalArrayT& disp)
{
  // derivatives dNa/dX 
  ComputeGDNa_C(fLocInitX, fLDNa, fGDNa);

  // displacement gradient dU/dX
  Jacobian(disp, fGDNa, fGradU);

  // deformation gradient F = dx/dX
  return FDContinuumT::F(fGradU);
}