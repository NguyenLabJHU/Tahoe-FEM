/* $Id: LocalCrystalPlastFp_C.cpp,v 1.10 2004-10-14 20:24:51 paklein Exp $ */
#include "LocalCrystalPlastFp_C.h"
#include "LatticeOrient.h"
#include "CrystalElasticity.h"
#include "VoceHardening.h"
#include "ElementCardT.h"

#include "Utils.h"
#include "ContinuumElementT.h" // needed for initial coordinates

using namespace Tahoe;

/* spatial dimensions of the problem */
const int kNSD = 3; 
const double sqrt23 = sqrt(2.0/3.0);

/* element output data */
const int kNumOutput = 3;
static const char* Labels[kNumOutput] = {"VM_stress", "IterNewton", "IterState"};

LocalCrystalPlastFp_C::LocalCrystalPlastFp_C(ifstreamT& in, const FSMatSupportT& support) :
	ParameterInterfaceT("local_crystal_plasticity_Fp_C"),
  LocalCrystalPlastFp(in, support),
  fLocInitX (ContinuumElement().InitialCoordinates()),
  //fNNodes   (fLocInitX.NumberOfNodes()),  /* what is it wrong here? */
  //fLNa      (1, fNNodes),
  //fLDNa     (NumSD(), fNNodes),
  //fGDNa     (NumSD(), fNNodes),
  fGradU    (NumSD())
{
  // number of nodes per element
  fNNodes = fLocInitX.NumberOfNodes();

  // allocate arrays
  fLNa.Dimension(1, fNNodes);
  fLDNa.Dimension(NumSD(), fNNodes);
  fGDNa.Dimension(NumSD(), fNNodes);

  // set shape functions and their derivatives in parent domain (at center)
  SetLocalShape_C(fLNa, fLDNa);
}

LocalCrystalPlastFp_C::~LocalCrystalPlastFp_C() { } 

int LocalCrystalPlastFp_C::NumVariablesPerElement()
{
  int d_size = 0;
  int dim = dSymMatrixT::NumValues(kNSD);

  // variables per crystal per ip
  d_size += kNSD*kNSD;                       // fRotMat (const)
  d_size += kNSD*kNSD;                       // fFp_n   (t_n)
  d_size += kNSD*kNSD;                       // fFe_n   (t_n)
  d_size += kNSD*kNSD;                       // fFp     (t)
  d_size += kNSD*kNSD;                       // fFe     (t)
  d_size += fNumSlip;                        // fDGamma (t)
  d_size += dim;                             // fs_ij   (t)
  d_size += fHardening->NumberOfVariables(); // Hard Vars at t_n and t

  // total # crystal variables per element (at all IP's)
  d_size *= 1 * fNumGrain;

  // averaged (aggregate) stress and moduli (at all IP's)
  d_size += 1 * dim;                   // fsavg_ij   (t)
  d_size += 1 * dim * dim;             // fcavg_ijkl (t)

  return d_size;
}

const dSymMatrixT& LocalCrystalPlastFp_C::s_ij()
{
  // fetch current element
  ElementCardT& element = CurrentElement();
  int intpt = 0;

  // recover aggregate data
  LoadAggregateData(element, intpt);

  // compute state and stress at center of element
  if (fFSMatSupport->RunState() == GlobalT::kFormRHS && CurrIP() == 0)
    {
      // reset iteration counter to check NLCSolver
      fIterCount = 0;
      fIterState = 0;

      // initialize average stress and moduli at center
      fsavg_ij = 0.0;

      // shape function derivatives dNa/dX in physical domain (at center) 
      ComputeGDNa_C(fLocInitX, fLDNa, fGDNa);

      // deformation gradient at center of element
      DeformationGradient_C(fLocDisp, fFtot);
      DeformationGradient_C(fLocLastDisp, fFtot_n);

      // crystal state/stress
      for (int igrn = 0; igrn < fNumGrain; igrn++)
	{
	  // recover local data
	  LoadCrystalData(element, intpt, igrn);
	  
          // elasticity matrix in Bbar configuration
          if (!fElasticity->IsIsotropic())
             {
                fElasticity->ComputeModuli(fRank4);
                FFFFC_3D(fcBar_ijkl, fRank4, fRotMat);
             }
          else
                fElasticity->ComputeModuli(fcBar_ijkl);


          // compute crystal Cauchy stress (elastic predictor at first iteration)
          if (fFSMatSupport->IterationNumber() <= -1)
             {
               // defomation gradient
               fMatx1.SetToCombination(1., fFtot, -1., fFtot_n);
               fFt.SetToCombination(1., fFtot_n, 100.0 / 100.0, fMatx1);

               // elastic tensors
               fFpi.Inverse(fFp_n);
               fFe.MultAB(fFt, fFpi);
               fCeBar.MultATA(fFe);

               // 2nd Piola Kirchhoff stress
               fEeBar.SetToCombination(0.5, fCeBar, -0.5, fISym);
               fSBar.A_ijkl_B_kl(fcBar_ijkl, fEeBar);

               // Cauchy stress
               fs_ij.MultQBQT(fFe, fSBar);
               fs_ij /= fFe.Det();
             }
          else
             {
               // Schmidt tensor in Bbar configuration (sample axes)
               for (int i = 0; i < fNumSlip; i++)
                   fZ[i].MultQBQT(fRotMat, fZc[i]);

               // compute crystal state
               SolveCrystalState();

               // compute crystal Cauchy stress
               CrystalS_ij();
             }

	  // add stress and moduli to corresponding averaged quantities
	  fsavg_ij.AddScaled(1./fNumGrain, fs_ij);
	}
    }

  // return averaged stress
  return fsavg_ij;
}

const dMatrixT& LocalCrystalPlastFp_C::c_ijkl()
{
  // fetch current element
  ElementCardT& element = CurrentElement();
  int intpt = 0;

  // load aggregate data
  LoadAggregateData(element, intpt);

  // compute moduli at center of element
  if (CurrIP() == 0) 
    {
      // initialize averaged moduli
      fcavg_ijkl = 0.0;

      // compute averaged moduli
      for (int igrn = 0; igrn < fNumGrain; igrn++)
         {
            // recover local data
            LoadCrystalData(element, intpt, igrn);

            // elasticity matrix in Bbar configuration
            if (!fElasticity->IsIsotropic())
               {
                  fElasticity->ComputeModuli(fRank4);
                  FFFFC_3D(fcBar_ijkl, fRank4, fRotMat);
               }
            else
                  fElasticity->ComputeModuli(fcBar_ijkl);

            // compute consistent tangent (elastic predictor at fisrt iteration)
            if (fFSMatSupport->IterationNumber() <= 0)
                {
                   // elastic crystal stiffness
                   FFFFC_3D(fc_ijkl, fcBar_ijkl, fFe);
                }
            else
                {
                   // Schmidt tensor in Bbar configuration (sample axes)
                   for (int i = 0; i < fNumSlip; i++)
                       fZ[i].MultQBQT(fRotMat, fZc[i]);

                   // elasto-plastic crystal moduli (consistent tangent)
                   CrystalC_ijkl();
                }

            // add moduli to corresponding averaged quantity
            fcavg_ijkl.AddScaled(1./fNumGrain, fc_ijkl);
         }
    }

  // return averaged moduli
  return fcavg_ijkl;
}

void LocalCrystalPlastFp_C::UpdateHistory()
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
      fFp_n = fFp;
      fFe_n = fFe;
      fHardening->UpdateHistory();
    }
}

void LocalCrystalPlastFp_C::ResetHistory()
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
      fFp = fFp_n;
      fFe = fFe_n;
      fHardening->ResetHistory();
    }
}

int LocalCrystalPlastFp_C::NumOutputVariables() const {return kNumOutput;}

void LocalCrystalPlastFp_C::OutputLabels(ArrayT<StringT>& labels) const
{
  // allocate space for labels
  labels.Dimension(kNumOutput);
	  
  // copy labels
  for (int i = 0; i < kNumOutput; i++)
    labels[i] = Labels[i];
}

void LocalCrystalPlastFp_C::ComputeOutput(dArrayT& output)
{
  if (CurrIP() == 0)
    {
      // gather element information
      ElementCardT& element = CurrentElement();
      int group = ContinuumElement().ElementGroupNumber();
      int elem  = CurrElementNumber();
      int intpt = 0;

      // aggregate Cauchy stress
      LoadAggregateData(element, intpt);

      // aggregate Von Mises stress
      output[0] = sqrt(fSymMatx1.Deviatoric(fsavg_ij).ScalarProduct())/sqrt23;
      //cerr << " S_eq = " << output[0] << endl;
      //output[0] /= (fHardening->MaterialProperties())[1];

      // compute averaged equivalent stress
      if (elem == 0) fAvgStress = 0.0;
      fAvgStress.AddScaled(1./NumElements(), fsavg_ij);
      // cout << " elem = " << elem << "   intpt = " << intpt << endl;
      // cout << "    fsavg_ij = " << endl << fsavg_ij << endl;
      // cout << "    fAvgStress = " << endl << fAvgStress << endl;
      if (elem == (NumElements()-1))
         cerr << " step # " << fFSMatSupport->StepNumber()
              << "    S_eq_avg = "
              << sqrt(fSymMatx1.Deviatoric(fAvgStress).ScalarProduct())/sqrt23
              << "    Savg_12 = " << fAvgStress(0,1) << endl;

      // iteration counter for nlcsolver and state
      output[1] = fIterCount;
      output[2] = fIterState;

      // compute texture of aggregate, if requested
	int step = fFSMatSupport->StepNumber();
	int nsteps = fFSMatSupport->NumberOfSteps();

      if (step % fODFOutInc == 0 || step == nsteps)
	{
	  for (int igrn = 0; igrn < fNumGrain; igrn++)
	    {
	      // fetch crystal data 
	      LoadCrystalData(element, intpt, igrn);
	  
	      // texture: rotation tensor from fFe
              fSpecD.PolarDecomp(fFe, fRe, fUe, false);
	  
	      // texture: compute new crystal orientation matrix
	      fMatx1.MultAB(fRe, fRotMat);

	      // texture: compute Euler angles (in radians)
	      dArrayT& angles = fangles[igrn];
	      fLatticeOrient->RotMatrixToAngles(fMatx1, angles);
	    }
	  
	  // write texture at center of element
	  fLatticeOrient->WriteTexture(group, elem, intpt, fNumGrain, step, fangles);
	}
    }
}

/* PROTECTED MEMBER FUNCTIONS */

void LocalCrystalPlastFp_C::InitializeCrystalVariables(ElementCardT& element)
{
      int elem = CurrElementNumber();
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
	      
          // plastic/elastic deformation gradients
          fFp_n.Identity();
          fFe_n.Identity();
          fFp.Identity();
          fFe.Identity();

          // shear rates on slip systems
          fDGamma = 0.;

          // crystal Cauchy stress
          fs_ij = 0.;

          // hardening variables
          fHardening->InitializeHardVariables();
	}
}

void LocalCrystalPlastFp_C::DeformationGradient_C(const LocalArrayT& disp, dMatrixT& F_3D)
{
  // displacement gradient dU/dX
  Jacobian(disp, fGDNa, fGradU);

  // deformation gradient
  int nsd = NumSD();
  if (nsd == 3)
    {
       F_3D.Identity();
       F_3D += fGradU;
    }
  else if (nsd == 2)
    {
      dMatrixT matx(nsd);
      matx.Identity();
      matx += fGradU;

      // expand total deformation gradient: 2D -> 3D (plane strain)
      F_3D.Rank2ExpandFrom2D(matx);
      F_3D(2, 2) = 1.0;
    }
  else 
      throw ExceptionT::kGeneralFail;
}
