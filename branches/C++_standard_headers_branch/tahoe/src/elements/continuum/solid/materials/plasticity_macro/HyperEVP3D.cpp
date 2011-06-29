/* $Id: HyperEVP3D.cpp,v 1.11 2004-07-15 08:29:14 paklein Exp $ */
#include "HyperEVP3D.h"
#include "NLCSolver.h"
#include "ElementCardT.h"

#include "Utils.h"
#include "SimplePowerLaw.h"

using namespace Tahoe;

const double sqrt32 = sqrt(3.0/2.0);

/* printing for debugging */
const bool Hyper_MESSAGES = false;

/* spatial dimensions of the problem */
const int kNSD = 3;

/* number of internal scalar variables */
const int kNumInternal = 5;

/* element output data */
const int kNumOutput = 4;
static const char* Labels[kNumOutput] = {"EQP_strain","VM_stress","Pressure","Hardness"};

HyperEVP3D::HyperEVP3D(ifstreamT& in, const FSMatSupportT& support) :
	ParameterInterfaceT("HyperEVP_3D"),
	EVPFDBaseT(in, support),  

  // elastic def gradients
  fFeTr (kNSD,kNSD),
  fFe   (kNSD,kNSD),

  // plastic def gradients 
  fFpi_n (kNSD,kNSD),
  fFpi   (kNSD,kNSD),
  fDFpi  (kNSD,kNSD),

  // sym tensors in interm config	
  fCeBarTr     (kNSD),
  fEeBarTr     (kNSD),
  fSigBarTr    (kNSD),
  fSigBarTrDev (kNSD),
  fSigBar      (kNSD),

  // 2nd order identity tensor
  fISym (kNSD),

  // spectral/polar decomposition
  fSpecD (kNSD),
  fEigs  (kNSD),
  fReTr  (kNSD,kNSD),
  fUeTr  (kNSD),
  
  // array for internal variables
  fInternal (kNumInternal),

  // incremental equivalent plastic strain
  fDEQP (1),

  // work spaces
  farray     (kNSD),
  fmatx1     (kNSD,kNSD),
  fRank4     (dSymMatrixT::NumValues(kNSD)),
  fIdentity4 (dSymMatrixT::NumValues(kNSD)),
  fsymmatx1  (kNSD),
  fsymmatx2  (kNSD)
{
  // set 2nd order unit tensor (sym matrix)
  fISym.Identity();

  // set 4th order identity tensor
  fIdentity4.ReducedIndexI();
}

HyperEVP3D::~HyperEVP3D() {} 

int HyperEVP3D::NumVariablesPerElement()
{
  int d_size = 0;
  int dim = dSymMatrixT::NumValues(kNSD);

  // variables per ip
  d_size += kNSD*kNSD;     // fFpi_n(t_n)
  d_size += kNSD*kNSD;     // fFpi(t)
  d_size += kNumInternal;  // fEQP_n(t_n), fEQP(t), fEQSig_n(t_n), fEQSig(t), fpressTr(t)
  d_size += dim;           // fs_ij(t)
  d_size += dim * dim;     // fc_ijkl(t)

  // total # variables per element
  d_size *= NumIP();

  return d_size;
}

const dSymMatrixT& HyperEVP3D::s_ij()
{
  // fetch current element and int point
  ElementCardT& element = CurrentElement();
  int intpt = CurrIP();

  // load element data
  LoadElementData(element, intpt);

  // compute state, stress and moduli 
  if (MaterialSupport().RunState() == GlobalT::kFormRHS)
    {
      // reset iteration counter to check NLCSolver
      if (intpt == 0) fIterCount = 0;

      // total deformation gradient
      // fFtot = F();
      //fFtot = DeformationGradient(fLocDisp);
      // fFtot = fContinuumElement.FEManager().DeformationGradient();

      //compute 3D total deformation gradient
      Compute_Ftot_3D(fFtot);

      // time step
      fdt = fFSMatSupport->TimeStep();

      // compute state (stress and state variables)
      IntegrateConstitutiveEqns();
    }

  // return cauchy stress
  return fs_ij;
}

const dMatrixT& HyperEVP3D::c_ijkl()
{
  // fetch current element and int point
  ElementCardT& element = CurrentElement();
  int intpt = CurrIP();

  // load element data
  LoadElementData(element, intpt);

  // return tangent moduli
  return fc_ijkl;
}

void HyperEVP3D::FormRHS(const dArrayT& deqp, dArrayT& rhs)
{
  // equivalent plastic strain
  fInternal[kEQP] = fInternal[kEQP_n] + deqp[0];

  // equivalent stress
  fInternal[kEQSig] = fEQSigTr - 3.*fmu*deqp[0];

  // residual
  rhs[0] = deqp[0] - fdt * fKineticEqn->f(fInternal[kEQSig], fInternal[kEQP]); 
}

void HyperEVP3D::FormLHS(const dArrayT& deqp, dMatrixT& lhs)
{
  // equivalent plastic strain
  fInternal[kEQP] = fInternal[kEQP_n] + deqp[0];

  // equivalent stress
  fInternal[kEQSig] = fEQSigTr - 3.*fmu*deqp[0];

  // local jacobian
  lhs(0,0) = 1. - fdt * ( fKineticEqn->DfDs(fInternal[kEQSig], fInternal[kEQP]) - 
		  3.*fmu*fKineticEqn->DfDsigma(fInternal[kEQSig], fInternal[kEQP]) );
}

void HyperEVP3D::UpdateHistory()
{
  // get element
  ElementCardT& element = CurrentElement();

  // update state at each integration point
  for (int intpt = 0; intpt < NumIP(); intpt++)
    {
      // recover local data
      LoadElementData(element, intpt);

      // effective stress (use in fwd gradient approx)
      fInternal[kEQSig_n] = fInternal[kEQSig];

      // update state
      fFpi_n = fFpi;
      fInternal[kEQP_n] = fInternal[kEQP];
    }
}

void HyperEVP3D::ResetHistory()
{
  // get element
  ElementCardT& element = CurrentElement();

  // reset state at each integration point
  for (int intpt = 0; intpt < NumIP(); intpt++)
    {
      // recover local data
      LoadElementData(element, intpt);

      // effective stress (use in fwd gradient approx)
      fInternal[kEQSig] = fInternal[kEQSig_n];

      // state
      fFpi = fFpi_n;
      fInternal[kEQP] = fInternal[kEQP_n];
    }
}

int HyperEVP3D::NumOutputVariables() const {return kNumOutput;}

void HyperEVP3D::OutputLabels(ArrayT<StringT>& labels) const
{
  // allocate space for labels
  labels.Dimension(kNumOutput);

  // copy labels
  for (int i = 0; i < kNumOutput; i++)
    labels[i] = Labels[i];
}

void HyperEVP3D::ComputeOutput(dArrayT& output)
{
  // gather element/integ point information
  ElementCardT& element = CurrentElement();
  int intpt = CurrIP();

  // load element data
  LoadElementData(element, intpt);

  // equivalent plastic strain
  output[0] = fInternal[kEQP];
  
  // Von Mises stress
  output[1] = fInternal[kEQSig];

  // pressure
  output[2] = fInternal[kpressTr];

  // iter counter for nlcsolver
  output[3] = fIterCount;

  if (Hyper_MESSAGES && intpt == 0)
     cerr << " step # " << fFSMatSupport->StepNumber()
          << " EQP-strain  "  << output[0] 
          << " VM-stress  "   << output[1] 
          << " pressure  "    << output[2] 
          << " iterCounter "  << (int)output[3] << endl;
}

GlobalT::SystemTypeT HyperEVP3D::TangentType() const
{
  return GlobalT::kSymmetric;
}

/* PROTECTED MEMBER FUNCTIONS */

int HyperEVP3D::GetNumberOfEqns()
{
  return fDEQP.Length();
}

void HyperEVP3D::SetKineticEquation()
{
  // read kinetic eqn code
  fInput >> fKinEqnCode;

  // select kinetic eqn model
  switch(fKinEqnCode)
    {
    case KineticEqnBase::kSimplePowLaw:
      fKineticEqn = new SimplePowerLaw(*this);
      break;
      
    default:
      throwRunTimeError("HyperEVP3D::SetKineticEquation: Bad fKinEqnCode");
    }
  if (!fKineticEqn) throwMemoryError("HyperEVP3D::SetKineticEquation");
}

void HyperEVP3D::InitializeVariables(ElementCardT& element)
{
	// initialize state at each integration point
	for (int intpt = 0; intpt < NumIP(); intpt++)
	{
		// load element data
		LoadElementData(element, intpt);

		// plastic deformation gradients
		fFpi_n.Identity();
		fFpi.Identity();
	      
		// scalar internal variables
		fInternal = 0.;

		// Cauchy stress
		fs_ij = 0.;

		// tangent moduli
		fc_ijkl = 0.;
	}
}

void HyperEVP3D::LoadElementData(ElementCardT& element, int intpt)
{
  // fetch internal variable array
  dArrayT& d_array = element.DoubleData();
  
  // decode
  int dim = dSymMatrixT::NumValues(kNSD);
  int block = 2*kNSD*kNSD + kNumInternal + dim + dim*dim;
  int dex = intpt*block;

  fFpi_n.Set   (kNSD,kNSD,    &d_array[dex                ]);     
  fFpi.Set     (kNSD,kNSD,    &d_array[dex += kNSD*kNSD   ]); 
  fInternal.Set(kNumInternal, &d_array[dex += kNSD*kNSD   ]); 
  fs_ij.Set    (kNSD,         &d_array[dex += kNumInternal]);
  fc_ijkl.Set  (dim,dim,      &d_array[dex += dim         ]);
}

void HyperEVP3D::IntegrateConstitutiveEqns()
{
  // step 1. trial elastic deformation gradient
  fFeTr.MultAB(fFtot, fFpi_n);

  // step 2. polar decomposition of fFeTr
  PolarDecomposition();

  // step 3. trial logarithmic elastic strain
  LogarithmicStrain();

  // step 4. trial stresses
  ElasticTrialStress();

  // check for inelastic process
  if ( fEQSigTr > (1.+1.e-6)*fKineticEqn->g(fInternal[kEQP_n]) )
    {
      // step 5. solve for state variable(s)
      Solve();

      // step 6. compute Cauchy stress
      CauchyStress();

      // step 7. plastic deformation gradient
      FPInverse();
  
      // step 8. elasto-plastic moduli
      TangentModuli();
    }
  else  // elastic process
    {
      // step 5. reset values
      fInternal[kEQP]   = fInternal[kEQP_n];
      fInternal[kEQSig] = fEQSigTr;
      fFpi = fFpi_n;
      
      // step 6. compute Cauchy stress
      CauchyStress();

      // step 7. elastic moduli (fEta = 1)
      ElasticModuli();
    }
}

void HyperEVP3D::PolarDecomposition()
{
  // polar decomposition of fFeTr (using true does not work!!)
  fSpecD.PolarDecomp(fFeTr, fReTr, fUeTr, false);

  // recover eigenvalues of right stretch tensor
  fEigs = fSpecD.Eigenvalues();

/*
  // elastic right Cauchy-Green tensor
  fCeBarTr.MultATA(fFeTr);

  // spectral decomposition of CeBarTr
  //fCeBarTr.PrincipalValues(fEigs);
  //fSpecD.SpectralDecomp(fCeBarTr, fEigs, false);
  fSpecD.SpectralDecomp_new(fCeBarTr, false);
  fEigs = fSpecD.Eigenvalues();

  // principal values of right stretch tensor
  for (int i = 0; i < kNSD; i++)
    {
      if (fEigs[i] < 0.) 
        throwRunTimeError("HyperEVP3D::PolarDecomposition: eigs < 0");
      fEigs[i] = sqrt(fEigs[i]);
    }

  // temporaries
  double i1 = fEigs[0] + fEigs[1] + fEigs[2];
  double i2 = fEigs[0]*fEigs[1] + fEigs[0]*fEigs[2] + fEigs[1]*fEigs[2];
  double i3 = fEigs[0]*fEigs[1]*fEigs[2];
  double D = i1*i2 - i3;
  if (D < 0.0) throwRunTimeError("HyperEVP3D::PolarDecomposition: D < 0");
  
  // elastic right stretch tensor
  fsymmatx1.MultAB(fCeBarTr, fCeBarTr);
  fUeTr.SetToCombination(-1., fsymmatx1, (i1*i1-i2), fCeBarTr, i1*i3, fISym);
  fUeTr /= D;

  // inverse of fUeTr
  fsymmatx1.SetToCombination(1., fCeBarTr, -i1, fUeTr, i2, fISym);
  fsymmatx1 /= i3;
  fsymmatx1.ToMatrix(fmatx1);

  // rotation tensor
  fReTr.MultAB(fFeTr, fmatx1); 
*/
}

void HyperEVP3D::LogarithmicStrain()
{
  // eigenvalues of logarithmic strain
  for (int i = 0; i < kNSD; i++)
    farray[i] = log(fEigs[i]);
 
  // logarithmic strain (fEeBarTr = log(fUeTr))
  fEeBarTr = fSpecD.EigsToRank2(farray);
}

void HyperEVP3D::ElasticTrialStress()
{
  // trial rotated neutralized elastic stress
  fSigBarTr.SetToCombination(2.*fmu, fEeBarTr, flambda*fEeBarTr.Trace(), fISym);

  // deviatoric and volumetric parts of trial stress
  fSigBarTrDev.Deviatoric(fSigBarTr);
  fInternal[kpressTr] = -fSigBarTr.Trace() / 3.;

  // trial equivalent stress
  fEQSigTr = sqrt32*sqrt(fSigBarTrDev.ScalarProduct());
}

void HyperEVP3D::Solve()
{
  int ierr = 0;

  // forward gradient estimate for incremental equiv plastic strain
  ForwardGradientEstimate();

  // solve for incremental shear strain
  fSolver->Solve(fSolverPtr, fDEQP, ierr);

  if (ierr == 1) {
    writeMessage("HyperEVP3D::Solve: Convergence problems");
    writeWarning("HyperEVP3D::Solve: Will use unconverged DEQP");
    //    throw ExceptionT::kGeneralFail;	
  }

  // update iteration count from NLCSolver
  fIterCount = max(fIterCount, fSolver->GetIterationCount());
}

void HyperEVP3D::CauchyStress()
{
  // radial return factor
  fEta = fInternal[kEQSig] / fEQSigTr;

  // rotated neutralized stress
  fSigBar.SetToCombination(fEta, fSigBarTrDev, -fInternal[kpressTr], fISym);

  // Cauchy stress
  fs_ij.MultQBQT(fReTr, fSigBar);
  fs_ij *= exp(fInternal[kpressTr]/fbulk);
}

void HyperEVP3D::FPInverse()
{
  // incremental plastic deformation gradient 
  double d = fEigs[0]*fEigs[1]*fEigs[2];
  d = pow(d, 1./3.);
  for (int i = 0; i < kNSD; i++)
    farray[i] = pow(fEigs[i]/d, fEta-1.);
  fsymmatx1 = fSpecD.EigsToRank2(farray);
  fsymmatx1.ToMatrix(fDFpi);

  // current plastic deformation gradient
  fFpi.MultAB(fFpi_n, fDFpi);

  // check/normalize plastic deformation gradient
  //  if (fFpi.Det() <= 0.0) 
  //    throwRunTimeError("HyperEVP3D::FPInverse: det(Fpi) < 0");
  //  fFpi /= pow(fFpi.Det(), 1./3.);
}

/*void HyperEVP3D::FPInverse()
{
  // eigenvalues of fUe
  double d = fEigs[0]*fEigs[1]*fEigs[2];
  d = pow(d, (fEta-1.)/3.);
  for (int i = 0; i < kNSD; i++)
    farray[i] = pow(fEigs[i], fEta) / d;

  // current elastic deformation gradient 
  fsymmatx1 = fSpecD.EigsToRank2(farray);
  fFe.MultAB(fReTr, fsymmatx1);

  // current plastic deformation gradient
  fFpi.MultAB(fFtot.Inverse(), fFe);
}*/

void HyperEVP3D::ForwardGradientEstimate()
{
  // some references
  double& EQP_n   = fInternal[kEQP_n];
  double& EQSig_n = fInternal[kEQSig_n];

  // set effective stress at begining of plastic process
  if (EQP_n <= 1.e-8) EQSig_n = fKineticEqn->g(EQP_n);

  // kinetic equation functions at t_n
  double f_n      = fKineticEqn->f(EQSig_n, EQP_n);
  double dfdsig_n = fKineticEqn->DfDsigma(EQSig_n, EQP_n);
  double dfdeqp_n = fKineticEqn->DfDs(EQSig_n, EQP_n);

  // forward gradient estimate for DEQP
  fDEQP[0] = fdt * (f_n + dfdsig_n * (fEQSigTr - EQSig_n));
  fDEQP[0] /= (1. - fdt * (-3.*fmu*dfdsig_n + dfdeqp_n)) ;
}

/* consistent tangent based on updated Lagrangian procedure */
void HyperEVP3D::TangentModuli()
{
  // elastic-like terms
  ElasticModuli();

  // preliminary computations for inelastic contribution
  double dfdsig = fKineticEqn->DfDsigma(fInternal[kEQSig], fInternal[kEQP]);
  double dfdeqp = fKineticEqn->DfDs(fInternal[kEQSig], fInternal[kEQP]);
  double c = (1.-fdt*dfdeqp)/(1.-fdt*(dfdeqp-3.*fmu*dfdsig));

  // inelastic part: term N(x)N
  fsymmatx1.SetToScaled(sqrt32/fEQSigTr, fSigBarTrDev);
  fsymmatx2.MultQBQT(fReTr, fsymmatx1);
  fRank4.Outer(fsymmatx2, fsymmatx2);
  fc_ijkl.AddScaled(-2.*fmu*(fEta-c), fRank4); 
}

void HyperEVP3D::ElasticModuli()
{
  // first term I_ijkl 
  fc_ijkl.SetToScaled(2.*fmu*fEta, fIdentity4);

  // second term I(x)I
  fRank4.Outer(fISym, fISym);
  fc_ijkl.AddScaled((fbulk-2./3.*fmu*fEta), fRank4);
}
