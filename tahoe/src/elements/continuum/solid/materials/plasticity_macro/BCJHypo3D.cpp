/* $Id: BCJHypo3D.cpp,v 1.20 2004-07-15 08:29:14 paklein Exp $ */
#include "BCJHypo3D.h"
#include "NLCSolver.h"
#include "ElementCardT.h"

#include "Utils.h"
#include "BCJKineticEqn.h"
#include "MaterialSupportT.h"

using namespace Tahoe;

const double sqrt32 = sqrt(3.0/2.0);

/* printing for debugging */
const bool BCJ_MESSAGES = true;
const int IPprint = -1;

/* spatial dimensions of the problem */
const int kNSD = 3;

/* number of internal scalar variables */
const int kNumInternal = 3;
const int kNumEQValues = 5;

/* number of material properties */
const int kNumMatProp = 6;

/* element output data */
const int kNumOutput = 4;
static const char* Labels[kNumOutput] = {"EQP","VMISES","PRESS","KAPPA"};
//const int kNumOutput = 5;
//static const char* Labels[kNumOutput] = {"EQP","VMISES","PRESS","KAPPA","EQPDOT"};

BCJHypo3D::BCJHypo3D(ifstreamT& in, const FSMatSupportT& support) :
	ParameterInterfaceT("BCJHypo_3D"),
  EVPFDBaseT(in, support),  

  // some constants
  fNumInternal (kNumInternal),   // fDEQP, fALPH, fKAPP
  fNumEQValues (kNumEQValues),   // fEQP_n, fEQP, fEQXi_n, fEQXi, fPress

  // array for material parameters
  fMatProp (kNumMatProp),

  // deformation gradients
  fFtot_n (kNSD),
  fF      (kNSD),
  fFr     (kNSD),

  // tensors at t_n rotated to current configuration 
  fsigma_n (kNSD),
  falpha_n (kNSD),

  // trial stresses
  fSigTr    (kNSD),
  fSigTrDev (kNSD),
  fXiTr     (kNSD),

  // trial backstress
  falphaTr  (kNSD),

  // incremental strains
  fDEBar (kNSD),
  fDE    (kNSD),

  // other helpful tensors
  fA        (kNSD),
  fX        (kNSD),
  fUnitNorm (kNSD),
  fUnitM    (kNSD),

  // 2nd order identity tensor
  fISym (kNSD),
  fI    (kNSD, kNSD),

  // spectral/polar decomposition
  fSpecD (kNSD),
  fEigs  (kNSD),
  fC     (kNSD),
  fU     (kNSD),
  fR     (kNSD,kNSD),
  
  // internal variables
  fInternal_n (fNumInternal),
  fInternal   (fNumInternal),
  fInt_save   (fNumInternal),
  fEQValues   (fNumEQValues),
  fInternalTr (fNumInternal),

  // work spaces
  fRHS       (fNumInternal),
  fLHS       (fNumInternal),
  farray     (fNumInternal),
  fmatx1     (kNSD,kNSD),
  fmatx2     (kNSD,kNSD),
  fmatx3     (kNSD,kNSD),
  fmatx4     (kNSD,kNSD),
  fRank4     (dSymMatrixT::NumValues(kNSD)),
  fIdentity4 (dSymMatrixT::NumValues(kNSD)),
  fsymmatx1  (kNSD)
{
  // set 2nd order unit tensor (sym matrix)
  fISym.Identity();
  fI.Identity();
  
  // set 4th order identity tensor
  fIdentity4.ReducedIndexI();

  // read temperature
  fInput >> fTheta;

  // constants for kinematic hardening properties
  // (C1...C6 read in BCJKineticEqn class)
  fInput >> fC7 >> fC8 >> fC9 >> fC10 >> fC11 >> fC12;

  // constants for isotropic hardening properties
  fInput >> fC13 >> fC14 >> fC15 >> fC16 >> fC17 >> fC18;

  // compute material parameters
  ComputeMaterialProperties(fTheta);
}

BCJHypo3D::~BCJHypo3D() {} 

int BCJHypo3D::NumVariablesPerElement()
{
  int d_size = 0;
  int dim = dSymMatrixT::NumValues(kNSD);

  // variables per ip
  d_size += dim;             // fs_ij_n
  d_size += dim;             // fs_ij
  d_size += dim;             // falph_ij_n
  d_size += dim;             // falph_ij
  d_size += fNumInternal;    // internal variables at t_n
  d_size += fNumInternal;    // internal variables at t
  d_size += fNumEQValues;    // equivalent quantities at t_n & t
  d_size += dim * dim;       // fc_ijkl(t)

  // total # variables per element
  d_size *= NumIP();

  return d_size;
}

const dSymMatrixT& BCJHypo3D::s_ij()
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

      //compute 3D total deformation gradient
      Compute_Ftot_3D(fFtot);
      Compute_Ftot_last_3D(fFtot_n);

      // compute state (stress and state variables)
      SolveState();
    }

  // return cauchy stress
  return fs_ij;
}

const dMatrixT& BCJHypo3D::c_ijkl()
{
  // fetch current element and int point
  ElementCardT& element = CurrentElement();
  int intpt = CurrIP();

  // load element data
  LoadElementData(element, intpt);

  // return tangent moduli
  return fc_ijkl;
}

void BCJHypo3D::FormRHS(const dArrayT& array, dArrayT& rhs)
{
  // preliminary computations
  ComputeInternalQntsRHS(array);

  // form residuals
  // 1. from evolution equation of cauchy stress: function G1
  rhs[0] = fKineticEqn->h(array[kDEQP]/fdt, array[kKAPP])
             - fEQXiTr + (3.*fmu + fMatProp[1]*(1.-fEta))*array[kDEQP];

  // 2. from evolution equation of alpha: function G2
  rhs[1] = array[kALPH] - (1.-fEta)*fXMag;

  // 3. from evolution equation of kappa: function G3
  rhs[2] = array[kKAPP] - fInternal_n[kKAPP]  
//             - (fMatProp[4]-fMatProp[3]*array[kKAPP])*array[kDEQP]
             - (fMatProp[4]-fMatProp[3]*array[kKAPP]*array[kKAPP])*array[kDEQP]
             + fMatProp[5]*fdt*array[kKAPP]*array[kKAPP];
}

void BCJHypo3D::FormLHS(const dArrayT& array, dMatrixT& lhs)
{
  // preliminary computations ("ComputeInternalQntsRHS" has been already called)
  ComputeInternalQntsLHS(array);

  // previous operations to evaluate local Jacobian
  fsymmatx1.SetToCombination(1., fX, -sqrt32*fMatProp[1]*(1.-fEta)*array[kDEQP], fA);
  fUnitNorm.ToMatrix(fmatx1);
  fX.ToMatrix(fmatx2);
  fUnitM.ToMatrix(fmatx3);
  fsymmatx1.ToMatrix(fmatx4);

  double cnx  = dMatrixT::Dot(fmatx1, fmatx2);
  double cmn  = dMatrixT::Dot(fmatx1, fmatx3);
  double cmxa = dMatrixT::Dot(fmatx3, fmatx4);

  // Jacobian
  // 1. d(G1)/d(DEQP), d(G1)/d(ALPH), d(G1)/(dKAPP)
  lhs(0,0) = fKineticEqn->DhDeqpdot(array[kDEQP]/fdt, array[kKAPP]) / fdt
             - 1./sqrt32*fDEtaDDEQP*cnx + (3.*fmu+fMatProp[1]*(1.-fEta));
  lhs(0,1) = - 1./sqrt32*fDEtaDALPH*cnx;
  lhs(0,2) = fKineticEqn->DhDs(array[kDEQP]/fdt, array[kKAPP]);

  // 2. d(G2)/d(DEQP), d(G2)/d(ALPH), d(G2)/(dKAPP)
  lhs(1,0) = fDEtaDDEQP*cmxa - sqrt32*(1.-fEta)*fMatProp[1]*cmn;
  lhs(1,1) = 1. + cmxa*fDEtaDALPH;
  lhs(1,2) = 0.0;

  // 3. d(G3)/d(DEQP), d(G3)/d(ALPH), d(G3)/(dKAPP)
//  lhs(2,0) = -fMatProp[4] + fMatProp[3]*array[kKAPP];
  lhs(2,0) = -fMatProp[4] + fMatProp[3]*array[kKAPP]*array[kKAPP];
  lhs(2,1) = 0.0;
//  lhs(2,2) = 1. + fMatProp[3]*array[kDEQP] + 2.*fMatProp[5]*fdt * array[kKAPP];
  lhs(2,2) = 1. + 2.*(fMatProp[3]*array[kDEQP] + fMatProp[5]*fdt) * array[kKAPP];
}

void BCJHypo3D::UpdateHistory()
{
  // get element
  ElementCardT& element = CurrentElement();

  // update state at each integration point
  for (int intpt = 0; intpt < NumIP(); intpt++)
    {
      // recover local data
      LoadElementData(element, intpt);

      // update state
      fs_ij_n     = fs_ij;
      falph_ij_n  = falph_ij;
      fInternal_n = fInternal;

      // update useful scalar values
      fEQValues[kEQXi_n] = fEQValues[kEQXi];
      fEQValues[kEQP_n]  = fEQValues[kEQP];
    }
}

void BCJHypo3D::ResetHistory()
{
  // get element
  ElementCardT& element = CurrentElement();

  // reset state at each integration point
  for (int intpt = 0; intpt < NumIP(); intpt++)
    {
      // recover local data
      LoadElementData(element, intpt);

      // reset state
      fs_ij     = fs_ij_n;
      falph_ij  = falph_ij_n;
      fInternal = fInternal_n;

      // reset useful scalar values
      fEQValues[kEQXi] = fEQValues[kEQXi_n];
      fEQValues[kEQP]  = fEQValues[kEQP_n];
    }
}

int BCJHypo3D::NumOutputVariables() const {return kNumOutput;}

void BCJHypo3D::OutputLabels(ArrayT<StringT>& labels) const
{
  // allocate space for labels
  labels.Dimension(kNumOutput);

  // copy labels
  for (int i = 0; i < kNumOutput; i++)
    labels[i] = Labels[i];
}

void BCJHypo3D::ComputeOutput(dArrayT& output)
{
  // gather element/integ point information
  ElementCardT& element = CurrentElement();
  int intpt = CurrIP();

  // load element data
  LoadElementData(element, intpt);

  // equivalent plastic strain
  output[0] = fEQValues[kEQP];
  
  // mises stress
  output[1] = sqrt(fsymmatx1.Deviatoric(fs_ij).ScalarProduct())*sqrt32;

  // pressure
  output[2] = fEQValues[kPress];

  // isotropic hardening variable
  output[3] = fInternal[kKAPP];
  
  // equivalent plastic strain rate
  //output[4] = fInternal[kDEQP]/fdt;

  // effective overstress
  //output[4] = fEQValues[kEQXi];

  // iter counter
  //output[5] = fIterCount;

  if (BCJ_MESSAGES && intpt == 0 && CurrElementNumber() == 0)
     cerr << " step # " << MaterialSupport().StepNumber()
          << " EQP  "   << fEQValues[kEQP]
          << " EQXi "   << fEQValues[kEQXi]
          << " PRESS "  << fEQValues[kPress]
          << " STRESS33 " << fs_ij(2,2)
	  << " VMISES " << output[1]
	  << " DEQP "   << fInternal[kDEQP]
	  << " ALPHA "  << fInternal[kALPH]
	  << " KAPPA "  << fInternal[kKAPP]
          << " ITERS "  << fIterCount << endl;
}

#if 0
  // print temperature and material constants
  out << "    Material constants for isotropic/kinematic hardening laws\n";
  out << "       Temperature [K] . . . . . . . . . . . . . = " << fTheta << "\n";
  out << "       C7  [1/MPa] . . . . . . . . . . . . . . . = " << fC7    << "\n";
  out << "       C8  [K] . . . . . . . . . . . . . . . . . = " << fC8    << "\n";
  out << "       C9  [MPa] . . . . . . . . . . . . . . . . = " << fC9    << "\n";
  out << "       C10 [MPa/K] . . . . . . . . . . . . . . . = " << fC10   << "\n";
  out << "       C11 [1/(MPa-s)] . . . . . . . . . . . . . = " << fC11   << "\n";
  out << "       C12 [K] . . . . . . . . . . . . . . . . . = " << fC12   << "\n";
  out << "       C13 [1/MPa] . . . . . . . . . . . . . . . = " << fC13   << "\n";
  out << "       C14 [K] . . . . . . . . . . . . . . . . . = " << fC14   << "\n";
  out << "       C15 [MPa] . . . . . . . . . . . . . . . . = " << fC15   << "\n";
  out << "       C16 [MPa/K] . . . . . . . . . . . . . . . = " << fC16   << "\n";
  out << "       C17 [1/(MPa-s)] . . . . . . . . . . . . . = " << fC17   << "\n";
  out << "       C18 [K] . . . . . . . . . . . . . . . . . = " << fC18   << "\n";

  // print material properties
  out << "       rd  [1/Mpa] . . . . . . . . . . . . . . . = " << sqrt32*fMatProp[0] << "\n";
  out << "       h   [MPa] . . . . . . . . . . . . . . . . = " << fMatProp[1] << "\n";
  out << "       rs  [1/(MPa-s)] . . . . . . . . . . . . . = " << sqrt32*fMatProp[2] << "\n";
  out << "       Rd  [1/MPa] . . . . . . . . . . . . . . . = " << fMatProp[3] << "\n";
  out << "       H   [MPa] . . . . . . . . . . . . . . . . = " << fMatProp[4] << "\n";
  out << "       Rs  [1/(MPa-s)] . . . . . . . . . . . . . = " << fMatProp[5] << "\n";
#endif

GlobalT::SystemTypeT BCJHypo3D::TangentType() const
{
  return GlobalT::kNonSymmetric;
}

/* PROTECTED MEMBER FUNCTIONS */

void BCJHypo3D::ComputeMaterialProperties(double theta)
{
  // kinematic hardening matl properties
  fMatProp[0] = fC7*exp(-fC8/theta);     // rd
//  fMatProp[1] = fC9 - fC10*theta;        // h
  fMatProp[1] = max(0.0, fC9 - fC10*theta);        // h
  fMatProp[2] = fC11*exp(-fC12/theta);   // rs

  // isotropic hardening matl properties
  fMatProp[3] = fC13*exp(-fC14/theta);   // Rd
//  fMatProp[4] = fC15 - fC16*theta;       // H
  fMatProp[4] = max(0.0, fC15 - fC16*theta);       // H
  fMatProp[5] = fC17*exp(-fC18/theta);   // Rs

  // local fix for rd & rs
  fMatProp[0] *= 1./sqrt32;
  fMatProp[2] *= 1./sqrt32;
}

int BCJHypo3D::GetNumberOfEqns()
{
  return fInternal.Length();
}

void BCJHypo3D::SetKineticEquation()
{
  // read kinetic eqn code
  fInput >> fKinEqnCode;

  // select kinetic eqn model
  switch(fKinEqnCode)
    {
    case KineticEqnBase::kBCJKinEqn:
      fKineticEqn = new BCJKineticEqn(*this);
      break;
      
    default:
      throwRunTimeError("BCJHypo3D::SetKineticEquation: Bad fKinEqnCode");
    }

  if (!fKineticEqn) throwMemoryError("BCJHypo3D::SetKineticEquation");
}

void BCJHypo3D::InitializeVariables(ElementCardT& element)
{
	// initialize state at each integration point
	for (int intpt = 0; intpt < NumIP(); intpt++)
	{
		// load element data
		LoadElementData(element, intpt);

		// Cauchy stress
		fs_ij_n = 0.;
		fs_ij   = 0.;
	      
		// backstress
		falph_ij_n = 0.;
		falph_ij   = 0.;

		// scalar internal variables
		fInternal_n = 0.;
		fInternal   = 0.;

		// other scalar values
		fEQValues = 0.;

		// tangent moduli
		fc_ijkl = 0.;
	}
}

void BCJHypo3D::LoadElementData(ElementCardT& element, int intpt)
{
  // fetch internal variable array
  dArrayT& d_array = element.DoubleData();
  
  // decode
  int dim = dSymMatrixT::NumValues(kNSD);
  int block = 4*dim + 2*fNumInternal + fNumEQValues + dim*dim;
  int dex = intpt*block;

  fs_ij_n.Set    (kNSD,         &d_array[dex                ]);     
  fs_ij.Set      (kNSD,         &d_array[dex += dim         ]); 
  falph_ij_n.Set (kNSD,         &d_array[dex += dim         ]);     
  falph_ij.Set   (kNSD,         &d_array[dex += dim         ]); 
  fInternal_n.Set(fNumInternal, &d_array[dex += dim         ]); 
  fInternal.Set  (fNumInternal, &d_array[dex += fNumInternal]); 
  fEQValues.Set  (fNumEQValues, &d_array[dex += fNumInternal]);
  fc_ijkl.Set    (dim,dim,      &d_array[dex += fNumEQValues]);
}

/* driver to solve state using subincrementation */
void BCJHypo3D::SolveState()
{
  // flag to track convergence of state
  bool stateConverged;

  // counters for subincrementation
  int subIncr = 1;
  int totSubIncrs = 1;

  // time step and relative deformation gradient
  fdt = MaterialSupport().TimeStep();
  fmatx1.Inverse(fFtot_n);
  fF.MultAB(fFtot, fmatx1);
  fFr = fF;

  // integrate BCJ constitutive equations equations
  IntegrateConstitutiveEqns(stateConverged, subIncr, totSubIncrs);

  // if converged -> return; else -> do subincrementation
  if (stateConverged) return;

  // loop for subincrementation procedure
  for(;;)
    {
      // increase number of subincrements if not converged
      if (!stateConverged)
        {
          subIncr = 2 * subIncr - 1;
          totSubIncrs = 2 * totSubIncrs;
          if (totSubIncrs > 128) {
              if (BCJ_MESSAGES) {
                 cout << "BCJHypo3D::SolveState: totSubIncrs > 128 \n";
                 cout << "  **will throw 'EBadJacobianDet' to force dtime decrease** \n";
	      }
              throw ExceptionT::kBadJacobianDet;
	  }
          if (subIncr > 1) fInternal = fInt_save;
	  if (BCJ_MESSAGES) cout << " BCJHypo3D::SolveState: will try " 
		                 << subIncr << "/" << totSubIncrs << "\n";
        }

      // if converged, adjust subincrements
      else if (subIncr < totSubIncrs)
        {
          if ((subIncr/2*2) == subIncr)
            {
              subIncr = subIncr / 2 + 1;
              totSubIncrs = totSubIncrs / 2;
            }
          else
            subIncr = subIncr + 1;
        }

      // succesful return for subincrementation
      else
	{
          if (BCJ_MESSAGES) cout << " BCJHypo3D::SolveState: SubIncrs worked!!! " 
	                         << subIncr << "/" << totSubIncrs << "\n";
          break;
	}
 
      // time step
      double tmp = (float)subIncr/(float)totSubIncrs;
      fdt = MaterialSupport().TimeStep() * tmp;
      
      // relative deformation gradient for subincrement
      fFr.SetToCombination((1.-tmp), fI, tmp, fF);

      // save current converged solution before getting next solution
      if (subIncr > 1 && stateConverged) {
         fInt_save = fInternal;
      }

      // integrate BCJ equations
      IntegrateConstitutiveEqns(stateConverged, subIncr, totSubIncrs);
    }
}

void BCJHypo3D::IntegrateConstitutiveEqns(bool& converged, int subIncr, 
		                          int totSubIncrs)
{
  // initialize flag to track convergence
  converged = true; 
  
  // step 1. polar decomposition of relative deformation gradient
  PolarDecomposition();

  // step 2. incremental total strain
  IncrementalStrain();

  // step 3. rotate tensorial state variables at t_n to current configuration
  RotateTensorialVariables();

  // step 4. trial stresses
  ElasticTrialStress();

  // check for inelastic process
  if ( fEQXiTr > (1.+1.e-6)*fKineticEqn->h(((fabs(fdt) > kSmall) ? fInternal_n[kDEQP]/fdt : 0.0),fInternalTr[kKAPP]) )
    {
      // step 5. forward gradient estimate
      if (subIncr == 1) ForwardGradientEstimate();
      
      // step 6. solve for state variables
      Solve(converged);
      if (!converged) {
	 if (BCJ_MESSAGES)
            cout << " Did not converged at elem # " << CurrElementNumber()
	         << ";  IP # " << CurrIP() << ";  subIncr/totSunIncrs = " 
	         << subIncr << "/" << totSubIncrs << endl;
	 return;
      }
      if (subIncr < totSubIncrs) return;
      
      // step 7. compute Cauchy stress and backstress
      UpdateStresses();

      // step 8. elasto-plastic moduli
      TangentModuli();
    }
  else  // elastic process
    {
      // step 5. Elastic update
      UpdateElasticProcess(subIncr, totSubIncrs);

      // step 6. elastic moduli
      ElasticModuli(fmu, fbulk);
    }

  // state has converged
  //converged = true;
}

void BCJHypo3D::PolarDecomposition()
{
  // polar decomposition
  fSpecD.PolarDecomp(fFr, fR, fU, false);

  // recover eigenvalues
  fEigs = fSpecD.Eigenvalues();

/*
  // right Cauchy-Green tensor
  fC.MultATA(fFr);

  // spectral decomposition of fC
  // fC.PrincipalValues(fEigs);
  // fSpecD.SpectralDecomp(fC, fEigs);
  fSpecD.SpectralDecomp_new(fC, false);
  fEigs = fSpecD.Eigenvalues();

  // principal values of right stretch tensor
  for (int i = 0; i < kNSD; i++)
    {
      if (fEigs[i] < 0.) 
        throwRunTimeError("BCJHypo3D::PolarDecomposition: eigs < 0");
      fEigs[i] = sqrt(fEigs[i]);
    }

  // temporaries
  double i1 = fEigs[0] + fEigs[1] + fEigs[2];
  double i2 = fEigs[0]*fEigs[1] + fEigs[0]*fEigs[2] + fEigs[1]*fEigs[2];
  double i3 = fEigs[0]*fEigs[1]*fEigs[2];
  double D = i1*i2 - i3;
  if (D < 0.0) throwRunTimeError("BCJHypo3D::PolarDecomposition: D < 0");
  
  // right stretch tensor
  fsymmatx1.MultAB(fC, fC);
  fU.SetToCombination(-1., fsymmatx1, (i1*i1-i2), fC, i1*i3, fISym);
  fU /= D;

  // inverse of fU
  fsymmatx1.SetToCombination(1., fC, -i1, fU, i2, fISym);
  fsymmatx1 /= i3;
  fsymmatx1.ToMatrix(fmatx1);

  // rotation tensor
  fR.MultAB(fFr, fmatx1); 
*/
}

void BCJHypo3D::IncrementalStrain()
{
  // eigenvalues of logarithmic strain
  for (int i = 0; i < kNSD; i++)
    fEigs[i] = log(fEigs[i]);
 
  // logarithmic strain (fDEBar = log(fU))
  fDEBar = fSpecD.EigsToRank2(fEigs);
}

void BCJHypo3D::RotateTensorialVariables()
{
  // rotate Cauchy stress and backstress at t_n
  fsigma_n.MultQBQT(fR, fs_ij_n);
  falpha_n.MultQBQT(fR, falph_ij_n);

  // rotate incremental strain
  fDE.MultQBQT(fR, fDEBar);
}

void BCJHypo3D::ElasticTrialStress()
{
  // trial stress (elastic predictor)
  fSigTr.SetToCombination(1., fsigma_n, 2.*fmu, fDE, 
			  flambda*fDE.Trace(), fISym);

  // deviatoric and volumetric parts of trial stress
  fSigTrDev.Deviatoric(fSigTr);
  fEQValues[kPress] = -fSigTr.Trace() / 3.;

  // trial state variables
  fInternalTr[kDEQP] = 0.0;
  double a = fdt*fMatProp[2];
  double b = 1.0;
  double c = -fInternal_n[kALPH];
  if (fabs(a) < 1.0e-16)
     fInternalTr[kALPH] = -c;
  else
     fInternalTr[kALPH] = (-b+sqrt(b*b-4.0*a*c))/(2.0*a);
  a = fdt*fMatProp[5];
  b = 1.0;
  c = -fInternal_n[kKAPP];
  if (fabs(a) < 1.0e-16)
     fInternalTr[kKAPP] = -c;
  else
     fInternalTr[kKAPP] = (-b+sqrt(b*b-4.0*a*c))/(2.0*a);

  // trial back stress
  double factor = 1.0/(1.0 + fMatProp[2]*fdt*fInternalTr[kALPH]);
  falphaTr.SetToScaled(factor, falpha_n);

  // trial (elastic) equivalent stress
  fXiTr.SetToCombination(1., fSigTrDev, -2./3., falphaTr);
  fEQXiTr = sqrt32*sqrt(fXiTr.ScalarProduct());
}

void BCJHypo3D::Solve(bool& converged)
{
  int ierr = 0;

  // solve for primary unknowns (internal variables)
  try { fSolver->Solve(fSolverPtr, fInternal, ierr); }
  catch (ExceptionT::CodeT code)
     {
       if (BCJ_MESSAGES)
          writeWarning("BCJHypo3D::Solve: exception caught at Solve() -> subincrementation");
       converged = false;
       return;
     }

  // check for problems in NLCSolver
  if ( ierr != 0 ) {
     if (BCJ_MESSAGES) writeWarning("BCJHypo3D::Solve: ierr!= 0 -> subincrementation");
     converged = false;
     return;
  }
 
  // check for negative values of internal variables
  if ( IsSolnVariableNegative() ) {
     if (BCJ_MESSAGES) {
        writeWarning("BCJHypo3D::Solve: negative solution variable -> subincrementation");
        cout << "   fInternal: " << fInternal << endl;
     }
     converged = false;
     return;
  }

  // update iteration count from NLCSolver
  fIterCount = max(fIterCount, fSolver->GetIterationCount());
}

void BCJHypo3D::UpdateStresses()
{
  // preliminary computations
  // ComputeInternalQntsRHS(fInternal);

  // current equivalent stress
  // fEQValues[kEQXi] = fKineticEqn->h(fInternal[kDEQP]/fdt, fInternal[kKAPP]);
  fEQValues[kEQXi] = fEQXiTr - (3*fmu+fMatProp[1]*(1.-fEta))*fInternal[kDEQP];

  // current equivalent plastic strain
  fEQValues[kEQP] = fEQValues[kEQP_n] + fInternal[kDEQP];

  // radial factor
  fradial = fEQValues[kEQXi]/fEQXiTr;

  // backstress
  falph_ij.SetToScaled((1.-fEta), fX);

  // Cauchy stress
  fs_ij.SetToCombination(fradial, fXiTr, 2./3., falph_ij, -fEQValues[kPress], fISym);
}

/* Consisten Tangent based on Updated Lagrangian Procedure */
void BCJHypo3D::TangentModuli()
{
  // preliminary computations
  //ComputeInternalQntsRHS(fInternal);
  //ComputeInternalQntsLHS(fInternal);

  FormLHS(fInternal, fLHS);
  dMatrixT Jaci = MatrixInversion(fLHS);

  ffactor = fradial + fMatProp[1]*(1.-fEta)*fInternal[kDEQP]/fEQXiTr;
  
  double tmpa = sqrt32*(2.*fmu);
  if (fXMag <= 1.e-16) fXMag = 1.0;
  double tmpb = 1.5*(2.*fmu)/fXMag*(ffactor-fradial);

  double a1 = tmpa*Jaci(0,0);
  double a2 = tmpa*Jaci(1,0);
  double b1 = tmpb*Jaci(0,1);
  double b2 = tmpb*Jaci(1,1);

  double a4 = a1*fDEtaDDEQP + a2*fDEtaDALPH;
  double b4 = b1*fDEtaDDEQP + b2*fDEtaDALPH;

  // elastic-like terms: 2*mu*f*(I) +(k-2/3*mu*f)*(1(x)1)
  ElasticModuli(fmu*ffactor, fbulk);

  // inelastic part: term  -c1*(N(x)N)
  double c = 2.*fmu*ffactor - 2.*fmu*(1.-sqrt32*a1);
  fRank4.Outer(fUnitNorm, fUnitNorm);
  fc_ijkl.AddScaled(-c, fRank4);

  // inelastic part: term  -c2*(N(x)A)
  c = 3.*fmu*fEQXiTr*b1;
  fRank4.Outer(fUnitNorm, fA);
  fc_ijkl.AddScaled(-c, fRank4);

  // inelastic part: term  -c3*(A(x)N)
  c = 2.*sqrt32*fmu*fInternal[kDEQP]*a4;
  fRank4.Outer(fA, fUnitNorm);
  fc_ijkl.AddScaled(-c, fRank4);

  // inelastic part: term  -c4*(A(x)A)
  c = 3.*fmu*fInternal[kDEQP]*fEQXiTr*b4;
  fRank4.Outer(fA, fA);
  fc_ijkl.AddScaled(-c, fRank4);
}

void BCJHypo3D::ElasticModuli(double mu, double bulk)
{
  // first term I_ijkl 
  fc_ijkl.SetToScaled(2.*mu, fIdentity4);

  // second term I(x)I
  fRank4.Outer(fISym, fISym);
  fc_ijkl.AddScaled((bulk-2./3.*mu), fRank4);
}

void BCJHypo3D::ForwardGradientEstimate()
{
  // residual at t_n (they should be different from zero)
  FormRHS(fInternal_n, fRHS);

  // inverse of Jacobian at t_n (MatrixInversion changes fLHS)
  FormLHS(fInternal_n, fLHS);
  dMatrixT Jaci = MatrixInversion(fLHS);

  // forward gradient estimate
  Jaci.Multx(fRHS, farray);
  fInternal.SetToCombination(1., fInternal_n, -1., farray);
}

void BCJHypo3D::UpdateElasticProcess(int subIncr, int totSubIncrs)
{
   // update state variables
   fInternal = fInternalTr;
   if (subIncr < totSubIncrs) return;

   // update backstress
   falph_ij  = falphaTr;

   //update useful scalar qnts
   fEQValues[kEQXi] = fEQXiTr;
   fEQValues[kEQP]  = fEQValues[kEQP_n];

   // update Cauchy stress
   fs_ij = fSigTr;
}

bool BCJHypo3D::IsSolnVariableNegative()
{
  // test for negative values of internal variables
  int i = 0;
  bool isNegative = false;
  while (i < fNumInternal && !isNegative)
    {
      isNegative = (fInternal[i] < 0.0);
      i++;
    }

  // check for zero values of internal variables
  if (!isNegative) {
     if (fInternal[kDEQP] < 1.e-16) fInternal[kDEQP] = 1.e-16;
//     if (fInternal[kALPH] < 1.e-16) fInternal[kALPH] = 1.e-16;
//     if (fInternal[kKAPP] < 1.e-16) fInternal[kKAPP] = 1.e-16;
  }

  return isNegative;
}

/* PRIVATE MEMBER FUNCTIONS */

void BCJHypo3D::ComputeInternalQntsRHS(const dArrayT& array)
{
  // factor eta
  fEta = fMatProp[0]*array[kALPH]*array[kDEQP] + fMatProp[2]*array[kALPH]*fdt;
  fEta = fEta/(1.+fEta);

  // "new" trial equivalente stress and unit normal in stress space
  fXiTr.SetToCombination(1., fSigTrDev, -2./3.*(1.-fEta), falpha_n);
  fEQXiTr = sqrt32*sqrt(fXiTr.ScalarProduct());
  fUnitNorm.SetToScaled(sqrt32/fEQXiTr, fXiTr);

  // factor (tensor) in integrated eqn for backstress: alpha=(1-eta)*X
  fX.SetToCombination(1., falpha_n, sqrt32*fMatProp[1]*array[kDEQP], fUnitNorm);
  fXMag = sqrt(fX.ScalarProduct());
//  if (fXMag <= 1.e-16) fXMag = 1.;  // to avoid NaN in fUnitM in first incr
}

void BCJHypo3D::ComputeInternalQntsLHS(const dArrayT& array)
{
  // derivatives of factor Eta
  fDEtaDDEQP = (1.-fEta)*(1.-fEta)*(fMatProp[0]*array[kALPH]);
  fDEtaDALPH = (1.-fEta)*(1.-fEta)*(fMatProp[0]*array[kDEQP]+fMatProp[2]*fdt);

  // unit normal M: M=X/||X||
//  fUnitM.SetToScaled(1./fXMag, fX);
  double tmp;
  if (fXMag <= 1.e-16) 
     tmp = 0.0;
  else
     tmp = 1.0/fXMag;
  fUnitM.SetToScaled(tmp, fX);

  // derivative factor (tensor) for unit normal: dN/d()=A*dEta/d()
  //  fRank4.Outer(fUnitNorm, fUnitNorm);
  //  fc_ijkl.SetToCombination(1., fIdentitiy4, -1., fRank4);  // fc_ijkl is dummy
  //  fsymmatx1.A_ijkl_B_kl(fc_ijkl, falpha_n);
  falpha_n.ToMatrix(fmatx1);
  fUnitNorm.ToMatrix(fmatx2);
  fsymmatx1.SetToCombination(1., falpha_n, -dArrayT::Dot(fmatx1,fmatx2), fUnitNorm);
  fA.SetToScaled(1./(sqrt32*fEQXiTr), fsymmatx1);
}
