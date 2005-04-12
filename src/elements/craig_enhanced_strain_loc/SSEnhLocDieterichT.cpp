#include "SSEnhLocDieterichT.h"
#include "ShapeFunctionT.h"
#include "SSSolidMatT.h"

#include "math.h"

using namespace Tahoe;

/* constructor */
SSEnhLocDieterichT::SSEnhLocDieterichT(const ElementSupportT& support):
  SSEnhLocCraigT(support),
  fDieterichBand(NULL)
  //fDieterichBand(static_cast<DieterichBandT*> (fBand))
{
  //fDieterichBand =  dynamic_cast<DieterichBandT*> (fBand);

  SmallStrainT::SetName("small_strain_enh_loc_dieterich");
}

/* implementation of the ParameterInterfaceT interface */
void SSEnhLocDieterichT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SSEnhLocCraigT::DefineParameters(list);

	/*PARAMETERS FOR Dieterich model*/
	list.AddParameter(fMu_star, "ref_friction_coeff_mu_star");
	list.AddParameter(fTheta_star, "ref_isv__theta_star");
	list.AddParameter(fV_star, "ref_velocity_v_star");
	list.AddParameter(fFrictionA, "friction_parameter_A");
	list.AddParameter(fFrictionB, "friction_parameter_B");
	list.AddParameter(fD_c, "evolution_paramter_D_c");
	list.AddParameter(fTheta_0, "initial_isv_theta_0");
}

void SSEnhLocDieterichT::TakeParameterList(const ParameterListT& list)
{
  /* inherited */
  SSEnhLocCraigT::TakeParameterList(list);


  /*PARAMETERS FOR Dieterich Model*/
  fMu_star = list.GetParameter("ref_friction_coeff_mu_star");
  fTheta_star = list.GetParameter("ref_isv__theta_star");
  fV_star = list.GetParameter("ref_velocity_v_star");
  fFrictionA = list.GetParameter("friction_parameter_A");
  fFrictionB = list.GetParameter("friction_parameter_B");
  fD_c = list.GetParameter("evolution_paramter_D_c");
  fTheta_0 = list.GetParameter("initial_isv_theta_0");
}

void SSEnhLocDieterichT::FormStiffness(double constK)
{

  if (fStrainDispOpt == kMeanDilBbar)
    cout << "Warning SSEnhLocCraigT::FormStiffness, b-bar integration not implemented for enhanced strain element, postlocalization results may be garbage.";

  if (!IsElementTraced() || !fBand->IsActive() )
    {
    /* form stiffness in standard way */
    SmallStrainT::FormStiffness(constK);

    }
  else //if already localized, use localized stiffness routine
    {
      //cout << "SSEnhLocDieterichT::FormStiffness\n";


	      cout << "Dieterich FormStiffness, fBand->JumpIncrement() = "
		   << fBand->JumpIncrement() << endl;

	/* matrix format */
	dMatrixT::SymmetryFlagT format = dMatrixT::kWhole;

	/* integrate element stiffness */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	int ndof = NumDOF();
	int nen = NumElementNodes();
	int nedof = nen*ndof;//# of element dof
	double k_zeta_zeta = 0.0, area = 0.0;
	dArrayT k_d_zeta(nedof), k_zeta_d(nedof);
	dArrayT k_d_zeta_work(nedof), k_zeta_d_work(nedof);
	k_d_zeta = 0.0;
	k_zeta_d = 0.0;

	dSymMatrixT gradActiveTensorFlowDir(ndof), dGdSigmaPlus(ndof), nTensorn(ndof);
	dMatrixT fLHSWork(nedof),fDfB((fCurrMaterial->ce_ijkl().Rows()),nedof);

	dGdSigmaPlus = FormdGdSigma(ndof, fBand->JumpIncrement(), fDieterichBand->DeltaTheta());
	nTensorn.Outer(fBand->Normal());
	nTensorn.ScaleOffDiagonal(2.0);
	/*dGdSigmaPlus.AddScaled(BigConstant(fBand->JumpIncrement(),fDieterichBand->DeltaTheta()),
	nTensorn); */

	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		area += (*Det)*(*Weight);
	        double scale = constK*(*Det++)*(*Weight++);
	
		//form ke_dd - assume elastic
		/* strain displacement matrix */
		if (fStrainDispOpt == kMeanDilBbar)
			Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
		else
			Set_B(fShapes->Derivatives_U(), fB);

		/* get D matrix */
		fD.SetToScaled(scale, fCurrMaterial->ce_ijkl());

		/* multiply b(transpose) * db, taking account of symmetry, */
		/* and accumulate in elstif */
		fDfB.MultAB(fD, fB);
		fLHSWork.MultATB(fB, fDfB);
		fLHS +=fLHSWork;
		/* update nMatrixT to get rid of fLHSwork */

		//form k_d_zeta
		gradActiveTensorFlowDir =
		FormGradActiveTensorFlowDir(ndof, CurrIP());

		fDfB.MultTx(gradActiveTensorFlowDir, k_d_zeta_work);
		k_d_zeta += k_d_zeta_work;

		//form k_zeta_d
		fDfB.MultTx(dGdSigmaPlus, k_zeta_d_work);
		k_zeta_d += k_zeta_d_work;

	}
	k_zeta_d *= 1.0/area;

	k_zeta_zeta = DdeltaGdJump(fBand->JumpIncrement(),fDieterichBand->DeltaTheta());
	//	k_zeta_zeta = DdeltaGdJumpGlobal(fBand->JumpIncrement(),fDieterichBand->DeltaTheta());


  //k_zeta_zeta *= 1.0/area;
  //k_zeta_zeta += fBand->EffectiveSoftening();
  k_zeta_zeta *= -1.0;

  fLHS.Outer(k_d_zeta, k_zeta_d, -1.0/k_zeta_zeta, dMatrixT::kAccumulate);
    }
}

double SSEnhLocDieterichT::CalculateJumpIncrement()
{
  if(!IsBandActive())
    return 0.0;


  // cout << "fBand->SlipDir =\n" << fBand->SlipDir() << endl;

  double newtonTol = 1.0e-12;
  double jumpIncrement = fDieterichBand->JumpIncrLast();
  if (jumpIncrement <= 0.0)
    jumpIncrement = 3.013904e-03; //think of better estimate

  double deltaTheta = 0.0;
  int newtonCounter = 0;
  int maxIter = 20;

  // what about coming off an elastic step? 
  double deltaG = DeltaG(jumpIncrement, deltaTheta);
  cout << "deltaG = " << deltaG << endl;

      // make this a relative tolerance
  while (fabs(deltaG) > newtonTol)
    {
      /* if too many iterations, stop and cut load step */
      if(newtonCounter++ >maxIter)
	{
	  cout << "SSEnhLocDieterichT::CalculateJumpIncrement, Newton iteration did not converge, cutting load step\n";
	  throw ExceptionT::kGeneralFail;
	}

      /* update jump increment via Newton iteration */
      //cout << "DdeltaGdJump = " << DdeltaGdJump(jumpIncrement, deltaTheta) << endl;

      jumpIncrement -= deltaG/DdeltaGdJump(jumpIncrement, deltaTheta);

      /* can't allow non-positive jump increment - in log fn creates
	 non-real answer */
      if (jumpIncrement <= 0.0)
	jumpIncrement = newtonTol;
      //cout << "deltaG = " << deltaG << endl;
      //cout << "DdeltaGdJump = " << DdeltaGdJump(jumpIncrement, deltaTheta) << endl;
      //cout << "jumpIncrement = " << jumpIncrement << endl;

      /*update increment of ISV */
      deltaTheta = DeltaTheta(jumpIncrement, deltaTheta);
      //cout << "deltaTheta = " << deltaTheta << endl;
      /* update regular strains? */

      /*reform DeltaG*/
      deltaG = DeltaG(jumpIncrement, deltaTheta);
      cout << "deltaG = " << deltaG << endl;
    }

  fBand -> StoreJumpIncrement(jumpIncrement);
  fDieterichBand -> StoreDeltaTheta(deltaTheta);

  return jumpIncrement;
}      

bool SSEnhLocDieterichT::IsBandActive()
{
   /*calculate average stress assuming no jump*/
  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();
    
  double area = 0.0;
  //double normalStress = 0.0;
  double shearStress = 0.0;

  /* collect incremental heat */
  //bool need_heat = fElementHeat.Length() == fShapes->NumIP();
  
  fShapes->TopIP();
  while (fShapes->NextIP())
    {
      /* strain displacement matrix */
      if (fStrainDispOpt == kMeanDilBbar)
	Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
      else
	Set_B(fShapes->Derivatives_U(), fB);
      
      /* B^T * Cauchy stress */
      dSymMatrixT strainIncr = fStrain_List [CurrIP()];
      strainIncr -= fStrain_last_List [CurrIP()];
      //    cout << "strainIncr =\n" << strainIncr << endl;
      dSymMatrixT stressIncr(NumSD());
      stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);
      stressIncr += fBand->Stress_List(CurrIP());
      
      area += (*Det)*(*Weight);
      double scale = (*Det++)*(*Weight++);

      //normalStress += scale * stressIncr.MultmBn(fBand->Normal(),fBand-> Normal());
      shearStress += scale * stressIncr.MultmBn(fBand->PerpSlipDir(),fBand-> Normal()); 
    }

  if (shearStress < 0.0)
    {
      /* align slip direction with shear stress direction to get 
     correct yield surface */
      fBand-> FlipSlipDir();
      shearStress *= -1.0;
    }

  //normalStress/= area;
  shearStress = shearStress/area;

  double neededCohesion = shearStress; // + normalStress * fLocalizedFrictionCoeff;
  cout << "ResidualCohesion = " << fBand->ResidualCohesion() << endl; 
  cout << "neededCohesion = " << neededCohesion << endl; 

  if (fBand-> ResidualCohesion() < neededCohesion)
    {
      fBand-> SetActive(true);
      cout << "Band is active\n";
      return true;
    }
  else
    {
      fBand -> SetActive(false);
      cout << "Band is NOT active\n";
      return false;
    }

}

/*--------------------------------math functions------------------*/

double SSEnhLocDieterichT::DeltaTheta(double jumpIncrement, double deltaTheta)
{
  double dt = ElementSupport().TimeStep();
  return  (dt * fD_c - fDieterichBand->Theta()*jumpIncrement)/ (
  fD_c + jumpIncrement);

  //return 1.0 - (fDieterichBand->Theta() + deltaTheta ) *jumpIncrement/(dt*fD_c);
  //return 1.0 - (fDieterichBand->Theta()) *jumpIncrement/(dt*fD_c);
}

double SSEnhLocDieterichT::DeltaG(double jumpIncrement, double deltaTheta)
{
  dSymMatrixT stressIncr = StressIncrOnBand(jumpIncrement);
  dSymMatrixT currStress = LastStressOnBand();
  currStress += stressIncr;

  dSymMatrixT dGdSigma = FormdGdSigma(NumDOF(), jumpIncrement, deltaTheta);
  //dGdSigma.ScaleOffDiagonal(0.5);

  double deltaG = dGdSigma.Dot(dGdSigma, currStress); 

  //double normalStress = currStress.MultmBn(fBand->Normal(), fBand->Normal());

  //deltaG += BigConstant(jumpIncrement, deltaTheta) * normalStress;
  
  //cout << "BigConstant = " << BigConstant(jumpIncrement, deltaTheta) << endl;

  //address end of cohesion
  deltaG -= (fBand->ResidualCohesion() + fBand->H_delta()*jumpIncrement);

  return deltaG;
}

dSymMatrixT SSEnhLocDieterichT::StressIncrOnBand(double jumpIncrement)
{
  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();      
  double area = 0.0;
  int ndof = NumDOF();
  dSymMatrixT gradActiveTensorFlowDir(ndof);
  dSymMatrixT avgStressIncr(ndof);
  avgStressIncr = 0.0;

  fShapes->TopIP();
  while (fShapes->NextIP())
    {
      double scale = (*Det++)*(*Weight++);
      area += scale;
       
      dSymMatrixT strainIncr = fStrain_List [CurrIP()];
      strainIncr -= fStrain_last_List [CurrIP()];

      gradActiveTensorFlowDir = FormGradActiveTensorFlowDir(ndof, CurrIP());
      gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);
      strainIncr.AddScaled(-1.0*jumpIncrement, gradActiveTensorFlowDir);

      dSymMatrixT stressIncr(NumSD());
      stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);

      avgStressIncr.AddScaled(scale, stressIncr);
    }
  avgStressIncr/=area;

  return avgStressIncr;
}

dSymMatrixT SSEnhLocDieterichT::LastStressOnBand()
{
  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();      
  double area = 0.0;
  dSymMatrixT avgStress(NumDOF());
  avgStress = 0.0;

  fShapes->TopIP();
  while (fShapes->NextIP())
    {
      double scale = (*Det++)*(*Weight++);
      area += scale;

      avgStress.AddScaled(scale, fBand->Stress_List(CurrIP()));
    }
  avgStress/=area;

  return avgStress;
}

//rename
//this is currently in vector form (doubled off-diags)
// - may want to change to matrix
dSymMatrixT SSEnhLocDieterichT::AvgStrainRelaxation(double jumpIncrement)
{
  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();      
  double area = 0.0;
  int ndof = NumDOF();
  dSymMatrixT gradActiveTensorFlowDir(NumDOF());
  dSymMatrixT avgStrain(NumDOF());
  avgStrain = 0.0;

  fShapes->TopIP();
  while (fShapes->NextIP())
    {
      double scale = (*Det++)*(*Weight++);
      area += scale;
      //gradActiveTensorFlowDir = FormGradActiveTensorFlowDir(ndof, CurrIP());
      avgStrain.AddScaled(scale,FormGradActiveTensorFlowDir(ndof, CurrIP()));
    }
  //avgStrain*=jumpIncrement/area;
  avgStrain /= area;

  return avgStrain;
}


double SSEnhLocDieterichT::BigConstant(double jumpIncrement, double deltaTheta) 
{
  /*  
  cout << "bc fBand = " << fBand << flush << endl;  
  cout << "bc fDieterich = " << fDieterichBand  << flush << endl;
  cout << "bc Jump = " << fBand -> Jump() << flush << endl;
  cout << "bc Jump = " << fDieterichBand -> Jump() << flush << endl;
  cout << "bc Theta = " << fDieterichBand -> Theta() << flush << endl;
  */

  double bc = fFrictionA * (1 - fDieterichBand->JumpIncrLast()/jumpIncrement);
  bc += fFrictionB * ( 1 - fDieterichBand->Theta()/(fDieterichBand->Theta() +
  deltaTheta));

  //cout << "bc = " << bc << endl;

  return bc;
}

double SSEnhLocDieterichT::DdeltaGdJump(double jumpIncr, double deltaTheta)
{
  return DdeltaGdJumpAtConstTheta(jumpIncr, deltaTheta) +
  DdeltaGdTheta(jumpIncr, deltaTheta) * DThetaDJump(jumpIncr, deltaTheta);

}

double SSEnhLocDieterichT::DdeltaGdJumpGlobal(double jumpIncr, double deltaTheta)
{
  return DdeltaGdJumpAtConstTheta(jumpIncr, deltaTheta) +
  DdeltaGdThetaGlobal(jumpIncr, deltaTheta) * DThetaDJump(jumpIncr, deltaTheta);
}

double SSEnhLocDieterichT:: DdeltaGdJumpAtConstTheta(double jumpIncrement,
double deltaTheta) 
{
  
  dSymMatrixT avgStrainRelaxation = AvgStrainRelaxation(jumpIncrement);
  //currently in vector form - make into matrix
  avgStrainRelaxation.ScaleOffDiagonal(0.5);
  //avgStrainRelaxation *= -1.0;

  dSymMatrixT dGdSigmaPlus = FormdGdSigma(NumDOF(), jumpIncrement, deltaTheta);
  dGdSigmaPlus.ScaleOffDiagonal(0.5);

  
  dSymMatrixT nTensorn(NumDOF());
  nTensorn.Outer(fBand->Normal());

  /*
  dGdSigmaPlus.AddScaled(BigConstant(jumpIncrement, deltaTheta),
  nTensorn);
  */

  dSymMatrixT avgStressRelaxation(NumDOF());
  avgStressRelaxation.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(),
  avgStrainRelaxation);
  double answer = -1.0 * dGdSigmaPlus.ScalarProduct(avgStressRelaxation);

  answer -= fBand->H_delta();



  dSymMatrixT stressIncr = StressIncrOnBand(jumpIncrement);
  dSymMatrixT currStress = LastStressOnBand();
  currStress += stressIncr;
  double normalStress = currStress.MultmBn(fBand->Normal(),
  fBand->Normal());
  //double normalStressIncr = stressIncr.MultmBn(fBand->Normal(), fBand->Normal());

  answer += fFrictionA * normalStress/( jumpIncrement);
  /*
  answer += fFrictionA * normalStressIncr/jumpIncrement;
  */
  return answer;
}

double SSEnhLocDieterichT::DdeltaGdTheta(double jumpIncrement, double deltaTheta) 
{
    dSymMatrixT stressIncr = StressIncrOnBand(jumpIncrement);
  dSymMatrixT currStress = LastStressOnBand();
  currStress += stressIncr;
  double normalStress = currStress.MultmBn(fBand->Normal(), fBand->Normal());
  double normalStressIncr = stressIncr.MultmBn(fBand->Normal(), fBand->Normal());

  double currentTheta = fDieterichBand->Theta() + deltaTheta;

  //    return fFrictionB  *normalStress/(currentTheta); 

    
 return fFrictionB  * normalStress/currentTheta; 
}

double SSEnhLocDieterichT::DdeltaGdThetaGlobal(double jumpIncrement, double deltaTheta) 
{
    dSymMatrixT stressIncr = StressIncrOnBand(jumpIncrement);
  dSymMatrixT currStress = LastStressOnBand();
  currStress += stressIncr;
  double normalStress = currStress.MultmBn(fBand->Normal(), fBand->Normal());
  double normalStressIncr = stressIncr.MultmBn(fBand->Normal(), fBand->Normal());

  double currentTheta = fDieterichBand->Theta() + deltaTheta;

    return fFrictionB * ( 2*  fDieterichBand->Theta()) *
      normalStress/(currentTheta *currentTheta); 

  /* return fFrictionB * fDieterichBand->Theta()  *
      normalStress/(currentTheta *currentTheta) + fFrictionB *
      normalStressIncr/currentTheta; */
}

double SSEnhLocDieterichT::DThetaDJump(double jumpIncrement, double deltaTheta)
{
  double dt = ElementSupport().TimeStep();
    double work = (fD_c + jumpIncrement);

   return  -1.0 * fD_c * (fDieterichBand->Theta() + dt)/(work*work);

  //return  -1.0*(fDieterichBand->Theta() + deltaTheta)/(fD_c * dt);
  //return  -1.0*(fDieterichBand->Theta())/(fD_c * dt);
}

/*------end math functions----------------------------------------*/


BandT* SSEnhLocDieterichT::FormNewBand(dArrayT normal, dArrayT slipDir,
				   dArrayT perpSlipDir, dArrayT coords, double residCohesion, ArrayT<dSymMatrixT>
stressList)
{
return new DieterichBandT(normal, slipDir, perpSlipDir, coords,
				   fH_delta_0, residCohesion, stressList,
				   this, fTheta_0);
}

#if 0
void SSEnhLocDieterichT::CloseStep()
{

Top();
while (NextElement())
  {
    if (IsElementTraced())
      fDieterichBand->UpdateDieterichBand();
  }

 #if 0
  ArrayT<BandT*> locElements;

  /* loop over traced elements to update theta */
  TracedElements()->Ascending(locElements);
  int length = locElements.Length();

  for(int i = 0; i < length; i++)
    {
      //locElements[i]
      fDieterichBand->UpdateDieterichBand();
    }

  SSEnhLocCraigT::CloseStep();
#endif
}
#endif


void SSEnhLocDieterichT::LoadBand(int elementNumber)
{
  SSEnhLocCraigT::LoadBand(elementNumber);
  //change to static cast?
  fDieterichBand = dynamic_cast<DieterichBandT*> (fBand);
}

dSymMatrixT SSEnhLocDieterichT::FormdGdSigma(int ndof, double
jumpIncrement, double deltaTheta)
{
  dSymMatrixT dGdSigma(ndof), work(ndof);
  dMatrixT dGNonSym(ndof);

  dGNonSym.Outer(fBand->Normal(), fBand->PerpSlipDir());
  dGdSigma.Symmetrize(dGNonSym);

  double dt = ElementSupport().TimeStep();
  double frictionCoeff = fLocalizedFrictionCoeff;
  frictionCoeff += fFrictionA * log(jumpIncrement/(dt * fV_star));
  frictionCoeff += fFrictionB * log((fDieterichBand->Theta()+deltaTheta)/fTheta_star); 

  work.Outer(fBand->Normal());
  dGdSigma.AddScaled(frictionCoeff, work);

  /* use vector values for mult */
  dGdSigma.ScaleOffDiagonal(2.0);

  //cout << "dGdSigm =\n" << dGdSigma << endl; 
  //cout << "deltaTheta = " << deltaTheta << endl;

  return dGdSigma;
}
