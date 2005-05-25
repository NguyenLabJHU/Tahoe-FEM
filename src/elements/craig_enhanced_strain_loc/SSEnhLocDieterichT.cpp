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
	list.AddParameter(fBeta_zeta,"time_integration_parameter_beta_zeta");
	list.AddParameter(fBeta_theta,"time_integration_parameter_beta_theta");        
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
  fBeta_zeta = list.GetParameter("time_integration_parameter_beta_zeta");
  fBeta_theta = list.GetParameter("time_integration_parameter_beta_theta");
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

      // cout << "Dieterich FormStiffness, fBand->JumpIncrement() = "
      //	   << fBand->JumpIncrement() << endl;
      
      double slipRate = fDieterichBand->SlipRate();
      double thetaNew = fDieterichBand->Theta();
      double jumpIncrement = fBand -> JumpIncrement();

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

	dSymMatrixT gradActiveTensorFlowDir(ndof), dGdSigma(ndof), nTensorn(ndof);
	dMatrixT fLHSWork(nedof),fDfB((fCurrMaterial->ce_ijkl().Rows()),nedof);

	dGdSigma = FormdGdSigma(ndof, slipRate, thetaNew);
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
		fDfB.MultTx(dGdSigma, k_zeta_d_work);
		k_zeta_d += k_zeta_d_work;

	}
	k_d_zeta *= fBeta_zeta * ElementSupport().TimeStep();

	k_zeta_d *= 1.0/area;

	k_zeta_zeta = DPhidSlipRate(slipRate, jumpIncrement, thetaNew);
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
  	{
 	   double dt = ElementSupport().TimeStep();	
       double jumpAtZeroVel = dt * (1.0 - fBeta_zeta) * fDieterichBand -> SlipRateLast(); 
 
 	   fBand -> StoreJumpIncrement(jumpAtZeroVel);
  	   fDieterichBand -> StoreTheta(fDieterichBand -> ThetaLast() 
  	   	+ dt * ((1 - fBeta_theta) * fDieterichBand -> ThetaRateLast() + fBeta_theta));
  	   fDieterichBand -> StoreSlipRate(0.0);
  
    	return jumpAtZeroVel;
	}

  cout << "Element # " << CurrElementNumber() << endl;
  // cout << "fBand->SlipDir =\n" << fBand->SlipDir() << endl;

  double newtonTol = 1.0e-10;
  //double dt = ElementSupport().TimeStep();
  
  double slipRate = 0.0;
  double jumpIncrement = JumpIncrement(slipRate);
  double thetaNew = ThetaNew(slipRate);
  int newtonCounter = 0;
  int maxIter = 20;

  cout << "jumpIncrement = " << jumpIncrement << ", thetaNew = " << thetaNew;

  // what about coming off an elastic step? 
  double yieldFn = Phi(slipRate, jumpIncrement, thetaNew);
  double yieldFn0 = yieldFn;
  cout << ", yieldFn = " << yieldFn << endl;

      // make this a relative tolerance
  while (fabs(yieldFn/yieldFn0) > newtonTol)
    {
      /* if too many iterations, stop and cut load step */
      if(newtonCounter++ >maxIter)
	{
	  cout << "SSEnhLocDieterichT::CalculateJumpIncrement, Newton iteration did not converge, cutting load step\n";
	  throw ExceptionT::kGeneralFail;
	}

      /* update jump increment via Newton iteration */
      slipRate -= yieldFn/DPhidSlipRate(slipRate, jumpIncrement, thetaNew);

      cout << "dPhiSlipRate = " << DPhidSlipRate(slipRate, jumpIncrement,
      					 thetaNew);
      cout << " slipRate = " << slipRate;

      /*update increment of ISV */
      jumpIncrement = JumpIncrement(slipRate);
      thetaNew = ThetaNew(slipRate);

      //cout << "jumpIncrement = " << jumpIncrement << endl;
      //cout << "thetaNew = " << thetaNew << endl;

         /* update regular strains? */

      /*reform DeltaG*/
      yieldFn = Phi(slipRate, jumpIncrement, thetaNew);
      cout << ", yieldFn = " << yieldFn  << endl;
    }


  fBand -> StoreJumpIncrement(jumpIncrement);
  fDieterichBand -> StoreTheta(thetaNew);
  fDieterichBand -> StoreSlipRate(slipRate);

  //cout << endl;
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
      
      double jumpAtZeroVel = ElementSupport().TimeStep() * (1.0 - fBeta_zeta) * fDieterichBand -> SlipRateLast(); 
      
      dSymMatrixT gradActiveTensorFlowDir = FormGradActiveTensorFlowDir(NumDOF(), CurrIP());
      gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);
      strainIncr.AddScaled(-1.0*jumpAtZeroVel, gradActiveTensorFlowDir);
      
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
  //cout << "ResidualCohesion = " << fBand->ResidualCohesion() << endl; 
  //cout << "neededCohesion = " << neededCohesion << endl; 

  if (fBand-> ResidualCohesion() < neededCohesion)
    {
      fBand-> SetActive(true);
      // cout << "Band is active\n";
      return true;
    }
  else
    {
      fBand -> SetActive(false);
      //cout << "Band is NOT active\n";
      return false;
    }

}

/*--------------------------------math functions------------------*/

double SSEnhLocDieterichT::JumpIncrement(double slipRate)
{
  return ElementSupport().TimeStep() * ((1-fBeta_zeta)*
	 fDieterichBand->SlipRateLast() + fBeta_zeta * slipRate);
}

double SSEnhLocDieterichT::ThetaNew(double slipRate)
{
  double dt = ElementSupport().TimeStep();

  return fD_c * (fDieterichBand->ThetaLast() + dt * ((1 - fBeta_theta) *
		 fDieterichBand -> ThetaRateLast() + fBeta_theta))
    / (fD_c + fBeta_theta* slipRate * dt);
}

double SSEnhLocDieterichT::Phi(double slipRate, double jumpIncrement, double thetaNew)
{
  dSymMatrixT stressIncr = StressIncrOnBand(jumpIncrement);
  dSymMatrixT currStress = LastStressOnBand();

  //cout << "currStress =\n" << currStress << endl;

  currStress += stressIncr;

  dSymMatrixT dPhidSigma = FormdGdSigma(NumDOF(), slipRate, thetaNew);
  //dGdSigma.ScaleOffDiagonal(0.5);

  double phi = dPhidSigma.Dot(dPhidSigma, currStress); 

  //cout << "dPhidSigma =\n" << dPhidSigma << endl;
  //cout << "currStress =\n" << currStress << endl;

  //cout << "phi = " << phi; 

  //address end of cohesion
  double newCohesion = fBand->ResidualCohesion() +
    fBand->H_delta()*jumpIncrement;

  if (newCohesion > 0.0)
    phi -= newCohesion;
  //else newCohesion is really 0.0, no need to subtract anything

  cout << ", phi = " << phi << endl;

  //cout << ", newCohesion = " << newCohesion << ", ";

  return phi;
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


double SSEnhLocDieterichT::DPhidSlipRate(double slipRate, double jumpIncr, double thetaNew)
{

  dSymMatrixT dSigmadSlipRate = DSigmadSlipRate(jumpIncr);

  dSymMatrixT dPhidSigma = FormdGdSigma(NumDOF(), slipRate, thetaNew);
  //dGdSigma.ScaleOffDiagonal(0.5);

  double dPhi = dPhidSigma.Dot(dPhidSigma, dSigmadSlipRate); 

  //cout << "dPhi = " << dPhi << endl;

  dPhi += DmudSlipRate(slipRate, thetaNew) * NormalStress(jumpIncr);

  //cout << "dPhi = " << dPhi << endl;

  /* if there's still cohesion, take softening of band into account,
     otherwise, no effect */
  if (fBand->ResidualCohesion() + fBand->H_delta()*jumpIncr > 0.0)
    dPhi -= fBand -> H_delta() * DjumpdSlipRate();

  //cout << "dPhi = " << dPhi << endl;

  return dPhi;
}

dSymMatrixT SSEnhLocDieterichT::DSigmadSlipRate(double jumpIncrement)
{
  dSymMatrixT dSigmadSlipRate(NumSD());

  dSymMatrixT avgStrainRelaxation = AvgStrainRelaxation(jumpIncrement);
  avgStrainRelaxation.ScaleOffDiagonal(0.5);

  dSigmadSlipRate.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), avgStrainRelaxation);

  dSigmadSlipRate *= -1.0*DjumpdSlipRate();
  return dSigmadSlipRate;
}

double SSEnhLocDieterichT::DjumpdSlipRate()
{
  return fBeta_zeta * ElementSupport().TimeStep();
}

double SSEnhLocDieterichT::DmudSlipRate(double slipRate, double thetaNew)
{
  if (slipRate == 0.0)
    return fFrictionA/(2.0* fV_star) * exp((fMu_star + fFrictionB * log( thetaNew/fTheta_star))/fFrictionA);
  else
    {
      double arg = ArcSinhArg(slipRate, thetaNew);
      return arg * ( fFrictionA /slipRate + fFrictionB *
      DthetadSlipRate(slipRate) / thetaNew) / sqrt(1 + arg*arg);
    }
}

double SSEnhLocDieterichT::ArcSinhArg(double slipRate, double theta)
{
  return slipRate/(2.0* fV_star) * exp((fMu_star + fFrictionB * log( theta/fTheta_star))/fFrictionA);
}

double SSEnhLocDieterichT::DthetadSlipRate(double slipRate)
{
  double dt = ElementSupport().TimeStep();

  double numerator = fDieterichBand-> ThetaLast();
  numerator += dt * ( ( 1 - fBeta_theta) * fDieterichBand -> ThetaRateLast()
		      + fBeta_theta); 
 numerator *= -1.0 * fBeta_theta * dt * fD_c;

 double sqrtDenom = (fD_c + fBeta_theta * slipRate * dt);

 return numerator/(sqrtDenom * sqrtDenom);  
}

double SSEnhLocDieterichT::NormalStress(double jumpIncr)
{
  dSymMatrixT currStress = LastStressOnBand();
  currStress += StressIncrOnBand(jumpIncr);

  return currStress.MultmBn(fBand->Normal(), fBand->Normal());
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

void SSEnhLocDieterichT::CloseStep()
{
  /* inhertied */
  SSEnhLocCraigT::CloseStep();

  if (LocalizationHasBegun())
    {
      Top();
      while (NextElement())
	if (IsElementTraced())
	  {
	    fDieterichBand -> UpdateThetaRate(fD_c);
	  }
    }
}



void SSEnhLocDieterichT::LoadBand(int elementNumber)
{
  SSEnhLocCraigT::LoadBand(elementNumber);
  //change to static cast?
  fDieterichBand = dynamic_cast<DieterichBandT*> (fBand);
}

double SSEnhLocDieterichT::FrictionCoeff(double slipRate, double theta)
{
  return fFrictionA * arcsinh(ArcSinhArg(slipRate, theta));
}

double SSEnhLocDieterichT::arcsinh(double arg)
{
  return log(arg + sqrt(arg*arg + 1));
}

dSymMatrixT SSEnhLocDieterichT::FormdGdSigma(int ndof, double
slipRate, double thetaNew)
{
  dSymMatrixT dGdSigma(ndof), work(ndof);
  dMatrixT dGNonSym(ndof);

  dGNonSym.Outer(fBand->Normal(), fBand->PerpSlipDir());
  dGdSigma.Symmetrize(dGNonSym);

  //double dt = ElementSupport().TimeStep();
  double frictionCoeff = FrictionCoeff(slipRate, thetaNew);
 
  work.Outer(fBand->Normal());
  dGdSigma.AddScaled(frictionCoeff, work);

  /* use vector values for mult */
  dGdSigma.ScaleOffDiagonal(2.0);

  //cout << "dGdSigm =\n" << dGdSigma << endl; 
  //cout << "deltaTheta = " << deltaTheta << endl;

  return dGdSigma;
}
