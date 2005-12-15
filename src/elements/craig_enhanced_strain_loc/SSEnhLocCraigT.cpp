/* $Id: SSEnhLocCraigT.cpp,v 1.21 2005-12-15 18:39:07 cfoster Exp $ */
#include "SSEnhLocCraigT.h"
#include "ShapeFunctionT.h"
#include "SSSolidMatT.h"

/* materials lists */
#include "SSSolidMatList1DT.h"
#include "SSSolidMatList2DT.h"
#include "SSSolidMatList3DT.h"
#include "SSMatSupportT.h"
#include "ParameterContainerT.h"
#include "ModelManagerT.h"
#include "OutputSetT.h"

#include "DetCheckT.h"
#include <math.h>

using namespace Tahoe;

/*initialize static variables */
bool SSEnhLocCraigT::fLocalizationHasBegun = false;
double SSEnhLocCraigT::fDetAMin = 1.0e99;
int SSEnhLocCraigT::fLeastDetEle = -1;

/* constructor */
SSEnhLocCraigT::SSEnhLocCraigT(const ElementSupportT& support):
	SmallStrainT(support),
	fBand(NULL)
	//fLeastDetEle(-1)
{
	SmallStrainT::SetName("small_strain_enh_loc_craig");
	jump_out.open("jump.info");
	fEdgeOfBandElements.Free();
	fEdgeOfBandCoords.Free();
}

/* destructor */
/*
SSEnhLocCraigT::~SSEnhLocCraigT(void)
{}
*/

/* implementation of the ParameterInterfaceT interface */
void SSEnhLocCraigT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SmallStrainT::DefineParameters(list);

	/*PARAMETERS FOR ENHANCED STRAIN*/
	list.AddParameter(fH_delta_0, "Post-Localization_softening_parameter_H_Delta"); 
	list.AddParameter(fNoBandDilation, "Disallow_Dilation_on_Band");
	list.AddParameter(fLocalizedFrictionCoeff, "Localized_Friction_Coefficient");
	//list.AddParameter(fMultiBand, "Allow_Multiple_Bands");
	ParameterT multiBand(fMultiBand, "Allow_Multiple_Bands");
	multiBand.SetDefault(false);
	list.AddParameter(multiBand);
	
	
	list.AddParameter(fBVPType, "BVP_type");
	list.AddParameter(fFirstElementToLocalize, "First_element_to_localize");
	
	double dummy;
	
	list.AddParameter(dummy, "First_component_of_normal_if_prefailed");
	list.AddParameter(dummy, "Second_component_of_normal_if_prefailed");
	list.AddParameter(dummy, "First_component_of_slip_dir_if_prefailed");
	list.AddParameter(dummy, "Second_component_of_slip_dir_if_prefailed");

    //list.AddParameter(fPreFailedNormal[0], "First_component_of_normal_if_prefailed");
	//list.AddParameter(fPreFailedNormal[1], "Second_component_of_normal_if_prefailed");
	//list.AddParameter(fPreFailedSlipDir[0], "First_component_of_slip_dir_if_prefailed");
	//list.AddParameter(fPreFailedSlipDir[1], "Second_component_of_slip_dir_if_prefailed");
	
}

/* information about subordinate parameter lists */
void SSEnhLocCraigT::DefineSubs(SubListT& sub_list) const
{	
  /* inherited */
  SmallStrainT::DefineSubs(sub_list);
}

void SSEnhLocCraigT::TakeParameterList(const ParameterListT& list)
{
  const char caller[] = "SmallStrainT::TakeParameterList";
  /* strain displacement option before calling
     SolidElementT::TakeParameterList */
  int b = list.GetParameter("strain_displacement");
  fStrainDispOpt = (b == kStandardB) ? kStandardB : kMeanDilBbar;
  
  /* inherited */
  SolidElementT::TakeParameterList(list);

  /* dimension workspace */
  fGradU.Dimension(NumSD());        
  if (fStrainDispOpt == kMeanDilBbar) {
    fLocDispTranspose.Dimension(fLocDisp.Length());
    fMeanGradient.Dimension(NumSD(), NumElementNodes());
  }        
  
  /* offset to class needs flags */
  fNeedsOffset = fMaterialNeeds[0].Length();
  
  /* set material needs */
  for (int i = 0; i < fMaterialNeeds.Length(); i++)
    {      
      /* needs array */
      ArrayT<bool>& needs = fMaterialNeeds[i];
      
      /* resize array */
      needs.Resize(needs.Length() + 2, true);
      
      /* casts are safe since class contructs materials list */
      ContinuumMaterialT* pcont_mat = (*fMaterialList)[i];
      SSSolidMatT* mat = (SSSolidMatT*) pcont_mat;
      
      /* collect needs */
      needs[fNeedsOffset + kstrain     ] = true; //mat->Need_Strain();
      needs[fNeedsOffset + kstrain_last] = true; //mat->Need_Strain_last();

      /* consistency */
      needs[kNeedDisp] = needs[kNeedDisp] || needs[fNeedsOffset +
						   kstrain];
      needs[KNeedLastDisp] = needs[KNeedLastDisp] || needs[fNeedsOffset
							   + kstrain_last];
    }

  /* what's needed */
  bool need_strain = true;
  bool need_strain_last = true;
  for (int i = 0; i < fMaterialNeeds.Length(); i++) {
    const ArrayT<bool>& needs = fMaterialNeeds[i];
    need_strain = need_strain || needs[fNeedsOffset + kstrain];
    need_strain_last = need_strain_last || needs[fNeedsOffset +
						 kstrain_last];
  }
  
  /* allocate strain list */
  if (need_strain) {
    fStrain_List.Dimension(NumIP());
    for (int i = 0; i < NumIP(); i++)
      fStrain_List[i].Dimension(NumSD());
  }
  
  /* allocate "last" strain list */
  if (need_strain_last) {
    fStrain_last_List.Dimension(NumIP());
    for (int i = 0; i < NumIP(); i++)
      fStrain_last_List[i].Dimension(NumSD());
  }

  /*PARAMETERS FOR ENHANCED STRAIN*/
  fH_delta_0 = list.GetParameter("Post-Localization_softening_parameter_H_Delta"); 
  fNoBandDilation = list.GetParameter("Disallow_Dilation_on_Band");
  fLocalizedFrictionCoeff = list.GetParameter("Localized_Friction_Coefficient");
  fMultiBand = list.GetParameter("Allow_Multiple_Bands");
  fBVPType = list.GetParameter("BVP_type");
  fFirstElementToLocalize = list.GetParameter("First_element_to_localize");
  fFirstElementToLocalize -= 1;	

  if (fBVPType == kPreFailed)
    {
      fPreFailedNormal.Dimension(2);
      fPreFailedSlipDir.Dimension(2);

      fPreFailedNormal[0] = list.GetParameter("First_component_of_normal_if_prefailed");
      fPreFailedNormal[1] = list.GetParameter("Second_component_of_normal_if_prefailed");
      fPreFailedSlipDir[0] = list.GetParameter("First_component_of_slip_dir_if_prefailed");
      fPreFailedSlipDir[1] = list.GetParameter("Second_component_of_slip_dir_if_prefailed");
	  
	  fPreFailedNormal.UnitVector();
	  fPreFailedSlipDir.UnitVector();
	  
	  //PreFailElements();
    }
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct a new material support and return a pointer */
MaterialSupportT* SSEnhLocCraigT::NewMaterialSupport(MaterialSupportT* p) const
{
  return SmallStrainT::NewMaterialSupport(p);
}


/* return a pointer to a new material list */
MaterialListT* SSEnhLocCraigT::NewMaterialList(const StringT& name, int size)
{
  return SmallStrainT::NewMaterialList(name, size);
}


/* calculate the internal force contribution ("-k*d") */
void SSEnhLocCraigT::FormKd(double constK)
{

  if(!IsElementTraced())
    SmallStrainT::FormKd(constK);
  else 
      {
	//cout << "SSEnhLocCraigT::FormKd\n";
	  
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/* collect incremental heat */
	bool need_heat = fElementHeat.Length() == fShapes->NumIP();
	
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

	    
	    dSymMatrixT gradActiveTensorFlowDir =
	    FormGradActiveTensorFlowDir(NumSD(), CurrIP());
	    gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);

	    strainIncr.AddScaled(-1.0*fBand->JumpIncrement(), gradActiveTensorFlowDir);
	    dSymMatrixT stressIncr(NumSD());
	    stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);
	    stressIncr += fBand->Stress_List(CurrIP());
	    fB.MultTx(stressIncr, fNEEvec);
	    
	    /* accumulate */
	    fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);
	    
	    /* incremental heat generation */
	    if (need_heat) 
	      fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	  }     
      }
}

/* form the element stiffness matrix */
void SSEnhLocCraigT::FormStiffness(double constK)
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

	dSymMatrixT gradActiveTensorFlowDir(ndof), dGdSigma(ndof);
	dMatrixT fLHSWork(nedof),fDfB((fCurrMaterial->ce_ijkl().Rows()),nedof);

	dGdSigma = FormdGdSigma(ndof);

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
                //fLHS.MultATB(fB, fD, format, dMatrixT::kAccumulate);	

		//form k_d_zeta
		gradActiveTensorFlowDir =
		FormGradActiveTensorFlowDir(ndof, CurrIP());

		fDfB.MultTx(gradActiveTensorFlowDir, k_d_zeta_work);
		k_d_zeta += k_d_zeta_work;

		//form k_zeta_d
		//fDfB.MultTx(dGdSigma, k_zeta_d_work, dMatrixT::kOverwrite);
		fDfB.MultTx(dGdSigma, k_zeta_d_work);
		k_zeta_d += k_zeta_d_work;

		//form k_zeta_zeta
		k_zeta_zeta +=
		fD.MultmBn(dGdSigma,gradActiveTensorFlowDir);
	}
	k_zeta_d *= 1.0/area;

	k_zeta_zeta *= 1.0/area;
	k_zeta_zeta += fBand->EffectiveSoftening();
	//k_zeta_zeta *= -1.0;

	fLHS.Outer(k_d_zeta, k_zeta_d, -1.0/k_zeta_zeta, dMatrixT::kAccumulate);
    }
   //cout << "Element Number " << CurrElementNumber() << ". InitialCoodinates() =\n" << InitialCoordinates() << endl;

}

/* compute the measures of strain/deformation over the element */
void SSEnhLocCraigT::SetGlobalShape(void)
{
  SmallStrainT::SetGlobalShape();
  
  /* subtract band deformation */
  if (IsElementTraced())
    {
      int ndof = NumDOF();
      dSymMatrixT gradActiveTensorFlowDir(ndof);
      
      if (fStrainDispOpt == kMeanDilBbar)
	cout << "Warning - B-bar not implemented for localized element\n";
      else
	{
	  int material_number = CurrentElement().MaterialNumber();
	  const ArrayT<bool>& needs = fMaterialNeeds[material_number];
	  
	  double jumpIncrement = CalculateJumpIncrement();

#if 0
	  /* loop over integration points again */
	  for (int i = 0; i < NumIP(); i++)
	    {
	      gradActiveTensorFlowDir =
		FormGradActiveTensorFlowDir(ndof, i);
	      
	      /*change shear strains back to matrix values */
	      /* vector values are used when created */
	      gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);
				
	      /* deformation gradient */
	      if (needs[fNeedsOffset + kstrain])
		{			  
		  fStrain_List[i].AddScaled(-(fBand->Jump()+
		     jumpIncrement),gradActiveTensorFlowDir);
		}
		
	      /* "last" deformation gradient */ //is this right?
	      if (needs[fNeedsOffset + kstrain_last])
		{
		  fStrain_last_List[i].AddScaled(-(fBand->Jump()),
						 gradActiveTensorFlowDir);
		}
	    }
#endif
	}
    }
}

/***********************************************************************
 * Protected
 ***********************************************************************/



double SSEnhLocCraigT::CalculateJumpIncrement()
{
  if (!IsBandActive())
    {
      fBand -> StoreJumpIncrement(0.0);
      return 0.0;
    }
  
  int ndof = NumDOF();
  dSymMatrixT dGdSigma = FormdGdSigma(ndof);
  dSymMatrixT gradActiveTensorFlowDir(ndof);
    
  //fD.SetToScaled(1.0, HookeanMatT::Modulus());
  fD.SetToScaled(1.0, fCurrMaterial->ce_ijkl());
  dArrayT dGfD(fD.Rows());
  fD.MultTx(dGdSigma, dGfD);
	    
  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();	    
  double area = 0.0;
  double jumpIncrement = 0.0;
  double jumpWork = 0.0;

  /* loop over integration points */
  for (int i = 0; i < NumIP(); i++)
    {
      double scale = (*Det++)*(*Weight++);
      area += scale;
      dSymMatrixT strainIncr = fStrain_List [i];
      strainIncr -= fStrain_last_List [i];
      strainIncr.ScaleOffDiagonal(2.0);

      gradActiveTensorFlowDir = FormGradActiveTensorFlowDir(ndof, i);
      
      jumpIncrement += scale * strainIncr.Dot(dGfD, strainIncr);
      jumpWork += scale * gradActiveTensorFlowDir.Dot(dGfD,gradActiveTensorFlowDir);		  
    }

  jumpIncrement /= (jumpWork + area * fBand->H_delta());

  double trialDeltaResidCohesion = -1.0*fabs(jumpIncrement)*fBand->H_delta();
  /* check to see that residual cohsion does not drop below 0, adjust if nec */
  if (fBand->ResidualCohesion() < trialDeltaResidCohesion)
    {
      //full step with no softening
      jumpIncrement *= (jumpWork + area * fBand -> H_delta())/jumpWork;

      // *fraction of step over which no softening occurs
      jumpIncrement *= (trialDeltaResidCohesion -
			fBand->ResidualCohesion())/trialDeltaResidCohesion;
      double jumpIncrementSign = jumpIncrement/fabs(jumpIncrement); 

      //plus amount to get cohesion to zero
      jumpIncrement -= jumpIncrementSign * (fBand->ResidualCohesion())/(fBand-> H_delta());
      
      fBand->SetEffectiveSoftening(0.0);
    }
  else
    {
      fBand->SetEffectiveSoftening(fBand->H_delta()); 
    } 
  fBand -> StoreJumpIncrement(jumpIncrement);

  return jumpIncrement;
}

bool SSEnhLocCraigT::IsBandActive()
{
  /*calculate average stress assuming no jump*/
  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();
	
  double area = 0.0;
  double normalStress = 0.0;
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

      dSymMatrixT stressIncr(NumSD());
      stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);
      stressIncr += fBand->Stress_List(CurrIP());
      
      area += (*Det)*(*Weight);
      double scale = (*Det++)*(*Weight++);

      normalStress += scale * stressIncr.MultmBn(fBand->Normal(),fBand-> Normal());
      shearStress += scale * stressIncr.MultmBn(fBand->PerpSlipDir(),fBand-> Normal()); 
    }

  if (shearStress < 0.0)
    {
      /* align slip direction with shear stress direction to get 
	 correct yield surface */
      fBand-> FlipSlipDir();
      shearStress *= -1.0;
      cout << "Slip direction flipped\n";
    }

  normalStress/= area;
  shearStress = shearStress/area;

  double neededCohesion = shearStress + normalStress * fLocalizedFrictionCoeff;
 
  if (fBand-> ResidualCohesion() < neededCohesion)
    {
      fBand-> SetActive(true);
      //cout << "Band is active.\n";
      return true;
    }
  else
    {
      fBand -> SetActive(false);
      cout << "Band is traced but not active this step.\n";
      return false;
    }
}


void SSEnhLocCraigT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
               const iArrayT& e_codes, dArray2DT& e_values)
{

  cout << "ComputeOutput \n\n";


/* number of output values */
    int n_out = n_codes.Sum();
    int e_out = e_codes.Sum();
        
    int n_simo, n_extrap;
    if (qUseSimo) {
        n_simo = n_out - n_codes[iNodalDisp] - n_codes[iNodalCoord];
        n_extrap = n_codes[iNodalDisp] + n_codes[iNodalCoord];
    } else {
        n_simo = 0;
        n_extrap = n_out;
    }
            
    /* nothing to output */
    if (n_out == 0 && e_out == 0) return;

    /* dimensions */
    int nsd = NumSD();
    int ndof = NumDOF();
    int nen = NumElementNodes();
    int nnd = ElementSupport().NumNodes();
    int nstrs = fB.Rows();

    /* reset averaging workspace */
    ElementSupport().ResetAverage(n_extrap);

    /* allocate element results space */
    e_values.Dimension(NumElements(), e_out);

    /* nodal work arrays */
    dArray2DT nodal_space(nen, n_extrap);
    dArray2DT nodal_all(nen, n_extrap);
    dArray2DT coords, disp;
    dArray2DT nodalstress, princstress, matdat;
    dArray2DT energy, speed;
    dArray2DT Poynting;

    /* ip values */
    dSymMatrixT cauchy((nstrs != 4) ? nsd : dSymMatrixT::k3D_plane), nstr_tmp;
    dArrayT ipmat(n_codes[iMaterialData]), ipenergy(1);
    dArrayT ipspeed(nsd), ipprincipal(nsd), ipPoynting(nsd);

    /* set shallow copies */
    double* pall = nodal_space.Pointer();
    coords.Alias(nen, n_codes[iNodalCoord], pall); pall += coords.Length();
    disp.Alias(nen, n_codes[iNodalDisp], pall)   ; pall += disp.Length();

    /* work space for Poynting vector */
    dMatrixT stress_mat, F_inv, PbyJ;
    if (n_codes[iPoyntingVector]) {
        stress_mat.Dimension(nsd);
        F_inv.Dimension(nsd);
        PbyJ.Dimension(nsd);
    }

    /* workspaces for Simo */
    int simo_offset = coords.MinorDim() + disp.MinorDim();
    dArray2DT simo_space(nen,qUseSimo ? n_simo : 0);
    dArray2DT simo_all(nen,qUseSimo ? n_simo : 0);
    dArray2DT simoNa_bar(nen,qUseSimo ? 1 : 0);
    dArray2DT simo_force;
    dArray2DT simo_mass;
    iArrayT simo_counts;
        
    if (!qUseSimo) 
    {
        nodalstress.Alias(nen, n_codes[iNodalStress], pall); pall += nodalstress.Length();
        princstress.Alias(nen, n_codes[iPrincipal], pall)  ; pall += princstress.Length();
        energy.Alias(nen, n_codes[iEnergyDensity], pall)   ; pall += energy.Length();
        speed.Alias(nen, n_codes[iWaveSpeeds], pall)       ; pall += speed.Length();
        matdat.Alias(nen, n_codes[iMaterialData], pall)    ; pall += matdat.Length();
        Poynting.Alias(nen, n_codes[iPoyntingVector], pall);
    }
    else
    {
        simo_force.Dimension(ElementSupport().NumNodes(),qUseSimo ? n_simo : 0);
        simo_mass.Dimension(ElementSupport().NumNodes(),qUseSimo ? 1 : 0);
        simo_counts.Dimension(ElementSupport().NumNodes());

        pall = simo_space.Pointer();
        nodalstress.Alias(nen, n_codes[iNodalStress], pall); pall += nodalstress.Length();
        princstress.Alias(nen, n_codes[iPrincipal], pall)  ; pall += princstress.Length();
        energy.Alias(nen, n_codes[iEnergyDensity], pall)   ; pall += energy.Length();
        speed.Alias(nen, n_codes[iWaveSpeeds], pall)       ; pall += speed.Length();
        matdat.Alias(nen, n_codes[iMaterialData], pall)    ; pall += matdat.Length();
        Poynting.Alias(nen, n_codes[iPoyntingVector], pall);

        simo_mass = 0.;
        simo_force = 0.;
        simo_counts = 0;
    }
        
    /* element work arrays */
    dArrayT element_values(e_values.MinorDim());
    pall = element_values.Pointer();
    dArrayT centroid, ip_centroid, ip_mass;
    dArrayT ip_coords(nsd);
    if (e_codes[iCentroid])
    {
        centroid.Alias(nsd, pall); pall += nsd;
        ip_centroid.Dimension(nsd);
    }
    if (e_codes[iMass]) {
        ip_mass.Alias(NumIP(), pall); 
        pall += NumIP();
    }
    double w_tmp, ke_tmp;
    double mass;
    double& strain_energy = (e_codes[iStrainEnergy]) ? *pall++ : w_tmp;
    double& kinetic_energy = (e_codes[iKineticEnergy]) ? *pall++ : ke_tmp;
    dArrayT linear_momentum, ip_velocity;
    if (e_codes[iLinearMomentum])
    {
        linear_momentum.Alias(ndof, pall); pall += ndof;
        ip_velocity.Dimension(ndof);
    }
    else if (e_codes[iKineticEnergy]) ip_velocity.Dimension(ndof);

    dArray2DT ip_stress;
    if (e_codes[iIPStress])
    {
        ip_stress.Alias(NumIP(), e_codes[iIPStress]/NumIP(), pall);
        pall += ip_stress.Length();
    }
    dArray2DT ip_material_data;
    if (e_codes[iIPMaterialData])
    {
        ip_material_data.Alias(NumIP(), e_codes[iIPMaterialData]/NumIP(), pall);
        pall += ip_material_data.Length();
        ipmat.Dimension(ip_material_data.MinorDim());
    }

    /* check that degrees are displacements */
    int interpolant_DOF = InterpolantDOFs();

    bool is_axi = Axisymmetric();
    double Pi2 = 2.0*acos(-1.0);

    Top();
    while (NextElement())
        if (CurrentElement().Flag() != ElementCardT::kOFF)
        {
            /* initialize */
            nodal_space = 0.0;
            simo_space = 0.;
            simo_all = 0.;
            simoNa_bar = 0.;

            /* global shape function values */
            SetGlobalShape();
			
			/**********************************************************************/
	     
            /* output jump for element */
			jump_out.open_append("jump.info");
			if(IsElementTraced())
				jump_out << setw(16) << fBand->JumpIncrement() + fBand->Jump() << " ";
			else
			    jump_out << setw(16) << "0.0 ";	
	    


            /* collect nodal values */
            if (e_codes[iKineticEnergy] || e_codes[iLinearMomentum] || n_codes[iPoyntingVector]) {
                if (fLocVel.IsRegistered())
                    SetLocalU(fLocVel);
                else
                    fLocVel = 0.0;
            }

            /* coordinates and displacements all at once */
            if (n_codes[iNodalCoord]) fLocInitCoords.ReturnTranspose(coords);
            if (n_codes[ iNodalDisp]) {
                if (interpolant_DOF)
                    fLocDisp.ReturnTranspose(disp);
                else
                    NodalDOFs(CurrentElement().NodesX(), disp);
            }
 
            /* initialize element values */
            mass = strain_energy = kinetic_energy = 0;
            if (e_codes[iCentroid]) centroid = 0.0;
            if (e_codes[iLinearMomentum]) linear_momentum = 0.0;
            const double* j = fShapes->IPDets();
            const double* w = fShapes->IPWeights();

            /* integrate */
            dArray2DT Na_X_ip_w;
            fShapes->TopIP();
            while (fShapes->NextIP())
            {
                /* density may change with integration point */
                double density = fCurrMaterial->Density();

                /* element integration weight */
                double ip_w = (*j++)*(*w++);
                if (is_axi) {
                    fShapes->IPCoords(ip_coords);
                    ip_w *= Pi2*ip_coords[0]; /* radius is the x1 coordinate */
                }

                if (qUseSimo || qNoExtrap)
                {
                    Na_X_ip_w.Dimension(nen,1);
                    if (qUseSimo)
                    {
                        const double* Na_X = fShapes->IPShapeX();
                        Na_X_ip_w = ip_w;
                        for (int k = 0; k < nen; k++)
                            Na_X_ip_w(k,0) *= *Na_X++;
                        simoNa_bar += Na_X_ip_w;
                    }
                    else
                        for (int k = 0; k < nen; k++)
                            Na_X_ip_w(k,0) = 1.;
                }

                /* get Cauchy stress */
                //const dSymMatrixT stress(NumSD());                
                const dSymMatrixT& stress = fCurrMaterial->s_ij();
                //dSymMatrixT stress(NumSD());
                
        // #if 1
                if (IsElementTraced())
          {
            dSymMatrixT strainIncr = fStrain_List [CurrIP()];
            strainIncr -= fStrain_last_List [CurrIP()];
            
            dSymMatrixT gradActiveTensorFlowDir =
              FormGradActiveTensorFlowDir(NumSD(), CurrIP());
            gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);
            
            strainIncr.AddScaled(-1.0*fBand->JumpIncrement(), gradActiveTensorFlowDir);
            dSymMatrixT stressIncr(NumSD());
            stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);
            stressIncr += fBand->Stress_List(CurrIP());
            //stress = stressIncr;
            //cout << "stressIncr = " << stressIncr << endl;
            cauchy.Translate(stressIncr);
          }
                else
          {
            //const dSymMatrixT& stress = fCurrMaterial->s_ij();
            //cout << "stress = " << stress << endl;
            cauchy.Translate(stress);
            //cout << "cauchy = \n" << cauchy << endl;
          }    
                 
        //#endif
        //cout << "hi\n";    
        //cauchy.Translate(stress);
        //cout << "cauchy = \n" << cauchy << endl;


                /* stresses */
                if (n_codes[iNodalStress]) {        
                    if (qNoExtrap)
                        for (int k = 0; k < nen; k++)
                            nodalstress.AddToRowScaled(k,Na_X_ip_w(k,0),cauchy);
                    else
                        fShapes->Extrapolate(cauchy, nodalstress);
                }

                if (e_codes[iIPStress]) {
                    double* row = ip_stress(fShapes->CurrIP());
                    nstr_tmp.Set(nsd, row);
                    nstr_tmp = cauchy;
                    row += cauchy.Length();
                    nstr_tmp.Set(nsd, row);
                    fCurrMaterial->Strain(nstr_tmp);
                }

                /* wave speeds */
                if (n_codes[iWaveSpeeds])
                {
                    /* acoustic wave speeds */
                    fCurrMaterial->WaveSpeeds(fNormal, ipspeed);
                    if (qNoExtrap)
                        for (int k = 0; k < nen; k++)
                            speed.AddToRowScaled(k,Na_X_ip_w(k,0),ipspeed);
                        else
                            fShapes->Extrapolate(ipspeed, speed);
                }

                /* principal values - compute principal before smoothing */
                if (n_codes[iPrincipal])
                {
                    /* compute eigenvalues */
                    cauchy.PrincipalValues(ipprincipal);
                    if (qNoExtrap)
                        for (int k = 0; k < nen; k++)
                            princstress.AddToRowScaled(k,Na_X_ip_w(k,0),ipprincipal);
                    else
                        fShapes->Extrapolate(ipprincipal, princstress);        
                }

                /* strain energy density */
                if (n_codes[iEnergyDensity] || e_codes[iStrainEnergy] || n_codes[iPoyntingVector])
                {
                    double ip_strain_energy = fCurrMaterial->StrainEnergyDensity();

                    /* nodal average */
                    if (n_codes[iEnergyDensity])
                    {
                        ipenergy[0] = ip_strain_energy;
                        if (qNoExtrap)
                            for (int k = 0; k < nen; k++)
                                energy.AddToRowScaled(k,Na_X_ip_w(k,0),ipenergy);
                        else
                            fShapes->Extrapolate(ipenergy,energy);
                    }

                    /* integrate over element */
                    if (e_codes[iStrainEnergy])
                        strain_energy += ip_w*ip_strain_energy;
                        
                    /* Poynting vector */
                    if (n_codes[iPoyntingVector]) {
                        ipPoynting = 0.0;
                        ipPoynting[0] += fv_ss*ip_strain_energy;
                    }
                }

                /* material stuff */
                if (n_codes[iMaterialData] || e_codes[iIPMaterialData])
                {
                    /* compute material output */
                    fCurrMaterial->ComputeOutput(ipmat);

                    /* store nodal data */
                    if (n_codes[iMaterialData])
                    {
                        if (qNoExtrap)
                            for (int k = 0; k < nen; k++)
                                matdat.AddToRowScaled(k,Na_X_ip_w(k,0),ipmat);
                        else 
                            fShapes->Extrapolate(ipmat, matdat);
                    }

                    /* store element data */
                    if (e_codes[iIPMaterialData]) ip_material_data.SetRow(fShapes->CurrIP(), ipmat);
                }

                /* mass averaged centroid */
                if (e_codes[iCentroid] || e_codes[iMass])
                {
                    /* mass */
                    mass += ip_w*density;

                    /* integration point mass */
                    if (e_codes[iMass]) ip_mass[fShapes->CurrIP()] = ip_w*density;

                    /* moment */
                    if (e_codes[iCentroid]) {
                        fShapes->IPCoords(ip_centroid);
                        centroid.AddScaled(ip_w*density, ip_centroid);
                    }
                }

                /* kinetic energy/linear momentum */
                if (e_codes[iKineticEnergy] || e_codes[iLinearMomentum] || n_codes[iPoyntingVector])
                {
                    /* velocity at integration point */
                    fShapes->InterpolateU(fLocVel, ip_velocity);                                        
                    double ke_density = 0.5*density*dArrayT::Dot(ip_velocity, ip_velocity);

                    /* kinetic energy */
                    if (e_codes[iKineticEnergy])
                        kinetic_energy += ip_w*ke_density;
        
                    /* linear momentum */
                    if (e_codes[iLinearMomentum])
                        linear_momentum.AddScaled(ip_w*density, ip_velocity);
                
                    /* Poynting vector */
                    if (n_codes[iPoyntingVector]) {
                        ipPoynting[0] += fv_ss*ke_density;

                        /* is finite strain */
                        
                        #if 0
                        FSSolidMatT* fs_mat = TB_DYNAMIC_CAST(FSSolidMatT*, fCurrMaterial);
                        if (fs_mat)
                        {
                            /* compute PK1/J */
                            cauchy.ToMatrix(stress_mat);
                            F_inv.Inverse(fs_mat->F());
                            PbyJ.MultABT(stress_mat, F_inv);
                        
                            PbyJ.Multx(ip_velocity, ipPoynting, 1.0/F_inv.Det(), dMatrixT::kAccumulate);
                        } 
                        else { /* small strain */
                        
                        #endif
                        
                            cauchy.Multx(ip_velocity, ipPoynting, 1.0, dMatrixT::kAccumulate);
                        
                        //}
                        
                        /* extrapolate */
                        ipPoynting *= -1; /* positive toward crack tip */
                        fShapes->Extrapolate(ipPoynting, Poynting);
                    }
                }
            }

            /* copy in the cols */
            int colcount = 0;
            nodal_all.BlockColumnCopyAt(disp       , colcount); colcount += disp.MinorDim();
            nodal_all.BlockColumnCopyAt(coords     , colcount); colcount += coords.MinorDim();

            if (!qUseSimo)
            {
                if (qNoExtrap)
                {
                    double nip(fShapes->NumIP());
                    nodalstress /= nip;
                    princstress /= nip;
                    energy /= nip;
                    speed /= nip;
                    matdat /= nip;
                }
                nodal_all.BlockColumnCopyAt(nodalstress, colcount); colcount += nodalstress.MinorDim();
                nodal_all.BlockColumnCopyAt(princstress, colcount); colcount += princstress.MinorDim();
                nodal_all.BlockColumnCopyAt(energy     , colcount); colcount += energy.MinorDim();
                nodal_all.BlockColumnCopyAt(speed      , colcount); colcount += speed.MinorDim();
                nodal_all.BlockColumnCopyAt(matdat     , colcount); colcount += matdat.MinorDim();
                nodal_all.BlockColumnCopyAt(Poynting   , colcount); colcount += Poynting.MinorDim();
            }
            else
            {        
                colcount = 0;
                simo_all.BlockColumnCopyAt(nodalstress, colcount); colcount += nodalstress.MinorDim();
                simo_all.BlockColumnCopyAt(princstress, colcount); colcount += princstress.MinorDim();
                simo_all.BlockColumnCopyAt(energy     , colcount); colcount += energy.MinorDim();
                simo_all.BlockColumnCopyAt(speed      , colcount); colcount += speed.MinorDim();
                simo_all.BlockColumnCopyAt(matdat     , colcount); colcount += matdat.MinorDim();
                simo_all.BlockColumnCopyAt(Poynting   , colcount); colcount += Poynting.MinorDim();

                iArrayT currIndices = CurrentElement().NodesX();
                simo_force.Accumulate(currIndices,simo_all);
                simo_mass.Accumulate(currIndices,simoNa_bar);
                for (int i = 0; i < currIndices.Length(); i++)
                    simo_counts[currIndices[i]]++;
            }

            /* accumulate - extrapolation done from ip's to corners => X nodes */
            ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_all);

            /* element values */
            if (e_codes[iCentroid]) centroid /= mass;

            /* store results */
            e_values.SetRow(CurrElementNumber(), element_values);
        }

    /* get nodally averaged values */
    const OutputSetT& output_set = ElementSupport().OutputSet(fOutputID);
    const iArrayT& nodes_used = output_set.NodesUsed();
    dArray2DT extrap_values(nodes_used.Length(), n_extrap);
    extrap_values.RowCollect(nodes_used, ElementSupport().OutputAverage());

    int tmpDim = extrap_values.MajorDim();
    n_values.Dimension(tmpDim,n_out);
    n_values.BlockColumnCopyAt(extrap_values,0);
    if (qUseSimo)
    {        
        int rowNum = 0;
//      iArrayT nodes_used(tmpDim);
        dArray2DT tmp_simo(tmpDim, n_simo);
        for (int i = 0; i < simo_force.MajorDim(); i++)
            if (simo_counts[i] > 0) {
//              nodes_used[rowNum] = i;
                simo_force.ScaleRow(i, 1./simo_mass(i,0));
                tmp_simo.SetRow(rowNum, simo_force(i));
                rowNum++;
            }

        /* collect final values */
        n_values.BlockColumnCopyAt(tmp_simo, simo_offset);

        /* write final values back into the averaging workspace */
        if (extrap_values.MinorDim() != n_values.MinorDim()) {
            ElementSupport().ResetAverage(n_values.MinorDim());
            ElementSupport().AssembleAverage(nodes_used, n_values);
        }
    }

  /* inherited */
  //SmallStrainT::ComputeOutput(n_codes, n_values, e_codes, e_values);
  
  jump_out << endl;  

}


void SSEnhLocCraigT::CloseStep(void)
{

cout << "\n CloseStep \n \n";

/*update traced elements */ 
if (fLocalizationHasBegun)
	{
	  Top();
      while (NextElement())
	{
	  if (IsElementTraced())
	    {
		SetGlobalShape();	

	    	/* loop over integration points */
	    	fShapes->TopIP();
	      while (fShapes->NextIP())
			{    
	      	dSymMatrixT strainIncr = fStrain_List [CurrIP()];
	      	strainIncr -= fStrain_last_List [CurrIP()]; 
	      
	     	dSymMatrixT gradActiveTensorFlowDir =
			FormGradActiveTensorFlowDir(NumSD(), CurrIP());
	      	gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);      	

	      	strainIncr.AddScaled(-1.0*fBand->JumpIncrement(), gradActiveTensorFlowDir);
	      
	     	dSymMatrixT stressIncr(NumSD());
	      	stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);
	      	fBand -> IncrementStress(stressIncr, CurrIP());
	      	}
	      fBand -> CloseStep();
	      	
			cout << "residual cohesion = " << fBand->ResidualCohesion();
	      	cout << "fBand->JumpIncrement = " << fBand->JumpIncrement();
	      	cout << ", fBand->Jump() = " << fBand->Jump() << endl;	
	    }			
	}
	}

bool finishedTracing = false;
	
/* check for newly localized elements on existing bands */
//fJustChecking = false;
fEdgeOfBandElements.Top();
fEdgeOfBandCoords.Top();

while(!finishedTracing)
{
while(fEdgeOfBandElements.Next())
	{
	  fEdgeOfBandCoords.Next();
	  cout << "fEdgeOfBandElements.Current() = " << flush << fEdgeOfBandElements.Current() << endl << flush;
	  cout << "fEdgeOfBandCoords.Current() = " << fEdgeOfBandCoords.Current() << endl;	  
	  GetElement(fEdgeOfBandElements.Current());
	  if (IsElementLocalized())
		TraceElement();
	}

//reset to check for elements not on an existing band
fDetAMin = 1.0e99;
finishedTracing = true;
//fEdgeOfBandElements.Previous();
fEdgeOfBandElements.Current(fEdgeOfBandElements.Position() -1);

/* check for newly localized elements not on existing bands */
//fJustChecking = true;
	  cout << "hi hi\n";
    
if (!fLocalizationHasBegun || fMultiBand) 
    {
	  cout << "hi\n";
      //choose first element then let band progress
 
      Top();
      while (NextElement())
	{
		  //cout << "1\n";
		if (!IsElementTraced())
		{
		    cout << "CurrElementNumber = " << CurrElementNumber() << "\n";
			GetElement(CurrElementNumber());

			if (IsElementLocalized())
				finishedTracing = false;
		}
	}
			if (!finishedTracing)
			{
				fLocalizationHasBegun = true;
				//fEdgeOfBandElements.Current(0); //should be at
	  
	  			cout << "fLeastDetEle = " << fLeastDetEle << endl;
				GetElement(fLeastDetEle);
				//fEdgeOfBandCoords.Free();
				fEdgeOfBandElements.Append(fLeastDetEle);
				fEdgeOfBandCoords.Append(Centroid());	  
				//fEdgeOfBandElements.Top();
				//fEdgeOfBandCoords.Top();
			}
    }
}

 //cout << "hi\n" << flush; 

 if (fBVPType == kPreFailed && ElementSupport().Time() == 0.0)
   {
     //cout << "time = 0.0\n" << flush; 
     PreFailElements();
   }

  SmallStrainT::CloseStep();
 












#if 0
 if (fLocalizationHasBegun)
    {
      /*update traced elements */ 
      Top();
      while (NextElement())
	{
	  if (IsElementTraced())
	    {
	      //cout << "fBand->JumpIncrement = " << fBand->JumpIncrement() << endl;

	      //fBand -> CloseStep();

		SetGlobalShape();	


	    	/* loop over integration points */
	    	fShapes->TopIP();
	      while (fShapes->NextIP())
			{    
	      	dSymMatrixT strainIncr = fStrain_List [CurrIP()];
	      	strainIncr -= fStrain_last_List [CurrIP()]; 
	      
	     	dSymMatrixT gradActiveTensorFlowDir =
			FormGradActiveTensorFlowDir(NumSD(), CurrIP());
	      	gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);      	

	      	strainIncr.AddScaled(-1.0*fBand->JumpIncrement(), gradActiveTensorFlowDir);
	      
	     	dSymMatrixT stressIncr(NumSD());
	      	stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);
	      	fBand -> IncrementStress(stressIncr, CurrIP());
		if (CurrIP() == 0)
		  {
		    //cout << "strain = \n" << fStrain_List[CurrIP()] << endl;
		    //cout << "strain_last = \n" << fStrain_last_List[CurrIP()] << endl;			  
		    //cout << "strainIncr = \n" << strainIncr << endl; 
		    //cout << "stressIncr = \n" << stressIncr << endl; 
		  }
	      	}
	      	


	      fBand -> CloseStep();
	      	
	      	cout << "residual cohesion = " << fBand->ResidualCohesion();
	      	cout << "fBand->JumpIncrement = " << fBand->JumpIncrement();
	      	cout << ", fBand->Jump() = " << fBand->Jump() << endl;	
	    }
		
	
	}
      /* check for newly localized elements */
      fEdgeOfBandElements.Top();
      fEdgeOfBandCoords.Top();
      while(fEdgeOfBandElements.Next())
	{
	  //    cout << "checking for localization of element " <<
	  //    	fEdgeOfBandElements.Current() << endl;
	  fEdgeOfBandCoords.Next();
	  GetElement(fEdgeOfBandElements.Current());
	  IsElementLocalized();
	}
	  // if fMultiBand..
    }
  else
    {
      //choose first element then let band progress
      bool localizationHasBegun = false;
      Top();
      //cout << "1 " << flush;
      while (NextElement())
	{
	  GetElement(CurrElementNumber());
	  
	  //cout << "2 " << flush;

	  if (IsElementLocalized())
	    localizationHasBegun = true;
	}
      if (localizationHasBegun)
	{
          //cout << "3 " << flush;

	  fLocalizationHasBegun = true;
	  
	  fEdgeOfBandElements.Current(0);


	  cout << fEdgeOfBandElements.Current() << endl << flush;
	  
	  GetElement(fEdgeOfBandElements.Current());
	  fEdgeOfBandCoords.Free();
	  fEdgeOfBandCoords.Append(Centroid());	  
	  fEdgeOfBandElements.Top();
	  fEdgeOfBandCoords.Top();
	  	  //localize 1st element?
	  while(fEdgeOfBandElements.Next())
	    {
	      cout << "checking for localization of element " <<
	      	fEdgeOfBandElements.Current() << endl; 

		  fEdgeOfBandCoords.Next();
	      GetElement(fEdgeOfBandElements.Current());    
	      IsElementLocalized();
	    }
	}
    }

 //cout << "hi\n" << flush; 

 if (fBVPType == kPreFailed && ElementSupport().Time() == 0.0)
   {
     //cout << "time = 0.0\n" << flush; 
     PreFailElements();
   }

  SmallStrainT::CloseStep();
 
  /* 	
  Top();
  while (NextElement())
    {	
      SetGlobalShape();
      
      jump_out.open_append("jump.info");
      if(IsElementTraced())
	jump_out << setw(16) << fBand->JumpIncrement() + fBand->Jump() << " ";
      else
	jump_out << setw(16) << "0.0 ";
    }
  */
 
 #endif
 
}


/***********************************************************************
 * Protected
 ***********************************************************************/

/* current element operations */
void SSEnhLocCraigT::GetElement(int elementNumber)
{
  /* inherited */
  //bool result = ContinuumElementT::NextElement();
  fElementCards.Current(elementNumber);  

  /* get material pointer */
  ContinuumMaterialT* pcont_mat = (*fMaterialList)[CurrentElement().MaterialNumber()];
      
  /* cast is safe since class contructs materials list */
  fCurrMaterial = (SolidMaterialT*) pcont_mat;
  
  //SetLocalX(fLocInitCoords);
  SetGlobalShape();
 }




//move to surface mat model?
dSymMatrixT SSEnhLocCraigT::FormdGdSigma(int ndof)
{
  dSymMatrixT dGdSigma(ndof), work(ndof);
  dMatrixT dGNonSym(ndof);

  dGNonSym.Outer(fBand->Normal(), fBand->PerpSlipDir());
  dGdSigma.Symmetrize(dGNonSym);
  work.Outer(fBand->Normal());
  dGdSigma.AddScaled(fLocalizedFrictionCoeff, work);

  /* use vector values for mult */
  dGdSigma.ScaleOffDiagonal(2.0);

  return dGdSigma;
}

dSymMatrixT SSEnhLocCraigT::FormGradActiveTensorFlowDir(int ndof, int ip)
{
  dSymMatrixT G(ndof);
  dMatrixT G_NonSym(ndof);
  iAutoArrayT activeNodes = fBand->ActiveNodes();
  int A;
  dArrayT grad_f(ndof);
  grad_f = 0.0;

  fShapes->SetIP(ip);
  Set_B(fShapes->Derivatives_U(), fB);

  for (int i=0; i<ndof; i++)
    {
      activeNodes.Top();

      while(activeNodes.Next())      
	{
	  A = activeNodes.Current();  
	  grad_f[i] += fB(i, (A)*ndof +i);
	}
    }

  G_NonSym.Outer(grad_f, fBand->SlipDir());
  G.Symmetrize(G_NonSym);

  /* use vector values for mult */
  G.ScaleOffDiagonal(2.0);

  return G;
}


bool SSEnhLocCraigT::IsElementTraced()
{
 	 //int elementNumber = CurrElementNumber();
 	 return IsElementTraced(CurrElementNumber());
}
  
bool SSEnhLocCraigT::IsElementTraced(int elementNumber)
{
  bool isTraced = fTracedElements.HasKey(elementNumber);

  if (isTraced)
    LoadBand(elementNumber);

  return isTraced;
}

void SSEnhLocCraigT::LoadBand(int elementNumber)
{
  fBand = fTracedElements[elementNumber];
}

bool SSEnhLocCraigT::IsElementLocalized()
{
  /* newly localized elements only */
  if (IsElementTraced())
	return false;

  bool locCheck = false;
  double detA, detAMin  = 1.0e99;
  AutoArrayT <dArrayT> normals;
  AutoArrayT <dArrayT> slipDirs;
  //AutoArrayT <dArrayT> bestNormals;
  //AutoArrayT <dArrayT> bestSlipDirs;

  /* loop over integration points */
  fShapes->TopIP();
  while ( fShapes->NextIP() )
    {
      //if (CurrIP() == 0 && CurrElementNumber() == 0)
      //cout << "stress = " << fCurrMaterial->s_ij() << endl;

      //checker.SetfStructuralMatSupport(*fSSMatSupport);
      if (fCurrMaterial->IsLocalized(normals, slipDirs, detA))
	{
	  locCheck = true;
	  normals.Top();
	  while (normals.Next())

	  if (detA < detAMin)
	    detAMin = detA;

	  cout << "detAMin = " << flush << detAMin << flush << endl;
	  cout << "fDetAMin = " << fDetAMin << endl << endl;
	
      if ((fBVPType == kNonhomogeneous && detAMin < fDetAMin) || (fBVPType !=
		kNonhomogeneous && CurrElementNumber() == fFirstElementToLocalize))
      {
	  
		fDetAMin = detAMin;
		fLeastDetEle = CurrElementNumber();
		cout << "detAMin = " << detAMin << endl;
		cout << "fDetAMin = " << fDetAMin << endl;
		cout << "fLeastDetEle = " << fLeastDetEle << endl;
		//save best normals and slipdirs here?
      }
	}
	}

  return locCheck;
}

bool SSEnhLocCraigT::TraceElement()
{

  bool locCheck = false;
  double detA, detAMin  = 1.0e99;
  AutoArrayT <dArrayT> normals;
  AutoArrayT <dArrayT> slipDirs;
  AutoArrayT <dArrayT> bestNormals;
  AutoArrayT <dArrayT> bestSlipDirs;

  /* loop over integration points */
  fShapes->TopIP();
  while ( fShapes->NextIP() )
    {
      //if (CurrIP() == 0 && CurrElementNumber() == 0)
      //cout << "stress = " << fCurrMaterial->s_ij() << endl;

      //checker.SetfStructuralMatSupport(*fSSMatSupport);
      if (fCurrMaterial->IsLocalized(normals, slipDirs, detA))
	{
	  locCheck = true;
	  normals.Top();
	  while (normals.Next())

	  if (detA < detAMin)
	    {
	      detAMin = detA;
	      bestNormals.Dimension(normals);
	      normals.CopyInto(bestNormals);
	      bestSlipDirs.Dimension(slipDirs);
	      slipDirs.CopyInto(bestSlipDirs);
	    }
	}
    }

  if (locCheck)
    {
      ChooseNormals(bestNormals, bestSlipDirs);
      /* remove element from list of Edge elements*/
      fEdgeOfBandElements.DeleteAt(fEdgeOfBandElements.Position());
      fEdgeOfBandElements.Current(fEdgeOfBandElements.Position() -1);
      fEdgeOfBandCoords.DeleteAt(fEdgeOfBandCoords.Position());  
      fEdgeOfBandCoords.Current(fEdgeOfBandCoords.Position() - 1);
     }

  return locCheck;
}

void SSEnhLocCraigT::PreFailElements()
{
	fLocalizationHasBegun = true;
	
	AutoArrayT<dArrayT> normals;
	AutoArrayT<dArrayT> slipDirs;
	
	normals.Append(fPreFailedNormal);
	slipDirs.Append(fPreFailedSlipDir);

	GetElement(fFirstElementToLocalize);
	fEdgeOfBandCoords.Append(Centroid());
	fEdgeOfBandCoords.Current(0);		
	ChooseNormals(normals, slipDirs);
	fEdgeOfBandCoords.DeleteAt(0);
  
	fEdgeOfBandElements.Top();
	fEdgeOfBandCoords.Top();
	
	while(fEdgeOfBandElements.Next())
	    {
	      cout << "checking for localization of element " <<
	      	fEdgeOfBandElements.Current() << endl; 

		  fEdgeOfBandCoords.Next();
	      GetElement(fEdgeOfBandElements.Current());    
		  ChooseNormals(normals, slipDirs);	
	
      /* remove element from list of Edge elements*/
      fEdgeOfBandElements.DeleteAt(fEdgeOfBandElements.Position());
      fEdgeOfBandElements.Current(fEdgeOfBandElements.Position() -1);
      fEdgeOfBandCoords.DeleteAt(fEdgeOfBandCoords.Position());  
      fEdgeOfBandCoords.Current(fEdgeOfBandCoords.Position() - 1);	
	}
}

void SSEnhLocCraigT::ChooseNormals(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipDirs)
{
  normals.Top();
  slipDirs.Top();

  int ndof = NumDOF();

  dArrayT normal(ndof);
  dArrayT slipDir(ndof);
  dArrayT perpSlipDir(ndof);

  //determine strain on band
  //should we use elastic strain?
  dMatrixT avgGradU(ndof);
  avgGradU = 0.0;

  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();
  double area = 0.0;

  fShapes->TopIP();
  while (fShapes-> NextIP())
    {
      double scale = (*Det++)*(*Weight++);
      area += scale;
      fShapes->GradU(fLocDisp, fGradU, CurrIP());
      avgGradU.AddScaled(scale, fGradU);
      //cout << "fGradU = \n" << fGradU << endl;
    }
  avgGradU /= area;

  double prod, maxProd = -1.0;

  while(normals.Next())
    {
      cout << "current normal = \n" << normals.Current() << endl; 
      slipDirs.Next();
      prod = fabs( avgGradU.MultmBn(normals.Current(), slipDirs.Current()));
      //prod = (normals.Current() [1]) * (normals.Current() [2]);
      if (prod > maxProd)
	{
	  cout << "best normal = \n" << normals.Current() << endl; 
	  normal = normals.Current(); 
	  slipDir = slipDirs.Current();
          //maxProd = prod;
	}
    }

  //make sure slip direction is dilatant
  if (normal.Dot(normal, slipDir)<0.0)
    slipDir *= -1.0;

  normal *= -1.0;
  slipDir *= -1.0;

  perpSlipDir = slipDir;
  perpSlipDir.AddScaled(-1.0*slipDir.Dot(slipDir, normal), normal);
  perpSlipDir.UnitVector();

  if (fNoBandDilation)
    {
      slipDir = perpSlipDir;
    }


  fBand = FormNewBand(normal, slipDir, perpSlipDir, fEdgeOfBandCoords.Current(), area);

  //cout << "1 " << flush;

  fTracedElements.Insert(CurrElementNumber(), fBand);

  //cout << "2 " << flush;

  AddNewEdgeElements(CurrElementNumber());
  
  //cout << "3 " << flush;
}

BandT* SSEnhLocCraigT::FormNewBand(dArrayT normal, dArrayT slipDir,
dArrayT perpSlipDir, dArrayT coords, double area)
{
  ArrayT<dSymMatrixT> stressList;
  stressList.Dimension(NumIP());

  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();
  //double area = 0.0;

  /*for residual cohesion*/
  double normalStress = 0.0;
  double shearStress = 0.0;

  Det    = fShapes->IPDets();
  Weight = fShapes->IPWeights();

  fShapes->TopIP();
  for (int i = 0; i < NumIP(); i++)
    {
      fShapes -> NextIP();
      stressList[i].Dimension(NumSD());
      stressList[i] = fCurrMaterial -> s_ij();

      double scale = (*Det++)*(*Weight++);

      normalStress += scale * stressList[i].MultmBn(normal, normal);
      shearStress += scale * stressList[i].MultmBn(perpSlipDir, normal); 
    }

  normalStress/= area;
  shearStress = fabs(shearStress)/area;

  double residCohesion = shearStress + normalStress * fLocalizedFrictionCoeff;

  if (fBVPType == kPreFailed && ElementSupport().Time() == 0.0)
    residCohesion = 0.0;

  cout << "residCohesion = " << residCohesion << endl;


  //cout << "coords =\n" << coords << endl;
return new BandT(normal, slipDir, perpSlipDir, coords, fH_delta_0, residCohesion, stressList, this);
}

void SSEnhLocCraigT::AddNewEdgeElements(int elementNumber)
{
  ModelManagerT& model = ElementSupport().ModelManager();
  iArray2DT neighbors;
  ArrayT<StringT> ids;
  
  ElementBlockIDs(ids);  
  //cout << "hi " << flush;

  model.ElementNeighbors(ids, neighbors);

  //cout << "neighbors =\n" << neighbors << endl;

  //2D
  int numSides;
  int numSidesFound = 0;
  //int numElementSides = 0;
  iAutoArrayT activeNodes = fBand->ActiveNodes();

  switch (GeometryCode())
{
 case GeometryT::kQuadrilateral:
   {
     numSides = 4;
     break;
   }
 case GeometryT::kTriangle: 
   {
     numSides = 3;
     break;
   }
 default:
   {
     cout << "SSEnhLocCraigT::AddNewEdgeElements, geometry not implemented. \n" << flush;
     throw ExceptionT::kGeneralFail;
   }
}


  LocalArrayT nodalCoords = InitialCoordinates();
  dArrayT nodalCoord1(NumSD()), nodalCoord2(NumSD()); //coords a particular node

  for(int i = 0; i < numSides; i++)
    if ((activeNodes.HasValue((i+1) % numSides) && !activeNodes.HasValue(i))
	|| (!activeNodes.HasValue((i+1) % numSides) && activeNodes.HasValue(i)))
      {
	//cout << "i = " << i << ", neighbor = " << neighbors(elementNumber,
	//						    i) << endl;	
		if (!(neighbors(elementNumber,i) == -1 || IsElementTraced(neighbors(elementNumber ,i))))
	  	{
	    	//get coords
	    	dArrayT localizedEleCoord = fBand -> Coords();

		    for (int j = 0; j < nodalCoords.MinorDim(); j++)
		      {		
				nodalCoord1 [j] = nodalCoords(i,j);
				nodalCoord2 [j] = nodalCoords((i+1) % numSides, j);
	      		}
		    /*
	      	cout << CurrElementNumber() << endl;
	      	cout << "nodalCoords =\n" << nodalCoords;
	      	cout << "neighbors =\n" << neighbors << endl;	
	      	cout << "nodalCoord1 = \n" << nodalCoord1 << endl;	
	      	cout << "nodalCoord2 = \n" << nodalCoord2 << endl;
		    */


	    	dArrayT interceptCoords = InterceptCoords(localizedEleCoord,
						      nodalCoord1, nodalCoord2);
	    	//stick element in fEdgeOfBandElements
	    	//EdgeOfBandElement* newElement = new EdgeOfBandElement;
	    	//newElement->elementNumber = CurrElementNumber();
	    	//newElement -> coords = interceptCoords;
	    	if(fEdgeOfBandElements.AppendUnique(neighbors(elementNumber,i)))
		  {
		    cout << "Successfully appended element " <<
		      neighbors(elementNumber,i) << endl << "coords =\n" << interceptCoords << endl;
		    fEdgeOfBandCoords.Append(interceptCoords);
		  }
	  	} 
		if (++numSidesFound > 1) 
	  		break;
      }
  cout << "numSidesFound = " << numSidesFound <<endl;

}

dArrayT SSEnhLocCraigT::InterceptCoords(dArrayT& localizedEleCoord,
dArrayT& nodalCoord1, dArrayT& nodalCoord2)
{
  //assumes staright sides
  dArrayT sideVector = nodalCoord2;
  sideVector -= nodalCoord1;

  dArrayT perpSlipDir = fBand -> PerpSlipDir();

  double alpha = sideVector[0] * (localizedEleCoord[1] - nodalCoord1[1]) -
  sideVector[1] * (localizedEleCoord[0] - nodalCoord1[0]);
  alpha /= sideVector[1] * perpSlipDir[0] - sideVector[0] *
  perpSlipDir[1];

  dArrayT interceptCoord = localizedEleCoord;
  interceptCoord.AddScaled(alpha, perpSlipDir);

  return interceptCoord;
}


dArrayT SSEnhLocCraigT::Centroid()
{
  dArrayT centroid(NumDOF());

  double area = 0.0;
  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();
  dArrayT coords(NumDOF());

  centroid = 0.0;

  fShapes->TopIP();
  while ( fShapes->NextIP() )
    {
      double scale = (*Det++)*(*Weight++);
      area += scale;
      fShapes->IPCoords(coords);
      centroid.AddScaled(scale,coords);
    }
  centroid *= 1.0/area;
  
  return centroid;
 
}
