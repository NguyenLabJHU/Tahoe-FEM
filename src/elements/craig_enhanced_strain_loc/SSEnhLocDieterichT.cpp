#include "SSEnhLocDieterichT.h"

using namespace Tahoe;

/* constructor */
SSEnhLocDieterichT::SSEnhLocDieterichT(const ElementSupportT& support):
	SSEnhLocCraigT(support)
{
	SmallStrainT::SetName("small_strain_enh_loc_dieterich");
}

/* implementation of the ParameterInterfaceT interface */
void SSEnhLocDieterichT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SSEnhLocCraigT::DefineParameters(list);


	/*PARAMETERS FOR Dieterich model*/
	

	/*
	list.AddParameter(fH_delta_0, "Post-Localization_softening_parameter_H_Delta"); 
	list.AddParameter(fNoBandDilation, "Disallow_Dilation_on_Band");
	list.AddParameter(fLocalizedFrictionCoeff,
	"Localized_Friction_Coefficient");
	*/
}

void SSEnhLocDieterichT::TakeParameterList(const ParameterListT& list)
{
  /* inherited */
  SSEnhLocCraigT::TakeParameterList(list);


  /*PARAMETERS FOR Dieterich Model*/


  /*
  fH_delta_0 = list.GetParameter("Post-Localization_softening_parameter_H_Delta"); 
  fNoBandDilation = list.GetParameter("Disallow_Dilation_on_Band");
  fLocalizedFrictionCoeff =
  list.GetParameter("Localized_Friction_Coefficient");
  */
}

void SSEnhLocDieterichT::FormStiffness(double constK)
{
  /*inherited */
  SSEnhLocCraigT::FormStiffness(constK);
}

double SSEnhLocDieterichT::CalculateJumpIncrement()
{
  return SSEnhLocCraigT::CalculateJumpIncrement();
}      

bool SSEnhLocDieterichT::IsBandActive()
{
  return SSEnhLocCraigT::IsBandActive();
}

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
  /* loop over traced elements to update theta */

  SSEnhLocCraigT::CloseStep();
}
