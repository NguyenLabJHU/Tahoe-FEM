/* $Header: /home/regueiro/tahoe_cloudforge_repo_snapshots/development/src/elements/fluid_element/materials/FluidMaterialT.cpp,v 1.2 2006-07-18 01:21:20 a-kopacz Exp $ */
/* created: tdnguye (07/12/2006) */
#include "FluidMaterialT.h"
#include "FluidMatSupportT.h"

#include "StringT.h"
#include "dArrayT.h"
#include "dSymMatrixT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<FluidMaterialT>::fByteCopy = false;
DEFINE_TEMPLATE_STATIC const bool ArrayT<FluidMaterialT*>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
FluidMaterialT::FluidMaterialT(void):
	ParameterInterfaceT("linear_fluid_material"),
	fFluidMatSupport(NULL),
	fDensity(0.0),
	fMu(0.0)
{

}

/* set support */
void FluidMaterialT::SetFluidMatSupport(const FluidMatSupportT* support)
{
  WriteCallLocation("SetFluidMatSupport"); //DEBUG
	/* inherited */
	SetMaterialSupport(support);
	fFluidMatSupport = support;

	/* dimension */
	int nsd = NumSD();
	fStrainRate.Dimension(nsd);
	fStress.Dimension(nsd);
	fModulus.Dimension(dSymMatrixT::NumValues(nsd));
	/* initialize */
	fStress = 0.0;
	fModulus = 0.0;
}

/* viscosity */
const dMatrixT& FluidMaterialT::c_ijkl(void)
{
  WriteCallLocation("c_ijkl"); //DEBUG

	double third = 1.0/3.0;
	
	if (NumSD() ==1)
		fModulus = 2.0*fMu;
	else if (NumSD() ==2 ) {
		/*"plane strain assumption for now.  i.e. no flow in out of plane direction*/
		fModulus(0,0) = fModulus(1,1) =  2.0*fMu*(1.0 - third);
		fModulus(1,0) = fModulus(0,1) = -2.0*fMu*third;
	}
	else if (NumSD() == 3) {
		fModulus(0,0) = fModulus(1,1) =  fModulus(2,2) = 2.0*fMu*(1.0 - third);
		fModulus(3,3) = fModulus(4,4) =  fModulus(5,5) = fMu;
		fModulus(0,1) = fModulus(0,2) =  fModulus(1,2) = -2.0*fMu*third;
		fModulus(1,0) = fModulus(2,0) =  fModulus(2,1) = -2.0*fMu*third;
	}
return fModulus; 
}

/* fluid stress */
const dSymMatrixT& FluidMaterialT::s_ij(void)
{
  WriteCallLocation("s_ij"); //DEBUG

  /* should be 1 row */
	fStrainRate.Symmetrize(fFluidMatSupport->VelGrad());
	fStress = fStrainRate;
	fStress *= 2.0*fMu;
	double pressure = fFluidMatSupport->Pressure();
	for (int i = 0; i < NumSD(); i++) {
		fStress[i] -= pressure;
	}
	return fStress;
}

/* describe the parameters needed by the interface */
void FluidMaterialT::DefineParameters(ParameterListT& list) const
{
  WriteCallLocation("DefineParameters"); //DEBUG

  /* inherited */
	ContinuumMaterialT::DefineParameters(list);

	/* define parameters */
	list.AddParameter(fDensity, "density");
	list.AddParameter(ParameterT::Double, "viscosity");
}

/* accept parameter list */
void FluidMaterialT::TakeParameterList(const ParameterListT& list)
{
  WriteCallLocation("TakeParameterList"); //DEBUG

  /* inherited */
	ContinuumMaterialT::TakeParameterList(list);

	/* get parameters */
	fDensity = list.GetParameter("density");
	fMu = list.GetParameter("viscosity");
}

/** FOR DEBUGGING PURPOSES ONLY */
void FluidMaterialT::WriteCallLocation( char* loc ) const {
cout << "Inside of FluidMaterialT::" << loc << endl;
}