/* $Header: /home/regueiro/tahoe_cloudforge_repo_snapshots/development/src/elements/fluid_element/materials/FluidMatSupportT.cpp,v 1.2 2006-07-18 01:21:21 a-kopacz Exp $ */
/* created: tdnguye (07/12/2006) */

#include "FluidMatSupportT.h"
#include "ElementsConfig.h"

#ifdef CONTINUUM_ELEMENT
#include "FluidElementT.h"
#endif

using namespace Tahoe;

/* constructor */
FluidMatSupportT::FluidMatSupportT(int ndof, int nip):
	MaterialSupportT(ndof, nip),
	fVel_list(NULL),
	fPres_list(NULL),
	fGradVel_list(NULL),
	fGradPres_list(NULL),
	fFluid(NULL)
{
  WriteCallLocation("FluidMatSupportT"); //DEBUG
}

/* set the source for the gradient information */
void FluidMatSupportT::SetField(const ArrayT<dArrayT>* fvel_list, const dArrayT* fpres_list)
{
  WriteCallLocation("SetField"); //DEBUG

  fVel_list = fvel_list;
	fPres_list = fpres_list;
}

/* set the source for the gradient information */
void FluidMatSupportT::SetGradient(const ArrayT<dMatrixT>* gradient_vel_list, const ArrayT<dArrayT>* gradient_pres_list)
{
  WriteCallLocation("SetGradient"); //DEBUG

//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
#if 0
	if (!gradient_list ||
	     gradient_list->Length() != NumIP())
	{
		cout << "\n FluidMatSupportT::SetGradient: inconsistent gradient source" 
		     << endl;
		throw ExceptionT::kGeneralFail;
	}
#endif

	fGradVel_list = gradient_vel_list;
	fGradPres_list = gradient_pres_list;
}

/* set the element group pointer */
void FluidMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
  WriteCallLocation("SetContinuumElement"); //DEBUG

  /* inherited */
	MaterialSupportT::SetContinuumElement(p);

#ifdef CONTINUUM_ELEMENT
	/* cast to small strain pointer */
	fFluid = TB_DYNAMIC_CAST(const FluidElementT*, p);
#endif
}

/** FOR DEBUGGING PURPOSES ONLY */
void FluidMatSupportT::WriteCallLocation( char* loc ) const {
cout << "Inside of FluidMatSupportT::" << loc << endl;
}