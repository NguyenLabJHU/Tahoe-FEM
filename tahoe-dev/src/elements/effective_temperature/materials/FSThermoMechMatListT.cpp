/* created: rxiao (01/01/2014) */
#include "FSThermoMechMatListT.h"
#include "FSThermoMechSupportT.h"
#include "DevelopmentElementsConfig.h"
/*  materials */
#include "SMP_coupled.h"

using namespace Tahoe;

/* constructors */
FSThermoMechMatListT::	FSThermoMechMatListT(int length, const FSThermoMechSupportT& support):
	SolidMatListT(length,support),
	fFSThermoMechSupport(&support)
{
	SetName("fsthermomech_material");
}

FSThermoMechMatListT::	FSThermoMechMatListT(void):
	fFSThermoMechSupport(NULL)
{
	SetName("fsthermomech_material");
}

/* information about subordinate parameter lists */
void FSThermoMechMatListT::DefineSubs(SubListT& sub_list) const
{
	//WriteCallLocation("DefineSubs"); //DEBUG

	/* inherited */
	SolidMatListT::DefineSubs(sub_list);

	/* an array of choices */
	sub_list.AddSub("fsthermomech_material_list", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void FSThermoMechMatListT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	//WriteCallLocation("DefineInlineSub"); //DEBUG

	/* list of choice of materials */
	if (name == "fsthermomech_material_list")
	{
		order = ParameterListT::Choice;
		/* fluid materials */
		sub_lists.AddSub("SMP_coupled");
	}	
	else /* inherited */
		SolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSThermoMechMatListT::NewSub(const StringT& name) const
{
	//WriteCallLocation("NewSub"); //DEBUG

	/* try to construct material */
	FSThermoMechMatT* material = NewFSThermoMechMat(name);
	if (material)
		return material;
	else /* inherited */
		return SolidMatListT::NewSub(name);
}

/* accept parameter list */
void FSThermoMechMatListT::TakeParameterList(const ParameterListT& list)
{
	//WriteCallLocation("TakeParameterList"); //DEBUG

	/* inherited */
	SolidMatListT::TakeParameterList(list);

	/* construct materials - NOTE: subs have been defined as a choice, but
	* here we construct as many materials as are passed in */
//	AutoArrayT<FluidMaterialT*> materials;
	const ArrayT<ParameterListT>& subs = list.Lists();
	for (int i = 0; i < subs.Length(); i++) 
	{
		const ParameterListT& sub = subs[i];
		FSThermoMechMatT* mat = NewFSThermoMechMat(sub.Name());
        int count = 0;
		if (mat) 
		{
            /* store pointer */
			(*this)[count++] = mat;
            //materials.Append(mat);
			mat->TakeParameterList(sub);
            if (mat->HasHistory()) fHasHistory = true;
			if (mat->HasThermalStrain()) fHasThermal = true;
			if (mat->HasLocalization()) fHasLocalizers = true;
		}
	}

	/* transfer */
//	Dimension(materials.Length());
//	for (int i = 0; i < materials.Length(); i++)
//		fArray[i] = materials[i];

}

/* construct the specified material or NULL if the request cannot be completed */
FSThermoMechMatT* FSThermoMechMatListT::NewFSThermoMechMat(const StringT& name) const
{
	//WriteCallLocation("NewFluidMaterial"); //DEBUG

	FSThermoMechMatT* mat = NULL;

	if (name == "SMP_coupled")
		mat = new SMP_coupled;

	/* set support */
	if (mat) mat->SetFSThermoMechSupport(fFSThermoMechSupport);

	return mat;
}

