/* $Id: FSSolidMixtureT.cpp,v 1.3 2004-11-05 22:53:49 paklein Exp $ */
#include "FSSolidMixtureT.h"
#include "ParameterContainerT.h"
//#include "FSSolidMixtureSupportT.h"
#include "FieldT.h"

/* stress functions */
#include "FDHookeanMatT.h"
#include "FDKStV.h"
#include "SimoIso3D.h"

using namespace Tahoe;

/* constructor */
FSSolidMixtureT::FSSolidMixtureT(void):
	ParameterInterfaceT("large_strain_solid_mixture"),
//	fFSSolidMixtureSupport(NULL),
	fStressSupport(NULL),
	fConc(LocalArrayT::kDisp)
{

}

/* destructor */
FSSolidMixtureT::~FSSolidMixtureT(void)
{
	/* free stress function for each species */
	for (int i = 0; i < fStressFunctions.Length(); i++)
		delete fStressFunctions[i];
	delete fStressSupport;
}

#if 0
/* set the material support or pass NULL to clear */
void FSSolidMixtureT::SetFSSolidMixtureSupport(const FSSolidMixtureSupportT* support)
{
	/* inherited */
	FSSolidMatT::SetFSMatSupport(support);
	
	fFSSolidMixtureSupport = support;
}
#endif

/* get all nodal concentrations over the current element */
void FSSolidMixtureT::UpdateConcentrations(void)
{
	/* current element */
	const ElementCardT& element = CurrentElement();

	/* collect nodal values */
	const iArrayT& nodes_u = element.NodesU();
	for (int i = 0; i < fFields.Length(); i++)
	{
		const FieldT& field = *(fFields[i]);
		const dArray2DT& c = field[0];
		double* pc = fConc(i);
		for (int j = 0; j < nodes_u.Length(); j++)
			*pc++ = c[nodes_u[j]];
	}
}

/* update given nodal concentrations over the current element */
void FSSolidMixtureT::UpdateConcentrations(int i)
{
	/* current element */
	const ElementCardT& element = CurrentElement();

	/* collect nodal values */
	const iArrayT& nodes_u = element.NodesU();
	const FieldT& field = *(fFields[i]);
	const dArray2DT& c = field[0];
	double* pc = fConc(i);
	for (int j = 0; j < nodes_u.Length(); j++)
		*pc++ = c[nodes_u[j]];
}

/** return the index of the species associated with the given field name */
int FSSolidMixtureT::SpeciesIndex(const StringT& field_name) const
{
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->FieldName() == field_name)
			return i;
	
	/* not found */
	return -1;
}

/* mass density */
double FSSolidMixtureT::Density(void)
{
	/* update concentrations */
	if (CurrIP() == 0) UpdateConcentrations();
	fFSMatSupport->Interpolate(fConc, fIPConc);

	fDensity = fIPConc.Sum();
	return fDensity;
}

/* strain energy density */
double FSSolidMixtureT::StrainEnergyDensity(void)
{
	/* current element information */
	const ElementCardT& element = CurrentElement();
	const dArrayT& conc_0 = element.DoubleData();

	/* update concentrations */
	if (CurrIP() == 0) UpdateConcentrations();
	fFSMatSupport->Interpolate(fConc, fIPConc);

	/* sum over species */
	double u = 0.0;
	for (int i = 0; i < fConc.Length(); i++)
	{
		/* compute mechanical strain */
		double rel_conc = fConc[i]/conc_0[i];
		fF_growth_inv.Identity(1.0/rel_conc);
		fF_species[0].MultAB(fFSMatSupport->DeformationGradient(), fF_growth_inv);

		/* compute modulus */
		u += fConc[i]*fStressFunctions[i]->StrainEnergyDensity();
	}

	return u;
}

/* total material tangent modulus */
const dMatrixT& FSSolidMixtureT::c_ijkl(void)
{
	/* current element information */
	const ElementCardT& element = CurrentElement();
	const dArrayT& conc_0 = element.DoubleData();

	/* update concentrations */
	if (CurrIP() == 0) UpdateConcentrations();
	fFSMatSupport->Interpolate(fConc, fIPConc);
	const dArrayT& conc = fIPConc;
	//const dArrayT& conc = fFSSolidMixtureSupport->Concentration();

	/* sum over species */
	fModulus = 0.0;
	for (int i = 0; i < conc.Length(); i++)
	{
		/* compute mechanical strain */
		double rel_conc = conc[i]/conc_0[i];
		fF_growth_inv.Identity(1.0/rel_conc);
		fF_species[0].MultAB(fFSMatSupport->DeformationGradient(), fF_growth_inv);

		/* compute modulus */
		fModulus.AddScaled(conc[i], fStressFunctions[i]->c_ijkl());
	}

	return fModulus;
}

/* partial material tangent modulus */
const dMatrixT& FSSolidMixtureT::c_ijkl(int i)
{
	/* current element information */
	const ElementCardT& element = CurrentElement();
	const dArrayT& conc_0 = element.DoubleData();

	/* concentrations */
	fFSMatSupport->Interpolate(fConc, fIPConc);
	const dArrayT& conc = fIPConc;	
//	const dArrayT& conc = fFSSolidMixtureSupport->Concentration();

	/* compute mechanical strain */
	double rel_conc = conc[i]/conc_0[i];
	fF_growth_inv.Identity(1.0/rel_conc);
	fF_species[0].MultAB(fFSMatSupport->DeformationGradient(), fF_growth_inv);

	/* compute stress */
	fModulus.SetToScaled(conc[i], fStressFunctions[i]->c_ijkl());

	return fModulus;
}

/* total Cauchy stress */
const dSymMatrixT& FSSolidMixtureT::s_ij(void)
{
	/* current element information */
	const ElementCardT& element = CurrentElement();
	const dArrayT& conc_0 = element.DoubleData();

	/* update concentrations */
	if (CurrIP() == 0) UpdateConcentrations();
	fFSMatSupport->Interpolate(fConc, fIPConc);
	const dArrayT& conc = fIPConc;
	//const dArrayT& conc = fFSSolidMixtureSupport->Concentration();

	/* sum over species */
	fStress = 0.0;
	for (int i = 0; i < conc.Length(); i++)
	{
		/* compute mechanical strain */
		double rel_conc = conc[i]/conc_0[i];
		fF_growth_inv.Identity(1.0/rel_conc);
		fF_species[0].MultAB(fFSMatSupport->DeformationGradient(), fF_growth_inv);

		/* compute stress */
		fStress.AddScaled(conc[i], fStressFunctions[i]->s_ij());
	}

	return fStress;
}

/* partial Cauchy stress */
const dSymMatrixT& FSSolidMixtureT::s_ij(int i)
{
	/* current element information */
	const ElementCardT& element = CurrentElement();
	const dArrayT& conc_0 = element.DoubleData();

	/* concentrations */
	fFSMatSupport->Interpolate(fConc, fIPConc);
	const dArrayT& conc = fIPConc;	
//	const dArrayT& conc = fFSSolidMixtureSupport->Concentration();

	/* compute mechanical strain */
	double rel_conc = conc[i]/conc_0[i];
	fF_growth_inv.Identity(1.0/rel_conc);
	fF_species[0].MultAB(fFSMatSupport->DeformationGradient(), fF_growth_inv);

	/* compute stress */
	fStress.SetToScaled(conc[i], fStressFunctions[i]->s_ij());

	return fStress;
}

/* initialization */
void FSSolidMixtureT::PointInitialize(void)
{
	/* info */
	int cip = CurrIP();
	int nip = NumIP();
	ElementCardT& element = CurrentElement();

	/* allocate element storage */
	if (cip == 0)
	{
		/* storage size (per integration point) */
		int i_size = 0;
		int d_size = 0;
		d_size += fStressFunctions.Length(); /* reference concentrations */
	
		/* allocate storage */
		element.Dimension(i_size*nip, d_size*nip);

		/* collect element concentrations */
		UpdateConcentrations();		
	}
	
	/* initialize concentrations */
	fFSMatSupport->Interpolate(fConc, fIPConc);
	int nc = fIPConc.Length();
//	element.DoubleData() = fIPConc;	
	element.DoubleData().CopyPart(cip*nc, fIPConc, 0, nc);
//	element.DoubleData() = fFSSolidMixtureSupport->Concentration();

	/* store results as last converged */
	if (cip == nip - 1) UpdateHistory();
}

/* describe the parameters needed by the interface */
void FSSolidMixtureT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FSSolidMatT::DefineParameters(list);
}

/* information about subordinate parameter lists */
void FSSolidMixtureT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);

	/* data for each species */
	sub_list.AddSub("mixture_species", ParameterListT::OnePlus);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSSolidMixtureT::NewSub(const StringT& name) const
{
	/* stress function */
	ParameterInterfaceT* solid = New(name);
	if (solid)
		return solid;

	if (name == "mixture_species")
	{
		ParameterContainerT* species = new ParameterContainerT(name);
		species->SetSubSource(this);

		/* name of the associated concentration field */
		species->AddParameter(ParameterT::Word, "concentration_field");
		
		/* choice of stress functions */
		species->AddSub("species_stress_function_choice", ParameterListT::Once, true);
	
		//other species properties ?????????
	
		return species;
	}
	else if (name == "species_stress_function_choice")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
		
		/* stress functions */
		choice->AddSub("large_strain_Hookean");
		choice->AddSub("large_strain_StVenant");
		choice->AddSub("Simo_isotropic");
	
		return choice;
	}
	else /* inherited */
		return FSSolidMatT::NewSub(name);
}

/* accept parameter list */
void FSSolidMixtureT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "FSSolidMixtureT::TakeParameterList";

	/* inherited */
	FSSolidMatT::TakeParameterList(list);

	/* construct support for stress functions */
	int nsd = NumSD();
	fStressSupport = new FSMatSupportT(nsd, 1); /* single integration point */
	fF_species.Dimension(1);
	fF_species[0].Dimension(nsd);
	fF_growth_inv.Dimension(nsd);
	fStressSupport->SetContinuumElement(MaterialSupport().ContinuumElement());
	fStressSupport->SetDeformationGradient(&fF_species);

	/* species */
	int num_species = list.NumLists("mixture_species");
	fStressFunctions.Dimension(num_species);
	fStressFunctions = NULL;
	fFields.Dimension(num_species);
	fFields = NULL;
	for (int i = 0; i < num_species; i++)
	{
		/* species */
		const ParameterListT& species = list.GetList("mixture_species", i);

		/* concentration field */
		const StringT& conc_field_name = list.GetParameter("concentration_field");
		const FieldT* conc_field = fFSMatSupport->Field(conc_field_name);
		if (!conc_field)
			ExceptionT::GeneralFail(caller, "could not resolve field \"%s\"", 
				conc_field_name.Pointer());
		if (conc_field->NumDOF() != 1)
			ExceptionT::GeneralFail(caller, "field \"%s\" has dimension %d not 1", 
				conc_field_name.Pointer(), conc_field->NumDOF());
		fFields[i] = conc_field;
	
		/* resolve stress function */
		const ParameterListT& solid_params = species.GetListChoice(*this, "species_stress_function_choice");

		/* construct material */	
		FSSolidMatT* solid = New(solid_params.Name());
		if (!solid)
			ExceptionT::GeneralFail(caller, "could not construct \"%s\"", 
				solid_params.Name().Pointer());
	
		/* initialize */
		solid->SetFSMatSupport(fStressSupport);
		solid->TakeParameterList(solid_params);
		
		/* store */
		fStressFunctions[i] = solid;
	}

	/* work space */
	fConc.Dimension(NumElementNodes(), num_species);
	fIPConc.Dimension(num_species);
}

/* return the specified stress function or NULL */
FSSolidMatT* FSSolidMixtureT::New(const StringT& name) const
{
	if (name == "large_strain_Hookean")
		return new FDHookeanMatT;
	else if (name == "large_strain_StVenant")
		return new FDKStV;
	else if (name == "Simo_isotropic")
		return new SimoIso3D;
	else
		return NULL;
}
