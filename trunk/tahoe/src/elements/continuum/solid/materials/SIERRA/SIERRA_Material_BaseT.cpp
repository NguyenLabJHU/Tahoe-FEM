/* $Id: SIERRA_Material_BaseT.cpp,v 1.15 2004-08-01 01:00:59 paklein Exp $ */
#include "SIERRA_Material_BaseT.h"
#include "SIERRA_Material_DB.h"
#include "SIERRA_Material_Data.h"
#include "SpectralDecompT.h"
#include "ParameterListT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"

using namespace Tahoe;

/* initialize static variables */
int SIERRA_Material_BaseT::sSIERRA_Material_count = 0;
const int kSIERRA_stress_dim = 6;

/* constructor */
SIERRA_Material_BaseT::SIERRA_Material_BaseT(void):
	ParameterInterfaceT("SIERRA_material"),
	fTangentType(GlobalT::kSymmetric),
	fSIERRA_Material_Data(NULL),
	fPressure(0.0),
	fdstran(kSIERRA_stress_dim),
	fDecomp(NULL),
	fDebug(false)
{
	const char caller[] = "SIERRA_Material_BaseT::SIERRA_Material_BaseT";

	/* instantiate materials database */
	if (++sSIERRA_Material_count == 1)
		SIERRA_Material_DB::Create();	
}

/* destructor */
SIERRA_Material_BaseT::~SIERRA_Material_BaseT(void)
{
	/* free spectral decomposition object */
	delete fDecomp;
	fDecomp = NULL;

	/* free materials database */
	if (--sSIERRA_Material_count == 0)
		SIERRA_Material_DB::Delete();
}

/* materials initialization */
bool SIERRA_Material_BaseT::NeedsPointInitialization(void) const { return true; }
void SIERRA_Material_BaseT::PointInitialize(void)
{
	/* allocate element storage */
	if (CurrIP() == 0)
	{
		ElementCardT& element = CurrentElement();
		element.Dimension(0, fBlockSize*NumIP());
	
		/* initialize */
		element.DoubleData() = 0.0;
	}

	/* load stored data */
	Load(CurrentElement(), CurrIP());
	
	/* parameters */
	int nelem = 1;
	double dt = 0.0; //OK?
	int nsv = fstate_old.Length();
	int ncd = 0;
	int matvals = fSIERRA_Material_Data->ID();

	/* call the initialization function */
	Sierra_function_material_init init_func = fSIERRA_Material_Data->InitFunction();
	init_func(&nelem, &dt, &nsv, fstate_old.Pointer(), fstate_new.Pointer(), &matvals, &ncd);

	/* write to storage */
	Store(CurrentElement(), CurrIP());

	/* store results as last converged */
	if (CurrIP() == NumIP() - 1) UpdateHistory();
}

/* update/reset internal variables */
void SIERRA_Material_BaseT::UpdateHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load stored data */
		Load(element, ip);
	
		/* assign "current" to "old" */	
		fstress_old = fstress_new;
		fstate_old = fstate_new;

		/* write to storage */
		Store(element, ip);
	}
}

void SIERRA_Material_BaseT::ResetHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load stored data */
		Load(element, ip);

		/* assign "old" to "current" */	
		fstress_new = fstress_old;
		fstate_new = fstate_old;

		/* write to storage */
		Store(element, ip);
	}
}

const dSymMatrixT& SIERRA_Material_BaseT::s_ij(void)
{
	/* call calc function */
	if (MaterialSupport().RunState() == GlobalT::kFormRHS ||
		MaterialSupport().RunState() == GlobalT::kFormLHS)
	{
		/* load stored data */
		Load(CurrentElement(), CurrIP());

		/* compute strains and other calc functions */
		Set_Calc_Arguments();
	
		/* parameters */
		int nelem = 1;
		double dt = fFSMatSupport->TimeStep();
		int nsv = fstate_old.Length();
		int ncd = 0;
		int matvals = fSIERRA_Material_Data->ID();
		int ivars_size = fSIERRA_Material_Data->InputVariables().Length();

		/* debug information */
		if (fDebug) {
			cout << "\n SIERRA_Material_BaseT::s_ij: IN\n"
				 << " element: " << CurrElementNumber()+1 << '\n'
				 << "      ip: " << CurrIP()+1 << '\n';
			
			cout << " rot strain inc = " << fdstran.no_wrap() << '\n';
			cout << " old stress = " << fstress_old.no_wrap() << '\n';
			cout << " old state =\n" << fstate_old.wrap(5) << '\n';
		}

		/* call the calc function */
		Sierra_function_material_calc calc_func = fSIERRA_Material_Data->CalcFunction();
		calc_func(&nelem, &dt, fdstran.Pointer(), &ivars_size,
			fstress_old_rotated.Pointer(), fstress_new.Pointer(), 
			&nsv, fstate_old.Pointer(), fstate_new.Pointer(), 
			&matvals, &ncd);

		/* debug information */
		if (fDebug) {
			cout << "\n SIERRA_Material_BaseT::s_ij: OUT\n";
			cout << " new stress = " << fstress_new.no_wrap() << '\n';
			cout << " new state =\n" << fstate_new.wrap(5) << '\n';
		}

		/* write to storage */
		Store(CurrentElement(), CurrIP());
	}
	else
		/* load stored data */
		Load(CurrentElement(), CurrIP());

	/* copy/convert stress */
	SIERRA_to_dSymMatrixT(fstress_new.Pointer(), fStress);
	fPressure = fStress.Trace()/3.0;
	return fStress;
}

/* return the pressure associated with the last call to SolidMaterialT::s_ij */
double SIERRA_Material_BaseT::Pressure(void) const
{
	/* load stored data (cast const-ness) */
	SIERRA_Material_BaseT* non_const_this = (SIERRA_Material_BaseT*) this;
	non_const_this->Load(CurrentElement(), CurrIP());
	return (fstress_new[3] + fstress_new[4] + fstress_new[5])/3.0;
}

/* returns the strain energy density for the specified strain */
double SIERRA_Material_BaseT::StrainEnergyDensity(void)
{
	return 0.0; /* not part of the Sierra materials interface */
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int SIERRA_Material_BaseT::NumOutputVariables(void) const
{
	/* set material output variables/labels */
	if (fOutputIndex.Length() == 0)
	{
		//TEMP - better place for this?
		SIERRA_Material_BaseT* tmp = (SIERRA_Material_BaseT*) this;
		tmp->SetOutputVariables(tmp->fOutputIndex, tmp->fOutputLabels);
	}
	return fOutputIndex.Length();
}

void SIERRA_Material_BaseT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(fOutputLabels.Length());
	for (int i = 0; i < labels.Length(); i++)
		labels[i] = fOutputLabels[i];
}

void SIERRA_Material_BaseT::ComputeOutput(dArrayT& output)
{
	/* check */
	if (output.Length() != fOutputIndex.Length())
		ExceptionT::SizeMismatch("SIERRA_Material_BaseT::ComputeOutput", 
			"output array should be length %d not %d", fOutputIndex.Length(), output.Length());

	/* load stored data */
	Load(CurrentElement(), CurrIP());

	/* collect variables */
	for (int i = 0; i < fOutputIndex.Length(); i++)
		output[i] = double(fstate_new[fOutputIndex[i]]);
}

/* describe the parameters needed by the interface */
void SIERRA_Material_BaseT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FSIsotropicMatT::DefineParameters(list);

	/* debugging flag */
	ParameterT debug(fDebug, "debug");
	debug.SetDefault(fDebug);
	list.AddParameter(debug);
	
	/* file with Sierra materials parameters */
	list.AddParameter(ParameterT::Word, "SIERRA_parameter_file");
}
	
/* accept parameter list */
void SIERRA_Material_BaseT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SIERRA_Material_BaseT::TakeParameterList";

	/* inherited */
	FSIsotropicMatT::TakeParameterList(list);

	/* 3D only */
	int nsd = NumSD();
	if (nsd != 3) ExceptionT::GeneralFail(caller, "3D only");

	/* dimension work space */
	fstress_old_rotated.Dimension(kSIERRA_stress_dim);
	fF_rel.Dimension(nsd);
	fA_nsd.Dimension(nsd);
	fU1.Dimension(nsd);
	fU2.Dimension(nsd); 
	fU1U2.Dimension(nsd);

	/* spectral decomp */
	fDecomp = new SpectralDecompT(nsd);

	/* call SIERRA registration function */
	Register_SIERRA_Material();

	/* extract parameters */
	fDebug = list.GetParameter("debug");
	
	/* open Sierra parameters file */
	StringT path;
	path.FilePath(MaterialSupport().InputFile());
	StringT params = list.GetParameter("SIERRA_parameter_file");
	params.ToNativePathName();
	params.Prepend(path);
	ifstreamT in('#', params);
	if (!in.is_open())
		ExceptionT::GeneralFail(caller, "could not open file \"%s\"",
			params.Pointer());

	/* read SIERRA-format input */
	StringT line;
	line.GetLineFromStream(in);
	StringT word;
	int count;
	word.FirstWord(line, count, true);
	if (word != "begin")
		ExceptionT::BadInputValue(caller, "expecting \"begin\":\n%s", line.Pointer());
	line.Tail(' ', word);
	word.ToUpper();
	ParameterListT param_list(word);
	Read_SIERRA_Input(in, param_list);
	fSIERRA_Material_Data = Process_SIERRA_Input(param_list);

	/* check parameters */
	Sierra_function_param_check param_check = fSIERRA_Material_Data->CheckFunction();
	int material_ID = fSIERRA_Material_Data->ID();
	param_check(&material_ID);

	/* set density */
	fDensity = fSIERRA_Material_Data->Property("DENSITY");

	/* set the modulus */
	double kappa  = fSIERRA_Material_Data->Property("BULK_MODULUS");
	double mu = fSIERRA_Material_Data->Property("TWO_MU")/2.0;
	IsotropicT::Set_mu_kappa(mu, kappa);

	/* storage block size (per ip) */
	int nsv = fSIERRA_Material_Data->NumStateVariables();
	fBlockSize = 0;
	fBlockSize += kSIERRA_stress_dim; // fstress_old
	fBlockSize += kSIERRA_stress_dim; // fstress_new
	fBlockSize += nsv;  // fstate_old
	fBlockSize += nsv;  // fstate_new
	
	/* argument array */
	fArgsArray.Dimension(fBlockSize);

	/* assign pointers */
	double* parg = fArgsArray.Pointer();
	
	fstress_old.Set(kSIERRA_stress_dim, parg); parg += kSIERRA_stress_dim;
	fstress_new.Set(kSIERRA_stress_dim, parg); parg += kSIERRA_stress_dim;
	fstate_old.Set(nsv, parg); parg += nsv;
	fstate_new.Set(nsv, parg);	

	/* notify */
	if (fThermal->IsActive())
		cout << "\n SIERRA_Material_BaseT::Initialize: thermal strains must\n"
		     <<   "    be handled within the UMAT\n" << endl;
	
	/* write properties array */
	ofstreamT& out = MaterialSupport().Output();
	out << " Material name . . . . . . . . . . . . . . . . . = " << fMaterialName << '\n';
	out << " Material model name . . . . . . . . . . . . . . = " << fSIERRA_Material_Data->Name() << '\n';
	out << " Number of state variables . . . . . . . . . . . = " << fSIERRA_Material_Data->NumStateVariables() << '\n';
	
	/* material properties */
	const ArrayT<StringT>& prop_names = fSIERRA_Material_Data->PropertyNames();
	const ArrayT<double>&  prop_values  = fSIERRA_Material_Data->PropertyValues();
	int d_width = OutputWidth(out, prop_values.Pointer());
	out << " Number of material properties . . . . . . . . . = " << prop_names.Length() << '\n';
	for (int i = 0; i < prop_names.Length(); i++)
		out << setw(d_width) << prop_values[i] << " : " << prop_names[i] << '\n';
	out.flush();

	out << "    SIERRA material: " << fSIERRA_Material_Data->Name() << '\n';	
}

/***********************************************************************
 * Protected
 ***********************************************************************/

void SIERRA_Material_BaseT::SIERRA_to_dSymMatrixT(const double* pA,
	dSymMatrixT& B) const
{
	double* pB = B.Pointer();
	*pB++ = pA[0]; // 11
	*pB++ = pA[1]; // 22
	*pB++ = pA[2]; // 33
	*pB++ = pA[4]; // 23
	*pB++ = pA[5]; // 13
	*pB   = pA[3]; // 12
}

void SIERRA_Material_BaseT::dSymMatrixT_to_SIERRA(const dSymMatrixT& A,
	double* pB) const
{
	const double* pA = A.Pointer();	
	*pB++ = pA[0]; // 11
	*pB++ = pA[1]; // 22
	*pB++ = pA[2]; // 33
	*pB++ = pA[5]; // 12
	*pB++ = pA[3]; // 23
	*pB   = pA[4]; // 31
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* read parameters from input stream */
void SIERRA_Material_BaseT::Read_SIERRA_Input(ifstreamT& in, 
	ParameterListT& param_list) const
{
	const char caller[] = "SIERRA_Material_BaseT::Read_SIERRA_Input";

	/* set comment character */
	char old_comment_marker = in.comment_marker();
	in.set_marker('#');

	bool C_word_only = true;
	int char_count;
	
	StringT line;
	line.GetLineFromStream(in);
	bool done = false;
	while (!done)
	{
		StringT word;
		word.FirstWord(line, char_count, C_word_only);
		word.ToUpper();
		
		/* start of new list */
		if (word == "BEGIN")
		{
			/* get list name */
			StringT name;
			line.Tail(' ', name);
			name.ToUpper();
			if (name.StringLength() == 0)
				ExceptionT::BadInputValue(caller, "could not list name from line:\n%s",
					line.Pointer());
					
			/* recursively construct list */
			ParameterListT param_sub_list(name);
			Read_SIERRA_Input(in, param_sub_list);
			
			/* add list to parameter list */
			if (!param_list.AddList(param_sub_list))
				ExceptionT::BadInputValue(caller, "list is duplicate: \"%s\"",
					param_sub_list.Name().Pointer());
		}
		else if (word == "END") /* end of this list */
		{
			/* get list name */
			StringT name;
			line.Tail(' ', name);
			name.ToUpper();
			if (name != param_list.Name())
				ExceptionT::BadInputValue(caller, "expecting end for \"%s\" not \"%s\"",
					param_list.Name().Pointer(), name.Pointer());
			done = true;
		}
		else /* split parameter name and value */
		{
			/* find equals sign */
			int equal_dex = line.FirstPositionOf('=');
			if (equal_dex < 1)
				ExceptionT::BadInputValue(caller, "line does not contain \"=\":\n%s",
					line.Pointer());
			StringT param_name;
			param_name.Take(line, equal_dex);
			param_name.DropLeadingSpace();
			param_name.DropTrailingSpace();
			param_name.Replace(' ', '_');
			param_name.ToUpper();
		
			/* get value */
			double value;
			if (!line.Tail('=', value))
				ExceptionT::BadInputValue(caller, "could not extract value from line:\n%s",
					line.Pointer());
			
			/* new parameter */
			ParameterT param(value, param_name);
			if (!param_list.AddParameter(param))
				ExceptionT::BadInputValue(caller, "parameter is duplicate: \"%s\"",
					param_name.Pointer());
		}
		
		/* get next line */
		if (!done) line.GetLineFromStream(in);
	}

	/* restore comment character */
	in.set_marker(old_comment_marker);
}

SIERRA_Material_Data* SIERRA_Material_BaseT::Process_SIERRA_Input(ParameterListT& param_list)
{
	const char caller[] = "SIERRA_Material_BaseT::Process_SIERRA_Input";
	fMaterialName = param_list.Name();

	/* sublist - gives model name */
	const ArrayT<ParameterListT>& sub_list = param_list.Lists();
	if (sub_list.Length() != 1)
		ExceptionT::BadInputValue(caller, "expecting 1 parameter sub-list: %d", sub_list.Length());

	const ParameterListT& model_param_list = sub_list[0];
	const StringT& model_name = model_param_list.Name();
	SIERRA_Material_Data* mat_data = SIERRA_Material_DB::Material(model_name);

	/* first add top level parameters */
	const ArrayT<ParameterT>& params = param_list.Parameters();
	for (int i = 0; i < params.Length(); i++)
		mat_data->AddProperty(params[i].Name(), params[i]);

	/* add model parameters */
	const ArrayT<ParameterT>& model_params = model_param_list.Parameters();
	for (int i = 0; i < model_params.Length(); i++)
		mat_data->AddProperty(model_params[i].Name(), model_params[i]);

	return mat_data;
}

/* load element data for the specified integration point */
void SIERRA_Material_BaseT::Load(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy values */
	fArgsArray.CopyPart(0, d_array, fBlockSize*ip, fBlockSize);
}

void SIERRA_Material_BaseT::Store(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* write back */
	d_array.CopyPart(fBlockSize*ip, fArgsArray, 0, fBlockSize);
}

/* set stress/strain arguments */
void SIERRA_Material_BaseT::Set_Calc_Arguments(void)
{
	const char caller[] = "SIERRA_Material_BaseT::Set_Calc_Arguments";

	/* determine material input */
	const ArrayT<StringT>& input = fSIERRA_Material_Data->InputVariables();
	
//TEMP
if (input.Length() != 1)
	ExceptionT::GeneralFail(caller, "expecting just 1 material input not %d", input.Length());

	/* relative deformation gradient */
	fA_nsd = F_total_last();
	const dMatrixT& F_n = F();
	fA_nsd.Inverse();
	fF_rel.MultAB(F_n, fA_nsd);

	/* polar decomposition */
	bool perturb_repeated_roots = false;
	fDecomp->PolarDecomp(fF_rel, fA_nsd, fU1, perturb_repeated_roots);
	fU2 = fU1;
	fU1.PlusIdentity(-1.0);
	fU2.PlusIdentity( 1.0);
	fU2.Inverse();
	fU1U2.MultAB(fU1, fU2);
	
	/* incremental strain */
	if (input[0] == "rot_strain_increment")
	{
		fdstran[0] = 2.0*fU1U2[0]; // 11
		fdstran[1] = 2.0*fU1U2[1]; // 22
		fdstran[2] = 2.0*fU1U2[2]; // 33
		fdstran[3] = 2.0*fU1U2[5]; // 12
		fdstran[4] = 2.0*fU1U2[3]; // 23
		fdstran[5] = 2.0*fU1U2[4]; // 31
	}
	/* incremental strain rate */
	else if (input[0] == "rot_strain_inc")
	{
		double dt = fFSMatSupport->TimeStep();
		double k = (fabs(dt) > kSmall) ? 2.0/dt : 0.0;
		fdstran[0] = k*fU1U2[0]; // 11
		fdstran[1] = k*fU1U2[1]; // 22
		fdstran[2] = k*fU1U2[2]; // 33
		fdstran[3] = k*fU1U2[5]; // 12
		fdstran[4] = k*fU1U2[3]; // 23
		fdstran[5] = k*fU1U2[4]; // 31
	}
	else 
		ExceptionT::GeneralFail(caller, "unrecognized input \"%s\"", 
			input[0].Pointer());

	/* rotate old stress to current configuration */
	SIERRA_to_dSymMatrixT(fstress_old.Pointer(), fU1);
	fU2.MultQBQT(fA_nsd, fU1);
	dSymMatrixT_to_SIERRA(fU2, fstress_old_rotated.Pointer());
}
