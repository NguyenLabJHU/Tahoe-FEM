/* $Id: SolidElementT.cpp,v 1.78 2007-10-09 23:17:43 r-jones Exp $ */
#include "SolidElementT.h"

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "ElementCardT.h"
#include "ShapeFunctionT.h"
#include "eIntegratorT.h"
#include "iAutoArrayT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "OutputSetT.h"

/* materials */
#include "SolidMaterialT.h"
#include "SolidMatSupportT.h"
#include "SolidMatListT.h"

/* exception codes */
#include "ExceptionT.h"

/* eigenvalue estimation */
#include "CCSMatrixT.h"
#include "NodeManagerT.h"

/* needed for Poynting vector output */
#include "FSSolidMatT.h"

using namespace Tahoe;

/* initialize static data */
const int SolidElementT::NumNodalOutputCodes = 8;
static const char* NodalOutputNames[] = {
	"coordinates",
	"displacements",
	"stress",
	"principal_stress",
	"strain_energy_density",
	"wave_speeds",
	"material_output",
	"Poynting_vector"};

const int SolidElementT::NumElementOutputCodes = 7;
static const char* ElementOutputNames[] = {
	"centroid",
	"mass",
	"strain_energy",
	"kinetic_energy",
	"linear_momentum",
	"stress",
	"material_output"};

/* constructor */
SolidElementT::SolidElementT(const ElementSupportT& support):
	ContinuumElementT(support),
	fLocLastDisp(LocalArrayT::kLastDisp),
	fLocVel(LocalArrayT::kVel),
	fLocAcc(LocalArrayT::kAcc),
	fLocTemp(NULL),
	fLocTemp_last(NULL),
	fStoreInternalForce(false),
	fMassType(kAutomaticMass),
	fEigenvalueInc(0),
	fv_ss(0)
{
	SetName("solid_element");
}

/* destructor */
SolidElementT::~SolidElementT(void)
{
	delete fLocTemp;
	delete fLocTemp_last;
}

/* close current time increment. Called if the integration over the
 * current time increment was successful. */
void SolidElementT::CloseStep(void)
{
	/* inherited */
	ContinuumElementT::CloseStep();
	
	/* check eigenvalue */
	if (fEigenvalueInc > 0 && ElementSupport().StepNumber() % fEigenvalueInc == 0)
	{
		/* get estimated max eigenvalue */
		double max_eig = MaxEigenvalue();
		ElementSupport().Output() << ElementSupport().StepNumber() << ": max_eig = " << max_eig << '\n';
	}
}

/* solution calls */
void SolidElementT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
	/* not my field */
	if (&field != &(Field())) return;

	/* quick exit */
	bool hasnode = false;
	for (int i=0; i < fBlockData.Length() && !hasnode; i++)
		if (fConnectivities[i]->HasValue(node)) hasnode = true;
	if (!hasnode) return;

	/* set components and weights */
	double constMa = 0.0;
	double constKd = 0.0;
	
	/* components dicated by the algorithm */
	int formMa = fIntegrator->FormMa(constMa);
	int formKd = fIntegrator->FormKd(constKd);

	/* body forces */
	int formBody = 0;
	if (fMassType != kNoMass &&
	   (fBodySchedule && fBody.Magnitude() > kSmall))
	{	
		formBody = 1;
		if (!formMa) constMa = 1.0; /* override */
	}

	/* override controller */
	if (fMassType == kNoMass) formMa = 0;
	
	/* temp for nodal force */
	dArrayT nodalforce;
	
	bool axisymmetric = Axisymmetric();
	Top();
	while (NextElement())
	{
		int nodeposition;
		const iArrayT& nodes_u = CurrentElement().NodesU();
		if (nodes_u.HasValue(node, nodeposition))
		{
			/* initialize */
			fRHS = 0.0;

			/* global shape function values */
			SetGlobalShape();

			/* internal force contribution */	
			if (formKd) FormKd(constKd);

			/* inertia forces */
			if (formMa || formBody)
			{
				/* nodal accelerations */
				if (formMa)
					SetLocalU(fLocAcc);
				else 
					fLocAcc = 0.0;
			
				/* body force contribution */
				if (formBody) AddBodyForce(fLocAcc);
				
				/* calculate inertial forces */
				if (!fCurrMaterial->HasChangingDensity())
					FormMa(fMassType, constMa*fCurrMaterial->Density(), axisymmetric, &fLocAcc, NULL, NULL);
				else /* need to compute density */
				{
					/* collect densities */
					fShapes->TopIP();
					while (fShapes->NextIP())
						fDensity[fShapes->CurrIP()] = fCurrMaterial->Density();

					FormMa(fMassType, constMa, axisymmetric, &fLocAcc, NULL, fDensity.Pointer());
				}
			}

			/* loop over nodes (double-noding OK) */
			int dex = 0;
			for (int i = 0; i < nodes_u.Length(); i++)
			{
				if (nodes_u[i] == node)
				{
					/* components for node */
					nodalforce.Set(NumDOF(), fRHS.Pointer(dex));
	
					/* accumulate */
					force += nodalforce;
				}
				dex += NumDOF();
			}
		}
	}
}

void SolidElementT::AddLinearMomentum(dArrayT& momentum)
{
	/* check */
	if (momentum.Length() != NumDOF()) throw ExceptionT::kSizeMismatch;
		
	/* loop over elements */
	Top();
	while (NextElement())
	{
		/* global shape function derivatives, jacobians, local coords */
		SetGlobalShape();
		
		/* get velocities */
		SetLocalU(fLocVel);

		/* integration */
		const double* Det    = fShapes->IPDets();
		const double* Weight = fShapes->IPWeights();
	
		fShapes->TopIP();
		while ( fShapes->NextIP() )
		{					
			/* density could change with position */
			double density = fCurrMaterial->Density();

			double temp  = density*(*Det++)*(*Weight++);

			/* integration point velocities */
			fShapes->InterpolateU(fLocVel, fDOFvec);

			double* p    = momentum.Pointer();
			double* pvel = fDOFvec.Pointer();					
			for (int dof = 0; dof < NumDOF(); dof++)			
				*p++ += temp*(*pvel++);
		}
	}
}

/* returns the energy as defined by the derived class types */
double SolidElementT::InternalEnergy(void)
{
	double energy = 0.0;

	Top();
	while ( NextElement() )
	{
		/* global shape function derivatives, jacobians, local coords */
		SetGlobalShape();

		/* integration */
		const double* Det    = fShapes->IPDets();
		const double* Weight = fShapes->IPWeights();

		fShapes->TopIP();
		while ( fShapes->NextIP() )
			energy += fCurrMaterial->StrainEnergyDensity()*(*Det++)*(*Weight++);
	}
	
	return energy;
}

void SolidElementT::SendOutput(int kincode)
{
	/* output flags */
	iArrayT flags(fNodalOutputCodes.Length());

	/* set flags to get desired output */
	flags = IOBaseT::kAtNever;
	switch (kincode)
	{
		case iNodalDisp:
		    flags[iNodalDisp] = 1;
			break;
		case iNodalStress:
		    flags[iNodalStress] = 1;
			break;
		case iEnergyDensity:
		    flags[iEnergyDensity] = 1;
			break;
		case iPrincipal:
			flags[iPrincipal] = 1;
			break;
		case iMaterialData:
			flags[iMaterialData] = 1;
			break;
		default:
			cout << "\n SolidElementT::SendOutput: invalid output code: " 
			     << kincode << endl;
	}

	/* number of output values */
	iArrayT n_counts;
	SetNodalOutputCodes(IOBaseT::kAtInc, flags, n_counts);

	/* reset averaging workspace */
	ElementSupport().ResetAverage(n_counts.Sum());
      
	/* no element output */
	iArrayT e_counts(fElementOutputCodes.Length());
	e_counts = 0;

	/* generate output */
	dArray2DT n_values, e_values;
	ComputeOutput(n_counts, n_values, e_counts, e_values);
}

/* contribution to the nodal residual forces */
const dArray2DT& SolidElementT::InternalForce(int group)
{
	const char caller[] = "SolidElementT::InternalForce";

	/* check */
	if (group != Group())
		ExceptionT::GeneralFail(caller, "expecting solver group %d not %d", 
			Group(), group);
			
	/* must be storing force */
	if (!fStoreInternalForce)
		ExceptionT::GeneralFail(caller, "internal force not being stored");
			
	return fForce;
}

/* describe the parameters needed by the interface */
void SolidElementT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ContinuumElementT::DefineParameters(list);

	/* mass type */
	ParameterT mass_type(ParameterT::Enumeration, "mass_type");
	mass_type.AddEnumeration("automatic", kAutomaticMass);
	mass_type.AddEnumeration("no_mass", kNoMass);
    mass_type.AddEnumeration("consistent_mass", kConsistentMass);
    mass_type.AddEnumeration("lumped_mass", kLumpedMass);
    mass_type.SetDefault(fMassType);
	list.AddParameter(mass_type);
	
	/* eigenvalue estimation increment */
	ParameterT eig_inc(fEigenvalueInc, "eigenvalue_inc");
	eig_inc.SetDefault(fEigenvalueInc);
	eig_inc.AddLimit(0, LimitT::LowerInclusive);
	list.AddParameter(eig_inc);
}

/* information about subordinate parameter lists */
void SolidElementT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ContinuumElementT::DefineSubs(sub_list);

	/* nodal output codes (optional) */
	sub_list.AddSub("solid_element_nodal_output", ParameterListT::ZeroOrOnce);
	sub_list.AddSub("solid_element_element_output", ParameterListT::ZeroOrOnce);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SolidElementT::NewSub(const StringT& name) const
{
	if (name == "solid_element_nodal_output")
	{
		ParameterContainerT* node_output = new ParameterContainerT(name);
		
		/* wave speed sampling direction */
		ParameterContainerT wave_direction("wave_direction");
		wave_direction.SetListOrder(ParameterListT::Choice);
		wave_direction.AddSub("Vector_2");
		wave_direction.AddSub("Vector_3");
		node_output->AddSub(wave_direction, ParameterListT::ZeroOrOnce);

		/* steady-state speed for calculating Poynting vector */
		node_output->AddParameter(ParameterT::Double, "v1_ss", ParameterListT::ZeroOrOnce);
		
		/* all false by default */
		for (int i = 0; i < NumNodalOutputCodes; i++) {
			ParameterT output(ParameterT::Integer, NodalOutputNames[i]);
			output.SetDefault(1);
			node_output->AddParameter(output, ParameterListT::ZeroOrOnce);
		}

		return node_output;
	}
	else if (name == "solid_element_element_output")
	{
		ParameterContainerT* element_output = new ParameterContainerT(name);
		
		/* all false by default */
		for (int i = 0; i < NumElementOutputCodes; i++) {
			ParameterT output(ParameterT::Integer, ElementOutputNames[i]);
			output.SetDefault(1);
			element_output->AddParameter(output, ParameterListT::ZeroOrOnce);
		}

		return element_output;	
	}
	else /* inherited */
		return ContinuumElementT::NewSub(name);
}

/* accept parameter list */
void SolidElementT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SolidElementT::TakeParameterList";

	/* set mass type before calling ContinuumElementT::TakeParameterList because
	 * it's needed for SolidElementT::TangentType. If the mass type is kAutomaticMass,
	 * this needs to be resolved; however, it cannot be resolved before ContinuumElementT::TakeParameterList
	 * has been called because fIntegrator won't be set. Therefore, kAutomaticMass is resolved
	 * both here and in SolidElementT::TangentType. */
	int mass_type = list.GetParameter("mass_type");
	fMassType = int2MassTypeT(mass_type);

	/* inherited */
	ContinuumElementT::TakeParameterList(list);

	/* check order of the time integrator */
	int order = fIntegrator->Order();
	if (order != 0 && order != 2)
		ExceptionT::GeneralFail(caller, "expecting time integrator order 0 or 2 not %d", order);

	/* resolve mass type */
	if (fMassType == kAutomaticMass) {
		if (fIntegrator->ImplicitExplicit() == IntegratorT::kImplicit)
			fMassType = kConsistentMass;
		else
			fMassType = kLumpedMass;
	}
	
	/* eigenvalue output increment */
	fEigenvalueInc = list.GetParameter("eigenvalue_inc");

	/* allocate work space */
	fB.Dimension(dSymMatrixT::NumValues(NumSD()), NumSD()*NumElementNodes());
	fD.Dimension(dSymMatrixT::NumValues(NumSD()));
	
	/* nodal output codes */
	fNodalOutputCodes.Dimension(NumNodalOutputCodes);
	fNodalOutputCodes = IOBaseT::kAtNever;
	qUseSimo = qNoExtrap = false;	
	const ParameterListT* node_output = list.List("solid_element_nodal_output");
	if (node_output) 
	{
		/* set flags */
		for (int i = 0; i < NumNodalOutputCodes; i++)
		{
			/* look for entry */
			const ParameterT* nodal_value = node_output->Parameter(NodalOutputNames[i]);
			if (nodal_value) {
				int do_write = *nodal_value;

				/* Additional smoothing flags */
				if (!qUseSimo && do_write == 3) {
	    			qUseSimo = qNoExtrap = true;
	    			fNodalOutputCodes[i] = IOBaseT::kAtInc;
	    		}
	    		else if (!qNoExtrap && do_write == 2) {
	    			qNoExtrap = true;
	    			fNodalOutputCodes[i] = IOBaseT::kAtInc;
	    		}
	    		else if (do_write == 1)
	    			fNodalOutputCodes[i] = IOBaseT::kAtInc;
			}
		}
	
		/* wave speed sampling direction */
		if (fNodalOutputCodes[iWaveSpeeds])
		{
			/* resolve list choice */
			ParameterInterfaceT* n_output = NewSub(node_output->Name());
			const ParameterListT& vec_params = node_output->GetListChoice(*n_output, "wave_direction");
			delete n_output;
			
			/* extract unit vector */
			VectorParameterT::Extract(vec_params, fNormal);
			if (fNormal.Length() != NumSD()) ExceptionT::GeneralFail(caller, "expecting normal length %d not %d",
				NumSD(), fNormal.Length());
			fNormal.UnitVector();
		}
		
		/* steady state crack speed */
		if (fNodalOutputCodes[iPoyntingVector])
		{
			const ParameterT* v_ss = node_output->Parameter("v1_ss");
			if (!v_ss)
				ExceptionT::GeneralFail(caller, "Poynting vector requires \"v1_ss\"");
			else
				fv_ss = *v_ss;
		}
	}

	/* element output codes */
	fElementOutputCodes.Dimension(NumElementOutputCodes);
	fElementOutputCodes = IOBaseT::kAtNever;
	const ParameterListT* element_output = list.List("solid_element_element_output");
	if (element_output)
		for (int i = 0; i < NumElementOutputCodes; i++)
		{
			/* look for entry */
			const ParameterT* element_value = element_output->Parameter(ElementOutputNames[i]);
			if (element_value) {
				int do_write = *element_value;
				if (do_write == 1)
					fElementOutputCodes[i] = IOBaseT::kAtInc;
			}
		}

	/* generate list of material needs */
	bool changing_density = false;
	bool constant_density = false;
	fMaterialNeeds.Dimension(fMaterialList->Length());
	for (int i = 0; i < fMaterialNeeds.Length(); i++)
	{
		/* allocate */
		ArrayT<bool>& needs = fMaterialNeeds[i];
		needs.Dimension(3);

		/* casts are safe since class contructs materials list */
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[i];
		SolidMaterialT* mat = (SolidMaterialT*) pcont_mat;

		/* collect needs */
		needs[kNeedDisp] = mat->NeedDisp();
		needs[kNeedVel] = mat->NeedVel();
		needs[KNeedLastDisp] = mat->NeedLastDisp();
		
		/* changing density */
		if (mat->HasChangingDensity()) 
			changing_density = true;
		else /* constant density */
			constant_density = true;
	}
	
	/* all must be the same */
	if (changing_density == constant_density)
		ExceptionT::GeneralFail(caller, "cannot mix materials with constant/changing density");

	/* work space for calculating integration point densities */
	if (changing_density)
		fDensity.Dimension(NumIP());
}

/* estimate the largest eigenvalue */
double SolidElementT::MaxEigenvalue(void)
{
	/* check */
	if (fMassType != kLumpedMass)
		ExceptionT::GeneralFail("SolidElementT::MaxEigenvalue", "mass matrix must be lumped");

	/* set up K */
	CCSMatrixT K(ElementSupport().Output(), CCSMatrixT::kNoCheck, ElementSupport().Communicator());
	
	/* collection equation numbers */
	AutoArrayT<const iArray2DT*> eq_1;
	AutoArrayT<const RaggedArray2DT<int>*> eq_2;
	Equations(eq_1, eq_2);
	for (int i = 0; i < eq_1.Length(); i++)
		K.AddEquationSet(*eq_1[i]);
	for (int i = 0; i < eq_2.Length(); i++)
		K.AddEquationSet(*eq_2[i]);
	int num_eq = ElementSupport().NodeManager().NumEquations(Group());
	K.Initialize(num_eq, num_eq, 1);
	K.Clear();
	
	/* mass matrix */
	dArrayT M(num_eq);
	M = 0.0;

	/* compute K, M */
	bool axisymmetric = Axisymmetric();
	Top();
	while (NextElement())
		if (CurrentElement().Flag() != ElementCardT::kOFF)
		{
			/* element equations */
			const iArrayT& eqs = CurrentElement().Equations();
		
			/* set shape function derivatives */
			SetGlobalShape();

			/* assemble element stiffness */
			fLHS = 0.0;
			FormStiffness(1.0);	
			K.Assemble(fLHS, eqs);

			/* assemble element mass */
			fLHS = 0.0;
			if (!fCurrMaterial->HasChangingDensity())
				FormMass(fMassType, fCurrMaterial->Density(), axisymmetric, NULL);
			else
			{
				/* collect densities */
				fShapes->TopIP();
				while (fShapes->NextIP())
					fDensity[fShapes->CurrIP()] = fCurrMaterial->Density();			
				FormMass(fMassType, 1.0, axisymmetric, fDensity.Pointer());
			}
			for (int i = 0; i < eqs.Length(); i++)
			{
				int eq = eqs[i];
				if (eq > 0) /* active equation */
					M[eq-1] += fLHS(i,i);
			}
		}

	/* estimate largest eigenvalue */
	double max_eig = 0.0;
	for (int i = 0; i < M.Length(); i++)
		if (M[i] > kSmall) {
			double abs_row_sum = K.AbsRowSum(i);
			double eig = abs_row_sum/M[i];
			max_eig = (eig > max_eig) ? eig : max_eig;
		}

	return max_eig;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct output labels array */
void SolidElementT::SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
	/* initialize */
	counts.Dimension(flags.Length());
	counts = 0;

	/* set output flags */
	if (flags[iNodalCoord] == mode)
		counts[iNodalCoord] = NumSD();
	if (flags[iNodalDisp] == mode)
		counts[iNodalDisp] = NumDOF();
	if (flags[iNodalStress] == mode)
		counts[iNodalStress] = fB.Rows();
	if (flags[iPrincipal] == mode)
		counts[iPrincipal] = NumSD();
	if (flags[iEnergyDensity] == mode)
		counts[iEnergyDensity] = 1;
	if (flags[iWaveSpeeds] == mode)
		counts[iWaveSpeeds] = NumSD();
	if (flags[iMaterialData] == mode)
		counts[iMaterialData] = (*fMaterialList)[0]->NumOutputVariables();
	if (flags[iPoyntingVector] == mode)
		counts[iPoyntingVector] = NumSD();
}

void SolidElementT::SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
	/* initialize */
	counts.Dimension(flags.Length());
	counts = 0;

	/* set output flags */
	if (fElementOutputCodes[iCentroid] == mode) counts[iCentroid] = NumSD();
	if (fElementOutputCodes[iMass] == mode) counts[iMass] = NumIP();
	if (fElementOutputCodes[iStrainEnergy] == mode) counts[iStrainEnergy] = 1;
	if (fElementOutputCodes[iKineticEnergy] == mode) counts[iKineticEnergy] = 1;
	if (fElementOutputCodes[iLinearMomentum] == mode) counts[iLinearMomentum] = NumDOF();
	if (fElementOutputCodes[iIPStress] == mode) counts[iIPStress] = 2*fB.Rows()*NumIP();
	if (fElementOutputCodes[iIPMaterialData] == mode) 
		counts[iIPMaterialData] = (*fMaterialList)[0]->NumOutputVariables()*NumIP();
}

/* initialize local arrays */
void SolidElementT::SetLocalArrays(void)
{
	/* inherited */
	ContinuumElementT::SetLocalArrays();

	/* allocate */
	int nen = NumElementNodes();
	fLocLastDisp.Dimension(nen, NumDOF());
	fLocAcc.Dimension(nen, NumDOF());
	fLocVel.Dimension(nen, NumDOF());

	/* register */
	Field().RegisterLocal(fLocLastDisp);
	if (fIntegrator->Order() == 2)
	{
		Field().RegisterLocal(fLocVel);
		Field().RegisterLocal(fLocAcc);
	}

	/* look for a temperature field */
	const FieldT* temperature = ElementSupport().Field("temperature");
	if (temperature) {
	
		/* construct */
		fLocTemp = new LocalArrayT(LocalArrayT::kDisp, nen, temperature->NumDOF());
		fLocTemp_last = new LocalArrayT(LocalArrayT::kLastDisp, nen, temperature->NumDOF());

		/* register */
		temperature->RegisterLocal(*fLocTemp);
		temperature->RegisterLocal(*fLocTemp_last);
	}
}

/* set the correct shape functions */
void SolidElementT::SetShape(void)
{
	/* construct shape functions */
	fShapes = new ShapeFunctionT(GeometryCode(), NumIP(), fLocInitCoords);
	if (!fShapes) throw ExceptionT::kOutOfMemory;

	/* initialize */
	fShapes->Initialize();
}

/* form shape functions and derivatives */
void SolidElementT::SetGlobalShape(void)
{
	/* inherited */
	ContinuumElementT::SetGlobalShape();

	/* material needs */
	const ArrayT<bool>& needs = fMaterialNeeds[CurrentElement().MaterialNumber()];

	/* material dependent local arrays */
	if (needs[kNeedDisp])     SetLocalU(fLocDisp);	
	if (needs[KNeedLastDisp]) SetLocalU(fLocLastDisp);	
	if (needs[kNeedVel])
	{
		/* have velocity */
		if (fLocVel.IsRegistered())
			SetLocalU(fLocVel);
		else /* finite difference approximation */
		{
			double onebydt = 1.0/ElementSupport().TimeStep();
			fLocVel.SetToCombination(onebydt, fLocDisp, -onebydt, fLocLastDisp);
		}
	}
	
	/* get nodal temperatures if available */
	if (fLocTemp)SetLocalU(*fLocTemp);
	if (fLocTemp_last) SetLocalU(*fLocTemp_last);
}

/* construct a new material support and return a pointer */
MaterialSupportT* SolidElementT::NewMaterialSupport(MaterialSupportT* p) const
{
	/* allocate new */
	if (!p) p = new SolidMatSupportT(NumDOF(), NumIP());

	/* inherited initializations */
	ContinuumElementT::NewMaterialSupport(p);

	/* set SolidMatSupportT fields */
	SolidMatSupportT* ps = TB_DYNAMIC_CAST(SolidMatSupportT*, p);
	if (ps) {
		/* set pointers to local arrays */
		ps->SetLocalArray(fLocLastDisp);
		ps->SetLocalArray(fLocVel);
		ps->SetLocalArray(fLocAcc);

		/* temperatures if available */
		if (fLocTemp) ps->SetTemperatures(*fLocTemp);

		if (fLocTemp_last) ps->SetLastTemperatures(*fLocTemp_last);
	}

	return p;
}

/* set the \e B matrix using the given shape function derivatives */
void SolidElementT::Set_B(const dArray2DT& DNa, dMatrixT& B) const
{
#if __option(extended_errorcheck)
	if (B.Rows() != dSymMatrixT::NumValues(DNa.MajorDim()) ||
	    B.Cols() != DNa.Length())
	    throw ExceptionT::kSizeMismatch;
#endif

	int nnd = DNa.MinorDim();
	double* pB = B.Pointer();

	/* 1D */
	if (DNa.MajorDim() == 1)
	{
		const double* pNax = DNa(0);
		for (int i = 0; i < nnd; i++)
			*pB++ = *pNax++;
	}
	/* 2D */
	else if (DNa.MajorDim() == 2)
	{
		const double* pNax = DNa(0);
		const double* pNay = DNa(1);
		for (int i = 0; i < nnd; i++)
		{
			/* see Hughes (2.8.20) */
			*pB++ = *pNax;
			*pB++ = 0.0;
			*pB++ = *pNay;

			*pB++ = 0.0;
			*pB++ = *pNay++;
			*pB++ = *pNax++;
		}
	}
	/* 3D */
	else		
	{
		const double* pNax = DNa(0);
		const double* pNay = DNa(1);
		const double* pNaz = DNa(2);
		for (int i = 0; i < nnd; i++)
		{
			/* see Hughes (2.8.21) */
			*pB++ = *pNax;
			*pB++ = 0.0;
			*pB++ = 0.0;
			*pB++ = 0.0;
			*pB++ = *pNaz;
			*pB++ = *pNay;

			*pB++ = 0.0;
			*pB++ = *pNay;
			*pB++ = 0.0;
			*pB++ = *pNaz;
			*pB++ = 0.0;
			*pB++ = *pNax;

			*pB++ = 0.0;
			*pB++ = 0.0;
			*pB++ = *pNaz++;
			*pB++ = *pNay++;
			*pB++ = *pNax++;
			*pB++ = 0.0;
		}
	}
}

/* set the \e B matrix for axisymmetric deformations */
void SolidElementT::Set_B_axi(const dArrayT& Na, const dArray2DT& DNa, 
	double r, dMatrixT& B) const
{
#if __option(extended_errorcheck)
	if (B.Rows() != 4 || /* (number of stress 2D components) + 1 */
	    B.Cols() != DNa.Length() ||
    DNa.MajorDim() != 2 ||
    DNa.MinorDim() != Na.Length())
			ExceptionT::SizeMismatch("SolidElementT::Set_B_axi");
#endif

	int nnd = DNa.MinorDim();
	double* pB = B.Pointer();

	const double* pNax = DNa(0);
	const double* pNay = DNa(1);
	const double* pNa  = Na.Pointer();
	for (int i = 0; i < nnd; i++)
	{
		/* see Hughes (2.12.8) */
		*pB++ = *pNax;
		*pB++ = 0.0;
		*pB++ = *pNay;
		*pB++ = *pNa++/r; /* about y-axis: u_r = u_x */

		*pB++ = 0.0;
		*pB++ = *pNay++;
		*pB++ = *pNax++;
		*pB++ = 0.0;
	}
}

/* set B-bar as given by Hughes (4.5.11-16) */
void SolidElementT::Set_B_bar(const dArray2DT& DNa, const dArray2DT& mean_gradient, 
	dMatrixT& B) const
{
#if __option(extended_errorcheck)
	if (B.Rows() != dSymMatrixT::NumValues(DNa.MajorDim()) ||
	    B.Cols() != DNa.Length() ||
	    mean_gradient.MinorDim() != DNa.MinorDim() ||
	    mean_gradient.MajorDim() != DNa.MajorDim())
	    throw ExceptionT::kSizeMismatch;
#endif

	int nnd = DNa.MinorDim();
	double* pB = B.Pointer();

	/* 1D */
	if (DNa.MajorDim() == 1)
	{
		cout << "\n SolidElementT::Set_B_bar: not implemented yet for 1D B-bar" << endl;
		throw ExceptionT::kGeneralFail;
	}
	/* 2D */
	else if (DNa.MajorDim() == 2)
	{
		const double* pNax = DNa(0);
		const double* pNay = DNa(1);
			
		const double* pBmx = mean_gradient(0);
		const double* pBmy = mean_gradient(1);
			
		for (int i = 0; i < nnd; i++)
		{
			double factx = ((*pBmx++) - (*pNax))/3.0;
			double facty = ((*pBmy++) - (*pNay))/3.0;
			
			/* Hughes (4.5.11-16) */
			*pB++ = *pNax + factx;
			*pB++ = factx;
			*pB++ = *pNay;
	
			*pB++ = facty;
			*pB++ = *pNay + facty;
			*pB++ = *pNax;
				
			pNax++; pNay++;
		}
	}
	/* 3D */
	else		
	{
		const double* pNax = DNa(0);
		const double* pNay = DNa(1);
		const double* pNaz = DNa(2);

		const double* pBmx = mean_gradient(0);
		const double* pBmy = mean_gradient(1);
		const double* pBmz = mean_gradient(2);
			
		for (int i = 0; i < nnd; i++)
		{
			double factx = ((*pBmx++) - (*pNax))/3.0;
			double facty = ((*pBmy++) - (*pNay))/3.0;
			double factz = ((*pBmz++) - (*pNaz))/3.0;

			/* Hughes (4.5.11-16) */
			*pB++ = *pNax + factx;
			*pB++ = factx;
			*pB++ = factx;
			*pB++ = 0.0;
			*pB++ = *pNaz;
			*pB++ = *pNay;

			*pB++ = facty;
			*pB++ = *pNay + facty;
			*pB++ = facty;
			*pB++ = *pNaz;
			*pB++ = 0.0;
			*pB++ = *pNax;
	
			*pB++ = factz;
			*pB++ = factz;
			*pB++ = *pNaz + factz;
			*pB++ = *pNay;
			*pB++ = *pNax;
			*pB++ = 0.0;
				
			pNax++; pNay++; pNaz++;
		}
	}
}

void SolidElementT::Set_B_bar_axi(const dArrayT& Na, const dArray2DT& DNa, const dArray2DT& mean_gradient, 
	double r, dMatrixT& B) const
{
	const char caller[] = "SolidElementT::Set_B_bar_axi";
	
#if __option(extended_errorcheck)
	if (B.Rows() != 4 || /* (number of stress 2D components) + 1 */
	    B.Cols() != DNa.Length() ||
	    DNa.MajorDim() != 2 ||
	    DNa.MinorDim() != Na.Length() ||
	    mean_gradient.MinorDim() != DNa.MinorDim() ||
	    mean_gradient.MajorDim() != DNa.MajorDim())
		ExceptionT::SizeMismatch("SolidElementT::Set_B_axi");
#endif

	int nnd = DNa.MinorDim();
	double* pB = B.Pointer();

	const double* pNax = DNa(0);
	const double* pNay = DNa(1);			
	const double* pNa  = Na.Pointer();

	const double* pBmx = mean_gradient(0);
	const double* pBmy = mean_gradient(1);			
	for (int i = 0; i < nnd; i++)
	{
		/* exchange volumetric part: b_bar_vol - b_vol */
		double factx = ((*pBmx++) - (*pNax + *pNa/r))/3.0;
		double facty = ((*pBmy++) - (*pNay))/3.0;
			
		*pB++ = *pNax + factx;
		*pB++ = factx;
		*pB++ = *pNay; /* shear */
		*pB++ = *pNa/r + factx; /* about y-axis: u_r = u_x */
	
		*pB++ = facty;
		*pB++ = *pNay + facty;
		*pB++ = *pNax; /* shear */
		*pB++ = facty;
				
		pNax++; pNay++; pNa++;
	}
}

/* construct the effective mass matrix */
void SolidElementT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
	/* inherited */
	ContinuumElementT::LHSDriver(sys_type);

	/* element contribution */
	ElementLHSDriver();
}

void SolidElementT::ElementLHSDriver(void)
{
	/* set components and weights */
	double constM = 0.0;
	double constK = 0.0;
	
	int formM = fIntegrator->FormM(constM);
	int formK = fIntegrator->FormK(constK);

	/* override algorithm */
	if (fMassType == kNoMass) formM = 0;

	/* quick exit */
	if ((formM == 0 && formK == 0) ||
	    (fabs(constM) < kSmall &&
	     fabs(constK) < kSmall)) return;

	/* loop over elements */
	bool axisymmetric = Axisymmetric();
	Top();
	while (NextElement())
		if (CurrentElement().Flag() != ElementCardT::kOFF)
		{
			double constKe = constK;
			double constMe = constM;
	
			/* initialize */
			fLHS = 0.0;
		
			/* set shape function derivatives */
			SetGlobalShape();

			/* element mass */
			if (fabs(constMe) > kSmall) {
				if (!fCurrMaterial->HasChangingDensity())
					FormMass(fMassType, constMe*(fCurrMaterial->Density()), axisymmetric, NULL);
				else
				{
					/* collect densities */
					fShapes->TopIP();
					while (fShapes->NextIP())
						fDensity[fShapes->CurrIP()] = fCurrMaterial->Density();
				
					FormMass(fMassType, constMe, axisymmetric, fDensity.Pointer());
				}
			}

			/* element stiffness */
			if (fabs(constKe) > kSmall)
				FormStiffness(constKe);

			/* add to global equations */
			AssembleLHS();
		}
}

void SolidElementT::RHSDriver(void)
{
	/* inherited */
	ContinuumElementT::RHSDriver();

	/* element contribution */
	ElementRHSDriver();
}

/* form the residual force vector */
void SolidElementT::ElementRHSDriver(void)
{
	/* heat source if needed */
	const FieldT* temperature = ElementSupport().Field("temperature");

	/* initialize sources */
	if (temperature && fIncrementalHeat.Length() == 0) {
	
		/* allocate the element heat */
		fElementHeat.Dimension(fShapes->NumIP());
			
		/* initialize heat source arrays */
		fIncrementalHeat.Dimension(fBlockData.Length());
		for (int i = 0; i < fIncrementalHeat.Length(); i++)
		{
			/* dimension */
			fIncrementalHeat[i].Dimension(fBlockData[i].Dimension(), NumIP());

			/* register */
			temperature->RegisterSource(fBlockData[i].ID(), fIncrementalHeat[i]);
		}
	}
	
	/* storage for entire internal force */
	if (fStoreInternalForce) {
		fForce.Dimension(ElementSupport().NumNodes(), NumDOF());
		fForce = 0.0;	
	} else 
		fForce.Dimension(0, NumDOF());

	/* set components and weights */
	double constMa = 0.0;
	double constKd = 0.0;
	
	/* components dicated by the algorithm */
	int formMa = fIntegrator->FormMa(constMa);
	int formKd = fIntegrator->FormKd(constKd);

	/* body forces */
	int formBody = 0;
	if (fMassType != kNoMass &&
	   (fBodySchedule && fBody.Magnitude() > kSmall))
	{	
		formBody = 1;
		if (!formMa) constMa = 1.0; /* override */
	}

	/* override controller */
	if (fMassType == kNoMass) formMa = 0;

	bool axisymmetric = Axisymmetric();
	int block_count = 0, block_dex = 0;
	Top();
	while (NextElement())
	{
		/* advance to block [skip empty blocks] */		
		while (block_count == fBlockData[block_dex].Dimension()) {
			block_count = 0;
			block_dex++;
		}

		/* current element */
		const ElementCardT& element = CurrentElement();

		if (element.Flag() != ElementCardT::kOFF)
		{
			/* initialize */
			fRHS = 0.0;
			fElementHeat = 0.0;
		
			/* global shape function values */
			SetGlobalShape();
		
			/* internal force contribution */	
			if (formKd) FormKd(-constKd);
				
			/* inertia forces */
			if (formMa || formBody)
			{
				/* nodal accelerations */
				if (formMa)
					SetLocalU(fLocAcc);
				else 
					fLocAcc = 0.0;
			
				/* body force contribution */
				if (formBody) AddBodyForce(fLocAcc);
		
				if (!fCurrMaterial->HasChangingDensity())
					FormMa(fMassType, -constMa*fCurrMaterial->Density(), axisymmetric, &fLocAcc, NULL, NULL);
				else
				{
					/* collect densities */
					fShapes->TopIP();
					while (fShapes->NextIP())
						fDensity[fShapes->CurrIP()] = fCurrMaterial->Density();
				
					FormMa(fMassType, -constMa, axisymmetric, &fLocAcc, NULL, fDensity.Pointer());
				}
			}
		
			/* store incremental heat */
			if (temperature)
				fIncrementalHeat[block_dex].SetRow(block_count, fElementHeat);

			/* global assembly */
			if (fStoreInternalForce) 
			{
				double* pRHS = fRHS.Pointer();
				int ndof = NumDOF();
				const iArrayT& nodes_u = CurrentElement().NodesU();
				for (int i = 0; i < nodes_u.Length(); i++) 
				{
					double* pForce = fForce(nodes_u[i]);
					for (int j = 0; j < ndof; j++)
						*pForce++ += *pRHS++;
				}
			}
			else /* assemble element-by-element */
				AssembleRHS();
		}
		
		/* next block */
		block_count++;
	}
	
	/* assemble all at once */
	if (fStoreInternalForce)
		ElementSupport().AssembleRHS(Group(), fForce, Field().Equations());
}

/* current element operations */
bool SolidElementT::NextElement(void)
{
	/* inherited */
	bool result = ContinuumElementT::NextElement();
	
	/* get material pointer */
	if (result)
	{
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[CurrentElement().MaterialNumber()];
	
		/* cast is safe since class contructs materials list */
		fCurrMaterial = (SolidMaterialT*) pcont_mat;
	}
	
	return result;
}

/* form of tangent matrix */
GlobalT::SystemTypeT SolidElementT::TangentType(void) const
{
	/* resolve mass type - see note in SolidElementT::TakeParameterList
	 * regarding when the mass type is set */
	MassTypeT mass_type = fMassType;
	if (mass_type == kAutomaticMass) {
		if (fIntegrator->ImplicitExplicit() == IntegratorT::kImplicit)
			mass_type = kConsistentMass;
		else
			mass_type = kLumpedMass;
	}	

	/* special case */
	if (fIntegrator->Order() > 0 &&
	    fIntegrator->ImplicitExplicit() ==  eIntegratorT::kExplicit &&
	    (mass_type == kNoMass ||
	     mass_type == kLumpedMass))
		return GlobalT::kDiagonal;
	else
		/* inherited */
		return ContinuumElementT::TangentType();
}

/* return the materials list */
const SolidMatListT& SolidElementT::StructuralMaterialList(void) const
{
	return TB_DYNAMIC_CAST(const SolidMatListT&, MaterialsList());
}

/* extrapolate the integration point stresses and strains and extrapolate */
void SolidElementT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	const iArrayT& e_codes, dArray2DT& e_values)
{
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
	dMatrixT ippvector(nsd);

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
				const dSymMatrixT& stress = fCurrMaterial->s_ij();
				cauchy.Translate(stress);

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
//        cauchy.Eigensystem(ipprincipal,ippvector,true);
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
							cauchy.Multx(ip_velocity, ipPoynting, 1.0, dMatrixT::kAccumulate);
						}
						
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

#if 0
				/* nodal reactions */
				if (n_codes[iReaction])
				{
					AddNodalForce(FieldT(), node, dArrayT& force)
				}
#endif


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
//		iArrayT nodes_used(tmpDim);
		dArray2DT tmp_simo(tmpDim, n_simo);
		for (int i = 0; i < simo_force.MajorDim(); i++)
			if (simo_counts[i] > 0) {
//				nodes_used[rowNum] = i;
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
}

/***********************************************************************
* Private
***********************************************************************/

/* construct output labels array */
void SolidElementT::GenerateOutputLabels(const iArrayT& n_codes, ArrayT<StringT>& n_labels, 
	const iArrayT& e_codes, ArrayT<StringT>& e_labels) const
{
	const char caller[] = "SolidElementT::GenerateOutputLabels";

	/* allocate */
	n_labels.Dimension(n_codes.Sum());

	int count = 0;
	if (n_codes[iNodalDisp])
	{
		/* labels from the field */
		const ArrayT<StringT>& labels = Field().Labels();
		for (int i = 0; i < labels.Length(); i++)
			n_labels[count++] = labels[i];
	}

	if (n_codes[iNodalCoord])
	{
		const char* xlabels[] = {"x1", "x2", "x3"};
		for (int i = 0; i < NumSD(); i++)
			n_labels[count++] = xlabels[i];
	}

	if (n_codes[iNodalStress])
	{
		const char* slabels1D[] = {"s11"};
		const char* slabels2D[] = {"s11", "s22", "s12"};
		const char* slabels2D_axi[] = {"srr", "szz", "srz", "stt"};
		const char* slabels3D[] = {"s11", "s22", "s33", "s23", "s13", "s12"};
		int nstrs = fB.Rows();
		const char** slabels = NULL;
		if (nstrs == 1)
			slabels = slabels1D;
		else if (nstrs == 3)
			slabels = slabels2D;
		else if (nstrs == 4)
			slabels = slabels2D_axi;
		else if (nstrs == 6)
			slabels = slabels3D;
		else
			ExceptionT::GeneralFail(caller);

		for (int i = 0; i < nstrs; i++)
			n_labels[count++] = slabels[i];
	}
		
	if (n_codes[iPrincipal])
	{
		const char* plabels[] = {"s1", "s2", "s3"};
		for (int i = 0; i < NumSD(); i++)
			n_labels[count++] = plabels[i];
	}
		
	if (n_codes[iEnergyDensity]) n_labels[count++] = "phi";
	if (n_codes[iWaveSpeeds])
	{
		const char* clabels2D[] = {"cd", "cs"};
		const char* clabels3D[] = {"cd", "cs_min", "cs_max"};
		const char**    clabels = (NumSD() == 2) ? clabels2D : clabels3D;
		for (int i = 0; i < NumSD(); i++)
			n_labels[count++] = clabels[i];		
	}

	/* material output labels */
	if (n_codes[iMaterialData]) {
		ArrayT<StringT> matlabels;
		(*fMaterialList)[0]->OutputLabels(matlabels);	
		
		for (int i = 0; i < matlabels.Length(); i++)
			n_labels[count++] = matlabels[i];
	}
	
	/* Poynting vector */
	if (n_codes[iPoyntingVector]) {
		const char* fp_labels[] = {"Fp_X", "Fp_Y", "Fp_Z"};
		for (int i = 0; i < NumSD(); i++)
			n_labels[count++] = fp_labels[i];		
	}

	/* allocate */
	e_labels.Dimension(e_codes.Sum());
	count = 0;
	if (e_codes[iCentroid])
	{
		const char* xlabels[] = {"xc_1", "xc_2", "xc_3"};
		for (int i = 0; i < NumSD(); i++)
			e_labels[count++] = xlabels[i];
	}
	if (e_codes[iMass])
	{
		/* over integration points */
		for (int j = 0; j < NumIP(); j++)
		{
			StringT ip_label;
			ip_label.Append("ip", j+1);
			ip_label.Append(".mass");	
			e_labels[count++] = ip_label;
		}
	}
	if (e_codes[iStrainEnergy]) e_labels[count++] = "U";
	if (e_codes[iKineticEnergy]) e_labels[count++] = "T";
	if (e_codes[iLinearMomentum])
	{
		const char* plabels[] = {"L_X", "L_Y", "L_Z"};
		for (int i = 0; i < NumDOF(); i++)
			e_labels[count++] = plabels[i];
	}
	if (e_codes[iIPStress])
	{
		const char* slabels1D[] = {"s11", "e11"};
		const char* slabels2D[] = {"s11", "s22", "s12","e11", "e22", "e12"};
		const char* slabels2D_axi[] = {"srr", "szz", "srz", "stt", "err", "ezz", "erz", "ett"};
		const char* slabels3D[] = {"s11", "s22", "s33", "s23", "s13", "s12", "e11", "e22", "e33", "e23", "e13", "e12"};
		int nstrs = fB.Rows();
		const char** slabels = NULL;
		if (nstrs == 1)
			slabels = slabels1D;
		else if (nstrs == 3)
			slabels = slabels2D;
		else if (nstrs == 4)
			slabels = slabels2D_axi;
		else if (nstrs == 6)
			slabels = slabels3D;
		else
			ExceptionT::GeneralFail(caller);

		/* over integration points */
		for (int j = 0; j < NumIP(); j++)
		{
			StringT ip_label;
			ip_label.Append("ip", j+1);
			
			/* over stress/strain components */
			for (int i = 0; i < 2*nstrs; i++)
			{
				e_labels[count].Clear();
				e_labels[count].Append(ip_label, ".", slabels[i]);
				count++;
			}
		}		
	}

	/* material output labels */
	if (e_codes[iIPMaterialData])
	{
		ArrayT<StringT> matlabels;
		(*fMaterialList)[0]->OutputLabels(matlabels);	

		/* over integration points */
		for (int j = 0; j < NumIP(); j++)
		{
			StringT ip_label;
			ip_label.Append("ip", j+1);
			
			/* over stress components */
			for (int i = 0; i < matlabels.Length(); i++)
			{
				e_labels[count].Clear();
				e_labels[count].Append(ip_label, ".", matlabels[i]);
				count++;
			}
		}		
	}
}
