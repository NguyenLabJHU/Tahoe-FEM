/* $Id: SolidElementT.cpp,v 1.23 2002-05-14 23:05:56 cjkimme Exp $ */
/* created: paklein (05/28/1996) */

#include "SolidElementT.h"

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include "Constants.h"

#include "fstreamT.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "ElementCardT.h"
#include "ShapeFunctionT.h"
#include "eControllerT.h"
#include "StructuralMaterialT.h"
#include "MaterialListT.h"
#include "iAutoArrayT.h"

/* materials lists */
#include "MaterialList2DT.h"
#include "MaterialList3DT.h"

/* exception codes */
#include "ExceptionCodes.h"

/* initialize static data */
const int SolidElementT::NumNodalOutputCodes = 7;
const int SolidElementT::NumElementOutputCodes = 7;

/* constructor */
SolidElementT::SolidElementT(FEManagerT& fe_manager):
	ContinuumElementT(fe_manager),
	fLocLastDisp(LocalArrayT::kLastDisp),
	fLocVel(LocalArrayT::kVel),
	fLocAcc(LocalArrayT::kAcc),
	fStress(fNumSD),
	fD(dSymMatrixT::NumValues(fNumSD))
{
	/* check base class initializations */
	if (fNumDOF != fNumSD) throw eGeneralFail;

	ifstreamT& in  = fFEManager.Input();
	ostream&    out = fFEManager.Output();
	
	/* control parameters */
	in >> fMassType;		
	in >> fStrainDispOpt;

	/* checks */
	if (fMassType < kNoMass ||
	    fMassType > kLumpedMass) throw eBadInputValue;
	
	if (fStrainDispOpt != ShapeFunctionT::kStandardB &&
	    fStrainDispOpt != ShapeFunctionT::kMeanDilBbar) throw eBadInputValue;
}

/* data initialization */
void SolidElementT::Initialize(void)
{
	/* inherited */
	ContinuumElementT::Initialize();

	/* allocate strain-displacement matrix */
	fB.Allocate(dSymMatrixT::NumValues(fNumSD), fNumSD*fNumElemNodes);

	/* override */
	if (fController->ImplicitExplicit() == eControllerT::kExplicit &&
	    fMassType != kNoMass &&
	    fMassType != kLumpedMass) fMassType = kLumpedMass;

	/* setup for material output */
	if (fNodalOutputCodes[iMaterialData] || fElementOutputCodes[iIPMaterialData])
	{
		/* check compatibility of output */
		if (!CheckMaterialOutput())
		{
			cout << "\n SolidElementT::Initialize: error with material output" << endl;
			throw eBadInputValue;
		}
		/* no material output variables */
		else if ((*fMaterialList)[0]->NumOutputVariables() == 0)
		{
			cout << "\n SolidElementT::Initialize: there are no material outputs" << endl;
			fNodalOutputCodes[iMaterialData] = IOBaseT::kAtNever;
			fElementOutputCodes[iIPMaterialData] = IOBaseT::kAtNever;
		}
	}	
}

/* set the controller */
void SolidElementT::SetController(eControllerT* controller)
{
	/* inherited */
	ContinuumElementT::SetController(controller);
	
	/* check consistency */
	int is_explicit = (fController->ImplicitExplicit() == eControllerT::kExplicit);
	int is_dynamic  = (fController->Order() > 0);

	/* warning about no mass */
	if (is_dynamic && is_explicit)
	{
		if (fMassType == kNoMass)
		{
			cout << "\n SolidElementT::SetController: WARNING: element group ";
			cout << fFEManager.ElementGroupNumber(this) + 1 << " has mass type " << kNoMass << endl;
		}
		else if (fMassType != kLumpedMass) /* override all others */
//{
//cout << "\n SolidElementT::SetController: SKIPPING lumped mass override for explicit dynamics" << endl;		
		fMassType = kLumpedMass;
//}
	}
}

/* solution calls */
void SolidElementT::AddNodalForce(int node, dArrayT& force)
{
	/* quick exit */
        bool hasnode = false;
	for (int i=0; i < fBlockData.Length() && !hasnode; i++)
	  if (fConnectivities[i]->HasValue(node)) hasnode = true;
	if (!hasnode) return;

	/* set components and weights */
	double constMa = 0.0;
	double constCv = 0.0;
	double constKd = 0.0;
	
	/* components dicated by the algorithm */
	int formMa = fController->FormMa(constMa);
	int formCv = fController->FormCv(constCv);
	int formKd = fController->FormKd(constKd);

	/* body forces */
	int formBody = 0;
	if (fMassType != kNoMass &&
	   (fBodyForceLTf > -1 && fBody.Magnitude() > kSmall))
	{	
		formBody = 1;
		if (!formMa) constMa = 1.0; /* override */
	}

	/* override controller */
	if (fMassType == kNoMass) formMa = 0;
	
	/* temp for nodal force */
	dArrayT nodalforce;
	
	Top();
	while (NextElement())
	{
		int nodeposition;
		if (CurrentElement().NodesU().HasValue(node, nodeposition))
		{
			/* initialize */
			fRHS = 0.0;
	
			/* effective accelerations and displacements */
			//ComputeEffectiveDVA(formBody, formMa, constMa, formCv, constCv, formKd, constKd);
			//DEV - Rayleigh damping is poorly formulated

			/* global shape function values */
			SetGlobalShape();
	
			/* internal force contribution */	
			if (formKd) FormKd(1.0);
				
			/* damping */
			//if (formCv) FormCv(1.0);
			//DEV - computed at constitutive level

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
				FormMa(fMassType, fCurrMaterial->Density(), fLocAcc);
			}

			/* components for node */
			nodalforce.Set(fNumDOF, &fRHS[fNumDOF*nodeposition]);
	
			/* accumulate */
			force += nodalforce;
		}
	}
}

void SolidElementT::AddLinearMomentum(dArrayT& momentum)
{
	/* check */
	if (momentum.Length() != fNumDOF) throw eSizeMismatch;
		
	/* loop over elements */
	Top();
	while (NextElement())
	{
		/* global shape function derivatives, jacobians, local coords */
		SetGlobalShape();
		
		/* get velocities */
		SetLocalU(fLocVel);

		/* material density */
		double density = fCurrMaterial->Density();

		/* integration */
		const double* Det    = fShapes->IPDets();
		const double* Weight = fShapes->IPWeights();
	
		fShapes->TopIP();
		while ( fShapes->NextIP() )
		{					
			double temp  = density*(*Det++)*(*Weight++);

			/* integration point velocities */
			fShapes->InterpolateU(fLocVel, fDOFvec);

			double* p    = momentum.Pointer();
			double* pvel = fDOFvec.Pointer();					
			for (int dof = 0; dof < fNumDOF; dof++)			
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
		default:
			cout << "\n SolidElementT::SendKinematic: invalid output code: ";
			cout << kincode << endl;
	}

	/* number of output values */
	iArrayT n_counts;
	SetNodalOutputCodes(IOBaseT::kAtInc, flags, n_counts);

	/* reset averaging workspace */
	fNodes->ResetAverage(n_counts.Sum());
      
	/* no element output */
	iArrayT e_counts(fElementOutputCodes.Length());
	e_counts = 0;

	/* generate output */
	dArray2DT n_values, e_values;
	ComputeOutput(n_counts, n_values, e_counts, e_values);
}


/***********************************************************************
* Protected
***********************************************************************/

/* construct list of materials from the input stream */
void SolidElementT::ReadMaterialData(ifstreamT& in)
{
	/* inherited */
	ContinuumElementT::ReadMaterialData(in);
	
	/* generate list of material needs */
	fMaterialNeeds.Allocate(fMaterialList->Length());
	for (int i = 0; i < fMaterialNeeds.Length(); i++)
	{
		/* allocate */
		ArrayT<bool>& needs = fMaterialNeeds[i];
		needs.Allocate(3);

		/* casts are safe since class contructs materials list */
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[i];
		StructuralMaterialT* mat = (StructuralMaterialT*) pcont_mat;

		/* collect needs */
		needs[kNeedDisp] = mat->NeedDisp();
		needs[kNeedVel] = mat->NeedVel();
		needs[KNeedLastDisp] = mat->NeedLastDisp();
	}
}

/* print element group data */
void SolidElementT::PrintControlData(ostream& out) const
{
	/* inherited */
	ContinuumElementT::PrintControlData(out);

	/* control parameters */
	out << " Mass type code. . . . . . . . . . . . . . . . . = " << fMassType << '\n';
	out << "    eq." << kNoMass			<< ", no mass matrix\n";
	out << "    eq." << kConsistentMass	<< ", consistent mass matrix\n";
	out << "    eq." << kLumpedMass		<< ", lumped mass matrix\n";
	out << " Strain-displacement option. . . . . . . . . . . = " << fStrainDispOpt << '\n';
	out << "    eq.0, standard\n";
	out << "    eq.1, B-bar (mean dilatation)\n";
}

void SolidElementT::EchoOutputCodes(ifstreamT& in, ostream& out)
{
	/* allocate nodal output codes */
	fNodalOutputCodes.Allocate(NumNodalOutputCodes);

	/*Initialize smoothing flags */
	qUseSimo = qNoExtrap = false;

	/* read in at a time to allow comments */
	for (int i = 0; i < fNodalOutputCodes.Length(); i++)
	{
		in >> fNodalOutputCodes[i];
		
		/*Set stress smoothing flags before it's too late */
		if (fNodalOutputCodes[i] == 3)
	 	{
		  qUseSimo = qNoExtrap = true;
		}
		else if (fNodalOutputCodes[i] == 2)
		  qNoExtrap = true;
		  
		
		/*that's all */

		/* convert all to "at print increment" */
		if (fNodalOutputCodes[i] != IOBaseT::kAtNever)
			fNodalOutputCodes[i] = IOBaseT::kAtInc;
	
		if (i == iWaveSpeeds && fNodalOutputCodes[iWaveSpeeds] != IOBaseT::kAtNever)
		{
			fNormal.Allocate(fNumSD);
			in >> fNormal;
			fNormal.UnitVector();
		}
	}
		
	/* checks */
	if (fNodalOutputCodes.Min() < IOBaseT::kAtFail ||
	    fNodalOutputCodes.Max() > IOBaseT::kAtInc) throw eBadInputValue;

	/* echo */
	out << " Number of nodal output codes. . . . . . . . . . = " << fNodalOutputCodes.Length() << '\n';
	out << "    [" << fNodalOutputCodes[iNodalCoord   ] << "]: initial nodal coordinates\n";
	out << "    [" << fNodalOutputCodes[iNodalDisp    ] << "]: nodal displacements\n";
	out << "    [" << fNodalOutputCodes[iNodalStress  ] << "]: nodal stresses\n";
	out << "    [" << fNodalOutputCodes[iPrincipal    ] << "]: nodal principal stresses\n";
	out << "    [" << fNodalOutputCodes[iEnergyDensity] << "]: nodal strain energy density\n";
	out << "    [" << fNodalOutputCodes[iWaveSpeeds   ] << "]: wave speeds\n";
	out << "    [" << fNodalOutputCodes[iMaterialData ] << "]: nodal material output parameters\n";
	
	if (fNodalOutputCodes[iWaveSpeeds] == 1)
	{
		out << " Wave speed sampling direction:\n";
		for (int i = 0; i < fNumSD; i++)
			out << "   N[" << i+1 << "] = " << fNormal[i] << '\n';
	}

	/* allocate nodal output codes */
	fElementOutputCodes.Allocate(NumElementOutputCodes);
	fElementOutputCodes = IOBaseT::kAtNever;

//TEMP - backward compatibility
	if (StringT::versioncmp(fFEManager.Version(), "v3.01") < 1)
	{
		/* message */
		cout << "\n SolidElementT::EchoOutputCodes: use input file version newer than v3.01\n" 
		     <<   "     to enable element output control" << endl;
		out << "\n SolidElementT::EchoOutputCodes: use input file version newer than v3.01\n" 
		    <<   "     to enable element output control" << endl;	
	}
	else
	{
		int num_codes = (StringT::versioncmp(fFEManager.Version(), "v3.4.1") < 0) ? 5 : 7;
	
		/* read in at a time to allow comments */
		for (int j = 0; j < num_codes; j++)
		{
			in >> fElementOutputCodes[j];
		
			/* convert all to "at print increment" */
			if (fElementOutputCodes[j] != IOBaseT::kAtNever)
				fElementOutputCodes[j] = IOBaseT::kAtInc;
		}	

		/* defaults */
		if (fController->Order() == 0)
		{
			fElementOutputCodes[iKineticEnergy] = IOBaseT::kAtNever;
			fElementOutputCodes[iLinearMomentum] = IOBaseT::kAtNever;
		}

		/* checks */
		if (fElementOutputCodes.Min() < IOBaseT::kAtFail ||
		    fElementOutputCodes.Max() > IOBaseT::kAtInc) throw eBadInputValue;
	}

	/* echo */
	out << " Number of element output codes. . . . . . . . . = " << fElementOutputCodes.Length() << '\n';
	out << "    [" << fElementOutputCodes[iCentroid      ] << "]: reference centroid\n";
	out << "    [" << fElementOutputCodes[iMass          ] << "]: ip mass\n";
	out << "    [" << fElementOutputCodes[iStrainEnergy  ] << "]: strain energy\n";
	out << "    [" << fElementOutputCodes[iKineticEnergy ] << "]: kinetic energy\n";
	out << "    [" << fElementOutputCodes[iLinearMomentum] << "]: linear momentum\n";
	out << "    [" << fElementOutputCodes[iIPStress      ] << "]: ip stresses\n";
	out << "    [" << fElementOutputCodes[iIPMaterialData] << "]: ip material output parameters\n";
}

/* construct output labels array */
void SolidElementT::SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
	/* initialize */
	counts.Allocate(flags.Length());
	counts = 0;

	/* set output flags */
	if (flags[iNodalCoord] == mode)
		counts[iNodalCoord] = fNumSD;
	if (flags[iNodalDisp] == mode)
		counts[iNodalDisp] = fNumDOF;
	if (flags[iNodalStress] == mode)
		counts[iNodalStress] = dSymMatrixT::NumValues(fNumSD);
	if (flags[iPrincipal] == mode)
		counts[iPrincipal] = fNumSD;
	if (flags[iEnergyDensity] == mode)
		counts[iEnergyDensity] = 1;
	if (flags[iWaveSpeeds] == mode)
		counts[iWaveSpeeds] = fNumSD;
	if (flags[iMaterialData] == mode)
		counts[iMaterialData] = (*fMaterialList)[0]->NumOutputVariables();
}

void SolidElementT::SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
	/* initialize */
	counts.Allocate(flags.Length());
	counts = 0;

	/* set output flags */
	if (fElementOutputCodes[iCentroid] == mode) counts[iCentroid] = fNumSD;
	if (fElementOutputCodes[iMass] == mode) counts[iMass] = fNumIP;
	if (fElementOutputCodes[iStrainEnergy] == mode) counts[iStrainEnergy] = 1;
	if (fElementOutputCodes[iKineticEnergy] == mode) counts[iKineticEnergy] = 1;
	if (fElementOutputCodes[iLinearMomentum] == mode) counts[iLinearMomentum] = fNumDOF;
	if (fElementOutputCodes[iIPStress] == mode) counts[iIPStress] = dSymMatrixT::NumValues(fNumSD)*fNumIP;
	if (fElementOutputCodes[iIPMaterialData] == mode) 
		counts[iIPMaterialData] = (*fMaterialList)[0]->NumOutputVariables()*fNumIP;
}

/* initialize local arrays */
void SolidElementT::SetLocalArrays(void)
{
	/* inherited */
	ContinuumElementT::SetLocalArrays();

	/* dimension */
	fLocLastDisp.Allocate(fNumElemNodes, fNumDOF);
	fLocVel.Allocate(fNumElemNodes, fNumDOF);
	fLocAcc.Allocate(fNumElemNodes, fNumDOF);

	/* set source */
	fFEManager.RegisterLocal(fLocLastDisp);
	fFEManager.RegisterLocal(fLocVel);
	fFEManager.RegisterLocal(fLocAcc);
}

/* set the correct shape functions */
void SolidElementT::SetShape(void)
{
	/* construct shape functions */
	ShapeFunctionT::StrainOptionT strain_opt = (fStrainDispOpt == 0) ? 
		ShapeFunctionT::kStandardB : ShapeFunctionT::kMeanDilBbar;
	fShapes = new ShapeFunctionT(fGeometryCode, fNumIP,
		fLocInitCoords, strain_opt);
	if (!fShapes) throw eOutOfMemory;

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
			double onebydt = 1.0/fFEManager.TimeStep();
			fLocVel.SetToCombination(onebydt, fLocDisp, -onebydt, fLocLastDisp);
		}
	}	
}

/* construct the effective mass matrix */
void SolidElementT::LHSDriver(void)
{
	/* inherited */
	ContinuumElementT::LHSDriver();

	/* element contribution */
	ElementLHSDriver();
}

void SolidElementT::ElementLHSDriver(void)
{
	/* set components and weights */
	double constM = 0.0;
	double constC = 0.0;
	double constK = 0.0;
	
	int formM = fController->FormM(constM);
	int formC = fController->FormC(constC);
	int formK = fController->FormK(constK);

	/* override algorithm */
	if (fMassType == kNoMass) formM = 0;

	/* quick exit */
	if ((formM == 0 && formC == 0 && formK == 0) ||
	    (fabs(constM) < kSmall &&
	     fabs(constC) < kSmall &&
	     fabs(constK) < kSmall)) return;

	/* loop over elements */
	Top();
	while (NextElement())
	{
		double constKe = constK;
		double constMe = constM;
	
		/* initialize */
		fLHS = 0.0;
		
		/* set shape function derivatives */
		SetGlobalShape();
	
		/* Rayleigh damping */
//		if (formC)
//		{
//			constKe += constC*(fCurrMaterial->StiffnessDamping());
//			constMe += constC*(fCurrMaterial->MassDamping());
//		}
//DEV - Rayleigh damping is too ugly to keep, could add some
//      phenomenological damping to the constitutive models
		
		/* element mass */
		if (fabs(constMe) > kSmall)
			FormMass(fMassType, constMe*(fCurrMaterial->Density()));

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
	/* set components and weights */
	double constMa = 0.0;
	double constCv = 0.0;
	double constKd = 0.0;
	
	/* components dicated by the algorithm */
	int formMa = fController->FormMa(constMa);
	int formCv = fController->FormCv(constCv);
	int formKd = fController->FormKd(constKd);

	/* body forces */
	int formBody = 0;
	if (fMassType != kNoMass &&
	   (fBodyForceLTf > -1 && fBody.Magnitude() > kSmall))
	{	
		formBody = 1;
		if (!formMa) constMa = 1.0; /* override */
	}

	/* override controller */
	if (fMassType == kNoMass) formMa = 0;

	Top();
	while (NextElement())
	{
		/* effective accelerations and displacements */
		//ComputeEffectiveDVA(formBody, formMa, constMa, formCv, constCv, formKd, constKd);
		//DEV - Rayleigh damping should be moved to the constitutive level
	
		/* last check w/ effective a and d - override controller */
		//int eformMa = fLocAcc.AbsMax() > 0.0;
		//int eformCv = fLocVel.AbsMax() > 0.0;
		//int eformKd = (fLocDisp.AbsMax() > 0.0 ||
		//              fCurrMaterial->HasInternalStrain());
		//DEV - although this may speed things up in some cases
		//      it's a pain in general -> these could be moved
		//      within FormKd, FormCv, or FormMa since the state
		//      of the local nodal vectors should be fixed by then

		/* initialize */
		fRHS = 0.0;
		
		/* global shape function values */
		SetGlobalShape();
			
		/* internal force contribution */	
		if (formKd) FormKd(-1.0);
				
		/* damping */
		//if (formCv) FormCv(-1.0);
		//DEV - computed at constitutive level

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
		
			FormMa(fMassType, -(fCurrMaterial->Density()), fLocAcc);			  		
		}
								
		/* assemble */
		AssembleRHS();
	}
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
		fCurrMaterial = (StructuralMaterialT*) pcont_mat;
	}
	
	return result;
}

/* form the element stiffness matrix */
void SolidElementT::FormStiffness(double constK)
{
	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	/* integrate element stiffness */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		double scale = constK*(*Det++)*(*Weight++);
	
		/* strain displacement matrix */
		fShapes->B(fB);

		/* get D matrix */
		fD.SetToScaled(scale, fCurrMaterial->c_ijkl());
							
		/* multiply b(transpose) * db, taking account of symmetry, */
		/* and accumulate in elstif */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}
}

//DEV - Rayleigh damping should be added to the constitutive level
#if 0		
/* compute the effective acceleration and velocities based
* on the algorithmic flags formXx and the given constants
* constXx.
*
*		acc_eff  = constMa acc  + constCv a vel
*      vel_eff  = 0
*      disp_eff = constKd disp + constCv b vel
*
* where a and b are the Rayleigh damping coefficients.  No
* effective velocity since it's accounted for in the effective
* a and d.
*
* Note: In the process, the function collects the required
*       local arrays */
void SolidElementT::ComputeEffectiveDVA(int formBody,
	int formMa, double constMa, int formCv, double constCv,
	int formKd, double constKd)
{
	/* acceleration */
	if (formMa || formBody)
	{
		if (formMa)
			SetLocalU(fLocAcc);
		else
			fLocAcc = 0.0;
		
		if (formBody) AddBodyForce(fLocAcc);

		fLocAcc *= constMa;	
	}
	else
		fLocAcc = 0.0;
	
	/* displacement */
	if (formKd)
	{
		SetLocalU(fLocDisp);
		fLocDisp *= constKd;
	}
	else
		fLocDisp = 0.0;
	
	/* Rayleigh damping */
	if (formCv)
	{
		SetLocalU(fLocVel);
		fLocVel *= constCv;
		
		/* effective a and d */
		fLocAcc.AddScaled(fCurrMaterial->MassDamping(), fLocVel);
		fLocDisp.AddScaled(fCurrMaterial->StiffnessDamping(), fLocVel);
		
		/* effective v */
		fLocVel = 0.0;
	}
	else
		fLocVel = 0.0;
}	
#endif

/* form of tangent matrix */
GlobalT::SystemTypeT SolidElementT::TangentType(void) const
{
	/* special case */
	if (fController->Order() > 0 &&
	    fController->ImplicitExplicit() ==  eControllerT::kExplicit &&
	    (fMassType == kNoMass ||
	     fMassType == kLumpedMass))
		return GlobalT::kDiagonal;
	else
		/* inherited */
		return ContinuumElementT::TangentType();
}

/* fast approximation to Rayleigh damping for explicit
* dynamics. approximate mass as equally distributed
* among the nodes */
void SolidElementT::FormRayleighMassDamping(double constM)
{
	/* compute total mass */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	double totmas = 0.0;
	int    numint = fShapes ->NumIP();
	for (int i = 0; i < numint; i++)
		totmas += constM*(*Det++)*(*Weight++);
		
	/* equally distributed */
	double nodalmass = totmas/(fLocAcc.NumberOfNodes());
	
	/* assemble into RHS vector */
	fRHS.AddScaled(nodalmass,fLocAcc);		
}

/* calculate the damping force contribution ("-c*v") */
void SolidElementT::FormCv(double constC)
{
#pragma unused(constC)
}

/* calculate the internal force contribution ("-k*d") */
void SolidElementT::FormKd(double constK)
{
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		/* get strain-displacement matrix */
		fShapes->B(fB);

		/* B^T * Cauchy stress */
		fB.MultTx(fCurrMaterial->s_ij(), fNEEvec);
		
		/* accumulate */
		fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);
	}	
}

/* return a pointer to a new material list */
MaterialListT* SolidElementT::NewMaterialList(int size) const
{
	/* allocate */
	if (fNumSD == 2)
		return new MaterialList2DT(size, *this);
	else if (fNumSD == 3)
		return new MaterialList3DT(size, *this);
	else
		return NULL;		
}

/* extrapolate the integration point stresses and strains and extrapolate */
void SolidElementT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	const iArrayT& e_codes, dArray2DT& e_values)
{
	/* number of output values */
	int n_out = n_codes.Sum();
	int e_out = e_codes.Sum();
	
	int n_simo, n_extrap;
	if (qUseSimo)
	{
		n_simo = n_out - n_codes[iNodalDisp] - n_codes[iNodalCoord];
		n_extrap = n_codes[iNodalDisp] + n_codes[iNodalCoord];
	}
	else
	{
		n_simo = 0;
		n_extrap = n_out;
	}
		
	/* nothing to output */
	if (n_out == 0 && e_out == 0) return;
	
	/* dimensions */
	int nsd = NumSD();
	int ndof = NumDOF();
	int nen = fNumElemNodes;

	/* reset averaging workspace */
	fNodes->ResetAverage(n_extrap);
	
	/* allocate element results space */
	e_values.Allocate(NumElements(), e_out);

	/* nodal work arrays */
	dArray2DT nodal_space(nen, n_extrap);
	dArray2DT nodal_all(nen, n_extrap);
	dArray2DT coords, disp;
	dArray2DT nodalstress, princstress, matdat;
	dArray2DT energy, speed;

	/* ip values */
	dSymMatrixT cauchy(nsd);
	dArrayT ipmat(n_codes[iMaterialData]), ipenergy(1);
	dArrayT ipspeed(nsd), ipprincipal(nsd);

	/* set shallow copies */
	double* pall = nodal_space.Pointer();
	coords.Set(nen, n_codes[iNodalCoord], pall)      ; pall += coords.Length();
	disp.Set(nen, n_codes[iNodalDisp], pall)         ; pall += disp.Length();
	
	/* workspaces for Simo */
	int simo_offset = coords.MinorDim() + disp.MinorDim();
	dArray2DT simo_space(nen,qUseSimo ? n_simo : 0);
	dArray2DT simo_all(nen,qUseSimo ? n_simo : 0);
	dArray2DT simoNa_bar(nen,qUseSimo ? 1 : 0);
	dArray2DT simo_force;
	dArray2DT simo_mass;
	
	if (!qUseSimo) 
	{
		nodalstress.Set(nen, n_codes[iNodalStress], pall); pall += nodalstress.Length();
		princstress.Set(nen, n_codes[iPrincipal], pall)  ; pall += princstress.Length();
		energy.Set(nen, n_codes[iEnergyDensity], pall)   ; pall += energy.Length();
		speed.Set(nen, n_codes[iWaveSpeeds], pall)       ; pall += speed.Length();
		matdat.Set(nen, n_codes[iMaterialData], pall);  
	}
	else
	{
	        simo_force.Allocate(fNodes->NumNodes(),qUseSimo ? n_simo : 0);
		simo_mass.Allocate(fNodes->NumNodes(),qUseSimo ? 1 : 0);
	
		pall = simo_space.Pointer();
		nodalstress.Set(nen, n_codes[iNodalStress], pall); pall += nodalstress.Length();
		princstress.Set(nen, n_codes[iPrincipal], pall)  ; pall += princstress.Length();
		energy.Set(nen, n_codes[iEnergyDensity], pall)   ; pall += energy.Length();
		speed.Set(nen, n_codes[iWaveSpeeds], pall)       ; pall += speed.Length();
		matdat.Set(nen, n_codes[iMaterialData], pall); 
		
		simo_mass = 0.;
		simo_force = 0.;
	}
	
	/* element work arrays */
	dArrayT element_values(e_values.MinorDim());
	pall = element_values.Pointer();
	dArrayT centroid, ip_centroid, ip_mass;
	if (e_codes[iCentroid])
	{
		centroid.Set(nsd, pall); pall += nsd;
		ip_centroid.Allocate(nsd);
	}
	if (e_codes[iMass]) {
		ip_mass.Set(NumIP(), pall); 
		pall += NumIP();
	}
	double w_tmp, ke_tmp;
	double mass;
	double& strain_energy = (e_codes[iStrainEnergy]) ? *pall++ : w_tmp;
	double& kinetic_energy = (e_codes[iKineticEnergy]) ? *pall++ : ke_tmp;
	dArrayT linear_momentum, ip_velocity;
	if (e_codes[iLinearMomentum])
	{
		linear_momentum.Set(ndof, pall); pall += ndof;
		ip_velocity.Allocate(ndof);
	}
	dArray2DT ip_stress;
	if (e_codes[iIPStress])
	{
		ip_stress.Set(NumIP(), e_codes[iIPStress]/NumIP(), pall);
		pall += ip_stress.Length();
	}
	dArray2DT ip_material_data;
	if (e_codes[iIPMaterialData])
	{
		ip_material_data.Set(NumIP(), e_codes[iIPMaterialData]/NumIP(), pall);
		pall += ip_material_data.Length();
		ipmat.Allocate(ip_material_data.MinorDim());
	}

	/* check that degrees are displacements */
	int interpolant_DOF = InterpolantDOFs();

	Top();
	while (NextElement())
	{
		/* initialize */
	    nodal_space = 0.0;
	    simo_space = 0.;
	    simo_all = 0.;
	    simoNa_bar = 0.;

		/* global shape function values */
		SetGlobalShape();
		
		/* collect nodal values */
		if (e_codes[iKineticEnergy] || e_codes[iLinearMomentum])
			SetLocalU(fLocVel);
		
		/* coordinates and displacements all at once */
		if (n_codes[iNodalCoord]) fLocInitCoords.ReturnTranspose(coords);
		if (n_codes[ iNodalDisp])
		{
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
		double density = fCurrMaterial->Density();

		/* integrate */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/* element integration weight */
			double ip_w = (*j++)*(*w++);
			dArray2DT Na_X_ip_w;
			if (qUseSimo || qNoExtrap)
			{
			  if (qUseSimo)
			    {
				const double* Na_X = fShapes->IPShapeX();
				Na_X_ip_w.Allocate(fNumElemNodes,1);
				Na_X_ip_w = ip_w;
				for (int k = 0;k<fNumElemNodes;k++)
				{
					Na_X_ip_w(k,0) *= *Na_X++;
				}
				simoNa_bar += Na_X_ip_w;
			    }
			  else
			    for (int k = 0; k < fNumElemNodes;k++)
				Na_X_ip_w(k,0) = 1.;
			}
		
			/* get Cauchy stress */
			cauchy = fCurrMaterial->s_ij();

			/* stresses */
			if (n_codes[iNodalStress])
			{	
				if (qNoExtrap)
					for (int k = 0;k<fNumElemNodes;k++)
						nodalstress.AddToRowScaled(k,Na_X_ip_w(k,0),cauchy);
				else
					fShapes->Extrapolate(cauchy, nodalstress);
			}
			
			if (e_codes[iIPStress]) ip_stress.SetRow(fShapes->CurrIP(), cauchy);

			/* wave speeds */
			if (n_codes[iWaveSpeeds])
			{
				/* acoustic wave speeds */
				fCurrMaterial->WaveSpeeds(fNormal, ipspeed);
				if (qNoExtrap)
					for (int k = 0;k<fNumElemNodes;k++)
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
					for (int k = 0;k<fNumElemNodes;k++)
						princstress.AddToRowScaled(k,Na_X_ip_w(k,0),ipprincipal);
				else
					fShapes->Extrapolate(ipprincipal, princstress);	
			}

			/* strain energy density */
			if (n_codes[iEnergyDensity] || e_codes[iStrainEnergy])
			{
				double ip_strain_energy = fCurrMaterial->StrainEnergyDensity();
			
				/* nodal average */
				if (n_codes[iEnergyDensity])
				{
					ipenergy[0] = ip_strain_energy;
					if (qNoExtrap)
					  for (int k = 0;k<fNumElemNodes;k++)
						energy.AddToRowScaled(k,Na_X_ip_w(k,0),ipenergy);
					else
						fShapes->Extrapolate(ipenergy,energy);
				}
				
				/* integrate over element */
				if (e_codes[iStrainEnergy])
					strain_energy += ip_w*ip_strain_energy;
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
						for (int k = 0;k<fNumElemNodes;k++)
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
				if (e_codes[iCentroid])
				{
					fShapes->IPCoords(ip_centroid);
					centroid.AddScaled(ip_w*density, ip_centroid);
				}
			}
			
			/* kinetic energy/linear momentum */
			if (e_codes[iKineticEnergy] || e_codes[iLinearMomentum])
			{
				/* velocity at integration point */
				fShapes->InterpolateU(fLocVel, ip_velocity);
				
				/* kinetic energy */
				if (e_codes[iKineticEnergy])
					kinetic_energy += 0.5*ip_w*density*dArrayT::Dot(ip_velocity, ip_velocity);
					
				/* linear momentum */
				if (e_codes[iLinearMomentum])
					linear_momentum.AddScaled(ip_w*density, ip_velocity);
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
		}
		else
		{	
			colcount = 0;
			simo_all.BlockColumnCopyAt(nodalstress, colcount); colcount += nodalstress.MinorDim();
			simo_all.BlockColumnCopyAt(princstress, colcount); colcount += princstress.MinorDim();
			simo_all.BlockColumnCopyAt(energy     , colcount); colcount += energy.MinorDim();
			simo_all.BlockColumnCopyAt(speed      , colcount); colcount += speed.MinorDim();
			simo_all.BlockColumnCopyAt(matdat     , colcount); colcount += matdat.MinorDim();
			simo_force.Accumulate(CurrentElement().NodesX(),simo_all);
			simo_mass.Accumulate(CurrentElement().NodesX(),simoNa_bar);
		}
	    
		/* accumulate - extrapolation done from ip's to corners => X nodes */
		fNodes->AssembleAverage(CurrentElement().NodesX(), nodal_all);
		
		/* element values */
		if (e_codes[iCentroid]) centroid /= mass;
		
		/* store results */
		e_values.SetRow(CurrElementNumber(), element_values);
	}

	/* get nodally averaged values */
	dArray2DT extrap_values;
	fNodes->OutputUsedAverage(extrap_values);

	n_values.Allocate(extrap_values.MajorDim(),n_out);
	n_values.BlockColumnCopyAt(extrap_values,0);
	if (qUseSimo)
	{	
		for (int i = 0; i < simo_force.MajorDim();i++)
			simo_force.ScaleRow(i,1./simo_mass(i,0));	
		n_values.BlockColumnCopyAt(simo_force,simo_offset);
	}
}

/***********************************************************************
* Private
***********************************************************************/

/* construct output labels array */
void SolidElementT::GenerateOutputLabels(const iArrayT& n_codes, ArrayT<StringT>& n_labels, 
	const iArrayT& e_codes, ArrayT<StringT>& e_labels) const
{
	/* allocate */
	n_labels.Allocate(n_codes.Sum());

	int count = 0;
	if (n_codes[iNodalDisp])
	{
		if (fNumDOF > 3) throw eGeneralFail;
		const char* dlabels[3] = {"D_X", "D_Y", "D_Z"};
		for (int i = 0; i < fNumDOF; i++)
			n_labels[count++] = dlabels[i];
	}

	if (n_codes[iNodalCoord])
	{
		const char* xlabels[] = {"x1", "x2", "x3"};
		for (int i = 0; i < fNumSD; i++)
			n_labels[count++] = xlabels[i];
	}

	if (n_codes[iNodalStress])
	{
		const char* slabels2D[] = {"s11", "s22", "s12"};
		const char* slabels3D[] = {"s11", "s22", "s33", "s23", "s13", "s12"};
		const char**    slabels = (fNumSD == 2) ? slabels2D : slabels3D;
		for (int i = 0; i < dSymMatrixT::NumValues(fNumSD); i++)
			n_labels[count++] = slabels[i];
	}
		
	if (n_codes[iPrincipal])
	{
		const char* plabels[] = {"s1", "s2", "s3"};
		for (int i = 0; i < fNumSD; i++)
			n_labels[count++] = plabels[i];
	}
		
	if (n_codes[iEnergyDensity]) n_labels[count++] = "phi";
	if (n_codes[iWaveSpeeds])
	{
		const char* clabels2D[] = {"cd", "cs"};
		const char* clabels3D[] = {"cd", "cs_min", "cs_max"};
		const char**    clabels = (fNumSD == 2) ? clabels2D : clabels3D;
		for (int i = 0; i < fNumSD; i++)
			n_labels[count++] = clabels[i];		
	}

	/* material output labels */
	if (n_codes[iMaterialData])
	{
		ArrayT<StringT> matlabels;
		(*fMaterialList)[0]->OutputLabels(matlabels);	
		
		for (int i = 0; i < matlabels.Length(); i++)
			n_labels[count++] = matlabels[i];
	}

	/* allocate */
	e_labels.Allocate(e_codes.Sum());
	count = 0;
	if (e_codes[iCentroid])
	{
		const char* xlabels[] = {"xc_1", "xc_2", "xc_3"};
		for (int i = 0; i < fNumSD; i++)
			e_labels[count++] = xlabels[i];
	}
	if (e_codes[iMass])
	{
		/* over integration points */
		for (int j = 0; j < fNumIP; j++)
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
		for (int i = 0; i < fNumDOF; i++)
			e_labels[count++] = plabels[i];
	}
	if (e_codes[iIPStress])
	{
		const char* slabels2D[] = {"s11", "s22", "s12"};
		const char* slabels3D[] = {"s11", "s22", "s33", "s23", "s13", "s12"};
		const char**    slabels = (fNumSD == 2) ? slabels2D : slabels3D;

		/* over integration points */
		for (int j = 0; j < fNumIP; j++)
		{
			StringT ip_label;
			ip_label.Append("ip", j+1);
			
			/* over stress components */
			for (int i = 0; i < dSymMatrixT::NumValues(fNumSD); i++)
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
		for (int j = 0; j < fNumIP; j++)
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
