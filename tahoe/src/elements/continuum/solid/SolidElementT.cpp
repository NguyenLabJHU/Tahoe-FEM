/* $Id: SolidElementT.cpp,v 1.1.1.1 2001-01-29 08:20:39 paklein Exp $ */
/* created: paklein (05/28/1996)                                          */

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

//TEMP
#include "FDStructMatT.h"

/* exception codes */
#include "ExceptionCodes.h"

/* initialize static data */
const int SolidElementT::NumOutputCodes = 7;

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
	fB.Allocate((fNumSD == 2) ? 3 : 6, fNumSD*fNumElemNodes);

	/* override */
	if (fController->ImplicitExplicit() == eControllerT::kExplicit &&
	    fMassType != kNoMass &&
	    fMassType != kLumpedMass)
//{	
//cout << "\n SolidElementT::Initialize: SKIPPING lumped mass override for explicit dynamics" << endl;
		fMassType = kLumpedMass;
//}

	/* setup for material output */
	if (fOutputCodes[iMaterialData])
	{
		/* only 1 material allowed for the group */
		if (fMaterialList->Length() > 1)
		{
			cout << "\n SolidElementT::ReadMaterialData: if the material output flag is set,\n";
			cout <<   " there can be only 1 material defined for the element group."<< endl;

			throw eGeneralFail;
		}
		/* no material output variables */
	 	else if ((*fMaterialList)[0]->NumOutputVariables() == 0)
		{
			cout << "\n SolidElementT::ReadMaterialData: there are no material outputs. ";
			cout << endl;

			fOutputCodes[iMaterialData] = IOBaseT::kAtNever;
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
	int is_dynamic  = (fController->StaticDynamic() == eControllerT::kDynamic);

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
	if ( !fConnectivities.HasValue(node) ) return;

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
		if (!formMa) constMa = 1.0; // correct value ??
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
			ComputeEffectiveDVA(formBody, formMa, constMa, formCv, constCv, formKd, constKd);

			/* global shape function values */
			SetGlobalShape();
	
			/* internal force contribution */	
			if (formKd) FormKd(1.0);
				
			/* damping */
			if (formCv) FormCv(1.0);

			/* inertia forces */
			if (formMa) FormMa(fMassType, fCurrMaterial->Density(), fLocAcc);

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
		
		/* get displacements */
		SetLocalU(fLocDisp);
		
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
iArrayT flags(fOutputCodes.Length());

	/* set flags to get desired output */
	flags = IOBaseT::kAtNever;
	switch (kincode)
	{
		case iNodalDisp:
		    flags[iNodalDisp] = fNumDOF;
			break;
		case iNodalStress:
		    flags[iNodalStress] = dSymMatrixT::NumValues(fNumSD);
			break;
		case iEnergyDensity:
		    flags[iEnergyDensity] = 1;
			break;
		case iPrincipal:
			flags[iPrincipal] = fNumSD;
			break;
		default:
			cout << "\n SolidElementT::SendKinematic: invalid output code: ";
			cout << kincode << endl;
	}

	/* number of output values */
	iArrayT counts;
	SetOutputCodes(IOBaseT::kAtInc, flags, counts);
	int num_out = counts.Sum();

	/* reset averaging workspace */
	fNodes->ResetAverage(num_out);

	/* generate nodal values */
	ComputeNodalValues(counts);
}

/***********************************************************************
* Protected
***********************************************************************/

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
	/* allocate */
	fOutputCodes.Allocate(NumOutputCodes);

	/* read in at a time to allow comments */
	for (int i = 0; i < fOutputCodes.Length(); i++)
	{
		in >> fOutputCodes[i];
		
		/* convert all to "at print increment" */
		if (fOutputCodes[i] != IOBaseT::kAtNever)
			fOutputCodes[i] = IOBaseT::kAtInc;
	
		if (i == iWaveSpeeds && fOutputCodes[iWaveSpeeds] != IOBaseT::kAtNever)
		{
			fNormal.Allocate(fNumSD);
			in >> fNormal;
			fNormal.UnitVector();
		}
	}
		
	/* checks */
	if (fOutputCodes.Min() < IOBaseT::kAtFail ||
	    fOutputCodes.Max() > IOBaseT::kAtInc) throw eBadInputValue;

	/* default behavior with output formats */
//	fOutputCodes[iNodalCoord] = IOBaseT::kAtNever;
//	fOutputCodes[iNodalDisp ] = IOBaseT::kAtNever;
// what to do about the defaults.

	/* echo */
	out << " Number of nodal output codes. . . . . . . . . . = " << NumOutputCodes << '\n';
	out << "    [" << fOutputCodes[iNodalCoord   ] << "]: initial nodal coordinates\n";
	out << "    [" << fOutputCodes[iNodalDisp    ] << "]: nodal displacements\n";
	out << "    [" << fOutputCodes[iNodalStress  ] << "]: nodal stresses\n";
	out << "    [" << fOutputCodes[iPrincipal    ] << "]: nodal principal stresses\n";
	out << "    [" << fOutputCodes[iEnergyDensity] << "]: nodal strain energy density\n";
	out << "    [" << fOutputCodes[iWaveSpeeds   ] << "]: wave speeds\n";
	out << "    [" << fOutputCodes[iMaterialData ] << "]: nodal material output parameters\n";
	
	if (fOutputCodes[iWaveSpeeds] == 1)
	{
		out << " Wave speed sampling direction:\n";
		for (int i = 0; i < fNumSD; i++)
			out << "   N[" << i+1 << "] = " << fNormal[i] << '\n';
	}
}

/* construct output labels array */
void SolidElementT::SetOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts)
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

/* initialize local arrays */
void SolidElementT::SetLocalArrays(void)
{
	/* inherited */
	ContinuumElementT::SetLocalArrays();

	/* dimension */
	fLocVel.Allocate(fNumElemNodes, fNumDOF);
	fLocAcc.Allocate(fNumElemNodes, fNumDOF);

	/* set source */
	fFEManager.RegisterLocal(fLocVel);
	fFEManager.RegisterLocal(fLocAcc);
}

/* set the correct shape functions */
void SolidElementT::SetShape(void)
{
	fShapes = new ShapeFunctionT(fGeometryCode, fNumIP,
		fLocInitCoords, fStrainDispOpt);
	if (!fShapes) throw eOutOfMemory;

	fShapes->Initialize();
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
		if (formC)
		{
			constKe += constC*(fCurrMaterial->StiffnessDamping());
			constMe += constC*(fCurrMaterial->MassDamping());
		}
		
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
		if (!formMa) constMa = 1.0; //override controller value??
	}

	/* override controller */
	if (fMassType == kNoMass) formMa = 0;

	Top();
	while (NextElement())
	{
		/* effective accelerations and displacements */
		ComputeEffectiveDVA(formBody, formMa, constMa, formCv, constCv, formKd, constKd);
	
		/* last check w/ effective a and d - override controller */
		int eformMa = fLocAcc.AbsMax() > 0.0;
		int eformCv = fLocVel.AbsMax() > 0.0;
		int eformKd = (fLocDisp.AbsMax() > 0.0 ||
		               fCurrMaterial->HasInternalStrain());

		if (eformMa || eformCv || eformKd)
		{
			/* initialize */
			fRHS = 0.0;
		
			/* global shape function values */
			SetGlobalShape();
			
			/* internal force contribution */	
			if (eformKd) FormKd(-1.0);
				
			/* damping */
			if (eformCv) FormCv(-1.0);

			/* inertia forces */
			if (eformMa) FormMa(fMassType, -(fCurrMaterial->Density()), fLocAcc);			  		
								
			/* assemble */
			AssembleRHS();
		}
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

/* form of tangent matrix */
GlobalT::SystemTypeT SolidElementT::TangentType(void) const
{
	/* special case */
	if (fController->StaticDynamic() == eControllerT::kDynamic &&
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
void SolidElementT::ComputeNodalValues(const iArrayT& codes)
{
/* number of nodally smoothed values */
int num_out = codes.Sum();

	/* nothing to output */
	if (num_out == 0) return;

	/* work arrays */
	dArray2DT nodal_space(fNumElemNodes, num_out);
	dArray2DT nodal_all(fNumElemNodes, num_out);
	dArray2DT coords, disp;
	dArray2DT nodalstress, princstress, matdat;
	dArray2DT energy, speed;

	/* ip values */
	dSymMatrixT cauchy(fNumSD);
	dArrayT ipmat(codes[iMaterialData]), ipenergy(1);
	dArrayT ipspeed(fNumSD), ipprincipal(fNumSD);

	/* set shallow copies */
	double* pall = nodal_space.Pointer();
	coords.Set(fNumElemNodes, codes[iNodalCoord], pall);
	pall += coords.Length();
	disp.Set(fNumElemNodes, codes[iNodalDisp], pall);
	pall += disp.Length();
	nodalstress.Set(fNumElemNodes, codes[iNodalStress], pall);
	pall += nodalstress.Length();
	princstress.Set(fNumElemNodes, codes[iPrincipal], pall);
	pall += princstress.Length();
	energy.Set(fNumElemNodes, codes[iEnergyDensity], pall);
	pall += energy.Length();
	speed.Set(fNumElemNodes, codes[iWaveSpeeds], pall);
	pall += speed.Length();
	matdat.Set(fNumElemNodes, codes[iMaterialData], pall);

	/* check that degrees are displacements */
	int interpolant_DOF = InterpolantDOFs();

	Top();
	while ( NextElement() )
	{
//TEMP
//out << '\n' << fElementCards.Position() << '\n';

		/* initialize */
	    nodal_space = 0.0;

		/* global shape function values */
		SetGlobalShape();
		SetLocalU(fLocDisp);
		
		/* coordinates and displacements all at once */
		if (codes[iNodalCoord]) fLocInitCoords.ReturnTranspose(coords);
		if (codes[ iNodalDisp])
		{
			if (interpolant_DOF)
				fLocDisp.ReturnTranspose(disp);
			else
				NodalDOFs(CurrentElement().NodesX(), disp);
		}

		/* integrate */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
//TEMP - output ip strains
//const dSymMatrixT& strn = fCurrMaterial->Strain(fLocDisp);
//out << setw(kDoubleWidth) << strn(0,0);
//out << setw(kDoubleWidth) << strn(1,1);
//out << setw(kDoubleWidth) << strn(0,1) << '\n';

			/* get Cauchy stress */
			cauchy = fCurrMaterial->s_ij();

			/* nodal stress */
			if (codes[iNodalStress]) fShapes->Extrapolate(cauchy, nodalstress);

			/* wave speeds */
			if (codes[iWaveSpeeds])
			{
				/* acoustic wave speeds */
				fCurrMaterial->WaveSpeeds(fNormal, ipspeed);
				fShapes->Extrapolate(ipspeed, speed);
			}

			/* principal values - compute principal before smoothing */
			if (codes[iPrincipal])
			{
				/* compute eigenvalues */
				cauchy.PrincipalValues(ipprincipal);
				fShapes->Extrapolate(ipprincipal, princstress);	
			}

			/* strain energy density */
			if (codes[iEnergyDensity])
			{
				ipenergy[0] = fCurrMaterial->StrainEnergyDensity();
				fShapes->Extrapolate(ipenergy,energy);
			}

			/* material stuff */
			if (codes[iMaterialData])
			{
				fCurrMaterial->ComputeOutput(ipmat);
				fShapes->Extrapolate(ipmat, matdat);
			}
		}

		/* copy in the cols */
		int colcount = 0;
		nodal_all.BlockColumnCopyAt(disp       , colcount); colcount += disp.MinorDim();
		nodal_all.BlockColumnCopyAt(coords     , colcount); colcount += coords.MinorDim();
		nodal_all.BlockColumnCopyAt(nodalstress, colcount); colcount += nodalstress.MinorDim();
		nodal_all.BlockColumnCopyAt(princstress, colcount); colcount += princstress.MinorDim();
		nodal_all.BlockColumnCopyAt(energy     , colcount); colcount += energy.MinorDim();
		nodal_all.BlockColumnCopyAt(speed      , colcount); colcount += speed.MinorDim();
		nodal_all.BlockColumnCopyAt(matdat     , colcount);

		/* accumulate - extrapolation done from ip's to corners => X nodes */
		fNodes->AssembleAverage(CurrentElement().NodesX(), nodal_all);
	}
}

/***********************************************************************
* Private
***********************************************************************/

/* construct output labels array */
void SolidElementT::GenerateOutputLabels(const iArrayT& codes,
	ArrayT<StringT>& labels) const
{
	/* allocate */
	labels.Allocate(codes.Sum());

	int count = 0;
	if (codes[iNodalDisp])
	{
		if (fNumDOF > 3) throw eGeneralFail;
		const char* dlabels[3] = {"D_X", "D_Y", "D_Z"};
		for (int i = 0; i < fNumDOF; i++)
			labels[count++] = dlabels[i];
	}

	if (codes[iNodalCoord])
	{
		const char* xlabels[] = {"x1", "x2", "x3"};
		for (int i = 0; i < fNumSD; i++)
			labels[count++] = xlabels[i];
	}

	if (codes[iNodalStress])
	{
		const char* slabels2D[] = {"s11", "s22", "s12"};
		const char* slabels3D[] = {"s11", "s22", "s33", "s23", "s13", "s12"};
		const char**    slabels = (fNumSD == 2) ? slabels2D : slabels3D;
		for (int i = 0; i < dSymMatrixT::NumValues(fNumSD); i++)
			labels[count++] = slabels[i];
	}
		
	if (codes[iPrincipal])
	{
		const char* plabels[] = {"s1", "s2", "s3"};
		for (int i = 0; i < fNumSD; i++)
			labels[count++] = plabels[i];
	}
		
	if (codes[iEnergyDensity]) labels[count++] = "phi";
	if (codes[iWaveSpeeds   ])
	{
		const char* clabels2D[] = {"cd", "cs"};
		const char* clabels3D[] = {"cd", "cs_min", "cs_max"};
		const char**    clabels = (fNumSD == 2) ? clabels2D : clabels3D;
		for (int i = 0; i < fNumSD; i++)
			labels[count++] = clabels[i];		
	}

	/* material output labels */
	if (codes[iMaterialData])
	{
		ArrayT<StringT> matlabels;
		(*fMaterialList)[0]->OutputLabels(matlabels);	
		
		for (int i = 0; i < codes[iMaterialData]; i++)
			labels[count++] = matlabels[i];
	}
}
