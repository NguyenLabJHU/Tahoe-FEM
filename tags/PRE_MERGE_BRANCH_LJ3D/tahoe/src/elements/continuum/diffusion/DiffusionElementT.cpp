/* $Id: DiffusionElementT.cpp,v 1.20 2004-06-17 07:39:59 paklein Exp $ */
/* created: paklein (10/02/1999) */
#include "DiffusionElementT.h"

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include "toolboxConstants.h"

#include "ifstreamT.h"
#include "ElementCardT.h"
#include "ShapeFunctionT.h"
#include "eIntegratorT.h"
#include "iAutoArrayT.h"

/* materials */
#include "DiffusionMaterialT.h"
#include "DiffusionMatSupportT.h"
#include "DiffusionMatListT.h"

using namespace Tahoe;

/* initialize static data */
const int DiffusionElementT::NumOutputCodes = 3;

/* parameters */
const int kDiffusionNDOF = 1;

/* constructor */
DiffusionElementT::DiffusionElementT(const ElementSupportT& support, const FieldT& field):
	ContinuumElementT(support, field),
	fLocVel(LocalArrayT::kVel),
	fD(NumSD()),
	fq(NumSD()),
	fDiffusionMatSupport(NULL)
{
	SetName("diffusion");

	/* check base class initializations */
	if (NumDOF() != kDiffusionNDOF) {
		cout << "\n DiffusionElementT::DiffusionElementT: expecting field with " << kDiffusionNDOF << " dof/node not " 
		     << NumDOF() << endl;
		throw ExceptionT::kBadInputValue;
	}
}

DiffusionElementT::DiffusionElementT(const ElementSupportT& support):
	ContinuumElementT(support),
	fLocVel(LocalArrayT::kVel),
	fDiffusionMatSupport(NULL)
{
	SetName("diffusion");
}

/* destructor */
DiffusionElementT::~DiffusionElementT(void)
{
	delete fDiffusionMatSupport;
}

/* data initialization */
void DiffusionElementT::Initialize(void)
{
	/* inherited */
	ContinuumElementT::Initialize();

	/* allocate */
	fB.Dimension(NumSD(), NumElementNodes());
	fGradient_list.Dimension(NumIP());
	for (int i = 0; i < fGradient_list.Length(); i++)
		fGradient_list[i].Dimension(NumSD());

	/* setup for material output */
	if (fNodalOutputCodes[iMaterialData])
	{
		/* check compatibility of output */
		if (!CheckMaterialOutput())
		{
			cout << "\n DiffusionElementT::Initialize: error with material output" << endl;
			throw ExceptionT::kBadInputValue;
		}
		/* no material output variables */
	 	else if ((*fMaterialList)[0]->NumOutputVariables() == 0)
		{
			cout << "\n DiffusionElementT::ReadMaterialData: there are no material outputs"
			     << endl;
			fNodalOutputCodes[iMaterialData] = 0;
		}
	}
}

/* TEMPORARY */
void DiffusionElementT::InitialCondition(void)
{
	/* inherited */
	ContinuumElementT::InitialCondition();
	
	/* set the source for the iteration number */
	fDiffusionMatSupport->SetIterationNumber(ElementSupport().IterationNumber(Group()));
}

/* compute nodal force */
void DiffusionElementT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
	//not implemented
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
}
	
/* returns the energy as defined by the derived class types */
double DiffusionElementT::InternalEnergy(void)
{
	double energy = 0.0;

	Top();
	while ( NextElement() )
	{
		/* shape function derivatives, jacobians, local coords */
		SetGlobalShape();
		
		/* get displacements */
		SetLocalU(fLocDisp);
		
		/* integration */
		const double* Det    = fShapes->IPDets();
		const double* Weight = fShapes->IPWeights();

		/* material properties */
		double heat_capacity = fCurrMaterial->Capacity();

		fShapes->TopIP();
		while ( fShapes->NextIP() )
		{
			/* ip value */
			fShapes->InterpolateU(fLocDisp, fDOFvec);
		
			/* accumulate */
			energy += heat_capacity*(*Det++)*(*Weight++)*fDOFvec[0];
		}
	}
	return energy;
}

void DiffusionElementT::SendOutput(int kincode)
{
	/* output flags */
	iArrayT flags(fNodalOutputCodes.Length());

	/* set flags to get desired output */
	flags = IOBaseT::kAtNever;
	switch (kincode)
	{
		case iNodalDisp:
		    flags[iNodalDisp] = NumDOF();
			break;
		default:
			cout << "\n DiffusionElementT::SendKinematic: invalid output code: ";
			cout << kincode << endl;
	}

	/* number of output values */
	iArrayT n_counts;
	SetNodalOutputCodes(IOBaseT::kAtInc, flags, n_counts);

	/* reset averaging workspace */
	ElementSupport().ResetAverage(n_counts.Sum());

	/* no element output */
	iArrayT e_counts(fElementOutputCodes.Length());
	e_counts = 0;

	/* generate nodal values */
	dArray2DT e_values, n_values;
	ComputeOutput(n_counts, n_values, e_counts, e_values);
}

/***********************************************************************
* Protected
***********************************************************************/

/* print element group data */
void DiffusionElementT::PrintControlData(ostream& out) const
{
	/* inherited */
	ContinuumElementT::PrintControlData(out);

	// anything else
}

void DiffusionElementT::EchoOutputCodes(ifstreamT& in, ostream& out)
{
	/* allocate */
	fNodalOutputCodes.Dimension(NumOutputCodes);

	/* read in at a time to allow comments */
	for (int i = 0; i < fNodalOutputCodes.Length(); i++)
	{
		in >> fNodalOutputCodes[i];
		
		/* convert all to "at print increment" */
		if (fNodalOutputCodes[i] != IOBaseT::kAtNever)
			fNodalOutputCodes[i] = IOBaseT::kAtInc;
	}		

	/* checks */
	if (fNodalOutputCodes.Min() < IOBaseT::kAtFail ||
	    fNodalOutputCodes.Max() > IOBaseT::kAtInc) throw ExceptionT::kBadInputValue;

	/* control parameters */
	out << " Number of nodal output codes. . . . . . . . . . = " << NumOutputCodes << '\n';
	out << "    [" << fNodalOutputCodes[iNodalCoord   ] << "]: initial nodal coordinates\n";
	out << "    [" << fNodalOutputCodes[iNodalDisp    ] << "]: nodal displacements\n";
	out << "    [" << fNodalOutputCodes[iMaterialData ] << "]: nodal material output parameters\n";
}

/* initialize local arrays */
void DiffusionElementT::SetLocalArrays(void)
{
	/* inherited */
	ContinuumElementT::SetLocalArrays();

	/* allocate */
	fLocVel.Dimension(NumElementNodes(), NumDOF());

	/* nodal velocities */
	if (fIntegrator->Order() == 1)
		Field().RegisterLocal(fLocVel);
}

/* construct output labels array */
void DiffusionElementT::SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
	/* initialize */
	counts.Dimension(flags.Length());
	counts = 0;

	if (flags[iNodalCoord] == mode)
		counts[iNodalCoord] = NumSD();
	if (flags[iNodalDisp] == mode)
		counts[iNodalDisp] = NumDOF();
	if (flags[iMaterialData] == mode)
		counts[iMaterialData ] = (*fMaterialList)[0]->NumOutputVariables();
}

void DiffusionElementT::SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
#pragma unused(mode)
#pragma unused(flags)
	if (counts.Sum() != 0)
	{
		cout << "\n DiffusionElementT::SetElementOutputCodes: not yet supported" << endl;
		throw ExceptionT::kBadInputValue;
	}
}

/* set the correct shape functions */
void DiffusionElementT::SetShape(void)
{
	fShapes = new ShapeFunctionT(GeometryCode(), NumIP(), fLocInitCoords);
	if (!fShapes ) throw ExceptionT::kOutOfMemory;
	fShapes->Initialize();
}

/* construct the effective mass matrix */
void DiffusionElementT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
	/* inherited */
	ContinuumElementT::LHSDriver(sys_type);

	/* set components and weights */
	double constC = 0.0;
	double constK = 0.0;
	
	int formC = fIntegrator->FormC(constC);
	int formK = fIntegrator->FormK(constK);

	/* loop over elements */
	Top();
	while (NextElement())
	{
		/* initialize */
		fLHS = 0.0;
		
		/* set shape function derivatives */
		SetGlobalShape();

		/* element mass */
		if (formC) FormMass(kConsistentMass, constC*(fCurrMaterial->Capacity()));

		/* element stiffness */
		if (formK) FormStiffness(constK);
	
		/* add to global equations */
		AssembleLHS();		
	}
}

void DiffusionElementT::RHSDriver(void)
{
	/* inherited */
	ContinuumElementT::RHSDriver();

	/* set components and weights */
	double constCv = 0.0;
	double constKd = 0.0;
	
	/* components dicated by the algorithm */
	int formCv = fIntegrator->FormCv(constCv);
	int formKd = fIntegrator->FormKd(constKd);

	/* body forces */
	int formBody = 0;
	if (fBodySchedule && fBody.Magnitude() > kSmall)
	{	
		formBody = 1;
		if (!formCv) constCv = 1.0; // correct value ??
	}

	/* block info - needed for source terms */
	int block_dex = 0;
	const ElementBlockDataT* block_data = fBlockData.Pointer(block_dex);
	const dArray2DT* block_source = Field().Source(block_data->ID());
	dArray2DT ip_source;
	if (block_source) ip_source.Dimension(NumIP(), 1);
	int block_count = 0;

	double dt = ElementSupport().TimeStep();
	Top();
	while (NextElement())
	{
		/* reset block info (skip empty) */
		while (block_count == block_data->Dimension()) {
			block_data = fBlockData.Pointer(++block_dex);
			block_source = Field().Source(block_data->ID());
			block_count = 0;
		}
		
		/* convert heat increment/volume to rate */
		if (block_source) {
			block_source->RowCopy(block_count, ip_source);
			if (fabs(dt) > kSmall)
				ip_source /= dt;
			else /* for dt -> 0 */
				ip_source = 0.0;
		}
		block_count++;
		
		/* initialize */
		fRHS = 0.0;

		/* global shape function values */
		SetGlobalShape();

		/* conduction term */
		if (formKd) 
		{
			SetLocalU(fLocDisp);
			FormKd(-constKd);
		}

		/* capacity term */
		if (formCv || formBody)
		{
			if (formCv) SetLocalU(fLocVel);
			else fLocVel = 0.0;
			if (formBody) AddBodyForce(fLocVel);

			FormMa(kConsistentMass, -constCv*fCurrMaterial->Capacity(), 
				&fLocVel,
				(block_source) ? &ip_source : NULL);			  		
		}
				
		/* assemble */
		AssembleRHS();
	}
}

/* set the \e B matrix at the specified integration point */
void DiffusionElementT::B(int ip, dMatrixT& B_matrix) const
{
	const dArray2DT& DNa = fShapes->Derivatives_U(ip);
	int nnd = DNa.MinorDim();
	double* pB = B_matrix.Pointer();

	/* 2D */
	if (DNa.MajorDim() == 2)
	{
		const double* pNax = DNa(0);
		const double* pNay = DNa(1);

		for (int i = 0; i < nnd; i++)
		{
			*pB++ = *pNax++;
			*pB++ = *pNay++;
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
			*pB++ = *pNax++;
			*pB++ = *pNay++;
			*pB++ = *pNaz++;
		}
	}
}

/* current element operations */
bool DiffusionElementT::NextElement(void)
{
	/* inherited */
	bool result = ContinuumElementT::NextElement();
	
	/* get material pointer */
	if (result)
	{
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[CurrentElement().MaterialNumber()];
	
		/* cast is safe since class contructs materials list */
		fCurrMaterial = (DiffusionMaterialT*) pcont_mat;
	}
	
	return result;
}

/* form the element stiffness matrix */
void DiffusionElementT::FormStiffness(double constK)
{
	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	/* integration parameters */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	
	/* integrate element stiffness */
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		double scale = constK*(*Det++)*(*Weight++);
	
		/* strain displacement matrix */
		B(fShapes->CurrIP(), fB);

		/* get D matrix */
		fD.SetToScaled(scale, fCurrMaterial->k_ij());
							
		/* multiply b(transpose) * db, taking account of symmetry, */
		/* and accumulate in elstif */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}
}

/* calculate the internal force contribution ("-k*d") */
void DiffusionElementT::FormKd(double constK)
{
	/* integration parameters */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	
	int nsd = NumSD();
	dMatrixT grad;
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		/* set field gradient */
		grad.Set(1, nsd, fGradient_list[CurrIP()].Pointer());
		IP_ComputeGradient(fLocDisp, grad);

		/* get strain-displacement matrix */
		B(fShapes->CurrIP(), fB);

		/* compute heat flow */
		fB.MultTx(fCurrMaterial->q_i(), fNEEvec);

		/* accumulate */
		fRHS.AddScaled(-constK*(*Weight++)*(*Det++), fNEEvec);
	}	
}

/* construct a new material support and return a pointer */
MaterialSupportT* DiffusionElementT::NewMaterialSupport(MaterialSupportT* p) const
{
	/* allocate */
	if (!p) p = new DiffusionMatSupportT(NumSD(), NumDOF(), NumIP());

	/* inherited initializations */
	ContinuumElementT::NewMaterialSupport(p);
	
	/* set DiffusionMatSupportT fields */
	DiffusionMatSupportT* ps = TB_DYNAMIC_CAST(DiffusionMatSupportT*, p);
	if (ps) {
		ps->SetContinuumElement(this);
		ps->SetGradient(&fGradient_list);
	}

	return p;
}

/* return a pointer to a new material list */
MaterialListT* DiffusionElementT::NewMaterialList(int nsd, int size)
{
#pragma unused(nsd)
	if (size > 0)
	{
		/* material support */
		if (!fDiffusionMatSupport) {
			fDiffusionMatSupport = TB_DYNAMIC_CAST(DiffusionMatSupportT*, NewMaterialSupport());
			if (!fDiffusionMatSupport) ExceptionT::GeneralFail("DiffusionElementT::NewMaterialList");
		}

		/* allocate */
		return new DiffusionMatListT(size, *fDiffusionMatSupport);
	}
	else
		return new DiffusionMatListT;
}

/* driver for calculating output values */
void DiffusionElementT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	const iArrayT& e_codes, dArray2DT& e_values)
{
	/* number of output values */
	int n_out = n_codes.Sum();
	int e_out = e_codes.Sum();

	/* nothing to output */
	if (n_out == 0 && e_out == 0) return;

//TEMP
#pragma unused(e_values)
if (e_out > 0)
	ExceptionT::GeneralFail("DiffusionElementT::ComputeOutput", "element output not supported");

	/* dimensions */
	int nen = NumElementNodes();
	int nsd = NumSD();

	/* reset averaging workspace */
	ElementSupport().ResetAverage(n_out);

	/* work arrays */
	dArray2DT nodal_space(nen, n_out);
	dArray2DT nodal_all(nen, n_out);
	dArray2DT coords, disp;
	dArray2DT nodalstress, princstress, matdat;
	dArray2DT energy, speed;

	/* ip values */
	dSymMatrixT cauchy(nsd);
	dArrayT ipmat(n_codes[iMaterialData]), ipenergy(1);
	dArrayT ipspeed(nsd), ipprincipal(nsd);

	/* set shallow copies */
	double* pall = nodal_space.Pointer();
	coords.Set(nen, n_codes[iNodalCoord], pall);
	pall += coords.Length();
	disp.Set(nen, n_codes[iNodalDisp], pall);
	pall += disp.Length();
	matdat.Set(nen, n_codes[iMaterialData], pall);

	Top();
	while (NextElement())
	{
		/* initialize */
	    nodal_space = 0.0;

		/* global shape function values */
		SetGlobalShape();
		SetLocalU(fLocDisp);
		
		/* coordinates and displacements all at once */
		if (n_codes[iNodalCoord]) fLocInitCoords.ReturnTranspose(coords);
		if (n_codes[iNodalDisp])  fLocDisp.ReturnTranspose(disp);

		/* integrate */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/* material stuff */
			if (n_codes[iMaterialData])
			{
				fCurrMaterial->ComputeOutput(ipmat);
				fShapes->Extrapolate(ipmat,matdat);
			}
		}

		/* copy in the cols (in sequence of output) */
		int colcount = 0;
		nodal_all.BlockColumnCopyAt(disp  , colcount); colcount += disp.MinorDim();
		nodal_all.BlockColumnCopyAt(coords, colcount); colcount += coords.MinorDim();
		nodal_all.BlockColumnCopyAt(matdat, colcount);

		/* accumulate - extrapolation done from ip's to corners => X nodes */
		ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_all);
	}
	
	/* get nodally averaged values */
	ElementSupport().OutputUsedAverage(n_values);
}

/* information about subordinate parameter lists */
void DiffusionElementT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ContinuumElementT::DefineSubs(sub_list);

	sub_list.AddSub("diffusion_materials");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* DiffusionElementT::NewSub(const StringT& list_name) const
{
	if (list_name == "diffusion_materials")
	{
		/* non-const this */
		DiffusionElementT* non_const_this = const_cast<DiffusionElementT*>(this);
		return non_const_this->NewMaterialList(0,0);
	}
	else /* inherited */
		return ContinuumElementT::NewSub(list_name);
}

/***********************************************************************
* Private
***********************************************************************/

/* construct output labels array */
void DiffusionElementT::GenerateOutputLabels(const iArrayT& n_codes,
	ArrayT<StringT>& n_labels, const iArrayT& e_codes, 
	ArrayT<StringT>& e_labels) const
{
//TEMP - no element labels for now
#pragma unused(e_labels)

	/* allocate node labels */
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

	/* material output labels */
	if (n_codes[iMaterialData])
	{
		ArrayT<StringT> matlabels;
		(*fMaterialList)[0]->OutputLabels(matlabels);	
		
		for (int i = 0; i < n_codes[iMaterialData]; i++)
			n_labels[count++] = matlabels[i];
	}
	
	if (e_codes.Sum() != 0)
		ExceptionT::GeneralFail("DiffusionElementT::GenerateOutputLabels", 
			"not expecting any element output codes");
}
