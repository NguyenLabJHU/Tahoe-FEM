/* $Id: ContinuumElementT.cpp,v 1.52 2006-03-28 18:57:19 regueiro Exp $ */
/* created: paklein (10/22/1996) */
#include "ContinuumElementT.h"

#include <iostream.h>
#include <iomanip.h>

#include "ifstreamT.h"
#include "ModelManagerT.h"
#include "SolidMaterialT.h"
#include "ShapeFunctionT.h"
#include "eIntegratorT.h"
#include "Traction_CardT.h"
#include "iAutoArrayT.h"
#include "OutputSetT.h"
#include "ScheduleT.h"
#include "ParameterContainerT.h"
#include "CommunicatorT.h"

//TEMP: all this for general traction BC implementation?
#include "VariArrayT.h"
#include "nVariArray2DT.h"
#include "VariLocalArrayT.h"

/* materials lists */
#include "MaterialSupportT.h"
#include "MaterialListT.h"

#ifdef __DEVELOPMENT_COLORADO__
#include "DevelopmentColoradoElementsConfig.h"
#endif

const double Pi = acos(-1.0);

using namespace Tahoe;

/* constructor */
ContinuumElementT::ContinuumElementT(const ElementSupportT& support):
	ElementBaseT(support),
	fGroupCommunicator(NULL),
	fMaterialList(NULL),
	fBodySchedule(NULL),
	fTractionBCSet(0),
	fShapes(NULL),
	fStoreShape(false),
	fLocInitCoords(LocalArrayT::kInitCoords),
	fLocDisp(LocalArrayT::kDisp),
	fNumIP(0),
	fGeometryCode(GeometryT::kNone)
{
	SetName("continuum_element");
}

/* destructor */
ContinuumElementT::~ContinuumElementT(void)
{
	delete fGroupCommunicator;
	delete fMaterialList;
	delete fShapes;
}

/* accessors */
const int& ContinuumElementT::CurrIP(void) const
{
	return ShapeFunction().CurrIP();
}

/* the coordinates of the current integration points */
void ContinuumElementT::IP_Coords(dArrayT& ip_coords) const
{
	/* computed by shape functions */
	ShapeFunction().IPCoords(ip_coords);
}

/* interpolate the nodal field values to the current integration point */
void ContinuumElementT::IP_Interpolate(const LocalArrayT& nodal_u, dArrayT& ip_u) const
{
    /* computed by shape functions */
    ShapeFunction().InterpolateU(nodal_u, ip_u);
}

void ContinuumElementT::IP_Interpolate(const LocalArrayT& nodal_u, dArrayT& ip_u, int ip) const
{
    /* computed by shape functions */
    ShapeFunction().InterpolateU(nodal_u, ip_u, ip);
}

/* field gradients */
void ContinuumElementT::IP_ComputeGradient(const LocalArrayT& field, 
	dMatrixT& gradient) const
{
	/* computed by shape functions */
	ShapeFunction().GradU(field, gradient);
}

void ContinuumElementT::IP_ComputeGradient(const LocalArrayT& field, 
	dMatrixT& gradient, int ip) const
{
	/* computed by shape functions */
	ShapeFunction().GradU(field, gradient, ip);
}

/* extrapolate all integration point values to the nodes */
void ContinuumElementT::IP_ExtrapolateAll(const dArrayT& ip_values, 
	dArrayT& nodal_values) const
{
	/* computed by shape functions */
	ShapeFunction().ExtrapolateAll(ip_values, nodal_values);
}

void ContinuumElementT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	/* inherited */
	ElementBaseT::Equations(eq_1, eq_2);

	/* mark traction BC data as old */
	fTractionBCSet = 0;
}

/* form of tangent matrix */
GlobalT::SystemTypeT ContinuumElementT::TangentType(void) const
{
	/* initialize to lowest precedence */
	GlobalT::SystemTypeT type = GlobalT::kDiagonal;

	if (fMaterialList) 
	{
		for (int i = 0; i < fMaterialList->Length(); i++)
		{
			GlobalT::SystemTypeT e_type = (*fMaterialList)[i]->TangentType();
	
			/* using type precedence */
			type = (e_type > type) ? e_type : type;
		}
	}
	
	return type;
}

/* initialize/finalize step */
void ContinuumElementT::InitStep(void)
{
	/* inherited */
	ElementBaseT::InitStep();

	/* set material variables */
	if (fMaterialList)  fMaterialList->InitStep();
}

/* initialize/finalize step */
void ContinuumElementT::CloseStep(void)
{
	/* inherited */
	ElementBaseT::CloseStep();

	if (fMaterialList) 
	{
		/* set material variables */
		fMaterialList->CloseStep();

		/* update element level internal variables */
		if (fMaterialList->HasHistoryMaterials())
		{
			Top();
			while (NextElement())
			{
				const ElementCardT& element = CurrentElement();
				if (element.IsAllocated())
				{
					ContinuumMaterialT* pmat = (*fMaterialList)[element.MaterialNumber()];

				/* material update function */
					pmat->UpdateHistory();
				}
			}
		}
	}
}

/* resets to the last converged solution */
GlobalT::RelaxCodeT ContinuumElementT::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::ResetStep();

	/* update material internal variables */
	if (fMaterialList && fMaterialList->HasHistoryMaterials())
	{
		Top();
		while (NextElement())
		{
			const ElementCardT& element = CurrentElement();		
			if (element.IsAllocated())
			{
				ContinuumMaterialT* pmat = (*fMaterialList)[element.MaterialNumber()];

				/* material reset function */
				pmat->ResetHistory();
			}
		}
	}

	return relax;
}

/* element level reconfiguration for the current time increment */
GlobalT::RelaxCodeT ContinuumElementT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	/* loop over materials */
	for (int i = 0; i < fMaterialList->Length(); i++)
		relax = GlobalT::MaxPrecedence((*fMaterialList)[i]->RelaxCode(), relax);

	return relax;
}

/* restart operations */
void ContinuumElementT::ReadRestart(istream& in)
{
	/* inherited */
	ElementBaseT::ReadRestart(in);

	/* update element level internal variables */
	if (fMaterialList && fMaterialList->HasHistoryMaterials())
	{
		for (int i = 0; i < fElementCards.Length(); i++)
		{
			int isallocated;
			in >> isallocated;
			if (isallocated) fElementCards[i].ReadRestart(in);
		}
	}
}

void ContinuumElementT::WriteRestart(ostream& out) const
{
	/* inherited */
	ElementBaseT::WriteRestart(out);

	/* update element level internal variables */
	if (fMaterialList && fMaterialList->HasHistoryMaterials())
	{
		for (int i = 0; i < fElementCards.Length(); i++)
		{
			const ElementCardT& element = fElementCards[i];
			out << element.IsAllocated() << '\n';
			if (element.IsAllocated()) element.WriteRestart(out);
		}
	}
}

/* writing output */
void ContinuumElementT::RegisterOutput(void)
{
//NOTE: could loop over each output mode and register
//      it with the output separately. for now just register
//      "kAtInc"
	
	/* nodal output */
	iArrayT n_counts;
	SetNodalOutputCodes(IOBaseT::kAtInc, fNodalOutputCodes, n_counts);

	/* element output */
	iArrayT e_counts;
	SetElementOutputCodes(IOBaseT::kAtInc, fElementOutputCodes, e_counts);

	/* inherited */
	if (n_counts.Length() == 0 && e_counts.Length() == 0)
		ElementBaseT::RegisterOutput();
	else
	{
		/* collect variable labels */
		ArrayT<StringT> n_labels(n_counts.Sum());
		ArrayT<StringT> e_labels(e_counts.Sum());
		GenerateOutputLabels(n_counts, n_labels, e_counts, e_labels);

		/* block ID's used by the group */
		ArrayT<StringT> block_ID(fBlockData.Length());
		for (int i = 0; i < block_ID.Length(); i++)
			block_ID[i] = fBlockData[i].ID();

		/* set output specifier */
		OutputSetT output_set(fGeometryCode, block_ID, fConnectivities, n_labels, e_labels, false);
		
		/* register and get output ID */
		fOutputID = ElementSupport().RegisterOutput(output_set);
	}
}

//NOTE - this function is/was identical to CSEBaseT::WriteOutput
void ContinuumElementT::WriteOutput(void)
{
	/* regular output */
	IOBaseT::OutputModeT mode = IOBaseT::kAtInc;

	/* map output flags to count of values */
	iArrayT n_counts;
	SetNodalOutputCodes(mode, fNodalOutputCodes, n_counts);
	iArrayT e_counts;
	SetElementOutputCodes(mode, fElementOutputCodes, e_counts);

	/* inherited */
	if (n_counts.Length() == 0 && e_counts.Length() == 0)
		ElementBaseT::WriteOutput();
	else
	{
		/* calculate output values */
		dArray2DT n_values;
		dArray2DT e_values;
		ComputeOutput(n_counts, n_values, e_counts, e_values);

		/* send to output */
		ElementSupport().WriteOutput(fOutputID, n_values, e_values);
	}
}

/* resolve the output variable label into the output code and offset within the output */
void ContinuumElementT::ResolveOutputVariable(const StringT& variable, int& code, int& offset)
{
	/* search output labels */
	code = -1;
	offset = -1;
	iArrayT e_counts(fElementOutputCodes.Length());
	e_counts = 0;
	iArrayT n_codes(fNodalOutputCodes.Length());
	for (int i = 0; code == -1 && i < n_codes.Length(); i++)
	{
		ArrayT<StringT> n_labels, e_labels;
		n_codes = 0;
		n_codes[i] = 1;
		
		iArrayT n_counts;
		SetNodalOutputCodes(IOBaseT::kAtInc, n_codes, n_counts);
		GenerateOutputLabels(n_counts, n_labels, e_counts, e_labels);
		
		for (int j = 0; offset == -1 && j < n_labels.Length(); j++)
			if (n_labels[j] == variable) /* found */ {
				code = i;
				offset = j;
			}
	}
	
	/* inherited */
	if (code == -1 || offset == -1)
		ElementBaseT::ResolveOutputVariable(variable, code, offset);
}

/* return geometry and number of nodes on each facet */
void ContinuumElementT::FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geometry, 
	iArrayT& num_facet_nodes) const
{
	/* from integration domain */
	ShapeFunction().FacetGeometry(facet_geometry, num_facet_nodes);
}

void ContinuumElementT::SetStatus(const ArrayT<ElementCardT::StatusT>& status)
{
  /* work space */
  dArrayT state;
  dArrayT t_in;

  /* loop over elements and initial state variables */
  int elem_num = 0;
  Top();
  while (NextElement())
    {
      /* current element */
      ElementCardT::StatusT& flag = CurrentElement().Flag();
      flag = status[elem_num++];
      /* material pointer */
      ContinuumMaterialT* pmat = (*fMaterialList)[CurrentElement().MaterialNumber()];

      if (flag == ElementCardT::kMarkON){
	if (pmat->NeedsPointInitialization()){
	  /* global shape function values */
	  SetGlobalShape();

	  fShapes->TopIP();
	  while (fShapes->NextIP())
	    pmat->PointInitialize();
	}
	flag = ElementCardT::kON;
      }
      else if (flag == ElementCardT::kMarkOFF)
	flag = ElementCardT::kOFF;
    }
}

/* initial condition/restart functions (per time sequence) */
void ContinuumElementT::InitialCondition(void)
{
	/* inherited */
	ElementBaseT::InitialCondition();

	/* check for initialization materials */
	bool need_init = false;
	if (fMaterialList)
		for (int i = 0; i < fMaterialList->Length() && !need_init; i++)
			need_init = (*fMaterialList)[i]->NeedsPointInitialization();

	/* need to run through elements */
	if (fStoreShape || need_init)
	{
		/* initialize storage */
		if (fStoreShape) fShapes->InitStore(NumElements(), &(fElementCards.Position()));

		/* loop over elements */
		Top();
		while (NextElement())
		{
			/* compute shape function derivarives */
			SetGlobalShape();

			/* store */
			if (fStoreShape) fShapes->Store();
		
			/* initialize material */
			if (need_init) 
			{
				/* material pointer */
				ContinuumMaterialT* pmat = (*fMaterialList)[CurrentElement().MaterialNumber()];
				if (pmat->NeedsPointInitialization())
				{
					/* loop over integration points */
					fShapes->TopIP();
					while (fShapes->NextIP())
						pmat->PointInitialize();
				}
			}
		}
		
		/* finalize storage */
		if (fStoreShape) fShapes->CloseStore();
	}
}

ContinuumElementT::MassTypeT ContinuumElementT::int2MassTypeT(int i)
{
	if (i == kNoMass)
		return kNoMass;
	else if (i == kConsistentMass)
		return kConsistentMass;
	else if (i == kLumpedMass)
		return kLumpedMass;
	else if (i == kAutomaticMass)
		return kAutomaticMass;
	else
		ExceptionT::GeneralFail("ContinuumElementT::int2MassTypeT", 
			"could not translate %d", i);	
	return kNoMass;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

void ContinuumElementT::SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
#pragma unused(mode)
#pragma unused(flags)
	counts.Dimension(0);
}

void ContinuumElementT::SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
#pragma unused(mode)
#pragma unused(flags)
	counts.Dimension(0);
}

void ContinuumElementT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	const iArrayT& e_codes, dArray2DT& e_values)
{
#pragma unused(n_codes)
#pragma unused(e_codes)
	n_values.Dimension(0, n_values.MinorDim());
	e_values.Dimension(0, e_values.MinorDim());
}

void ContinuumElementT::GenerateOutputLabels( const iArrayT& n_codes, 
	ArrayT<StringT>& n_labels, const iArrayT& e_codes, ArrayT<StringT>& e_labels) const
{
#pragma unused(n_codes)
#pragma unused(e_codes)
	n_labels.Dimension(0);
	e_labels.Dimension(0);
}

/* initialize local arrays */
void ContinuumElementT::SetLocalArrays(void)
{
	/* dimension */
	fLocInitCoords.Dimension(NumElementNodes(), NumSD());
	fLocDisp.Dimension(NumElementNodes(), NumDOF());

	/* set source */
	ElementSupport().RegisterCoordinates(fLocInitCoords);
	Field().RegisterLocal(fLocDisp);	
}

/* form the residual force vector */
void ContinuumElementT::RHSDriver(void)
{
	/* contribution from tractions */
	ApplyTractionBC();

	// should call derived class function here	
}

/* compute contribution to RHS from traction BC's */
void ContinuumElementT::ApplyTractionBC(void)
{
	if (fTractionList.Length() > 0)
	{		
		/* dimensions */
		int nsd = NumSD();
		int ndof = NumDOF();
		bool is_axi = Axisymmetric();
		if (is_axi && nsd != 2) ExceptionT::GeneralFail();
	
		/* update equation numbers */
		if (!fTractionBCSet) SetTractionBC();
	
		/* force vector */
		dArrayT rhs;
		VariArrayT<double> rhs_man(25, rhs);
		
		/* local coordinates */
		LocalArrayT coords(LocalArrayT::kInitCoords);
		VariLocalArrayT coord_man(25, coords, nsd);
		ElementSupport().RegisterCoordinates(coords);
		
		/* nodal tractions */
		LocalArrayT tract(LocalArrayT::kUnspecified);
		VariLocalArrayT tract_man(25, tract, ndof);

		/* integration point tractions */
		dArray2DT ip_tract;
		nVariArray2DT<double> ip_tract_man(25, ip_tract, ndof);
		dArrayT tract_loc, tract_glb(ndof);
		dMatrixT Q(ndof);
		dArrayT ip_coords(2);
		
		/* Jacobian of the surface mapping */
		dMatrixT jacobian(nsd, nsd-1);
		double Pi2 = 2.0*Pi;
		for (int i = 0; i < fTractionList.Length(); i++)
		{
			const Traction_CardT& BC_card = fTractionList[i];

			/* dimension */
			const iArrayT& nodes = BC_card.Nodes();
			int nnd = nodes.Length();
			rhs_man.SetLength(nnd*ndof, false);
			coord_man.SetNumberOfNodes(nnd);
			tract_man.SetNumberOfNodes(nnd);
			
			/* local coordinates */
			coords.SetLocal(nodes);

			/* nodal traction vectors: (ndof x nnd) */
			BC_card.CurrentValue(tract);
			
			/* BC destination */
			int elem, facet;
			BC_card.Destination(elem, facet);
			
			/* boundary shape functions */
			const ParentDomainT& surf_shape = ShapeFunction().FacetShapeFunction(facet);
			int nip = surf_shape.NumIP();
			
			/* all ip tractions: (nip x ndof) */
			ip_tract_man.SetMajorDimension(nip, false);
			surf_shape.Interpolate(tract, ip_tract);

			/* traction vector coordinate system */
			if (BC_card.CoordSystem() == Traction_CardT::kCartesian)
			{
				/* integrate */			
				rhs = 0.0;
				const double* w = surf_shape.Weight();
				for (int j = 0; j < nip; j++)
				{
					/* coordinate mapping */
					surf_shape.DomainJacobian(coords, j, jacobian);
					double detj = surf_shape.SurfaceJacobian(jacobian);
	
					/* ip weight */
					double jwt = detj*w[j];
					if (is_axi) {
						surf_shape.Interpolate(coords, ip_coords, j);
						jwt *= Pi2*ip_coords[0];
					}
					
					/* ip traction */
					const double* tj = ip_tract(j);
					
					/* accumulate */
					for (int l = 0; l < ndof; l++)
					{
						/* nodal shape function */
						const double* Na = surf_shape.Shape(j);
					
						double* prhs = rhs.Pointer(l);
						double  fact = jwt*(*tj++);
						for (int k = 0; k < nnd; k++)
						{
							*prhs += fact*(*Na++);
							prhs += ndof;
						}
					}				
				}
			}
			else if (BC_card.CoordSystem() == Traction_CardT::kLocal)
			{
				/* integrate */			
				rhs = 0.0;
				const double* w = surf_shape.Weight();
				for (int j = 0; j < nip; j++)
				{
					/* coordinate mapping */
					surf_shape.DomainJacobian(coords, j, jacobian);
					double detj = surf_shape.SurfaceJacobian(jacobian, Q);
	
					/* ip weight */
					double jwt = detj*w[j];
					if (is_axi) {
						surf_shape.Interpolate(coords, ip_coords, j);
						jwt *= Pi2*ip_coords[0];
					}
					
					/* transform ip traction out of local frame */
					ip_tract.RowAlias(j, tract_loc);
					Q.Multx(tract_loc, tract_glb);

					/* ip traction */
					const double* tj = tract_glb.Pointer();
					
					/* accumulate */
					for (int l = 0; l < ndof; l++)
					{
						/* nodal shape function */
						const double* Na = surf_shape.Shape(j);
					
						double* prhs = rhs.Pointer(l);
						double  fact = jwt*(*tj++);
						for (int k = 0; k < nnd; k++)
						{
							*prhs += fact*(*Na++);
							prhs += ndof;
						}
					}				
				}
			}
			else
				throw ExceptionT::kGeneralFail;

			/* assemble */
			ElementSupport().AssembleRHS(Group(), rhs, BC_card.Eqnos());
		}
	}
}

/* form global shape function derivatives */
void ContinuumElementT::SetGlobalShape(void)
{
	/* fetch (initial) coordinates */
	SetLocalX(fLocInitCoords);
	
	/* compute shape function derivatives */
	fShapes->SetDerivatives();
}

/* form the element mass matrix */
void ContinuumElementT::FormMass(MassTypeT mass_type, double constM, bool axisymmetric, const double* ip_weight)
{
#if __option(extended_errorcheck)
	if (fLocDisp.Length() != fLHS.Rows()) ExceptionT::SizeMismatch("ContinuumElementT::FormMass");
#endif

	switch (mass_type)
	{
		case kNoMass:			/* no mass matrix */
		
			break;
		
		case kConsistentMass:	/* consistent mass	*/
		{
			// integration of the element mass is done
			// in the reference configuration since density
			// is mass/(undeformed volume)
			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();
			
			int nen = NumElementNodes();
			int nun = fLocDisp.NumberOfNodes();
			int ndof = NumDOF();
			
			/* matrix form */
			int a = 0, zero = 0;
			int& b_start = (fLHS.Format() == ElementMatrixT::kSymmetricUpper) ? a : zero;
			
			if (axisymmetric)
			{
				const LocalArrayT& coords = fShapes->Coordinates();
				fShapes->TopIP();	
				while ( fShapes->NextIP() )
				{
					/* compute radius */
					const double* NaX = fShapes->IPShapeX();
					const double* x_r = coords(0); /* r is x-coordinate */
					double r = 0.0;
					for (a = 0; a < nen; a++)
						r += (*NaX++)*(*x_r++);

					/* integration factor */				
					double temp = 2.0*Pi*r*constM*(*Weight++)*(*Det++);
					if (ip_weight) temp *= *ip_weight++;

					const double* Na = fShapes->IPShapeU();		
					for (a = 0; a < nun; a++)
						for (int i = 0; i < ndof; i++)
						{
							int p = a*ndof + i;
							
							/* upper triangle only */
							for (int b = b_start; b < nun; b++) //TEMP - interpolate at the same time?
								for (int j = 0; j < ndof; j++)
									if(i == j) {
										int q = b*ndof + j;
										fLHS(p,q) += temp*Na[a]*Na[b];
									}
						}
				}
			}			
			else /* not axisymmetric */
			{
				fShapes->TopIP();	
				while ( fShapes->NextIP() )
				{
					/* integration factor */
					double temp = constM*(*Weight++)*(*Det++);
					if (ip_weight) temp *= *ip_weight++;

					const double* Na = fShapes->IPShapeU();		
					for (a = 0; a < nun; a++)
						for (int i = 0; i < ndof; i++)
						{
							int p = a*ndof + i;
							
							/* upper triangle only */
							for (int b = b_start; b < nun; b++)
								for (int j = 0; j < ndof; j++)
									if(i == j) {
										int q = b*ndof + j;
										fLHS(p,q) += temp*Na[a]*Na[b];
									}
						}
				}
			}
			break;
		}

		case kLumpedMass:	/* lumped mass */
		{
			int nen = NumElementNodes();
			int nun = fLocDisp.NumberOfNodes();
			int ndof = NumDOF();

		    double dsum   = 0.0;
		    double totmas = 0.0;
		    fNEEvec = 0.0;

			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();

			/* total mass and diagonal sum */
			if (axisymmetric)
			{
				const LocalArrayT& coords = fShapes->Coordinates();
				fShapes->TopIP();
				while (fShapes->NextIP()) {

					/* compute radius */
					const double* NaX = fShapes->IPShapeX();
					const double* x_r = coords(0); /* r is x-coordinate */
					double r = 0.0;
					for (int a = 0; a < nen; a++)
						r += (*NaX++)*(*x_r++);

					/* integration factor */
					double temp1 = 2.0*Pi*r*constM*(*Weight++)*(*Det++);
					if (ip_weight) temp1 *= *ip_weight++;

					const double* Na = fShapes->IPShapeU();
					totmas += temp1;
					for (int lnd = 0; lnd < nun; lnd++) {
						double temp2 = temp1*Na[lnd]*Na[lnd];
						dsum += temp2;
						fNEEvec[lnd] += temp2;
					}
				}
			}
			else /* not axisymmetric */
			{
				fShapes->TopIP();
				while (fShapes->NextIP()) {

					/* integration factor */
					double temp1 = constM*(*Weight++)*(*Det++);
					if (ip_weight) temp1 *= *ip_weight++;

					const double* Na = fShapes->IPShapeU();
					totmas += temp1;
					for (int lnd = 0; lnd < nun; lnd++) {
						double temp2 = temp1*Na[lnd]*Na[lnd];
						dsum += temp2;
						fNEEvec[lnd] += temp2;
					}
				}
			}	
				
			/* scale diagonal to conserve total mass */
			double diagmass = totmas/dsum;
			
			/* lump mass onto diagonal */
			double* pmass = fLHS.Pointer();
			int inc = fLHS.Rows() + 1;
			for (int lnd = 0; lnd < nun; lnd++)
			{
				double temp = diagmass*fNEEvec[lnd];
				for (int ed = 0; ed < ndof; ed++)
				{
					*pmass += temp;
					pmass += inc;	
				}
			}
			break;
		}			
		default:
			ExceptionT::BadInputValue("ContinuumElementT::FormMass", "unknown mass matrix code");
	}
}

/* add contribution from the body force */
void ContinuumElementT::AddBodyForce(LocalArrayT& body_force) const
{
	if (fBodySchedule)
	{
		int ndof = NumDOF();
		int nen = body_force.NumberOfNodes();
		double loadfactor = fBodySchedule->Value();
		double* p = body_force.Pointer();

		for (int i = 0; i < ndof; i++)
		{
			double temp = -fBody[i]*loadfactor;
			for (int j = 0; j < nen; j++)
				*p++ = temp;
		}
	}
}

/* calculate the body force contribution */
void ContinuumElementT::FormMa(MassTypeT mass_type, double constM, bool axisymmetric,
	const LocalArrayT* nodal_values,
	const dArray2DT* ip_values,
	const double* ip_weight)
{
	const char caller[] = "ContinuumElementT::FormMa";

	/* quick exit */
	if (!nodal_values && !ip_values) return;

#if __option(extended_errorcheck)
	/* dimension checks */
	if (nodal_values && 
		fRHS.Length() != nodal_values->Length()) 
			ExceptionT::SizeMismatch(caller);

	if (ip_values &&
		(ip_values->MajorDim() != fShapes->NumIP() ||
		 ip_values->MinorDim() != NumDOF()))
			ExceptionT::SizeMismatch(caller);
#endif

	switch (mass_type)
	{
		case kConsistentMass:
		{
			int ndof = NumDOF();
			int  nen = NumElementNodes();
			int  nun = nodal_values->NumberOfNodes();

			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();

			if (axisymmetric)
			{
				const LocalArrayT& coords = fShapes->Coordinates();
				fShapes->TopIP();
				while (fShapes->NextIP())
				{
					/* compute radius */
					const double* NaX = fShapes->IPShapeX();
					const double* x_r = coords(0); /* r is x-coordinate */
					double r = 0.0;
					for (int a = 0; a < nen; a++)
						r += (*NaX++)*(*x_r++);
				
					/* interpolate nodal values to ip */
					if (nodal_values)
						fShapes->InterpolateU(*nodal_values, fDOFvec);
					
					/* ip sources */
					if (ip_values)
						fDOFvec -= (*ip_values)(fShapes->CurrIP());

					/* accumulate in element residual force vector */				
					double*	res      = fRHS.Pointer();
					const double* Na = fShapes->IPShapeU();
				
					/* integration factor */
					double temp = 2.0*Pi*r*constM*(*Weight++)*(*Det++);
					if (ip_weight) temp *= *ip_weight++;

					for (int lnd = 0; lnd < nun; lnd++)
					{
						double temp2 = temp*(*Na++);
						double* pacc = fDOFvec.Pointer();
						for (int dof = 0; dof < ndof; dof++)			
							*res++ += temp2*(*pacc++);
					}
				}
			}
			else /* not axisymmetric */
			{
				fShapes->TopIP();
				while (fShapes->NextIP())
				{					
					/* interpolate nodal values to ip */
					if (nodal_values)
						fShapes->InterpolateU(*nodal_values, fDOFvec);
					
					/* ip sources */
					if (ip_values)
						fDOFvec -= (*ip_values)(fShapes->CurrIP());

					/* accumulate in element residual force vector */				
					double*	res      = fRHS.Pointer();
					const double* Na = fShapes->IPShapeU();
				
					/* integration factor */
					double temp = constM*(*Weight++)*(*Det++);
					if (ip_weight) temp *= *ip_weight++;
					
					for (int lnd = 0; lnd < nun; lnd++)
					{
						double temp2 = temp*(*Na++);
						double* pacc = fDOFvec.Pointer();
						for (int dof = 0; dof < ndof; dof++)			
							*res++ += temp2*(*pacc++);
					}
				}
			}
			break;
		}	
		case kLumpedMass:
		{
			fLHS = 0.0; //hope there's nothing in there!
			FormMass(kLumpedMass, constM, axisymmetric,ip_weight);

			/* init nodal values */
			if (nodal_values)
				nodal_values->ReturnTranspose(fNEEvec);
			else {
				ExceptionT::GeneralFail(caller, "expecting nodal values for lumped mass");
			}
				
//TEMP - what to do with ip values?
if (ip_values)
	ExceptionT::GeneralFail(caller, "lumped mass not implemented for ip sources");

			double* pAcc = fNEEvec.Pointer();
			double* pRes = fRHS.Pointer();
			int     massdex = 0;
			
			int nee = nodal_values->Length();
			for (int i = 0; i < nee; i++)
			{
				*pRes++ += (*pAcc++)*fLHS(massdex,massdex);
				massdex++;
			}
			
			break;
		}
	}
}

/* extract natural boundary condition information */
void ContinuumElementT::TakeNaturalBC(const ParameterListT& list)
{
	const char caller[] = "ContinuumElementT::TakeTractionBC";

	int num_natural_bc = list.NumLists("natural_bc");
	if (num_natural_bc > 0)
	{
		/* model manager */
		ModelManagerT& model = ElementSupport().ModelManager();
	
		/* temp space */
		ArrayT<StringT> block_ID(num_natural_bc);
	    ArrayT<iArray2DT> localsides(num_natural_bc);
	    iArrayT LTf(num_natural_bc);
	    ArrayT<Traction_CardT::CoordSystemT> coord_sys(num_natural_bc);
	    ArrayT<dArray2DT> values(num_natural_bc);

	    /* nodes on element facets */
	    iArrayT num_facet_nodes;
	    fShapes->NumNodesOnFacets(num_facet_nodes);
	    
	    /* loop over natural BC's */
	    int tot_num_sides = 0;
	    for (int i = 0; i < num_natural_bc; i++) 
	   	{
	    	const ParameterListT& natural_bc = list.GetList("natural_bc", i);
	    
	    	/* side set */
	    	const StringT& ss_ID = natural_bc.GetParameter("side_set_ID");
			localsides[i] = model.SideSet(ss_ID);
			int num_sides = localsides[i].MajorDim();
			tot_num_sides += num_sides;
			if (num_sides > 0)
			{
				block_ID[i] = model.SideSetGroupID(ss_ID);
				LTf[i] = natural_bc.GetParameter("schedule");
				coord_sys[i] = Traction_CardT::int2CoordSystemT(natural_bc.GetParameter("coordinate_system"));

				/* switch to elements numbering within the group */
				iArray2DT& side_set = localsides[i];
				iArrayT elems(num_sides);
				side_set.ColumnCopy(0, elems);
				BlockToGroupElementNumbers(elems, block_ID[i]);
				side_set.SetColumn(0, elems);

				/* all facets in set must have the same number of nodes */
				int num_nodes = num_facet_nodes[side_set(0,1)];
				for (int f = 0; f < num_sides; f++)
					if (num_facet_nodes[side_set(f,1)] != num_nodes)
						ExceptionT::BadInputValue(caller, "faces side set \"%s\" have different numbers of nodes",
							ss_ID.Pointer());

				/* read traction nodal values */
				dArray2DT& nodal_values = values[i];
				nodal_values.Dimension(num_nodes, NumDOF());
				int num_traction_vectors = natural_bc.NumLists("DoubleList");
				if (num_traction_vectors != 1 && num_traction_vectors != num_nodes)
					ExceptionT::GeneralFail(caller, "expecting 1 or %d vectors not %d",
						num_nodes, num_traction_vectors);
						
				/* constant over the face */
				if (num_traction_vectors == 1) {
					const ParameterListT& traction_vector = natural_bc.GetList("DoubleList");
					int dim = traction_vector.NumLists("Double");
					if (dim != NumDOF())
						ExceptionT::GeneralFail(caller, "expecting traction vector length %d not %d",
							NumDOF(), dim);

					/* same for all face nodes */
					for (int f = 0; f < NumDOF(); f++) {
						double t = traction_vector.GetList("Double", f).GetParameter("value");
						nodal_values.SetColumn(f, t);
					}
				}
				else
				{
					/* read separate vector for each face node */
					dArrayT t;
					for (int f = 0; f < num_nodes; f++) {
						const ParameterListT& traction_vector = natural_bc.GetList("DoubleList", f);
					int dim = traction_vector.NumLists("Double");
						if (dim != NumDOF())
							ExceptionT::GeneralFail(caller, "expecting traction vector length %d not %d",
								NumDOF(), dim);

						nodal_values.RowAlias(f, t);
						for (int j = 0; j < NumDOF(); j++)
							t[j] = traction_vector.GetList("Double", j).GetParameter("value");
					}
				}
			}
	    }
#pragma message("OK with empty side sets?")

		/* allocate all traction BC cards */
	    fTractionList.Dimension(tot_num_sides);

	    /* correct numbering offset */
	    LTf--;

		/* define traction cards */
		if (tot_num_sides > 0)
		{
			iArrayT loc_node_nums;
			int dex = 0;
			for (int i = 0; i < num_natural_bc; i++)
			{
				/* set traction BC cards */
				iArray2DT& side_set = localsides[i];
				int num_sides = side_set.MajorDim();
				for (int j = 0; j < num_sides; j++)
				{					
					/* get facet local node numbers */
					fShapes->NodesOnFacet(side_set(j, 1), loc_node_nums);
					
					/* set and echo */
					fTractionList[dex++].SetValues(ElementSupport(), side_set(j,0), side_set (j,1), LTf[i],
						 coord_sys[i], loc_node_nums, values[i]);
				}
			}
		}

		/* check coordinate system specifications */
		if (NumSD() != NumDOF())
			for (int i = 0; i < fTractionList.Length(); i++)
				if (fTractionList[i].CoordSystem() != Traction_CardT::kCartesian)
					ExceptionT::BadInputValue(caller, "coordinate system must be Cartesian if (nsd != ndof) for card %d", i+1);
	}
}

/* return a pointer to a new material list */
MaterialListT* ContinuumElementT::NewMaterialList(const StringT& name, int size)
{
#pragma unused(name)
#pragma unused(size)
	return NULL;
}

/* construct a new material support and return a pointer */
MaterialSupportT* ContinuumElementT::NewMaterialSupport(MaterialSupportT* p) const
{
	if (!p) p = new MaterialSupportT(NumDOF(), NumIP());

	/* ContinuumElementT sources */
	p->SetContinuumElement(this);
	p->SetElementCards(const_cast<AutoArrayT<ElementCardT>* >(&fElementCards));
	p->SetCurrIP(CurrIP());
	p->SetGroup(Group());

	/* ElementSupportT sources */
//	const ElementSupportT& e_support = ElementSupport();
//	p->SetRunState(e_support.RunState());
//	p->SetStepNumber(e_support.StepNumber());
//	p->SetIterationNumber(e_support.IterationNumber(Group()));
//TEMP - solvers not set up yet. For now, the source for the iteration number will
//       be set in the InitialCondition call for the subclass.
//	p->SetTime(e_support.Time());                              
//	p->SetTimeStep(e_support.TimeStep());
//	p->SetNumberOfSteps(e_support.NumberOfSteps());

	/* set pointer to local array */
	p->SetLocalArray(fLocDisp);

	return p;
}

/* write all current element information to the stream */
void ContinuumElementT::CurrElementInfo(ostream& out) const
{
	/* inherited */
	ElementBaseT::CurrElementInfo(out);
	dArray2DT temp;
	temp.Dimension(fLocInitCoords.NumberOfNodes(), fLocInitCoords.MinorDim());
	
	out <<   " initial coords:\n";
	temp.Dimension(fLocInitCoords.NumberOfNodes(), fLocInitCoords.MinorDim());
	fLocInitCoords.ReturnTranspose(temp);
	temp.WriteNumbered(out);

	out <<   " displacements:\n";
	temp.Dimension(fLocDisp.NumberOfNodes(), fLocDisp.MinorDim());
	fLocDisp.ReturnTranspose(temp);
	temp.WriteNumbered(out);
}

/* check material outputs - return true if OK */
bool ContinuumElementT::CheckMaterialOutput(void) const
{
	/* check compatibility of output */
	if (fMaterialList && fMaterialList->Length() > 1)
	{
		/* check compatibility of material outputs */
		bool OK = true;
		int i, j;
		for (i = 0; OK && i < fMaterialList->Length(); i++)
		{
			ContinuumMaterialT* m_i = (*fMaterialList)[i];
			for (j = i+1; OK && j < fMaterialList->Length(); j++)
			{
				ContinuumMaterialT* m_j = (*fMaterialList)[j];
				OK = ContinuumMaterialT::CompatibleOutput(*m_i, *m_j);
			}
		}
		i--; j--;
			
		/* output not compatible */
		if (!OK)	
		{
#pragma message("report names")
			cout << "\n ContinuumElementT::CheckMaterialOutput: incompatible output\n"
			    <<    "     between materials " << i+1 << " and " << j+1 << ":\n";
//			(*fMaterialList)[i]->PrintName(cout);
			cout << '\n';
//			(*fMaterialList)[j]->PrintName(cout);
			cout << endl;
			return false;
		}
	}
	
	/* no problems */
	return true;
}

/* describe the parameters needed by the interface */
void ContinuumElementT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ElementBaseT::DefineParameters(list);

	/* flag to store shape functions and derivatives */
	ParameterT store_shape(fStoreShape, "store_shapefunctions");
	store_shape.SetDefault(fStoreShape);
	list.AddParameter(store_shape, ParameterListT::ZeroOrOnce);
}

/* information about subordinate parameter lists */
void ContinuumElementT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ElementBaseT::DefineSubs(sub_list);

	/* geometry and integration rule (inline) */
	sub_list.AddSub("element_geometry", ParameterListT::Once, true);

	/* optional body force */
	sub_list.AddSub("body_force", ParameterListT::ZeroOrOnce);
	
	/* tractions */
	sub_list.AddSub("natural_bc", ParameterListT::Any);
}

/* return the description of the given inline subordinate parameter list */
void ContinuumElementT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	/* geometry and integration rule (inline) */
	if (name == "element_geometry")
	{
		/* choice */
		order = ParameterListT::Choice;
	
		/* element geometries */
		sub_lists.AddSub(GeometryT::ToString(GeometryT::kQuadrilateral));
		#ifdef AG_ELEMENT_DEV
		sub_lists.AddSub(GeometryT::ToString(GeometryT::kQuadrilateralAG));
		#endif
		sub_lists.AddSub(GeometryT::ToString(GeometryT::kTriangle));
		sub_lists.AddSub(GeometryT::ToString(GeometryT::kHexahedron));
		sub_lists.AddSub(GeometryT::ToString(GeometryT::kTetrahedron));
		sub_lists.AddSub(GeometryT::ToString(GeometryT::kLine));
	}
	else
		ElementBaseT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ContinuumElementT::NewSub(const StringT& name) const
{
	/* create non-const this */
	ContinuumElementT* non_const_this = const_cast<ContinuumElementT*>(this);

	/* try material list */
	MaterialListT* material_list = non_const_this->NewMaterialList(name, 0);
	if (material_list)
		return material_list;
		
	/* try geometry */
	ParameterInterfaceT* geometry = GeometryT::New(name);
	if (geometry)
		return geometry;

	/* body force */
	if (name == "body_force")
	{
		ParameterContainerT* body_force = new ParameterContainerT(name);
	
		/* schedule number */
		body_force->AddParameter(ParameterT::Integer, "schedule");
	
		/* body force vector */
		body_force->AddSub("Double", ParameterListT::OnePlus); 		
		
		return body_force;
	}
	else if (name == "natural_bc") /* traction bc */
	{
		ParameterContainerT* natural_bc = new ParameterContainerT(name);

		natural_bc->AddParameter(ParameterT::Word, "side_set_ID");
		natural_bc->AddParameter(ParameterT::Integer, "schedule");

		ParameterT coord_sys(ParameterT::Enumeration, "coordinate_system");
		coord_sys.AddEnumeration("global", Traction_CardT::kCartesian);
		coord_sys.AddEnumeration( "local", Traction_CardT::kLocal);
		coord_sys.SetDefault(Traction_CardT::kCartesian);
		natural_bc->AddParameter(coord_sys);

		natural_bc->AddSub("DoubleList", ParameterListT::OnePlus); 		
		
		return natural_bc;
	}
	else /* inherited */
		return ElementBaseT::NewSub(name);
}

/* accept parameter list */
void ContinuumElementT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "ContinuumElementT::TakeParameterList";

	/* resolve geometry before calling inherited method - geometry code
	 * may be needed while reading connectivities */
	const ParameterListT& integration_domain = list.GetListChoice(*this, "element_geometry");
	fGeometryCode = GeometryT::string2CodeT(integration_domain.Name());

	/* inherited */
	ElementBaseT::TakeParameterList(list);

	/* construct group communicator */
	const CommunicatorT& comm = ElementSupport().Communicator();
	int color = (NumElements() > 0) ? 1 : CommunicatorT::kNoColor;
	fGroupCommunicator = new CommunicatorT(comm, color, comm.Rank());

	/* allocate work space */
	fNEEvec.Dimension(NumElementNodes()*NumDOF());
	fDOFvec.Dimension(NumDOF());

	/* initialize local arrays */
	SetLocalArrays();

	/* construct shape functions */
	fNumIP = integration_domain.GetParameter("num_ip");
	SetShape();

	/* flag to compute and store shape function derivatives */
	const ParameterT* store_shape = list.Parameter("store_shapefunctions");
	if (store_shape) fStoreShape = *store_shape;

	/* construct material list */
	ParameterListT mat_params;
	CollectMaterialInfo(list, mat_params);
	if (mat_params.NumLists() > 0)
	{
		fMaterialList = NewMaterialList(mat_params.Name(), mat_params.NumLists());
		if (!fMaterialList) ExceptionT::GeneralFail(caller, "could not construct material list \"%s\"", 
			mat_params.Name().Pointer());
		fMaterialList->TakeParameterList(mat_params);
	}

	/* get form of tangent */
	GlobalT::SystemTypeT type = TangentType();
	
	/* set form of element stiffness matrix */
	if (type == GlobalT::kSymmetric)
		fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);
	else if (type == GlobalT::kNonSymmetric)
		fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
	else if (type == GlobalT::kDiagonal)
		fLHS.SetFormat(ElementMatrixT::kDiagonal);

	/* body force */
	const ParameterListT* body_force = list.List("body_force");
	if (body_force) {
		int schedule = body_force->GetParameter("schedule");
		fBodySchedule = ElementSupport().Schedule(--schedule);

		/* body force vector */
		const ArrayT<ParameterListT>& body_force_vector = body_force->Lists();
		if (body_force_vector.Length() != NumDOF())
			ExceptionT::BadInputValue(caller, "body force is length %d not %d",
				body_force_vector.Length(), NumDOF());
		fBody.Dimension(NumDOF());
		for (int i = 0; i < fBody.Length(); i++)
			fBody[i] = body_force_vector[i].GetParameter("value");
	}
	
	/* extract natural boundary conditions */
	TakeNaturalBC(list);
}

/* extract the list of material parameters */
void ContinuumElementT::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{
#pragma unused(all_params)
	mat_params.Clear();
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* update traction BC data */
void ContinuumElementT::SetTractionBC(void)
{
//NOTE: With the possibility of variable global node numbers and
//		and equations, we assume as little as possible here with
//      regard to the validity of the node/equation numbers, requiring
//      only that NodesX in the element cards has the correct global
//      node numbers.

	/* dimensions */
	int ndof = NumDOF();

	/* echo values */
	iArray2DT nd_tmp, eq_tmp;
	for (int i = 0; i < fTractionList.Length(); i++)
	{
		Traction_CardT& BC_card = fTractionList[i];
			
		/* traction element/facet */
		int elem, facet;
		BC_card.Destination(elem, facet);

		/* set global node numbers */
		const iArrayT& loc_nodes = BC_card.LocalNodeNumbers();
		int nnd = loc_nodes.Length();
		
		iArrayT& nodes = BC_card.Nodes();
		nodes.Dimension(nnd);
		nodes.Collect(loc_nodes, fElementCards[elem].NodesX());
		
		/* set global equation numbers */
		iArrayT& eqnos = BC_card.Eqnos();
		eqnos.Dimension(ndof*nnd);
		
		/* get from node manager */
		nd_tmp.Set(1, nnd, nodes.Pointer());
		eq_tmp.Set(1, ndof*nnd, eqnos.Pointer());
		Field().SetLocalEqnos(nd_tmp, eq_tmp);
	}

	/* set flag */
	fTractionBCSet = 1;
}

/* return the default number of element nodes */
int ContinuumElementT::DefaultNumElemNodes(void) const
{
	switch (fGeometryCode)
	{
		case GeometryT::kLine:
			return 2;
		case GeometryT::kQuadrilateral:
			return 4;
		#ifdef AG_ELEMENT_DEV	
		case GeometryT::kQuadrilateralAG:
			return 9;
		#endif	
		case GeometryT::kTriangle:
			return 3;
		case GeometryT::kHexahedron:
			return 8;
		case GeometryT::kTetrahedron:
			return 4;
		case GeometryT::kPentahedron:
			return 6;
		default:
			cout << "\n ContinuumElementT::DefaultNumElemNodes: unknown geometry code: "
			     << fGeometryCode << endl;
			return 0;
	}
}
//NOTE: needed because ExodusII does not store ANY information about
//      empty element groups, which causes trouble for parallel execution
//      when a partition contains no element from a group.
