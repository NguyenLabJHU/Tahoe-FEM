/* $Id: ContinuumElementT.cpp,v 1.18 2002-07-02 19:55:23 cjkimme Exp $ */
/* created: paklein (10/22/1996) */

#include "ContinuumElementT.h"

#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
//#include "FEManagerT.h"
#include "ModelManagerT.h"
#include "StructuralMaterialT.h"
#include "ShapeFunctionT.h"
#include "DomainIntegrationT.h"
#include "eControllerT.h"
#include "Traction_CardT.h"
#include "iAutoArrayT.h"
#include "OutputSetT.h"
#include "ScheduleT.h"

//TEMP: all this for general traction BC implementation?
#include "VariArrayT.h"
#include "nVariArray2DT.h"
#include "VariLocalArrayT.h"

/* services */
#include "EdgeFinderT.h"
#include "GraphT.h"

/* materials lists */
#include "MaterialListT.h"
#include "Material2DT.h"

/* constructor */

using namespace Tahoe;

ContinuumElementT::ContinuumElementT(const ElementSupportT& support, 
	const FieldT& field):
	ElementBaseT(support, field),
	fMaterialList(NULL),
	fBodySchedule(NULL),
	fBody(NumDOF()),
	fTractionBCSet(0),
	fShapes(NULL),
	fLocInitCoords(LocalArrayT::kInitCoords),
	fLocDisp(LocalArrayT::kDisp),
	fDOFvec(NumDOF())
{
	ifstreamT&  in = ElementSupport().Input();
	ostream&    out = ElementSupport().Output();
		
	/* control parameters */
	in >> fGeometryCode; //TEMP - should actually come from the geometry database
	in >> fNumIP;
}

/* destructor */
ContinuumElementT::~ContinuumElementT(void)
{	
	delete fShapes;
	delete fMaterialList;
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

/* allocates space and reads connectivity data */
void ContinuumElementT::Initialize(void)
{
	/* inherited */
	ElementBaseT::Initialize();
	
	/* allocate work space */
	fNEEvec.Allocate(NumElementNodes()*NumDOF());

	/* initialize local arrays */
	SetLocalArrays();

	/* construct shape functions */
	SetShape();

	/* streams */
	ifstreamT& in = ElementSupport().Input();
	ostream&  out = ElementSupport().Output();

	/* output print specifications */
	EchoOutputCodes(in, out);

	/* body force specification (non virtual) */
	EchoBodyForce(in, out);
	
	/* echo traction B.C.'s (non virtual) */
	EchoTractionBC(in, out);

	/* echo material properties */
	ReadMaterialData(in);	
	WriteMaterialData(out);

	/* get form of tangent */
	GlobalT::SystemTypeT type = TangentType();
	
	/* set form of element stiffness matrix */
	if (type == GlobalT::kSymmetric)
		fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);
	else if (type == GlobalT::kNonSymmetric)
		fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
	else if (type == GlobalT::kDiagonal)
		fLHS.SetFormat(ElementMatrixT::kDiagonal);
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

	for (int i = 0; i < fMaterialList->Length(); i++)
	{
		GlobalT::SystemTypeT e_type = (*fMaterialList)[i]->TangentType();
	
		/* using type precedence */
		type = (e_type > type) ? e_type : type;
	}
	
	return type;
}

/* initialize/finalize step */
void ContinuumElementT::InitStep(void)
{
	/* inherited */
	ElementBaseT::InitStep();

	/* set material variables */
	fMaterialList->InitStep();
}

/* initialize/finalize step */
void ContinuumElementT::CloseStep(void)
{
	/* inherited */
	ElementBaseT::CloseStep();

	/* set material variables */
	fMaterialList->CloseStep();

	/* update element level internal variables */
	if (fMaterialList->HasHistoryMaterials())
	{
		Top();
		while (NextElement())
		{
			ElementCardT& element = CurrentElement();
			if (element.IsAllocated())
			{
				ContinuumMaterialT* pmat = (*fMaterialList)[element.MaterialNumber()];

				/* material update function */
				pmat->UpdateHistory();
			}
		}
	}
}

/* resets to the last converged solution */
void ContinuumElementT::ResetStep(void)
{
	/* inherited */
	ElementBaseT::ResetStep();

	/* update material internal variables */
	if (fMaterialList->HasHistoryMaterials())
	{
		Top();
		while (NextElement())
		{
			ElementCardT& element = CurrentElement();		
			if (element.IsAllocated())
			{
				ContinuumMaterialT* pmat = (*fMaterialList)[element.MaterialNumber()];

				/* material reset function */
				pmat->ResetHistory();
			}
		}
	}
}

/* restart operations */
void ContinuumElementT::ReadRestart(istream& in)
{
	/* inherited */
	ElementBaseT::ReadRestart(in);

	/* update element level internal variables */
	if (fMaterialList->HasHistoryMaterials())
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
	if (fMaterialList->HasHistoryMaterials())
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
	ArrayT<StringT> block_ID(fBlockData.Length());
	for (int i = 0; i < block_ID.Length(); i++)
		block_ID[i] = fBlockData[i].ID();

	/* collect variable labels */
	ArrayT<StringT> n_labels(n_counts.Sum());
	ArrayT<StringT> e_labels(e_counts.Sum());
	GenerateOutputLabels(n_counts, n_labels, e_counts, e_labels);

	/* set output specifier */
	StringT set_ID;
	set_ID.Append(ElementSupport().ElementGroupNumber(this) + 1);
	OutputSetT output_set(set_ID, fGeometryCode, block_ID, fConnectivities,
		n_labels, e_labels, false);
		
	/* register and get output ID */
	fOutputID = ElementSupport().RegisterOutput(output_set);
}

//NOTE - this function is/was identical to CSEBaseT::WriteOutput
void ContinuumElementT::WriteOutput(IOBaseT::OutputModeT mode)
{
//TEMP - not handling general output modes yet
	if (mode != IOBaseT::kAtInc)
	{
		cout << "\n ContinuumElementT::WriteOutput: only handling \"at increment\"\n"
		     <<   "     print mode. SKIPPING." << endl;
		return;
	}

	/* map output flags to count of values */
	iArrayT n_counts;
	SetNodalOutputCodes(mode, fNodalOutputCodes, n_counts);
	iArrayT e_counts;
	SetElementOutputCodes(mode, fElementOutputCodes, e_counts);

	/* calculate output values */
	dArray2DT n_values;
	dArray2DT e_values;
	ComputeOutput(n_counts, n_values, e_counts, e_values);

	/* send to output */
	ElementSupport().WriteOutput(fOutputID, n_values, e_values);
}

/* side set to nodes on facets data */
void ContinuumElementT::SideSetToFacets(const StringT& block_ID, const iArray2DT& sideset,
	iArray2DT& facets) const
{
	/* checks */
	if (sideset.MinorDim() != 2) throw eGeneralFail;
	
	/* empty set */
	if (sideset.MajorDim() == 0)
	{
		facets.Allocate(0, 0);
		return;
	}

// NOTE: faster to get all nodes_on_facet data at once. also
//       would be easier to check dimensions of facets.

	/* get block data */
	const ElementBlockDataT& block_data = BlockData(block_ID);
	
	int offset = block_data.StartNumber();
	iArrayT facet_nodes, facet_tmp;
	int num_facets = sideset.MajorDim();
	for (int i = 0; i < num_facets; i++)
	{
		int nel = sideset(i,0) + offset;
		int nft = sideset(i,1);
	
		/* get facet node map */
		ShapeFunction().NodesOnFacet(nft, facet_nodes);
		
		/* dimension check */
		if (i == 0)
			facets.Allocate(sideset.MajorDim(), facet_nodes.Length());
		else if (facets.MinorDim() != facet_nodes.Length())
		{
			cout << "\n ContinuumElementT::SideSetToFacets: all sides in set must have\n"
			     <<   "     the same number of nodes in element block" << block_ID << endl;
			throw eGeneralFail;
		}
		
		/* shape check */

		/* get node numbers */
		facets.RowAlias(i, facet_tmp);
		facet_tmp.Collect(facet_nodes, fElementCards[nel].NodesX());
	}
}

/* return geometry and number of nodes on each facet */
void ContinuumElementT::FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geometry, 
	iArrayT& num_facet_nodes) const
{
	/* from integration domain */
	ShapeFunction().FacetGeometry(facet_geometry, num_facet_nodes);
}

/* initial condition/restart functions (per time sequence) */
void ContinuumElementT::InitialCondition(void)
{
	/* inherited */
	ElementBaseT::InitialCondition();
	
	/* check for initialization materials */
	bool need_init = false;
	for (int i = 0; i < fMaterialList->Length() && !need_init; i++)
		need_init = (*fMaterialList)[i]->NeedsPointInitialization();

	/* initialize materials */
	if (need_init)
	{
		/* loop over elements */
		Top();
		while (NextElement())
		{
			/* material pointer */
			ContinuumMaterialT* pmat = (*fMaterialList)[CurrentElement().MaterialNumber()];
		
			if (pmat->NeedsPointInitialization())
			{
				/* global shape function values */
				SetGlobalShape();
			
				/* loop over integration points */
				fShapes->TopIP();
				while (fShapes->NextIP())
					pmat->PointInitialize();
			}
		}
	}
}

/* surface facets */
void ContinuumElementT::SurfaceFacets(GeometryT::CodeT& geometry,
	iArray2DT& surface_facets, iArrayT& surface_nodes) const
{
	/* surface facets must all have same geometry */
	ArrayT<GeometryT::CodeT> facet_geom;
	iArrayT facet_nodes;
	ShapeFunction().FacetGeometry(facet_geom, facet_nodes);
	if (facet_nodes.Count(facet_nodes[0]) != facet_geom.Length())
	{
		cout << "\n ContinuumElementT::SurfaceFacets: only support identical\n";
		cout <<   "     facet shapes" << endl;
		throw eGeneralFail;
	}
	geometry = facet_geom[0];

	/* find bounding elements */
	AutoArrayT<int> border_nodes;
	iArrayT   border_elems;
	iArray2DT border_neighs;
	BoundingElements(border_elems, border_neighs);
	
	/* check */
	if (ShapeFunction().NumFacets() != border_neighs.MinorDim())
		throw eSizeMismatch;
		
	/* collect nodes on facets info */
	ArrayT<iArrayT> facetnodemap(ShapeFunction().NumFacets());
	for (int i2 = 0; i2 < facetnodemap.Length(); i2++)
		ShapeFunction().NodesOnFacet(i2, facetnodemap[i2]);	

	/* collect surface facets (with "outward" normal ordering) */
	int surf_count = 0;
	int num_facets = facetnodemap.Length();
	int num_facet_nodes = facet_nodes[0];
	border_nodes.Allocate(0);
	surface_facets.Allocate(border_neighs.Count(-1), num_facet_nodes);
	for (int i = 0; i < border_elems.Length(); i++)
	{
		/* element connectivity */
	        const iArrayT& elemnodes = fElementCards[border_elems[i]].NodesX();
		int* elem = elemnodes.Pointer();

		/* find open sides */
		int found_open = 0;
		int* pneigh = border_neighs(i);
		for (int j = 0; j < num_facets; j++)
		{
			/* open face */
			if (*pneigh == -1)
			{
				/* set flag */
				found_open = 1;
				
				/* collect facet nodes */
				int* pfacet = surface_facets(surf_count++);
				int* facet_nodes = facetnodemap[j].Pointer();
				for (int k = 0; k < num_facet_nodes; k++)
				{
					int node = elem[*facet_nodes++];
					*pfacet++ = node;
					border_nodes.AppendUnique(node);
					// better just to keep a "nodes used" map?
				}
			}	
			pneigh++;
		}
	
		/* no open facet */	
		if (!found_open)
		{
			cout << "\n ContinuumElementT::SurfaceFacets: error building surface facet list" << endl;
			throw eGeneralFail;
		}	
	}

	/* return value */
	surface_nodes.Allocate(border_nodes.Length());
	border_nodes.CopyInto(surface_nodes);
}

/* with surface facets sorted into connected sets */
void ContinuumElementT::SurfaceFacets(GeometryT::CodeT& geometry,
	ArrayT<iArray2DT>& surface_facet_sets,
	iArrayT& surface_nodes) const
{
	/* collect all surface facets */
	iArray2DT surface_facets;
	SurfaceFacets(geometry, surface_facets, surface_nodes);

	/* graph object */
	GraphT graph;
	graph.AddGroup(surface_facets);
	graph.MakeGraph();

	iArrayT branch_map;
	graph.LabelBranches(surface_nodes, branch_map);
	
	/* sort surfaces */
	int num_branches = branch_map.Max() + 1;
	surface_facet_sets.Allocate(num_branches);
	if (num_branches == 1)
		surface_facet_sets[0] = surface_facets;
	else
	{
		/* surfaces in each set */
		iArrayT count(num_branches);
		int size = surface_facets.MinorDim();
		count = 0;
		int* psurf = surface_facets(0);
		for (int i = 0; i < surface_facets.MajorDim(); i++)
		{
			count[branch_map[*psurf]]++;
			psurf += size;
		}

		/* set to surfaces map */
		RaggedArray2DT<int> set_data;
		set_data.Configure(count);

		count = 0;
		psurf = surface_facets(0);
		for (int j = 0; j < surface_facets.MajorDim(); j++)
		{
			int branch = branch_map[*psurf];
			*(set_data(branch) + count[branch]) = j;

			count[branch]++;
			psurf += size;
		}
			
		/* copy in */
		for (int k = 0; k < num_branches; k++)
		{
			surface_facet_sets[k].Allocate(set_data.MinorDim(k), size);
			surface_facet_sets[k].RowCollect(set_data(k), surface_facets);
		}
	}
}		
	
/* surface nodes */
void ContinuumElementT::SurfaceNodes(iArrayT& surface_nodes) const
{
	/* work space */
	AutoArrayT<int> border_nodes;
	iArrayT   border_elems;
	iArray2DT border_neighs;

	/* find bounding elements */
	BoundingElements(border_elems, border_neighs);
	
	/* check */
	if (ShapeFunction().NumFacets() != border_neighs.MinorDim())
		throw eSizeMismatch;
		
	/* collect nodes on facets map */
	ArrayT<iArrayT> facetnodemap(ShapeFunction().NumFacets());
	for (int i2 = 0; i2 < facetnodemap.Length(); i2++)
		ShapeFunction().NodesOnFacet(i2, facetnodemap[i2]);
			
	/* collect surface nodes from border elems */
	border_nodes.Allocate(0);
	int numfacets = facetnodemap.Length();
	for (int i = 0; i < border_elems.Length(); i++)
	{
		/* element connectivity */
	        const iArrayT& elemnodes = fElementCards[border_elems[i]].NodesX();
		int* elem = elemnodes.Pointer();

		/* find open sides */
		int found_open = 0;
		int* pneigh = border_neighs(i);
		for (int j = 0; j < numfacets; j++)
		{
			/* open face */
			if (*pneigh == -1)
			{
				/* set flag */
				found_open = 1;
				
				/* collect facet nodes */
				int  num_facet_nodes = facetnodemap[j].Length();
				int*     facet_nodes = facetnodemap[j].Pointer();
				for (int k = 0; k < num_facet_nodes; k++)
					border_nodes.AppendUnique(elem[*facet_nodes++]);
			}	
			pneigh++;
		}
	
		/* no open facet */	
		if (!found_open)
		{
			cout << "\n ContinuumElementT::SurfaceNodes: error building surface node list" << endl;
			throw eGeneralFail;
		}	
	}
	
	/* return value */
	surface_nodes.Allocate(border_nodes.Length());
	border_nodes.CopyInto(surface_nodes);
}

/***********************************************************************
* Protected
***********************************************************************/

namespace Tahoe {

/* stream extraction operator */
istream& operator>>(istream& in, ContinuumElementT::MassTypeT& type)
{
	int i_type;
	in >> i_type;
	switch (i_type)
	{
		case ContinuumElementT::kNoMass:
			type = ContinuumElementT::kNoMass;
			break;
		case ContinuumElementT::kConsistentMass:
			type = ContinuumElementT::kConsistentMass;
			break;
		case ContinuumElementT::kLumpedMass:
			type = ContinuumElementT::kLumpedMass;
			break;
		default:
			cout << "\n ContinuumElementT::MassTypeT: unknown type: "
			<< i_type<< endl;
			throw eBadInputValue;	
	}
	return in;
}

}

/* initialize local arrays */
void ContinuumElementT::SetLocalArrays(void)
{
	/* dimension */
	fLocInitCoords.Allocate(NumElementNodes(), NumSD());
	fLocDisp.Allocate(NumElementNodes(), NumDOF());

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
		
		/* Jacobian of the surface mapping */
		dMatrixT jacobian(nsd, nsd-1);
		
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
			
#ifdef __NO_RTTI__
			/* default thickness */
			double thick = 1.0;
#else
			/* use thickness for 2D solid deformation elements */
			double thick = 1.0;
			if (ndof == 2 && ndof == 2) //better to do this once elsewhere?
			{
				/* get material pointer */
				const ElementCardT& elem_card = fElementCards[elem];
				ContinuumMaterialT* pmat = (*fMaterialList)[elem_card.MaterialNumber()];
			
				/* thickness from 2D material */
				Material2DT* pmat2D = dynamic_cast<Material2DT*>(pmat);
				if (pmat2D) thick = pmat2D->Thickness();
			}
#endif
			
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
					double jwt = detj*w[j]*thick;
					
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
					double jwt = detj*w[j]*thick;
					
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
				throw eGeneralFail;

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
void ContinuumElementT::FormMass(int mass_type, double constM)
{
#if __option(extended_errorcheck)
	if (fLocDisp.Length() != fLHS.Rows()) throw eSizeMismatch;
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
			
			int nen = fLocDisp.NumberOfNodes();
			int ndof = NumDOF();
			
			/* matrix form */
			int a = 0, zero = 0;
			int& b_start = (fLHS.Format() == ElementMatrixT::kSymmetricUpper) ? a : zero;
			
			fShapes->TopIP();	
			while ( fShapes->NextIP() )
			{
				double temp = constM*(*Weight++)*(*Det++);
				const double* Na = fShapes->IPShapeU();
								
				for (a = 0; a < nen; a++)
					for (int i = 0; i < ndof; i++)
					{
						int p = a*ndof + i;
						
						/* upper triangle only */
						for (int b = b_start; b < nen; b++)
							for (int j = 0; j < ndof; j++)
								if(i == j)
								{									
									int q = b*ndof + j;
									fLHS(p,q) += temp*Na[a]*Na[b];
								}
					}
			}
			break;
		}

		case kLumpedMass:	/* lumped mass */
		{
			int nen = fLocDisp.NumberOfNodes();
			int ndof = NumDOF();

		    double dsum   = 0.0;
		    double totmas = 0.0;
		    fNEEvec = 0.0;

			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();

			/* total mass and diagonal sum */
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				double temp1     = constM*(*Weight++)*(*Det++);
				const double* Na = fShapes->IPShapeU();

				totmas += temp1;
				for (int lnd = 0; lnd < nen; lnd++)
				{
					double temp2 = temp1*Na[lnd]*Na[lnd];
					dsum += temp2;
					fNEEvec[lnd] += temp2;
				}
			}	
				
			/* scale diagonal to conserve total mass */
			double diagmass = totmas/dsum;
			
			/* lump mass onto diagonal */
			double* pmass = fLHS.Pointer();
			int inc = fLHS.Rows() + 1;
			for (int lnd = 0; lnd < nen; lnd++)
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
		
			cout << "\n Elastic::FormMass: unknown mass matrix code\n" << endl;
			throw eBadInputValue;
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
void ContinuumElementT::FormMa(MassTypeT mass_type, double constM, 
	const LocalArrayT* nodal_values,
	const dArray2DT* ip_values)
{
	/* quick exit */
	if (!nodal_values && !ip_values) return;

#if __option(extended_errorcheck)
	/* dimension checks */
	if (nodal_values && 
		fRHS.Length() != nodal_values->Length()) 
		throw eSizeMismatch;

	if (ip_values &&
		(ip_values->MajorDim() != fShapes->NumIP() ||
		 ip_values->MinorDim() != NumDOF()))
		throw eSizeMismatch;
#endif

	switch (mass_type)
	{
		case kConsistentMass:
		{
			int ndof = NumDOF();
			int  nen = NumElementNodes();

			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();

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
				
				double temp = constM*(*Weight++)*(*Det++);				
				for (int lnd = 0; lnd < nen; lnd++)
				{
					double temp2 = temp*(*Na++);
					double* pacc = fDOFvec.Pointer();

					for (int dof = 0; dof < ndof; dof++)			
						*res++ += temp2*(*pacc++);
				}
			}
			break;
		}	
		case kLumpedMass:
		{
			fLHS = 0.0; //hope there's nothing in there!
			FormMass(kLumpedMass, constM);

			/* init nodal values */
			if (nodal_values)
				nodal_values->ReturnTranspose(fNEEvec);
			else {
				cout << "\n ContinuumElementT::FormMa: expecting nodal values for lumped mass" << endl;
				throw eGeneralFail;
			}
				
//TEMP - what to do with ip values?
if (ip_values) {
	cout << "\n ContinuumElementT::FormMa: lumped mass not implemented for ip sources" << endl;
	throw eGeneralFail;
}

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

/* print element group data */
void ContinuumElementT::PrintControlData(ostream& out) const
{
	/* inherited */
	ElementBaseT::PrintControlData(out);

	out << " Element geometry code . . . . . . . . . . . . . = " << fGeometryCode << '\n';
	out << "    eq." << GeometryT::kPoint         << ", point\n";
	out << "    eq." << GeometryT::kLine          << ", line\n";
	out << "    eq." << GeometryT::kQuadrilateral << ", quadrilateral\n";
	out << "    eq." << GeometryT::kTriangle	  << ", triangle\n";
	out << "    eq." << GeometryT::kHexahedron	  << ", hexahedron\n";
	out << "    eq." << GeometryT::kTetrahedron   << ", tetrahedron\n";
	out << " Number of integration points. . . . . . . . . . = " << fNumIP    << '\n';
}

void ContinuumElementT::ReadMaterialData(ifstreamT& in)
{
	/* construct material list */
	int size;
	in >> size;
	fMaterialList = NewMaterialList(size);
	if (!fMaterialList) throw eOutOfMemory;

	/* read */
	fMaterialList->ReadMaterialData(in);
	
	/* check range */
	for (int i = 0; i < fBlockData.Length(); i++)
		if (fBlockData[i].MaterialID() < 0 ||
		    fBlockData[i].MaterialID() >= size)
		{
			cout << "\n ContinuumElementT::ReadMaterialData: material number "
			     << fBlockData[i].MaterialID() + 1 << '\n';
			cout<<    "     for element block " << i + 1 << " is out of range" << endl;
			throw eBadInputValue;
		}
}

/* use in conjunction with ReadMaterialData */
void ContinuumElementT::WriteMaterialData(ostream& out) const
{
	fMaterialList->WriteMaterialData(out);

	/* flush buffer */
	out.flush();
}

void ContinuumElementT::EchoBodyForce(ifstreamT& in, ostream& out)
{
	/* schedule number and body force vector */
	int n_sched;
	in >> n_sched >> fBody;		
	n_sched--;

	/* no LTf => no body force */
	if (n_sched < 0) 
		fBody = 0.0;
	else
	{
		fBodySchedule = ElementSupport().Schedule(n_sched);
		if (!fBodySchedule) {
			cout << "\n ContinuumElementT::EchoBodyForce: could not resolve schedule " 
			     << n_sched + 1 << endl;
			throw eBadInputValue;
		}	
	}
	
	out << "\n Body force vector:\n";
	out << " Body force load-time function number. . . . . . = " << n_sched + 1<< '\n';
	out << " Body force vector components:\n";
	for (int j = 0 ; j < NumDOF(); j++)
	{
		out << "   x[" << j+1 << "] direction. . . . . . . . . . . . . . . . = ";
		out << fBody[j] << '\n';
	}
	out.flush();   	   	
}

void ContinuumElementT::EchoTractionBC(ifstreamT& in, ostream& out)
{
	out << "\n Traction boundary conditions:\n";
	
	/* read data from parameter file */
	int numlines, numsets;
	ModelManagerT& model = ElementSupport().Model();
	model.ReadNumTractionLines (in, numlines, numsets);

	if (numlines > 0)
	  {
	    /* temp space */
	    ArrayT<StringT> block_ID(numlines);
	    ArrayT<iArray2DT> localsides (numlines);
	    iArrayT LTf (numlines);
	    ArrayT<Traction_CardT::CoordSystemT> coord_sys (numlines);
	    ArrayT<dArray2DT> values (numlines);

	    /* nodes on element facets */
	    iArrayT num_facet_nodes;
	    fShapes->NumNodesOnFacets (num_facet_nodes);
	    
	    /* read data by blocks */
	    int line = 0;
	    int count = 0;
	    for (int blockset = 0; blockset < numsets; blockset++)
	      {
		/* read num of cards in each block */
		int setsize = -1;
		StringT set_ID;
		model.ReadTractionSetData (in, set_ID, setsize);
		for (int card=0; card < setsize; card++)
		  {
		    /* read side set for that card */
		    block_ID[line] = set_ID;
		    model.ReadTractionSideSet (in, block_ID[line], localsides[line]);

		    /* increment count */
		    int num_sides = localsides[line].MajorDim();
		    count += num_sides;

		    /* read data for that card */
		    in >> LTf[line] >> coord_sys[line];

		    /* skip if empty */
		    int num_nodes;
		    if (num_sides > 0)
		      {
			iArray2DT& side_set = localsides[line];

			/* switch to group numbering */
			const ElementBlockDataT& block_data = BlockData (block_ID[line]);
			iArrayT elems (num_sides);
			side_set.ColumnCopy (0, elems);

			/* check */
			int min, max;
			elems.MinMax (min, max);
			if (min < 0 || max > block_data.Dimension())
			  {
			    cout << "\n ContinuumElementT::EchoTractionBC_TahoeII: node numbers\n";
			    cout <<   "     {"<< min << "," << max << "} are out of range in ";
			    cout << " dataline " << line << endl;
			    throw eBadInputValue;
			  }
			/* shift */
			elems += block_data.StartNumber();
			side_set.SetColumn (0, elems);

			/* all facets in set must have the same number of nodes */
			num_nodes = num_facet_nodes [side_set (0,1)];
			for (int f=0; f < num_sides; f++)
			  if (num_facet_nodes[side_set(f,1)] != num_nodes)
			    {
			      cout << "\n ContinuumElementT::EchoTractionBC_TahoeII: sides specified\n";
			      cout <<   "     in line " << line << " have differing numbers of nodes";
			      cout << endl;
			      throw eBadInputValue;
			    }
		      }
		    else
		      {
			/* still check numbef of facet nodes */
			int min, max;
			num_facet_nodes.MinMax (min, max);
			if (min != max)
			  {
			    cout << "\n ContinuumElementT::EchoTractionBC_TahoeII: cannot determine number of\n"
				 <<   "     facet nodes for empty side set at line " << line << endl;
			    throw eBadInputValue;
			  }
			else
			  num_nodes = min;
		      }

		    /* read traction values */
		    dArray2DT& valueT = values[line];
		    valueT.Allocate (num_nodes, NumDOF());
		    in >> valueT;
		    //NOTE - cannot simply clear to the end of the line with empty side sets
		    //       because the tractions may be on multiple lines

		    line++;
		  }
	      }

	    /* allocate all traction BC cards */
	    fTractionList.Allocate (count);

	    /* correct numbering offset */
	    LTf--;

	    if (count > 0)
	      {
		out << '\n';
		out << setw (kIntWidth) << "no.";
		fTractionList[0].WriteHeader (out, NumDOF());

		iArrayT loc_node_nums;
		int dex = 0;
		for (int ii=0; ii < numlines; ii++)
		  {
		    /* set traction BC cards */
		    iArray2DT& side_set = localsides[ii];
		    int numsides = side_set.MajorDim();
		    for (int j=0; j < numsides; j++)
		      {
			out << setw (kIntWidth) << j+1;

			/* get facet local node numbers */
			fShapes->NodesOnFacet (side_set (j, 1), loc_node_nums);

			/* set and echo */
			fTractionList[dex++].EchoValues (ElementSupport(), side_set(j,0), side_set (j,1), LTf[ii],
							 coord_sys[ii], loc_node_nums, values[ii], out);
		      }
		    out << endl;
		  }
	      }
	  }


	if (NumSD() != NumDOF())
	{
		/* check coordinate system specifications */
		for (int i = 0; i < fTractionList.Length(); i++)
			if (fTractionList[i].CoordSystem() != Traction_CardT::kCartesian)
			{
				cout << "\n ContinuumElementT::EchoTractionBC: coordinate system must be\n"
				     <<   "    Cartesian:" << Traction_CardT::kCartesian
				     << " if (spatial dimensions != degrees of freedom)\n"
				     <<   "    for card " << i+1 << endl;
				throw eBadInputValue;
			}
	}
}

/* return the "bounding" elements and the corresponding
* neighbors, both dimensioned internally */
void ContinuumElementT::BoundingElements(iArrayT& elements, iArray2DT& neighbors) const
{
	//TEMP - not parallelized
	if (ElementSupport().Size() > 1)
		cout << "\n ContinuumElementT::BoundingElements: not extended to parallel" << endl;

	/* build element neighbor list */
	iArray2DT nodefacetmap;
	fShapes->NeighborNodeMap(nodefacetmap);
	EdgeFinderT edger(fConnectivities, nodefacetmap);
	const iArray2DT& all_neighbors = edger.Neighbors();

	/* collect list of bounding elements */
	AutoArrayT<int> borders;
	iArrayT element;
	int nel = NumElements();
	for (int i = 0; i < nel; i++)
	{
		all_neighbors.RowAlias(i, element);
	
		/* has "free" edge */
		if (element.HasValue(-1)) borders.Append(i);
	}
	elements.Allocate(borders.Length());
	borders.CopyInto(elements);
	
	/* copy bounding element neighbor lists */
	neighbors.Allocate(elements.Length(), all_neighbors.MinorDim());
	neighbors.RowCollect(elements, all_neighbors);
}

/* write all current element information to the stream */
void ContinuumElementT::CurrElementInfo(ostream& out) const
{
	/* inherited */
	ElementBaseT::CurrElementInfo(out);
	dArray2DT temp;
	temp.Allocate(fLocInitCoords.NumberOfNodes(), fLocInitCoords.MinorDim());
	
	out <<   " initial coords:\n";
	temp.Allocate(fLocInitCoords.NumberOfNodes(), fLocInitCoords.MinorDim());
	fLocInitCoords.ReturnTranspose(temp);
	temp.WriteNumbered(out);

	out <<   " displacements:\n";
	temp.Allocate(fLocDisp.NumberOfNodes(), fLocDisp.MinorDim());
	fLocDisp.ReturnTranspose(temp);
	temp.WriteNumbered(out);
}

/* check material outputs - return true if OK */
bool ContinuumElementT::CheckMaterialOutput(void) const
{
	/* check compatibility of output */
	if (fMaterialList->Length() > 1)
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
			cout << "\n ContinuumElementT::CheckMaterialOutput: incompatible output\n"
			    <<    "     between materials " << i+1 << " and " << j+1 << ":\n";
			(*fMaterialList)[i]->PrintName(cout);
			cout << '\n';
			(*fMaterialList)[j]->PrintName(cout);
			cout << endl;
			return false;
		}
	}
	
	/* no problems */
	return true;
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
		nodes.Allocate(nnd);
		nodes.Collect(loc_nodes, fElementCards[elem].NodesX());
		
		/* set global equation numbers */
		iArrayT& eqnos = BC_card.Eqnos();
		eqnos.Allocate(ndof*nnd);
		
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
