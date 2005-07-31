/* $Id: ContinuumElementT.cpp,v 1.1.1.1 2001-01-29 08:20:39 paklein Exp $ */
/* created: paklein (10/22/1996)                                          */

#include "ContinuumElementT.h"

#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "StructuralMaterialT.h"
#include "ShapeFunctionT.h"
#include "eControllerT.h"
#include "Traction_CardT.h"
#include "ExodusT.h"
#include "ModelFileT.h"
#include "iAutoArrayT.h"
#include "OutputSetT.h"

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
ContinuumElementT::ContinuumElementT(FEManagerT& fe_manager):
	ElementBaseT(fe_manager),
	fMaterialList(NULL),
	fBody(fNumDOF),
	fTractionBCSet(0),
	fShapes(NULL),
	fLocInitCoords(LocalArrayT::kInitCoords),
	fLocDisp(LocalArrayT::kDisp),
	fDOFvec(fNumDOF),
	fNSDvec(fNumSD)
{
	ifstreamT&  in = fFEManager.Input();
	ostream&    out = fFEManager.Output();
		
	/* control parameters */
	in >> fGeometryCode;
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
#if __option(extended_errorcheck)
	if (!fShapes) throw eGeneralFail;
#endif
	return fShapes->CurrIP();
}

/* allocates space and reads connectivity data */
void ContinuumElementT::Initialize(void)
{
	/* inherited */
	ElementBaseT::Initialize();
	
	/* allocate work space */
	fNEEvec.Allocate(fNumElemNodes*fNumDOF);

	/* initialize local arrays */
	SetLocalArrays();

	/* construct shape functions */
	SetShape();

	/* streams */
	ifstreamT& in = fFEManager.Input();
	ostream&   out = fFEManager.Output();

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

/* set element group for new global equations numbers */
void ContinuumElementT::Reinitialize(void)
{
	/* inherited */
	ElementBaseT::Reinitialize();

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
		Top();
		while (NextElement())
		{
			int isallocated;
			in >> isallocated;
			if (isallocated) CurrentElement().ReadRestart(in);
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
			out << setw(kIntWidth) << element.IsAllocated() << '\n';
			if (element.IsAllocated()) element.WriteRestart(out);
		}
	}
}

/* writing output */
void ContinuumElementT::RegisterOutput(void)
{
	/* ID is just group number */
	int ID = fFEManager.ElementGroupNumber(this) + 1;

//NOTE: could loop over each output mode and register
//      it with the output separately. for now just register
//      "kAtInc"
	
	/* get output configuration */
	iArrayT counts;
	SetOutputCodes(IOBaseT::kAtInc, fOutputCodes, counts);
	int num_out = counts.Sum();

	/* variable labels */
	ArrayT<StringT> n_labels(num_out);
	ArrayT<StringT> e_labels;
	GenerateOutputLabels(counts, n_labels);

	OutputSetT output_set(ID, fGeometryCode, fConnectivities,
		n_labels, e_labels, false);
		
	/* register and get output ID */
	fOutputID = fFEManager.RegisterOutput(output_set);
}

//NOTE - this function is identical to CSEBaseT::WriteOutput
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
	iArrayT counts;
	SetOutputCodes(mode, fOutputCodes, counts);
	int num_out = counts.Sum();

	dArray2DT group_n_values;
	dArray2DT group_e_values(0,0);
	if (num_out > 0)
	{
		/* reset averaging workspace */
		fNodes->ResetAverage(num_out);

		/* compute nodal values */
		ComputeNodalValues(counts);

		/* get nodal values */
		const iArrayT& node_used = fFEManager.OutputSet(fOutputID).NodesUsed();
		fNodes->OutputAverage(node_used, group_n_values);
	}

	/* send out */
	fFEManager.WriteOutput(fOutputID, group_n_values, group_e_values);
}

/* side set to nodes on facets data */
void ContinuumElementT::SideSetToFacets(int block_ID, const iArray2DT& sideset,
	iArray2DT& facets)
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
	const int* block_data = BlockData(block_ID);
	
	int offset = block_data[kStartNum];
	iArrayT facet_nodes, facet_tmp;
	int num_facets = sideset.MajorDim();
	for (int i = 0; i < num_facets; i++)
	{
		int nel = sideset(i,0) + offset;
		int nft = sideset(i,1);
	
		/* get facet node map */
		fShapes->NodesOnFacet(nft, facet_nodes);
		
		/* dimension/check */
		if (i == 0)
			facets.Allocate(sideset.MajorDim(), facet_nodes.Length());
		else if (facets.MinorDim() != facet_nodes.Length())
		{
			cout << "\n ContinuumElementT::SideSetToFacets: all sides in set must have\n"
			     <<   "     the same number of nodes in element block" << block_ID << endl;
			throw eGeneralFail;
		}

		/* get node numbers */
		facets.RowAlias(i, facet_tmp);
		facet_tmp.Collect(facet_nodes, fElementCards[nel].NodesX());
	}
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
				SetLocalU(fLocDisp);
			
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
	fShapes->FacetGeometry(facet_geom, facet_nodes);
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
	if (fShapes->NumFacets() != border_neighs.MinorDim())
		throw eSizeMismatch;
		
	/* collect nodes on facets info */
	ArrayT<iArrayT> facetnodemap(fShapes->NumFacets());
	for (int i2 = 0; i2 < facetnodemap.Length(); i2++)
		fShapes->NodesOnFacet(i2, facetnodemap[i2]);	

	/* collect surface facets (with "outward" normal ordering) */
	int surf_count = 0;
	int num_facets = facetnodemap.Length();
	int num_facet_nodes = facet_nodes[0];
	border_nodes.Allocate(0);
	surface_facets.Allocate(border_neighs.Count(-1), num_facet_nodes);
	for (int i = 0; i < border_elems.Length(); i++)
	{
		/* element connectivity */
		int* elem = fConnectivities(border_elems[i]);

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
	if (fShapes->NumFacets() != border_neighs.MinorDim())
		throw eSizeMismatch;
		
	/* collect nodes on facets map */
	ArrayT<iArrayT> facetnodemap(fShapes->NumFacets());
	for (int i2 = 0; i2 < facetnodemap.Length(); i2++)
		fShapes->NodesOnFacet(i2, facetnodemap[i2]);
			
	/* collect surface nodes from border elems */
	border_nodes.Allocate(0);
	int numfacets = facetnodemap.Length();
	for (int i = 0; i < border_elems.Length(); i++)
	{
		/* element connectivity */
		int* elem = fConnectivities(border_elems[i]);

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

/* initialize local arrays */
void ContinuumElementT::SetLocalArrays(void)
{
	/* dimension */
	fLocInitCoords.Allocate(fNumElemNodes, fNumSD);
	fLocDisp.Allocate(fNumElemNodes, fNumDOF);

	/* set source */
	fFEManager.RegisterLocal(fLocInitCoords);
	fFEManager.RegisterLocal(fLocDisp);	
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
		/* update equation numbers */
		if (!fTractionBCSet) SetTractionBC();
	
		/* force vector */
		dArrayT rhs;
		VariArrayT<double> rhs_man(25, rhs);
		
		/* local coordinates */
		LocalArrayT coords(LocalArrayT::kInitCoords);
		VariLocalArrayT coord_man(25, coords, fNumSD);
		fFEManager.RegisterLocal(coords);
		
		/* nodal tractions */
		LocalArrayT tract(LocalArrayT::kUnspecified);
		VariLocalArrayT tract_man(25, tract, fNumDOF);

		/* integration point tractions */
		dArray2DT ip_tract;
		nVariArray2DT<double> ip_tract_man(25, ip_tract, fNumDOF);
		dArrayT tract_loc, tract_glb(fNumDOF);
		dMatrixT Q(fNumDOF);
		
		/* Jacobian of the surface mapping */
		dMatrixT jacobian(fNumSD,fNumSD-1);
		
		for (int i = 0; i < fTractionList.Length(); i++)
		{
			const Traction_CardT& BC_card = fTractionList[i];

			/* dimension */
			const iArrayT& nodes = BC_card.Nodes();
			int nnd = nodes.Length();
			rhs_man.SetLength(nnd*fNumDOF, false);
			coord_man.SetNumberOfNodes(nnd);
			tract_man.SetNumberOfNodes(nnd);
			
			/* local coordinates */
			coords.SetLocal(nodes);

			/* nodal traction vectors: (ndof x nnd) */
			BC_card.CurrentValue(tract);
			
			/* BC destination */
			int elem, facet;
			BC_card.Destination(elem, facet);
			double thick = 1.0;
			if (fNumSD == 2) //better to do this once elsewhere?
			{
				/* get material pointer */
				const ElementCardT& elem_card = fElementCards[elem];
				ContinuumMaterialT* pmat = (*fMaterialList)[elem_card.MaterialNumber()];
			
#ifdef __NO_RTTI__
				/* assume it's OK */
				Material2DT* pmat2D = (Material2DT*) pmat;
#else
				Material2DT* pmat2D = dynamic_cast<Material2DT*>(pmat);
				if (!pmat2D) throw eGeneralFail;
#endif				
				thick = pmat2D->Thickness();
			}
			
			/* boundary shape functions */
			const ParentDomainT& surf_shape = fShapes->FacetShapeFunction(facet);
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
					for (int l = 0; l < fNumDOF; l++)
					{
						/* nodal shape function */
						const double* Na = surf_shape.Shape(j);
					
						double* prhs = rhs.Pointer(l);
						double  fact = jwt*(*tj++);
						for (int k = 0; k < nnd; k++)
						{
							*prhs += fact*(*Na++);
							prhs += fNumDOF;
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
					for (int l = 0; l < fNumDOF; l++)
					{
						/* nodal shape function */
						const double* Na = surf_shape.Shape(j);
					
						double* prhs = rhs.Pointer(l);
						double  fact = jwt*(*tj++);
						for (int k = 0; k < nnd; k++)
						{
							*prhs += fact*(*Na++);
							prhs += fNumDOF;
						}
					}				
				}
			}
			else
				throw eGeneralFail;

			/* assemble */
			fFEManager.AssembleRHS(rhs, BC_card.Eqnos());
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
			
			/* matrix form */
			int a = 0, zero = 0;
			int& b_start = (fLHS.Format() == ElementMatrixT::kSymmetricUpper) ? a : zero;
			
			fShapes->TopIP();	
			while ( fShapes->NextIP() )
			{
				double temp = constM*(*Weight++)*(*Det++);
				const double* Na = fShapes->IPShapeU();
								
				for (a = 0; a < nen; a++)
					for (int i = 0; i < fNumDOF; i++)
					{
						int p = a*fNumDOF + i;
						
						/* upper triangle only */
						for (int b = b_start; b < nen; b++)
							for (int j = 0; j < fNumDOF; j++)
								if(i == j)
								{									
									int q = b*fNumDOF + j;
									fLHS(p,q) += temp*Na[a]*Na[b];
								}
					}
			}
			break;
		}

		case kLumpedMass:	/* lumped mass */
		{
			int nen = fLocDisp.NumberOfNodes();

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
				for (int ed = 0; ed < fNumDOF; ed++)
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
	if (fBodyForceLTf > -1)
	{
		int nen = body_force.NumberOfNodes();
		double loadfactor = fFEManager.LoadFactor(fBodyForceLTf);
		double* p = body_force.Pointer();

		for (int i = 0; i < fNumDOF; i++)
		{
			double temp = -fBody[i]*loadfactor;
		
			for (int j = 0; j < nen; j++)
				*p++ = temp;
		}
	}
}

/* calculate the body force contribution */
void ContinuumElementT::FormMa(int mass_type, double constM, const LocalArrayT& body_force)
{
	switch (mass_type)
	{
		case kConsistentMass:	
		{
#if __option(extended_errorcheck)
			if (fRHS.Length() != body_force.Length()) throw eSizeMismatch;
#endif
			int nen = body_force.NumberOfNodes();

			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();

			fShapes->TopIP();
			while ( fShapes->NextIP() )
			{					
				/* integration point accelerations */
				fShapes->InterpolateU(body_force, fDOFvec);

				/* accumulate in element residual force vector */				
				double*	res      = fRHS.Pointer();
				const double* Na = fShapes->IPShapeU();

				double temp = constM*(*Weight++)*(*Det++);				
				for (int lnd = 0; lnd < nen; lnd++)
				{
					double  temp2 = temp*(*Na++);
					double*  pacc = fDOFvec.Pointer();
					
					for (int dof = 0; dof < fNumDOF; dof++)			
						*res++ += temp2*(*pacc++);
				}
			}
			break;
		}	
		case kLumpedMass:
		{
			//cout << "\n ContinuumElementT::FormMa: inertial forces with lumped mass not supported";
			//cout << endl;
			//throw eGeneralFail;
			
			//for now, no inertial force for lumped mass
			//but should probably generalize the FormMass and
			//FormStiffness routines by passing in a target object
			//in which to place the data

#if __option(extended_errorcheck)
			if (fLHS.Rows() != body_force.Length()) throw eSizeMismatch;
#endif
			
			fLHS = 0.0; //hope there's nothing in there!
			FormMass(kLumpedMass, constM);
			body_force.ReturnTranspose(fNEEvec);

			double* pAcc = fNEEvec.Pointer();
			double* pRes = fRHS.Pointer();
			int     massdex = 0;
			
			int nee = body_force.Length();
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
	for (int i = 0; i < fBlockData.MajorDim(); i++)
		if (fBlockData(i, kBlockMat) < 0 ||
		    fBlockData(i, kBlockMat) >= size)
		{
			cout << "\n ContinuumElementT::ReadMaterialData: material number "
			     << fBlockData(i, kBlockMat) + 1 << '\n';
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
	/* read LTf and force vector */
	in >> fBodyForceLTf >> fBody;		
	if (fBodyForceLTf < 0 || fBodyForceLTf > fFEManager.NumberOfLTf())
		throw eBadInputValue;
	fBodyForceLTf--;

	/* no LTf => no body force */
	if (fBodyForceLTf < 0) fBody = 0.0;
	
	out << "\n Body force vector:\n";
	out << " Body force load-time function number. . . . . . = " << fBodyForceLTf + 1<< '\n';
	out << " Body force vector components:\n";
	for (int j = 0 ; j < fNumDOF; j++)
	{
		out << "   x[" << j+1 << "] direction. . . . . . . . . . . . . . . . = ";
		out << fBody[j] << '\n';
	}
	out.flush();   	   	
}

void ContinuumElementT::EchoTractionBC(ifstreamT& in, ostream& out)
{
	out << "\n Traction boundary conditions:\n";
	
	/* dispatch */
	switch (fFEManager.InputFormat())
	{
		case IOBaseT::kTahoe:
			EchoTractionBC_ASCII(in, out);
			break;

		case IOBaseT::kTahoeII:
			EchoTractionBC_TahoeII(in, out);
			break;

		case IOBaseT::kExodusII:
			EchoTractionBC_ExodusII(in, out);
			break;

		default:

			cout << "\n ContinuumElementT::EchoTractionBC: unsupported input format: ";
			cout << fFEManager.InputFormat() << endl;
			throw eGeneralFail;
	}
	
	if (fNumSD != fNumDOF)
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

void ContinuumElementT::EchoTractionBC_ASCII(ifstreamT& in, ostream& out)
{
	int num_BC;
	in >> num_BC;		
	out << " Number of traction BC's . . . . . . . . . . . . = " << num_BC << '\n';

	if (num_BC > 0)
	{
		int num_sets = 1;
		if (fBlockData.MajorDim() > 1) in >> num_sets;

		/* allocate */
		fTractionList.Allocate(num_BC);
		fTractionList[0].WriteHeader(out, fNumDOF);

		/* external file */
		ifstreamT tmp;
		ifstreamT& in2 = fFEManager.OpenExternal(in, tmp, out, true,
			"ContinuumElementT::EchoTractionBC_ASCII: could not open file");

		/* echo */
		int count = 0;
		for (int k = 0; k < num_sets; k++)
		{
			/* block to group element number */
			int offset = 0;
			int size = fNumElements;
			int num = num_BC;
			int ID;
			if (fBlockData.MajorDim() > 1)
			{
				in2 >> ID >> num;
				
				/* correct offset */
				ID--;
			
				/* check */
				if (count + num > fNumElements)
				{
					cout << "\n ContinuumElementT::EchoTractionBC_ASCII: block size " << num << '\n';
					cout <<   "     exceeds total BC count " << num_BC << endl;
					throw eBadInputValue;
				}
				
				/* get block data */
				const int* block_data = BlockData(ID);
				offset = block_data[kStartNum ];
				size   = block_data[kBlockDim];
			}
			
			for (int i = 0; i < num; i++)
			{
				/* resolve group element number */
				int elem;
				in2 >> elem;
				elem--;
				if (elem < 0 || elem >= size)
				{
					cout << "\n ContinuumElementT::EchoTractionBC_ASCII: element " << elem + 1;
					if (fBlockData.MajorDim() > 1)
						cout << " in block " << ID << "\n    ";
					cout << " is out of range " << size << endl;
					throw eBadInputValue;
				}
				elem += offset;
				
				/* echo traction specification */
				fTractionList[count].EchoValues(fFEManager, *fShapes, elem, fNumDOF, in, out);

				/* next */
				count++;
			}
		}
		
		/* check */
		if (count != num_BC)
		{
			cout << "\n ContinuumElementT::EchoTractionBC_ASCII: " << count << " cards read from\n";
			cout <<   "     file does not equal total " << num_BC << endl;
			throw eBadInputValue;		
		}
		else out.flush();
		
		// could check for duplicates
	}
}

void ContinuumElementT::EchoTractionBC_TahoeII(ifstreamT& in, ostream& out)
{
	/* number of kinematic BC node sets */
	int num_BC_sets;
	in >> num_BC_sets;
	if (num_BC_sets < 0) throw eBadInputValue;
	out << " Number of traction BC side sets . . . . . . . . = ";
	out << num_BC_sets << "\n\n";

	if (num_BC_sets > 0)
	{
		/* open database */
		ModelFileT model_file;
		model_file.OpenRead(fFEManager.ModelFile());

		/* nodes on element facets */
		iArrayT num_facet_nodes;
		fShapes->NumNodesOnFacets(num_facet_nodes);

		/* temp space */
		iArrayT set_ID(num_BC_sets);
		iArrayT LTf(num_BC_sets);
		ArrayT<Traction_CardT::CoordSystemT> coord_sys(num_BC_sets);
		ArrayT<iArray2DT> side_sets(num_BC_sets);
		ArrayT<dArray2DT> valuesT(num_BC_sets);

		/* echo set specifiers */
		int count = 0;
		for (int i = 0; i < num_BC_sets; i++)
		{
			in >> set_ID[i] >> LTf[i] >> coord_sys[i];
			int num_sides;
			if (model_file.GetSideSetDimensions(set_ID[i], num_sides) !=
			    ModelFileT::kOK) throw eBadInputValue;

			/* increment total count */
			count += num_sides;

			/* skip if empty */
			int num_nodes;
			if (num_sides > 0)
			{
				/* read side set */
				int block_ID;
				iArray2DT& side_set = side_sets[i];
				side_set.Allocate(num_sides, 2);
				if (model_file.GetSideSet(set_ID[i], block_ID, side_set) !=
				    ModelFileT::kOK) throw eBadInputValue;
	
				/* correct offset */
				side_set--;
	
				/* switch to group numbering */
				const int* block_data = BlockData(block_ID);
				iArrayT elems(side_set.MajorDim());
				side_set.ColumnCopy(0, elems);
	
				/* check */
				int min, max;
				elems.MinMax(min, max);
				if (min < 0 || max > block_data[kBlockDim])
				{
					cout << "\n ContinuumElementT::EchoTractionBC_TahoeII: node numbers\n";
					cout <<   "     {"<< min << "," << max << "} are out of range in side";
					cout << " set ID " << set_ID[i] << endl;
					throw eBadInputValue;
				}

				/* shift */
				elems += block_data[kStartNum];
				side_set.SetColumn(0, elems);

				/* all facets in set must have the same number of nodes */
				num_nodes = num_facet_nodes[side_set(0,1)];
				for (int j = 1; j < num_sides; j++)
					if (num_facet_nodes[side_set(j,1)] != num_nodes)
					{
						cout << "\n ContinuumElementT::EchoTractionBC_TahoeII: sides specified\n";
						cout <<   "     in ID " << set_ID[i] << " have differing numbers of nodes";
						cout << endl;
						throw eBadInputValue;
					}
			}
			else
			{
				/* still number of facet nodes */
				int min, max;
				num_facet_nodes.MinMax(min, max);
				if (min != max)
				{
					cout << "\n ContinuumElementT::EchoTractionBC_TahoeII: cannot determine number of\n"
					     <<   "     facet nodes for empty side set with ID " << set_ID[i] << endl;
					throw eBadInputValue;
				}
				else
					num_nodes = min;
			}
				
			/* read traction BC values */
			dArray2DT& valueT = valuesT[i];
			valueT.Allocate(num_nodes, fNumDOF);
			in >> valueT;
			//NOTE - cannot simply clear to the end of the line with empty side sets
			//       because the tractions may be on multiple lines
		}

		/* allocate all traction BC cards */
		fTractionList.Allocate(count);

		/* correct numbering offset */
		LTf--;

		/* echo node sets */
		iArrayT loc_node_nums;
		int dex = 0;
		for (int ii = 0; ii < num_BC_sets; ii++)
		{
			const iArray2DT& side_set = side_sets[ii];
			int num_sides = side_set.MajorDim();
	
			/* write header */
		  	out << " Exodus side set ID. . . . . . . . . . . . . . . = ";
			out << set_ID[ii] << '\n';
			out << " Number of traction BC cards . . . . . . . . . . = ";
			out << num_sides << endl;
			
			if (count > 0)
			{
				out << '\n';
				out << setw(kIntWidth) << "no.";
				fTractionList[0].WriteHeader(out, fNumDOF);
		
				/* set traction BC cards */
				for (int j = 0; j < num_sides; j++)
				{
					out << setw(kIntWidth) << j+1;			

					/* get facet local node numbers */
					fShapes->NodesOnFacet(side_set(j,1), loc_node_nums);

					/* set and echo */
					fTractionList[dex++].EchoValues(fFEManager, side_set(j,0), side_set(j,1), LTf[ii],
						coord_sys[ii], loc_node_nums, valuesT[ii], out);
				}			
				out << endl;
			}
		}
	}
}

void ContinuumElementT::EchoTractionBC_ExodusII(ifstreamT& in, ostream& out)
{
	/* number of kinematic BC node sets */
	int num_BC_sets;
	in >> num_BC_sets;
	if (num_BC_sets < 0) throw eBadInputValue;
	out << " Number of traction BC side sets . . . . . . . . = ";
	out << num_BC_sets << "\n\n";

	if (num_BC_sets > 0)
	{
		/* open database */
		ExodusT database(out);
		database.OpenRead(fFEManager.ModelFile());

		/* nodes on element facets */
		iArrayT num_facet_nodes;
		fShapes->NumNodesOnFacets(num_facet_nodes);

		/* temp space */
		iArrayT set_ID(num_BC_sets);
		iArrayT LTf(num_BC_sets);
		ArrayT<Traction_CardT::CoordSystemT> coord_sys(num_BC_sets);
		ArrayT<iArray2DT> side_sets(num_BC_sets);
		ArrayT<dArray2DT> valuesT(num_BC_sets);

		/* echo set specifiers */
		int count = 0;
		for (int i = 0; i < num_BC_sets; i++)
		{
			in >> set_ID[i] >> LTf[i] >> coord_sys[i];
			int num_sides = database.NumSidesInSet(set_ID[i]);

			/* increment total count */
			count += num_sides;
			
			/* skip if empty */
			int num_nodes;
			if (num_sides > 0)
			{
				/* read side set */
				int block_ID;
				iArray2DT& side_set = side_sets[i];
				side_set.Allocate(num_sides, 2);
				database.ReadSideSet(set_ID[i], block_ID, side_set);

				/* correct offset */
				side_set--;

				/* switch to group numbering */
				const int* block_data = BlockData(block_ID);
				iArrayT elems(side_set.MajorDim());
				side_set.ColumnCopy(0, elems);

				/* check */
				int min, max;
				elems.MinMax(min, max);
				if (min < 0 || max > block_data[kBlockDim])
				{
					cout << "\n ContinuumElementT::EchoTractionBC_Exodus: node numbers\n";
					cout <<   "     {"<< min << "," << max << "} are out of range in side";
					cout << " set ID " << set_ID[i] << endl;
					throw eBadInputValue;
				}

				/* shift */
				elems += block_data[kStartNum];
				side_set.SetColumn(0, elems);

				/* all facets in set must have the same number of nodes */
				num_nodes = num_facet_nodes[side_set(0,1)];
				for (int j = 1; j < num_sides; j++)
					if (num_facet_nodes[side_set(j,1)] != num_nodes)
					{
						cout << "\n ContinuumElementT::EchoTractionBC_Exodus: sides specified\n";
						cout <<   "     in ID " << set_ID[i] << " have differing numbers of nodes";
						cout << endl;
						throw eBadInputValue;
					}
			}
			else
			{
				/* still number of facet nodes */
				int min, max;
				num_facet_nodes.MinMax(min, max);
				if (min != max)
				{
					cout << "\n ContinuumElementT::EchoTractionBC_Exodus: cannot determine number of\n"
					     <<   "     facet nodes for empty side set with ID " << set_ID[i] << endl;
					throw eBadInputValue;
				}
				else
					num_nodes = min;
			}
				
			/* read traction BC values */
			dArray2DT& valueT = valuesT[i];
			valueT.Allocate(num_nodes, fNumDOF);
			in >> valueT;
			//NOTE - cannot simply clear to the end of the line with empty side sets
			//       because the tractions may be on multiple lines
		}

		/* allocate all traction BC cards */
		fTractionList.Allocate(count);

		/* correct numbering offset */
		LTf--;

		/* echo node sets */
		iArrayT loc_node_nums;
		int dex = 0;
		for (int ii = 0; ii < num_BC_sets; ii++)
		{
			const iArray2DT& side_set = side_sets[ii];
			int num_sides = side_set.MajorDim();
	
			/* write header */
		  	out << " Exodus side set ID. . . . . . . . . . . . . . . = ";
			out << set_ID[ii] << '\n';
			out << " Number of traction BC cards . . . . . . . . . . = ";
			out << num_sides << endl;
			
			if (count > 0)
			{
				out << '\n';
				out << setw(kIntWidth) << "no.";
				fTractionList[0].WriteHeader(out, fNumDOF);
		
				/* set traction BC cards */
				for (int j = 0; j < num_sides; j++)
				{
					out << setw(kIntWidth) << j+1;			

					/* get facet local node numbers */
					fShapes->NodesOnFacet(side_set(j,1), loc_node_nums);

					/* set and echo */
					fTractionList[dex++].EchoValues(fFEManager, side_set(j,0), side_set(j,1), LTf[ii],
						coord_sys[ii], loc_node_nums, valuesT[ii], out);
				}			
				out << endl;
			}
		}
	}
}

/* return the "bounding" elements and the corresponding
* neighbors, both dimensioned internally */
void ContinuumElementT::BoundingElements(iArrayT& elements, iArray2DT& neighbors) const
{
	//TEMP - not parallelized
	if (fFEManager.Size() > 1)
		cout << "\n ContinuumElementT::BoundingElements: not extended to parallel" << endl;

	/* build element neighbor list */
	iArray2DT nodefacetmap;
	fShapes->NeighborNodeMap(nodefacetmap);
	EdgeFinderT edger(fConnectivities, nodefacetmap);
	const iArray2DT& all_neighbors = edger.Neighbors();

	/* collect list of bounding elements */
	AutoArrayT<int> borders;
	iArrayT element;
	for (int i = 0; i < fConnectivities.MajorDim(); i++)
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
		eqnos.Allocate(fNumDOF*nnd);
		
		/* get from node manager */
		nd_tmp.Set(1, nnd, nodes.Pointer());
		eq_tmp.Set(1, fNumDOF*nnd, eqnos.Pointer());
		fNodes->SetLocalEqnos(nd_tmp, eq_tmp);
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
