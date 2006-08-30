/* $Id: CSEAnisoNodal.cpp,v 1.2 2006-08-30 17:49:21 thao Exp $ */
/* created: paklein (11/19/1997) */
#include "CSEAnisoNodal.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentElementsConfig.h"
#endif

#ifdef __DEVELOPMENT_COLORADO__
#include "DevelopmentColoradoElementsConfig.h"
#endif

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "SurfaceShapeT.h"
#include "SurfacePotentialT.h"
#ifndef _FRACTURE_INTERFACE_LIBRARY_
#include "eIntegratorT.h"
#include "NodeManagerT.h"
#endif
#include "ElementSupportT.h"
#include "dSymMatrixT.h"
#include "ParameterContainerT.h"

/* potential functions */
#ifndef _FRACTURE_INTERFACE_LIBRARY_
#include "TvergHutch2DT.h"
//#include "LinearDamage2DT.h"
#endif

using namespace Tahoe;

const double Pi = acos(-1.0);
/* vector functions */
inline static void CrossProduct(const double* A, const double* B, double* AxB)
{   
	AxB[0] = A[1]*B[2] - A[2]*B[1];
	AxB[1] = A[2]*B[0] - A[0]*B[2];
	AxB[2] = A[0]*B[1] - A[1]*B[0];

};

inline static double Dot(const double* A, const double* B)
{
		return(A[0]*B[0]+A[1]*B[1]+A[2]*B[2]);
};

const double perm3D[3][3][3] = { 0, 0, 0, // 1
	                              0, 0, 1,
	                              0,-1, 0,
                                      0, 0,-1, // 2
	                              0, 0, 0, 
	                              1, 0, 0,
                                      0, 1, 0, // 3
	                             -1, 0, 0,
	                              0, 0, 0};

const double perm2D[2][2] = { 0, 1, 
								-1, 0};
#ifndef _FRACTURE_INTERFACE_LIBRARY_
/* constructor */
CSEAnisoNodal::CSEAnisoNodal(const ElementSupportT& support):
	CSEAnisoT(support)
{
	SetName("CSE_rigid_aniso_nodal_int");
}
#else
CSEAnisoNodal::CSEAnisoNodal(ElementSupportT& support, bool rotate):
	CSEAnisoT(support, rotate)
{
	SetName("CSE_rigid_aniso_nodal_int");

	/* reset format for the element stiffness matrix */
	if (fRotate) fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
}
#endif


void CSEAnisoNodal::CloseStep(void)
{
	fStateVariables_n = fStateVariables;
}

/* resets to the last converged solution */
GlobalT::RelaxCodeT CSEAnisoNodal::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = CSEBaseT::ResetStep();
	fStateVariables = fStateVariables_n;
		
	return relax;
}

#ifndef _FRACTURE_INTERFACE_LIBRARY_
/* write restart data to the output stream. */
void CSEAnisoNodal::WriteRestart(ostream& out) const
{
	/* inherited */
	CSEBaseT::WriteRestart(out);
	
	out << fStateVariables.wrap_tight(5) << '\n';	
}

/* read restart data to the output stream */
void CSEAnisoNodal::ReadRestart(istream& in)
{
	/* inherited */
	CSEBaseT::ReadRestart(in);

	/* state variables */
	in >> fStateVariables;
	fStateVariables_n = fStateVariables;
}
#endif

/* describe the parameters needed by the interface */
void CSEAnisoNodal::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	CSEBaseT::DefineParameters(list);

	ParameterT rotate_frame(fRotate, "rotate_frame");
	rotate_frame.SetDefault(fRotate);
	list.AddParameter(rotate_frame);
}

/* information about subordinate parameter lists */
void CSEAnisoNodal::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	CSEBaseT::DefineSubs(sub_list);

	/* element block/material specification */
	sub_list.AddSub("rigid_aniso_CSE_element_block", ParameterListT::OnePlus);
}

/* return the description of the given inline subordinate parameter list */
void CSEAnisoNodal::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "cohesive_relations")
	{
		/* choice */
		order = ParameterListT::Choice;
		
		/* function types */
		sub_lists.AddSub("rigid_cohesive_relation_2D");
//		sub_lists.AddSub("rigid_cohesive_relation_3D");  /*add later*/
	}
	else /* inherited */
		CSEBaseT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT */
ParameterInterfaceT* CSEAnisoNodal::NewSub(const StringT& name) const
{
	/* try to construct cohesive relations */
	SurfacePotentialT* surf_pot = SurfacePotentialT::New(name);
	if (surf_pot)
		return surf_pot;

	if (name == "rigid_aniso_CSE_element_block")
	{
		ParameterContainerT* block = new ParameterContainerT(name);
		
		/* list of element block ID's (defined by ElementBaseT) */
		block->AddSub("block_ID_list", ParameterListT::Once);
	
		/* choice of materials lists (inline) */
		block->AddSub("cohesive_relations", ParameterListT::Once, true);
	
		/* set this as source of subs */
		block->SetSubSource(this);
		
		return block;
	}
	else if (name == "rigid_cohesive_relation_2D")
	{
		/* choice of 2D cohesive relations */
		ParameterContainerT* cz = new ParameterContainerT(name);
		cz->SetSubSource(this);
		cz->SetListOrder(ParameterListT::Choice);
	
		/* choices */
//		cz->AddSub("Linear_Damage_2D");
		cz->AddSub("Tvergaard-Hutchinson_2D");
		
		return cz;
	}
/*	else if (name == "rigid_cohesive_relation_3D")
	{
		ParameterContainerT* cz = new ParameterContainerT(name);
		cz->SetSubSource(this);
		cz->SetListOrder(ParameterListT::Choice);
	
		cz->AddSub("Linear_Damage_3D");
	
		return cz;	
	}*/ // add later
	else /* inherited */
		return CSEBaseT::NewSub(name);
}

/* accept parameter list */
void CSEAnisoNodal::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "CSEAnisoNodal::TakeParameterList";
	
	/* inherited */
	CSEBaseT::TakeParameterList(list);

	int nfn = NumFacetNodes();
	int nen = NumElementNodes();
	int nsd = NumSD();

	facet1.Dimension(nfn);
	facet2.Dimension(nfn);
	const iArray2DT& nodes_on_facets = fShapes->NodesOnFacets(); /*num_facet x num_facet_nodes  check*/
	nodes_on_facets.RowAlias(0,facet1);
	nodes_on_facets.RowAlias(1,facet2);

	/* dimension work space */
	fQ.Dimension(nsd);
	fdelta.Dimension(nsd);
	fT.Dimension(nsd);
	fddU.Dimension(nsd);

	/* rotating frame */
	fRotate = list.GetParameter("rotate_frame");
	if (fRotate) {
	
		/* reset format for the element stiffness matrix */
		fLHS.SetFormat(ElementMatrixT::kNonSymmetric);

		/* shape functions wrt. current coordinates (linked parent domains) */
		fCurrShapes = new SurfaceShapeT(*fShapes, fLocCurrCoords);
		if (!fCurrShapes) ExceptionT::OutOfMemory(caller);
		fCurrShapes->Initialize();
 		
		/* allocate work space */
		int nee = nen*nsd;
		fnsd_nee_1.Dimension(nsd, nee);
		fnsd_nee_2.Dimension(nsd, nee);
		fdQ.Dimension(nsd);
		for (int k = 0; k < nsd; k++)
			fdQ[k].Dimension(nsd, nee);
	}
	else
		fCurrShapes = fShapes;

	/* construct surface properties - one per block */
	int num_block = list.NumLists("rigid_aniso_CSE_element_block");
	fSurfPots.Dimension(num_block);
	
	for (int i = 0; i < fSurfPots.Length(); i++) {

		/* block information */
		const ParameterListT& block = list.GetList("rigid_aniso_CSE_element_block", i);
		
		/* resolve choices of properties choice by spatial dimension */
		const ParameterListT& mat_list_choice_choice = block.GetListChoice(*this, "cohesive_relations");

		/* resolve material choice */
		const ParameterListT& surf_pot_params = block.GetListChoice(*this, mat_list_choice_choice.Name());

		/* construct material */
		SurfacePotentialT* surf_pot = SurfacePotentialT::New(surf_pot_params.Name());
		if (!surf_pot) ExceptionT::BadInputValue(caller, "could not construct \"%s\"", surf_pot_params.Name().Pointer());
		surf_pot->SetTimeStep(ElementSupport().TimeStep());
		surf_pot->TakeParameterList(surf_pot_params);

		/* keep */
		fSurfPots[i] = surf_pot;
	}
	
	/*assume for now that num_state_vars are of equal length*/
	fNumStateVars = fSurfPots[0]->NumStateVariables();

#ifndef _FRACTURE_INTERFACE_LIBRARY_
	/* check compatibility of constitutive outputs */
	if (fSurfPots.Length() > 1 && fNodalOutputCodes[MaterialData])
		for (int k = 0; k < fSurfPots.Length(); k++)
		{
			const SurfacePotentialT* pot_k = fSurfPots[k];
			for (int i = k+1; i < fSurfPots.Length(); i++)
			{
				const SurfacePotentialT* pot_i = fSurfPots[i];
				if (!SurfacePotentialT::CompatibleOutput(*pot_k, *pot_i))
					ExceptionT::BadInputValue(caller, "incompatible output between potentials %d and %d",
						k+1, i+1);
			}
		}
#endif
	/* initialize state variable space */
	fNumIP = nfn;

	fNodePairMap.Dimension(ElementSupport().NumNodes());
	fNodePairMap = 0;
	int num_node_pairs = 0;
	if (fNumStateVars > 0)
	{	
		/* count cohesive zone nodes */
		int nel = NumElements();
		for (int i = 0; i < nel; i++) 
		{
			const ElementCardT& card = ElementCard(i);
			const iArrayT& nodes = card.NodesU();
			for (int j = 0; j < nfn; j++) 
			{
				int nd1 = nodes[facet1[j]];
				int nd2 = nodes[facet2[j]];
				if (nd1 != nd2 && fNodePairMap[nd1] == 0) /* first time and no self-self constraints */
				{
					fNodePairMap[nd1] = num_node_pairs;
					num_node_pairs++;
				}
			}
		}
			/*allocate space and initialize nodal state variables*/
		if (num_node_pairs <= 0)
			ExceptionT::GeneralFail(caller, "No node pairs found");
	}
		
	fStateVariables.Dimension(num_node_pairs, fNumStateVars); 
	fStateVariables_n.Dimension(num_node_pairs, fNumStateVars); 

#ifndef _FRACTURE_INTERFACE_LIBRARY_
	fSurfPots[0]->InitStateVariables(fStateVariables);
	fSurfPots[0]->InitStateVariables(fStateVariables_n);
#endif		

	/*set weights*/
	fNodalArea.Dimension(nfn);

	/*set shape functions*/
	fNa.Dimension(nfn);

	const dArray2DT& nd_coords = fShapes->ParentCoords();
	int parent_dim = nd_coords.MajorDim();
	fDNa.Dimension(parent_dim, nfn); 
	fDHa.Dimension(parent_dim, nen); 
	fCoords.Dimension(parent_dim);
	fJacobian.Dimension(nsd, parent_dim);

	/*global coordinates of the midplane in local surface ordering*/
	fFacetCoords.SetType(LocalArrayT::kUnspecified);
	fFacetCoords.Dimension(nfn,nsd);
		
	/*global gap (jump) vectors in local surface ordering*/
	fdelta_glob.Dimension(nsd);

	iArrayT jump_vector(nen);
	fShapes->SetJumpVector(jump_vector);
	
	/*operator to interpolate nodal displacements to nodal jump vector*/
	fjump.Dimension(nfn, nen);
	fjump = 0.0;
	
	/*operator to average nodal values; used to calculate midplane of cohesive surface*/
	fmid.Dimension(nfn, nen);
	fmid = 0.0;

	for (int j = 0; j < nfn; j++)
	{
		int node = nodes_on_facets(0,j);
		int pair = nodes_on_facets(1,j);
		fjump(j,node) = jump_vector[node];
		fjump(j,pair) = jump_vector[pair]; 
		
		fmid(j, node) = 0.5;
		fmid(j, pair) = 0.5;
	}
	fBa.Dimension(nen);

	/*workspaces*/
	/*derivatives of jacobian with respect to nodal values u*/
	fId.Dimension(nsd);
	fId.Identity(1.0);
	fperm_dx.Dimension(nsd);
	if (nsd == 3)
	{
		fperm_dx_dxi.Dimension(nsd);
		fperm_dx_deta.Dimension(nsd);
	}
}

/* set the active elements */
void CSEAnisoNodal::SetStatus(const ArrayT<ElementCardT::StatusT>& status)
{
	/* work space */
	dArrayT state;
	dArrayT t_in;
	iArrayT facet1;

	/* loop over elements and initial state variables */
	for (int i = 0; i < fElementCards.Length(); i++)
	{
		/* current element */
		ElementCardT& element = fElementCards[i];
		ElementCardT::StatusT& flag = element.Flag();
		flag = status[i];

		if (flag == ElementCardT::kMarkON)
			flag = ElementCardT::kON;
		else if (flag == ElementCardT::kMarkOFF)
			flag = ElementCardT::kOFF;
	}
}

void CSEAnisoNodal::InitializeStateVars(const iArrayT& elem_nodes, const dArray2DT& traction)
{
	/* get facet 1 nodes (global numbering) */
	fNodes1.Collect(facet1, elem_nodes);
	cout << "\nfNodes1: "<<fNodes1;	
	for (int nd = 0; nd < fNodes1.Length(); nd++)
	{
		int node_pair = fNodePairMap[fNodes1[nd]];
		cout << "\nnode_pair: "<<node_pair;
		const double* ptrac = traction(nd);
		double* pstate = fStateVariables(node_pair);
		if (NumSD() == 2)
		{
			pstate[0] = ptrac[0];
			pstate[1] = ptrac[1];
		}
		else if (NumSD() == 3)
		{
			pstate[0] = ptrac[0];
			pstate[1] = ptrac[1];
			pstate[2] = ptrac[2];
		} 
	}
	cout << "\nfStateVariables: "<<fStateVariables;
}

void CSEAnisoNodal::ResetStateVars(const iArrayT& elem_nodes, dArray2DT& traction)
{
	/* get facet 1 nodes (global numbering) */
	fNodes1.Collect(facet1, elem_nodes);	
	for (int nd = 0; nd < fNodes1.Length(); nd++)
	{
		int node_pair = fNodePairMap[fNodes1[nd]];
		double* ptrac = traction(nd);
		const double* pstate = fStateVariables(node_pair);
		if (NumSD() == 2)
		{
			ptrac[0] = pstate[0];
			ptrac[1] = pstate[1];
		}
		else if (NumSD() == 3)
		{
			ptrac[0] = pstate[0];
			ptrac[1] = pstate[1];
			ptrac[2] = pstate[2];
		} 
	}
}

int CSEAnisoNodal::ComputeTraction(const iArrayT& elem_nodes, const dArray2DT& nodal_stress, dArray2DT& traction)
{
	const char caller[] = "CSEAnisoNodal::ComputeTraction";
	
	int nfn = NumFacetNodes();
	int nsd = NumSD();
	
	if(elem_nodes.Length()/2 != nfn)
		ExceptionT::GeneralFail(caller, "Number of elem facet nodes, %d, exceed cse facet nodes %d", 
			elem_nodes.Length()/2, nfn);
			
//	const iArray2DT& nodes_on_facets = fShapes->NodesOnFacets(); /*num_facet x num_facet_nodes  check*/
//	nodes_on_facets.RowAlias(0,facet1);
//	nodes_on_facets.RowAlias(1,facet2);

	/* get ref geometry (1st facet only) */
	fNodes1.Collect(facet1, elem_nodes);
	fLocInitCoords1.SetLocal(fNodes1);
	fLocCurrCoords.SetLocal(elem_nodes);

		double tmax = 0.0;
		int nmax = 0;
	const dArray2DT& nd_coords = fShapes->ParentCoords();
	for (int nd = 0; nd < nfn; nd ++)
	{
		nd_coords.ColumnCopy(nd, fCoords);
		fShapes->EvaluateShapeFunctions(fCoords, fNa, fDNa);
		
		/*calculate coords at elem midplane*/
		if (fRotate)
			for (int i = 0; i < nfn; i++)
				for (int j = 0; j < nsd; j++)
					fFacetCoords(i,j) = fmid.DotRow(i, fLocCurrCoords(j));
			else
				fFacetCoords = fLocInitCoords1;
		/* coordinate transformations */
		fShapes->Jacobian(fFacetCoords, fDNa, fJacobian);
		double j = ComputeRotation(fJacobian, fQ);	

		double* ptrac = traction(nd);
		if (nsd == 2)
		{
			const double* normal = fQ(1); 
			const double* s1 = nodal_stress(facet1[nd]);
			const double* s2 = nodal_stress(facet2[nd]);
			double s11 = 0.5*(s1[0] + s2[0]);
			double s22 = 0.5*(s1[1] + s2[1]);
			double s12 = 0.5*(s1[2] + s2[2]);
				
			ptrac[0] = s11*normal[0] + s12*normal[1];
			ptrac[1] = s12*normal[0] + s22*normal[1];
			double sense = ptrac[0]*normal[0] + ptrac[1]*normal[1];
			if (sense < 0.0)
			{
				ptrac[0] = 0.0;
				ptrac[1] = 0.0;
			}
			
			double tmag = ptrac[0]*ptrac[0] + ptrac[1]*ptrac[1];
			if (tmag > tmax)
			{
				tmax = tmag;
				nmax = nd;
			}
		}
		else if (nsd == 3)
		{
			const double* normal = fQ(2);
				
			const double* s1 = nodal_stress(facet1[nd]);
			const double* s2 = nodal_stress(facet2[nd]);
			double s11 = 0.5*(s1[0] + s2[0]);
			double s22 = 0.5*(s1[1] + s2[1]);
			double s33 = 0.5*(s1[2] + s2[2]);
			double s23 = 0.5*(s1[3] + s2[3]);
			double s13 = 0.5*(s1[4] + s2[4]);
			double s12 = 0.5*(s1[5] + s2[5]);
				
			ptrac[0] = s11*normal[0] + s12*normal[1] + s13*normal[2];
			ptrac[1] = s12*normal[0] + s22*normal[1] + s23*normal[2];
			ptrac[2] = s13*normal[0] + s23*normal[1] + s33*normal[2];

			double sense = ptrac[0]*normal[0] + ptrac[1]*normal[1] +  ptrac[2]*normal[2];
			if (sense < 0.0)
			{
				ptrac[0] = 0.0;
				ptrac[1] = 0.0;
				ptrac[2] = 0.0;
			}

			double tmag = ptrac[0]*ptrac[0] + ptrac[1]*ptrac[1] + ptrac[2]*ptrac[2];
			if (tmag > tmax)
			{
				tmax = tmag;
				nmax = nd;
			}
		}				
	}
	return(nmax);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

void CSEAnisoNodal::LHSDriver(GlobalT::SystemTypeT)
{
	const char caller[] = "CSEAnisoNodal::LHSDriver";

	/* matrix format */
	dMatrixT::SymmetryFlagT format = (fRotate) ?
		dMatrixT::kWhole : dMatrixT::kUpperOnly;

	/* time-integration parameters */
	double constK = 1.0;
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;
#endif

	/* node map of facet 1 and 2*/
	int nfn = NumFacetNodes();
	int nsd = NumSD();
	int nen = NumElementNodes();
	
	const iArray2DT& nodes_on_facets = fShapes->NodesOnFacets(); /*num_facet x num_facet_nodes  check*/
	nodes_on_facets.RowAlias(0,facet1);
	
	/* work space for collecting element variables */
	LocalArrayT nodal_values(LocalArrayT::kUnspecified);
	dArray2DT elementVals;
	dArrayT localFrameIP;
	iArrayT ndIndices;

	/**/
	/*dimension work space*/
	const dArray2DT& nd_coords = fShapes->ParentCoords();
//	int parent_dim = nd_coords.MajorDim();
//	fDNa.Dimension(parent_dim, NumFacetNodes()); 
//	fDHa.Dimension(parent_dim,nen);
//	fCoords.Dimension(parent_dim);
//	fJacobian.Dimension(NumSD(), parent_dim);
	
//	AutoArrayT<double> state2;
	dArrayT state;
	Top();
	while (NextElement())
	{
		/* current element */
		const ElementCardT& element = CurrentElement();
		
		if (element.Flag() != ElementCardT::kOFF)
		{
	
		/* surface potential */
		SurfacePotentialT* surfpot = fSurfPots[element.MaterialNumber()];
		
		/* get ref geometry (1st facet only) */
		fNodes1.Collect(facet1, element.NodesX());
		fLocInitCoords1.SetLocal(fNodes1);

		/* get current geometry */
		SetLocalX(fLocCurrCoords); //EFFECTIVE_DVA

		/*calculate coords at elem midplane*/
		if (fRotate)
			for (int i = 0; i < nfn; i++)
				for (int j = 0; j < nsd; j++)
					fFacetCoords(i,j) = fmid.DotRow(i, fLocCurrCoords(j));
		else
			fFacetCoords = fLocInitCoords1;
		
		/*calculate contribution of area to element node*/
		fNodalArea = 0.0;
		/*integration occurs in undeformed configuration*/
		fShapes->TopIP();
		while(fShapes->NextIP())
		{
			double w = fShapes->IPWeight();	
			const double j0 = fShapes->Jacobian();
			if (j0 <= 0.0) ExceptionT::BadJacobianDet(caller);
			
			/* get shape functions */
			fShapes->Shapes(fNa);
			
			/* integrate */
			fNodalArea.AddScaled(w*j0, fNa);			
		}

		/* initialize */
		fLHS = 0.0;

		/* loop over integration points */
		for (int nd = 0; nd < nfn; nd++)
		{		
				/* set state variables */
				if (fNumStateVars > 0) 
				{
					int node_pair = fNodePairMap[fNodes1[nd]];
					double* pstate = fStateVariables(node_pair);
					state.Set(fNumStateVars, pstate);
//					cout << "\nLHS node_pair: "<<node_pair;
//					cout << "\nLHS state0: "<<state;
				}
			
			/*get shape functions and derivatives at nodes.  fNa are either 1 or 0.  fDNa depends on geometry*/
			nd_coords.ColumnCopy(nd, fCoords);
			fShapes->EvaluateShapeFunctions(fCoords, fNa, fDNa);
	
			/* coordinate transformations */
			fShapes->Jacobian(fFacetCoords, fDNa, fJacobian);
			
			if (fRotate)
				double j = ComputeRotation(fDNa, fJacobian, fQ, fdQ);
			else
				double j = ComputeRotation(fJacobian, fQ);
				
			if (fAxisymmetric) {
				ExceptionT::GeneralFail(caller, "Axisymmetric formulation not implemented");
			}
			
			/*calculate fBa*/
			fBa = 0.0;
			for (int i = 0; i < nfn; i++)
				for (int j = 0; j < nen; j++)
					fBa[j] += fjump(i,j)*fNa[i];

			/* gap vector from facet1 to facet2, i.e u2 - u1 */
			fdelta_glob = 0.0;
			for (int j = 0; j < nsd; j++)
				for (int b = 0; b < nen; b++)
					fdelta_glob[j] += fBa[b]*fLocCurrCoords(b,j);

			/* gap vector in local frame */
			fQ.MultTx(fdelta_glob, fdelta);			
			
			/* stiffness in local frame */
			double scale = constK*fNodalArea[nd];
			const dMatrixT& K = surfpot->Stiffness(fdelta, state, localFrameIP);
			
			/* rotation */
			if (fRotate)
			{
				/* traction in local frame */
				const dArrayT& T = surfpot->Traction(fdelta, state, localFrameIP, false);

				/* 1st term */
				u_i__Q_ijk(fdelta_glob, fdQ, fnsd_nee_1);	/*	d_m dQ_mlq								*/
				fnsd_nee_2.MultAB(K, fnsd_nee_1);			/*	K_kl dQ_mlq d_m							*/ 
				fnsd_nee_1.MultAB(fQ, fnsd_nee_2);			/*	fnsd_nee_1(i,q) = Q_ik K_kl dQ_mlq d_m	*/
				Q_ijk__u_j(fdQ, T, fnsd_nee_2);			/*	fnsd_nee_2(i,q) = dQ_ikq Tk				*/
				fddU.MultQBQT(fQ, K);						/*  Q_ik K_kl Q_jl*/

				/*add contributions*/
				for (int a = 0; a < nen; a++)
					for (int i = 0; i < nsd; i++)
					{
						int p = a*nsd + i;
						for (int b = 0; b < nen; b++)
							for (int j = 0; j < nsd; j++)
							{
								int q = b*nsd + j;
								fLHS(p,q) += scale*fBa[a]*(fnsd_nee_1(i,q) + fnsd_nee_2(i,q));
								fLHS(p,q) += scale* fBa[a]*fddU(i,j)*fBa[b];
							}
					}
			}
			else
			{
				fddU.MultQBQT(fQ, K);						/*  Q_ik K_kl Q_jl*/

				/*only one contribution*/
				for (int a = 0; a < nen; a++)
					for (int i = 0; i < nsd; i++)
					{
						int p = a*nsd + i;
						for (int b = 0; b < nen; b++)
							for (int j = 0; j < nsd; j++)
							{
								int q = b*nsd + j;
								fLHS(p,q) += scale*fBa[a]*fddU(i,j)*fBa[b];
							}
					}
			} /*rotation*/
	
/*debugging block*/
//				const dArrayT& traction = surfpot->Traction(fdelta, state, localFrameIP, false);
//				fQ.Multx(traction, fT);
//				cout << "\nLHS node: " << nd;
//				cout << "\nCurrCoords: "<<fLocCurrCoords;
//				cout << "\nfBa: "<<fBa;				
//				cout << "\nlocal delta: "<<fdelta;
//				cout << "\nglobal delta: "<<fdelta_glob;
//				cout << "\nlocal traction: "<< traction;				
//				cout << "\nglobal traction: "<< fT;
//				cout << "\nlocal stiffness: "<<K;
//				cout << "\nglobal STiffnes: "<<fddU;
//				cout << "\nfLHS: "<<fLHS;
//				cout <<"\nRotation Matrix: "<<fQ;
/*end debugging block */

		} /*sum over nfn*/

		/* assemble */
		AssembleLHS();
		
		} /* element.Flag() != kOFF */
	}	/*next element*/
}

void CSEAnisoNodal::RHSDriver(void)
{
	const char caller[] = "CSEAnisoNodal::RHSDriver";

	/* time-integration parameters */
	double constKd = 1.0;

#ifndef _FRACTURE_INTERFACE_LIBRARY_

	int formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* heat source if needed */
	const FieldT* temperature = ElementSupport().Field("temperature");

	/* initialize sources */
	if (temperature) 
		ExceptionT::GeneralFail(caller, "Temperature not enabled.");

#else // _FRACTURE_INTERFACE_LIBRARY_ defined

    /*Read in SIERRA's new state variables. We need their memory. */	
	/*need to fix this*/
	fStateVariables.Set(fStateVariables.MajorDim(),fStateVariables.MinorDim(),
		ElementSupport().StateVariableArray());

#endif

	fStateVariables = fStateVariables_n;

	int nfn = NumFacetNodes();
	int nen = NumElementNodes();
	int nsd = NumSD();

	/* node map of facet 1 and 2*/
	const iArray2DT& nodes_on_facets = fShapes->NodesOnFacets(); /*num_facet x num_facet_nodes  check*/
	nodes_on_facets.RowAlias(0,facet1);
	nodes_on_facets.RowAlias(1,facet2);
	/*dimension work space*/
	dArrayT localFrameIP;
	/* fracture surface area */
	fFractureArea = 0.0;

	const dArray2DT& nd_coords = fShapes->ParentCoords();
//	int parent_dim = nd_coords.MajorDim();
//	fDNa.Dimension(parent_dim, NumFacetNodes()); 
//	fCoords.Dimension(parent_dim);
//	fJacobian.Dimension(NumSD(), parent_dim);
	
	/* fracture surface area */
	fFractureArea = 0.0;

	int block_count = 0, block_dex = 0;
	dArrayT state;
/*
			double* p = fStateVariables.Pointer();
			for (int j = 0; j<fStateVariables.Length(); j++)
				cout << "\nRHS fStateVariables: "<<p[j];
*/
	Top();
	while (NextElement())
	{
		/* advance to block (skip empty blocks) */
		while (block_count == fBlockData[block_dex].Dimension()) {
			block_count = 0;
			block_dex++;		
		}

		/* current element */
		ElementCardT& element = CurrentElement();
	
		/* get ref geometry (1st facet only) */
		fNodes1.Collect(facet1, element.NodesX());
		fLocInitCoords1.SetLocal(fNodes1);
	  			
		fNodalArea = 0.0;
		/*integration occurs in undeformed configuration*/
		fShapes->TopIP();
		while(fShapes->NextIP())
		{
			double w = fShapes->IPWeight();	
			const double j0 = fShapes->Jacobian();
			if (j0 <= 0.0) ExceptionT::BadJacobianDet(caller);
			
			/* get shape functions */
			fShapes->Shapes(fNa);
			
			/* integrate */
			fNodalArea.AddScaled(w*j0, fNa);
			
		}

		if (element.Flag() != ElementCardT::kOFF)
		{
			/* surface potential */
			SurfacePotentialT* surfpot = fSurfPots[element.MaterialNumber()];
	
			/* get current geometry */
			SetLocalX(fLocCurrCoords);
				
			/*calculate coords at elem midplane*/
			if (fRotate)
				for (int i = 0; i < nfn; i++)
					for (int j = 0; j < nsd; j++)
						fFacetCoords(i,j) = fmid.DotRow(i, fLocCurrCoords(j));
			else
				fFacetCoords = fLocInitCoords1;

//			cout<< "\nelement nodes: "<<element.NodesX();
//			cout << "\nfacet1: "<<facet1;
//			cout << "\nfacet2: "<<facet2;
//			cout << "\nfacet_coords: "<<fFacetCoords;
//			cout << "\ndisp: "<<fLocCurrCoords;
		
			/*calculate contribution of area to element node*/
			/* initialize */
	  		fRHS = 0.0;
			
			/* loop over integration points */
			int all_failed = 1;
			for (int nd = 0; nd < nfn; nd ++)
			{
				/* set state variables */
				if (fNumStateVars > 0) 
				{
					int node_pair = fNodePairMap[fNodes1[nd]];
					double* pstate = fStateVariables(node_pair);
					state.Set(fNumStateVars, pstate);
//					cout << "\nnode_pair: "<<node_pair;
//					cout << "\nstate0: "<<state;
				}
				/*get shape functions and derivatives at nodes.  fNa are either 1 or 0.  fDNa depends on geometry*/
				nd_coords.ColumnCopy(nd, fCoords);
				fShapes->EvaluateShapeFunctions(fCoords, fNa, fDNa);

				/* coordinate transformations */
				fShapes->Jacobian(fFacetCoords, fDNa, fJacobian);

				double j = ComputeRotation(fJacobian, fQ);	
				
				if (fAxisymmetric) {
					ExceptionT::GeneralFail(caller, "Axisymmetric formulation not implemented");
				}
				
				/*calculate fBa*/
				fBa = 0.0;
				for (int i = 0; i < nfn; i++)
					for (int j = 0; j < nen; j++)
						fBa[j] += fjump(i,j)*fNa[i];

				/* gap vector from facet1 to facet2, i.e u2 - u1 */
				/* gap vector from facet1 to facet2, i.e u2 - u1 */
				fdelta_glob = 0.0;
				for (int j = 0; j < nsd; j++)
					for (int b = 0; b < nen; b++)
						fdelta_glob[j] += fBa[b]*fLocCurrCoords(b,j);

				/* rotate gap vector in local frame */
				fQ.MultTx(fdelta_glob, fdelta);
				
				/* calculate traction vector in/out of local frame */
				const dArrayT& traction = surfpot->Traction(fdelta, state, localFrameIP, true);
				fQ.Multx(traction, fT);

				double scale = -constKd*fNodalArea[nd];
				
				for (int a = 0; a < nen; a++)
				{
					for (int j = 0; j < nsd; j++)
					{
						int p = a*nsd + j;
						fRHS[p] += scale*fBa[a]*fT[j];
					}
				}
/*debugging block*/
//				cout << "\nRHS node: " << nd;
//				cout << "\nfCoords: "<<fCoords;
//				cout << "\nDna: "<<fDNa;
//				cout << "\njac: "<<fJacobian;
//				cout <<"\nrotation: "<<fQ;
//				cout << "\nCurrCoords: "<<fLocCurrCoords;
//				cout << "\nfBa: "<<fBa;				
//				cout << "\nlocal delta: "<<fdelta;
//				cout << "\nglobal delta: "<<fdelta_glob;
//				cout << "\nstate: "<<state;
//				cout << "\nlocal traction: "<< traction;				
//				cout << "\nglobal traction: "<< fT;
//				cout << "\nRHS: "<<fRHS;
//				cout <<"\nRotation Matrix: "<<fQ;
/*end debugging block*/
				/* check status */
				SurfacePotentialT::StatusT status = surfpot->Status(fdelta, state);
				if (status != SurfacePotentialT::Failed) all_failed = 0;
				
				/* fracture area */
				if (fOutputArea && status != SurfacePotentialT::Precritical)
					fFractureArea += fNodalArea[nd];
			}

			/* assemble */
			AssembleRHS();
			
			/* mark elements */
			if (all_failed)
			{
				ElementCardT::StatusT& flag = element.Flag();
				if (flag == ElementCardT::kON) flag = ElementCardT::kMarked;
			}
		}
		else if (fOutputArea)
		{
			/* integrate fracture area */
			for (int nd = 0; nd < fNumIP; nd ++)
			{
				fFractureArea += fNodalArea[nd];
				if (fAxisymmetric) {
					ExceptionT::GeneralFail(caller, "Axisymmetric formulation not implemented");
				}
			}
		}

		/* next in block */
		block_count++;
	}
	/*fix this*/
	fSurfPots[0]->UpdateStateVariables(fStateVariables);
//	cout<< "\nfStateVariables;
}


/******private*********/
double  CSEAnisoNodal::ComputeRotation(const dMatrixT& jacobian, dMatrixT& Q)
{
	const char caller[] = "CSEAnisoNodal::ComputeRotation";
	/*compute Q*/
	/* surface dimension */
	if (NumSD() == 2)
	{
		const double* t = jacobian.Pointer();
		double  j = sqrt(t[0]*t[0] + t[1]*t[1]);

		/* check */
		if (j <= 0.0) ExceptionT::BadJacobianDet(caller);

		/* column vectors */
		double* n1 = Q(0);
		double* n2 = Q(1);
		n1[0] = t[0]/j; // n1: tangent
		n1[1] = t[1]/j;

		n2[0] = n1[1];  // n2: normal (rotate -pi/2)
		n2[1] =-n1[0];
		
		return j;
	}
	else /*numsd = 3*/
	{
		/* column vectors */
		double* n1 = Q(0);
		double* n2 = Q(1);
		double* n3 = Q(2);
		
		const double* m1 = jacobian(0);
		const double* m2 = jacobian(1);
		CrossProduct(m1, m2, n3);
		
		double jn = sqrt(n3[0]*n3[0] + n3[1]*n3[1] + n3[2]*n3[2]);
		double j1 = sqrt(m1[0]*m1[0] + m1[1]*m1[1] + m1[2]*m1[2]);

		/* normalize */
		if (jn <= 0.0) ExceptionT::BadJacobianDet(caller);
		n3[0] /= jn;
		n3[1] /= jn;
		n3[2] /= jn;
		
		if (j1 <= 0.0) ExceptionT::BadJacobianDet(caller);
		n1[0] = m1[0]/j1;
		n1[1] = m1[1]/j1;
		n1[2] = m1[2]/j1;
		
		/* orthonormal, in-plane */
		CrossProduct(n3, n1, n2);
		return jn;
	}
} 

double CSEAnisoNodal::ComputeRotation(const dArray2DT& DNa, const dMatrixT& jacobian, dMatrixT& Q, ArrayT<dMatrixT>& dQ)
{
	const char caller[] = "CSEAnisoNodal::ComputeRotation";
#if __option(extended_errorcheck)
	
	if (dQ.Length() != NumSD()) throw ExceptionT::kSizeMismatch;
#endif

	/*compute Q*/
	/*compute Q*/
	/* surface dimension */
	int nnd = NumElementNodes();
	int nfn = NumFacetNodes();
	int nsd = NumSD();
	int npdim = jacobian.Rows();  /*spatial dimension of parent domain*/
	
	const iArray2DT& nodes_on_facets = fShapes->NodesOnFacets(); /*num_facet x num_facet_nodes  check*/
	int num_facets = nodes_on_facets.MajorDim();
	int num_facet_nodes = nodes_on_facets.MinorDim(); 
	double j;

	if (NumSD() == 2)
	{
		const double* t = jacobian.Pointer();
		j = sqrt(t[0]*t[0] + t[1]*t[1]);

		/* check */
		if (j <= 0.0) ExceptionT::BadJacobianDet(caller);

		/* column vectors */
		double* n1 = Q(0);
		double* n2 = Q(1);
		n1[0] = t[0]/j; // n1: tangent
		n1[1] = t[1]/j;

		n2[0] = n1[1];  // n2: normal (rotate -pi/2)
		n2[1] =-n1[0];
		
		
		/*calculate derivative of rotation matrix: dQ_ij/du_ak*/
		/* components of the rank 3 tensor */
		dMatrixT& dQ1 = dQ[0];
		dMatrixT& dQ2 = dQ[1];

		/*calculate fDHa*/
		fDHa = 0.0;
		double* p1 = fDHa(0);
		const double *pD1 = fDNa(0);
		
		for (int i = 0; i < nfn; i++)
			for (int j = 0; j < nnd; j++)
				p1[j] += fmid(i,j)*pD1[i];		
				
		/*calculate derivative of the jacobians with respect to current nodal coords x_ai*/
		const double* dx_dxi = jacobian(0);
		double dj_dxak;
		
		/*ep_mi dx_m/dxi*/
		fperm_dx[0] = dx_dxi[1];
		fperm_dx[1] = -dx_dxi[0];
		
		/*debugging block
		cout << "\nfDHa: "<<fDHa;
		cout << "\njacobian: "<<jacobian;
		cout << "\nfQ: "<<fQ;
		cout << "\ndx: "<< dx_dxi[0]<<"\t"<<dx_dxi[1];
		cout << "\nperm_dx: "<<fperm_dx;
		/*end debugging block*/
		
		double j_inv = 1.0/j;
	
		
		for (int a = 0; a < nnd; a++)
			for (int k = 0; k < nsd; k++)
			{
				int p = a*nsd + k;
				dj_dxak = j_inv*dx_dxi[k]*fDHa[a];
				
				dQ1(0,p) = j_inv*(fDHa[a]*fId(0,k) - j_inv*dj_dxak * dx_dxi[0]);			/*i = 1, j = 1*/
				dQ1(1,p) = j_inv*(fDHa[a]*perm2D[0][k] - j_inv*dj_dxak *fperm_dx[0]); /*i = 1, j = 2*/;

				dQ2(0,p) = j_inv*(fDHa[a]*fId(1,k) - j_inv*dj_dxak * dx_dxi[1]);			/*i = 2, j = 1*/
				dQ2(1,p) = j_inv*(fDHa[a]*perm2D[1][k] - j_inv*dj_dxak *fperm_dx[1]); /*i = 2, j = 2*/;
			}
	}
		
	else /*numsd = 3*/
	{
		/* column vectors */
		double* n1 = Q(0);
		double* n2 = Q(1);
		double* n3 = Q(2);
		
		const double* m1 = jacobian(0);
		const double* m2 = jacobian(1);
		CrossProduct(m1, m2, n3);
		
		double jn = sqrt(n3[0]*n3[0] + n3[1]*n3[1] + n3[2]*n3[2]);
		double j1 = sqrt(m1[0]*m1[0] + m1[1]*m1[1] + m1[2]*m1[2]);

		/* normalize */
		if (jn <= 0.0) ExceptionT::BadJacobianDet(caller);
		n3[0] /= jn;
		n3[1] /= jn;
		n3[2] /= jn;
		
		if (j1 <= 0.0) ExceptionT::BadJacobianDet(caller);
		n1[0] = m1[0]/j1;
		n1[1] = m1[1]/j1;
		n1[2] = m1[2]/j1;
		
		/* orthonormal, in-plane */
		CrossProduct(n3, n1, n2);
		j = jn;
	
		/*calculate derivative of rotation matrix dQ_ij/du_ak*/
		/* components of the rank 3 tensor */
		dMatrixT& dQ1 = dQ[0];
		dMatrixT& dQ2 = dQ[1];
		dMatrixT& dQ3 = dQ[2];

		/*calculate fDHa*/
		fDHa = 0.0;
		double* pDHa1 = fDHa(0);
		double* pDHa2 = fDHa(1);
		const double *pDNa1 = fDNa(0);
		const double *pDNa2 = fDNa(1);
		
		for (int i = 0; i < nfn; i++)
			for (int j = 0; j < nnd; j++)
			{
				pDHa1[j] += fmid(i,j)*pDNa1[i];
				pDHa2[j] += fmid(i,j)*pDNa2[i];
			}	

		const double* dx_dxi = jacobian(0);
		const double* dx_deta = jacobian(1);
		
		double djn_dxak;
		double dj1_dxak;
		
/*		fperm_dx[0] = perm[0][0][0]*dx_dxi[0]*dx_deta[0] + perm[0][0][1]*dx_dxi[0]*dx_deta[1] 
					+ perm[0][0][2]*dx_dxi[0]*dx_deta[2] + perm[0][1][0]*dx_dxi[1]*dx_deta[0] 
					+ perm[0][1][1]*dx_dxi[1]*dx_deta[1] + perm[0][1][2]*dx_dxi[1]*dx_deta[2]
					+ perm[0][2][0]*dx_dxi[2]*dx_deta[0] + perm[0][2][1]*dx_dxi[2]*dx_deta[1] 
					+ perm[0][2][2]*dx_dxi[2]*dx_deta[2];

		fperm_dx[1] = perm[1][0][0]*dx_dxi[0]*dx_deta[0] + perm[1][0][1]*dx_dxi[0]*dx_deta[1] 
					+ perm[1][0][2]*dx_dxi[0]*dx_deta[2] + perm[1][1][0]*dx_dxi[1]*dx_deta[0] 
					+ perm[1][1][1]*dx_dxi[1]*dx_deta[1] + perm[1][1][2]*dx_dxi[1]*dx_deta[2]
					+ perm[1][2][0]*dx_dxi[2]*dx_deta[0] + perm[1][2][1]*dx_dxi[2]*dx_deta[1] 
					+ perm[1][2][2]*dx_dxi[2]*dx_deta[2];

		fperm_dx[1] = perm[2][0][0]*dx_dxi[0]*dx_deta[0] + perm[2][0][1]*dx_dxi[0]*dx_deta[1] 
					+ perm[2][0][2]*dx_dxi[0]*dx_deta[2] + perm[2][1][0]*dx_dxi[1]*dx_deta[0] 
					+ perm[2][1][1]*dx_dxi[1]*dx_deta[1] + perm[2][1][2]*dx_dxi[1]*dx_deta[2]
					+ perm[2][2][0]*dx_dxi[2]*dx_deta[0] + perm[2][2][1]*dx_dxi[2]*dx_deta[1] 
					+ perm[2][2][2]*dx_dxi[2]*dx_deta[2];
*/

/*		fperm_dx_dxi(0,0) = fperm_dx_dxi(1,1) = fperm_dx_dxi(2,2) = 0.0;
		fperm_dx_dxi(0,1) = perm[2][1][0]*dx_dxi[2];
		fperm_dx_dxi(0,2) = perm[1][2][0]*dx_dxi[1];
		fperm_dx_dxi(1,0) = perm[2][0][1]*dx_dxi[2];
		fperm_dx_dxi(1,2) = perm[0][2][1]*dx_dxi[0];
		fperm_dx_dxi(2,0) = perm[1][0][2]*dx_dxi[1];
		fperm_dx_dxi(2,1) = perm[0][1][2]*dx_dxi[0];

		fperm_dx_deta(0,0) = fperm_dx_deta(1,1) = fperm_dx_deta(2,2) = 0.0;
		fperm_dx_deta(0,1) = perm[1][2][0]*dx_deta[2];
		fperm_dx_deta(0,2) = perm[2][1][0]*dx_deta[1];
		fperm_dx_deta(1,0) = perm[0][2][1]*dx_deta[2];
		fperm_dx_deta(1,2) = perm[2][0][1]*dx_deta[0];
		fperm_dx_deta(2,0) = perm[0][1][2]*dx_deta[1];
		fperm_dx_deta(2,1) = perm[1][0][2]*dx_deta[0];
*/
		/*ep_imn dx_m/dxi dx_n/deta*/
		fperm_dx[0] =  dx_dxi[1]*dx_deta[2] - dx_dxi[2]*dx_deta[1];
		fperm_dx[1] =  dx_dxi[2]*dx_deta[0] - dx_dxi[0]*dx_deta[2] ;
		fperm_dx[2] =  dx_dxi[0]*dx_deta[1] - dx_dxi[1]*dx_deta[0];

		for (int a = 0; a < nnd; a++)
			for (int k = 0; k < nsd; k++)
			{
				int p = a*nsd + k;
				double deta2 = Dot(dx_deta, dx_deta);
				double dxi2 = Dot(dx_dxi, dx_dxi);
				double dxideta = Dot(dx_dxi, dx_deta);
				
				dj1_dxak = 1.0/j1* dx_dxi[k]*pDHa1[a];
				djn_dxak = 1.0/jn*( dx_dxi[k]*(deta2*pDHa1[a] - dxideta*pDHa2[a]) + dx_deta[k]*(dxi2*pDHa2[a] 
						- dxideta*pDHa1[a]) );
				
				int n = nsd - k;
				dQ1(0,p) = 1.0/j1*(pDHa1[a]*fId(0,k) - 1.0/j1*dj1_dxak * dx_dxi[0]);					/*i = 1, j = 1*/
				dQ1(2,p) = 1.0/jn*(perm3D[0][k][n]*(pDHa1[a]*dx_deta[n] - pDHa2[a]*dx_dxi[n]) 
					- 1.0/j*djn_dxak *fperm_dx[0]);													/*i = 1, j = 3*/

				n = nsd - (k+1);
				dQ2(0,p) = 1.0/j1*(pDHa1[a]*fId(1,k) - 1.0/j1*dj1_dxak * dx_dxi[1]);					/*i = 2, j = 1*/
				dQ2(2,p) = 1.0/jn*(perm3D[1][k][n]*(pDHa1[a]*dx_deta[n] - pDHa2[a]*dx_dxi[n]) 
					- 1.0/j*djn_dxak *fperm_dx[1]);													/*i = 1, j = 3*/

				n = nsd - (k+2);
				dQ3(0,p) = 1.0/j1*(pDHa1[a]*fId(2,k) - 1.0/j1*dj1_dxak * dx_dxi[2]);					/*i = 2, j = 1*/
				dQ3(2,p) = 1.0/jn*(perm3D[2][k][n]*(pDHa1[a]*dx_deta[n] - pDHa2[a]*dx_dxi[n])
					- 1.0/j*djn_dxak *fperm_dx[2]);													/*i = 1, j = 3*/
			}
			
		
		const double* pt1 = Q(0);
		const double* pt3 = Q(2);
		/*= ep_imn * ( dt3_m/dx_ak t1_n + t3_m dt1_n/dx_ak)*/
		/*= ep_imn * ( dQ[m](3,p) Q(n,1) + Q(m,3) dQ[n](1,p); for p = (a-1)*nsd + k*/
		for (int p = 0; p < nnd*nsd; p++)   
		{
			dQ1(1,p) =  dQ2(2,p)*pt1[2] - dQ3(2,p)*pt1[1] + pt3[1]*dQ3(0,p) - pt3[2]*dQ2(0,p);      /*i = 1, j = 2*/
			dQ2(1,p) =  dQ3(2,p)*pt1[0] - dQ1(2,p)*pt1[2] + pt3[2]*dQ1(0,p) - pt3[0]*dQ3(0,p);		/*i = 2, j = 2*/	
			dQ3(1,p) =  dQ1(2,p)*pt1[1] - dQ2(2,p)*pt1[0] + pt3[0]*dQ2(0,p) - pt3[1]*dQ1(0,p);      /*i = 3, j = 2*/

		}
		
		return(j);
	}
}
