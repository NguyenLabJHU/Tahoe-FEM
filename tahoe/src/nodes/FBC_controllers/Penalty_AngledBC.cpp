/*  $Id:Penalty_AngledBC.h v 1.0 
 * 
 *
 *  Created by vicky on 3/5/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

/*header file*/
#include "Penalty_AngledBC.h"

#include "NodeManagerT.h"
#include "ModelManagerT.h"
#include "LocalArrayT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "iArray2DT.h"
#include "DomainIntegrationT.h"
#include "FieldT.h"
#include "ContinuumElementT.h"
#include "eIntegratorT.h"

using namespace Tahoe;

Penalty_AngledBC::Penalty_AngledBC():
	fLHS(ElementMatrixT::kSymmetric)
{
	SetName("angled_bc");
}

/*define parameters*/
void Penalty_AngledBC::DefineSubs(SubListT& sub_list) const
{
	/*inherited*/
	FBC_ControllerT::DefineSubs(sub_list);
	
	/*input displacement boundary conditions normal to displacement field */
	sub_list.AddSub("side_set_ID_list");
}

void Penalty_AngledBC::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FBC_ControllerT::DefineParameters(list);
	LimitT zero_bound(0.0, LimitT::Lower);

	/*boundary conditions parameters*/
	ParameterT BC_type(ParameterT::Enumeration, "type");
	BC_type.AddEnumeration("fixed", -1);
	BC_type.AddEnumeration("u", 0);
	BC_type.SetDefault(-1);
	list.AddParameter(BC_type);

	ParameterT schedule(ParameterT::Integer, "schedule");
	schedule.SetDefault(0);
	list.AddParameter(schedule);

	ParameterT value(ParameterT::Double, "value");
	value.SetDefault(0.0);
	list.AddParameter(value);


	/*penalty parameter*/
	ParameterT penalty(fK, "penalty_parameter");
	penalty.AddLimit(zero_bound);
	list.AddParameter(penalty);
}


/*read parameters*/
void Penalty_AngledBC::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "Penalty_AngledBC::TakeParameterList";
	
	/*inherited*/
	FBC_ControllerT::TakeParameterList(list);
	
	/*read penalty parameter*/
	fK = list.GetParameter("penalty_parameter");

	/*read BC parameters*/
	int type = list.GetParameter("type");							 /*bc type fixed or applied */
	KBC_CardT::CodeT code = KBC_CardT::int2CodeT(type + 1);			 /*set KBC code based on type*/
	int schedule_no = list.GetParameter("schedule"); schedule_no--;   /*read schedule number*/
	fValue = list.GetParameter("value");                        /*read applied displacement*/
	/* get the schedule */

	if (schedule_no > -1)
	{
		fSchedule = fFieldSupport->Schedule(schedule_no);
	}
	else {
		fSchedule = NULL;
		fValue = 0;
	}	
	/*read the BC cards and count the number of BC*/
	ReadAngledBC(list);
	
	int nsd = fFieldSupport->NumSD();
	int neq = fNumFacetNodes*nsd;
	fLHS.Dimension(neq);
	fRHS.Dimension(neq);
	
	fQ.Dimension(nsd);
	fnormal.Dimension(nsd);
	fjacobian.Dimension(nsd, nsd-1);
}

/*Collects nodes on facets both in local and global numbering and orders them such the that normal always points outwards*/
void Penalty_AngledBC::ReadAngledBC(const ParameterListT& list)
{	
	const char caller[] = "Penalty_AngledBC::ReadAngledBC";

	/*dimension workspace*/
	StringListT::Extract(list.GetList("side_set_ID_list"),  fside_set_IDs);
	int numlist = fside_set_IDs.Length();
	fBC_sides.Dimension(numlist); 
	fBC_global_nodes.Dimension(numlist);
	fBC_local_nodes.Dimension(numlist);
	fBC_eqnos.Dimension(numlist);
	
	fDomain.Dimension(numlist);
	
	/*count the number of  facets and facet nodes*/
	int nsd = fFieldSupport->NumSD();
	int countnodes=0;
	int countfacets=0;
	GeometryT::CodeT geom0;
	ModelManagerT& model_manager = fFieldSupport->ModelManager();
	
	for (int i = 0; i< numlist; i++)
	{
		const StringT& side_ID = fside_set_IDs[i];
		
		/*get element geometry*/
		StringT elemID;   /*element group ID*/
		elemID = model_manager.SideSetGroupID(side_ID);

		/*read and store sideset*/
		iArray2DT& side_set = fBC_sides[i];
		side_set = model_manager.SideSet(side_ID);

		/*get connectivities*/
		const iArray2DT& connectivities = model_manager.ElementGroup(elemID);
		/*number of nodes*/
		int nen = connectivities.MinorDim();
		/*get the geometry of the element */
		GeometryT::CodeT geometry_code = model_manager.ElementGroupGeometry(elemID);
		GeometryBaseT* geometry = GeometryT::New(geometry_code, nen);

		/* get geometry of the facet*/
		ArrayT<GeometryT::CodeT> facet_geom;   /*1 xnumfacets*/
		iArrayT facet_nodes;  /*number of facet nodes 1 x numfacets*/
		geometry->FacetGeometry(facet_geom, facet_nodes);

		if (i==0)
		{	
			fNumFacetNodes = facet_nodes[0];
			geom0 = geometry_code;
			fNumElementNodes = nen;
		}
		else {
			if (geom0 != geometry_code && nen != fNumElementNodes)
				ExceptionT::GeneralFail(caller, " all elements must be same type");
		}

		
		/*count number of nodes and facets*/
		int numfacets  = side_set.MajorDim();
		fBC_global_nodes[i].Dimension(numfacets, fNumFacetNodes);
		fBC_local_nodes[i].Dimension(numfacets, fNumFacetNodes);
		fBC_eqnos[i].Dimension(numfacets, fNumFacetNodes*nsd);
		
		countfacets += numfacets;
		int numfacetnodes = facet_nodes[0];  
//		countnodes += numfacetnodes*numfacets;
	
	}
	fNumFacets_Tot = countfacets;

	/*Dimension and register element work space*/
	fInitCoords.SetType(LocalArrayT::kInitCoords);
	fInitCoords.Dimension(fNumFacetNodes,nsd);
	fFieldSupport->RegisterCoordinates(fInitCoords);

	fDisp.SetType(LocalArrayT::kDisp);
	fDisp.Dimension(fNumFacetNodes,nsd);
	
	fnodes.Dimension(fNumFacetNodes);
	feqnos.Dimension(fNumFacetNodes*nsd);
	fip_disp.Dimension(nsd);
	fnormal.Dimension(nsd);
}

void Penalty_AngledBC::InitialCondition(void)
{
	/*create integration domains for all of the element*/
	int numlist = fside_set_IDs.Length();
	int nsd = fFieldSupport->NumSD();
	
	ModelManagerT& model_manager = fFieldSupport->ModelManager();
	iArray2DT nd_tmp, eq_tmp;
	fField->RegisterLocal(fDisp);


	for (int i = 0; i< numlist; i++)
	{
		const StringT& side_ID = fside_set_IDs[i];
		iArray2DT& side_set = fBC_sides[i];
		iArray2DT& facets_global = fBC_global_nodes[i];    /* number facets x num_facet nodes (global numbering) */
		iArray2DT& facets_local = fBC_local_nodes[i];     /* number facets x num_facet nodes (global numbering) */
		iArray2DT& facets_eqnos = fBC_eqnos[i];			/* number of facets x num facet nodes*nsd */
		
		int numfacets = side_set.MajorDim();
		/*get element geometry*/
		StringT elemID;   /*element group ID*/
		elemID = model_manager.SideSetGroupID(side_ID);

		/*get connectivities*/
		const iArray2DT& connectivities = model_manager.ElementGroup(elemID);

		/*get the geometry of the element */
		GeometryT::CodeT geometry_code = model_manager.ElementGroupGeometry(elemID);
		GeometryBaseT* geometry = GeometryT::New(geometry_code, fNumElementNodes);

		/*number integration points*/
		int group = model_manager.ElementGroupIndex(elemID);
		group--;
		ElementBaseT& element = fFieldSupport->ElementGroup(group);
		ContinuumElementT* continuum = (ContinuumElementT*) &element;
		int nip = continuum->NumIP();

		/*create ParentDomain for the element*/
		fDomain[i] = new DomainIntegrationT(geometry_code, nip, fNumElementNodes);
		
		iArrayT local_nodes(fNumFacetNodes);  /*for each facet*/
		iArrayT global_nodes(fNumFacetNodes); /*for each facet*/
		nd_tmp.Set(1,fNumFacetNodes, global_nodes.Pointer());
		for (int j = 0; j < numfacets; j++)
		{
			int elem = side_set(j,0);
			int side = side_set(j,1);
			/* gets nodes on faces in local numbering*/
			geometry->NodesOnFacet(side, local_nodes);
			facets_local.SetRow(j,local_nodes);

			/*gets nodes on faces in global numbering*/
			global_nodes.Collect(local_nodes, connectivities(elem));
			facets_global.SetRow(j,global_nodes);

			eq_tmp.Set(1,nsd*fNumFacetNodes, facets_eqnos(j));
			fField->SetLocalEqnos(nd_tmp, eq_tmp);
		}
/*
		cout << "\nsideset: "<<side_set;
		cout << "\nlocal nodes: "<<facets_local;
		cout << "\nglobal nodes: "<<facets_global;
		cout << "\neqnos: "<<facets_eqnos;
*/
	}
		
}
void Penalty_AngledBC::ApplyLHS(GlobalT::SystemTypeT sys_type)
{
	const char caller[] = "Penalty_AngledBC::FormLHS";
	ModelManagerT& model_manager = fFieldSupport->ModelManager();

	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	/* dimensions */
	int nsd = fFieldSupport->NumSD();
	int ndof = fField->NumDOF();	
	if (nsd != ndof) 				
		ExceptionT::GeneralFail(caller, " Can only be applied to displacements");
	
	StringT elemID;   /*element group ID*/

	int numlist = fBC_sides.Length();
	for (int i = 0; i < numlist; i++)
	{
		/*get element group number*/
		elemID = model_manager.SideSetGroupID(fside_set_IDs[i]);
		int group = model_manager.ElementGroupIndex(elemID);
		group--;
		const iArray2DT& side_set = fBC_sides[i];
		const iArray2DT& global_nodes = fBC_global_nodes[i];
		const iArray2DT& local_nodes = fBC_local_nodes[i];
		const iArray2DT& local_eqnos = fBC_eqnos[i];
		const DomainIntegrationT* shape = fDomain[i];
		int numfacets = side_set.MajorDim();

		for (int j = 0; j<numfacets; j++)
		{
			int facet = side_set(j,1);
			int elem = side_set(j,0);
			
			/*retrieve local nodes, equation numbers, and coords*/
			global_nodes.RowCopy(j, fnodes);
			local_eqnos.RowCopy(j, feqnos);
			fInitCoords.SetLocal(fnodes);
	
			const ParentDomainT& surf_shape = shape->FacetShapeFunction(facet);
			int facet_nip = surf_shape.NumIP();
			
			const double* w = surf_shape.Weight(); /*integration weights*/
			fLHS = 0.0;
			for (int k = 0; k < facet_nip; k++)
			{
				/*calculate surface normal using domain jacobian*/
				surf_shape.DomainJacobian(fInitCoords, k, fjacobian);
				double detj = surf_shape.SurfaceJacobian(fjacobian,fQ);
				double fe = constK*detj*w[k]*fK;
				const double* normal;
				if (nsd ==2){
					normal = fQ(1);
				}
				else if (nsd ==3){
					normal = fQ(2);
				}
				else {
					ExceptionT::GeneralFail(caller, "The normal is not defined for this geometry");
				}
				/* accumulate */
				const double* Na = surf_shape.Shape(k);
				for (int a = 0; a < fNumFacetNodes; a++)
				{
					/* nodal shape function of ip k*/
					for (int m = 0; m < nsd; m ++)
					{
						int p = a*nsd + m;
						for (int b = 0; b < fNumFacetNodes; b++)
						{
							for (int n = 0; n < nsd; n++)
							{
								int q = b*nsd + n;
								fLHS(p,q) += fe*normal[m]*normal[n]*Na[a]*Na[b];
							}/*loop n over nsd*/
						}/*loop b over nen*/
					}/*loop m over nsd*/
				}/*loop a over nen*/
				
			} /*loop over ip*/						
			/*assemble*/
			fFieldSupport->AssembleLHS(group, fLHS, feqnos);
		} /*loop over facets in sideset*/
		
	} /*loop over sidesets*/
}

void Penalty_AngledBC::ApplyRHS(void)
{
	const char caller[] = "Penalty_AngledBC::FormRHS";

	double constK = 0.0;
	int formK = fIntegrator->FormKd(constK);
	if (!formK) return;

	ModelManagerT& model_manager = fFieldSupport->ModelManager();
	/*retrieve applied normal displacement*/
	double time = fFieldSupport->Time();
	
	/*applied normal displacement*/
	double gn;
	if (fSchedule)
		gn = fValue*fSchedule->Value(time);
	else
		gn = 0;
		
	/* dimensions */
	int nsd = fFieldSupport->NumSD();
	int ndof = fField->NumDOF();	
	if (nsd != ndof) 				
		ExceptionT::GeneralFail(caller, " Can only be applied to displacements");
	
	/*allocate and dimension workspaces*/
	
	StringT elemID;   /*element group ID*/

	int numlist = fBC_sides.Length();
	for (int i = 0; i < numlist; i++)
	{
		/*get element group number*/
		elemID = model_manager.SideSetGroupID(fside_set_IDs[i]);
		int group = model_manager.ElementGroupIndex(elemID);
		group--;
		const iArray2DT& side_set = fBC_sides[i];
		const iArray2DT& global_nodes = fBC_global_nodes[i];
		const iArray2DT& local_nodes = fBC_local_nodes[i];
		const iArray2DT& local_eqnos = fBC_eqnos[i];
		const DomainIntegrationT* shape = fDomain[i];
		int numfacets = side_set.MajorDim();

		for (int j = 0; j<numfacets; j++)
		{
			int facet = side_set(j,1);
			int elem = side_set(j,0);
			
			/*retrieve local nodes, equation numbers, and coords*/

/*			cout << "\nelement: "<<elem;
			cout << "\nfacet: "<<facet;
			cout << "\nfnodes: "<<fnodes<<endl;
			cout << "\ninit coords: "<<fInitCoords<<endl;
			cout << "\ndisp: "<<fDisp<<endl;
*/
			global_nodes.RowCopy(j, fnodes);
			local_eqnos.RowCopy(j, feqnos);
			fInitCoords.SetLocal(fnodes);
			fDisp.SetLocal(fnodes);
			
			const ParentDomainT& surf_shape = shape->FacetShapeFunction(facet);
			int facet_nip = surf_shape.NumIP();
			
			const double* w = surf_shape.Weight(); /*integration weights*/
			fRHS = 0.0;
			for (int k = 0; k < facet_nip; k++)
			{
				/*calculate surface jacobian*/
				surf_shape.DomainJacobian(fInitCoords, k, fjacobian);

				/*integration weight*/
				double detj = surf_shape.SurfaceJacobian(fjacobian,fQ);

				/* interpolate displacements to ip*/
				surf_shape.Interpolate(fDisp,fip_disp,k);

				/*obtain normal from surface jacobian*/	
				if ((nsd <2) || (nsd > 3)){
					ExceptionT::GeneralFail(caller, "The normal is not defined for this geometry");
				}
				fQ.CopyColumn(nsd-1,fnormal);

				/* accumulate */
				/*compute normal compnent of disp*/
				double un = dArrayT::Dot(fip_disp, fnormal);

/*				cout << "\nnormal disp: "<< un;
				cout << "\napplied: "<<gn;
*/
				double fe = (-constK)*detj*w[k]*fK;
				
				for (int m = 0; m < nsd; m++)
				{
					/* nodal shape function */
					const double* Na = surf_shape.Shape(k);
					
					double* prhs = fRHS.Pointer(m);
					for (int a = 0; a < fNumFacetNodes; a++)
					{
						*prhs += fe*(un-gn)*fnormal[m]*(*Na++);
						prhs += nsd;
					}
				}				
				
			} /*loop over ip*/
			/*assemble*/
			fFieldSupport->AssembleRHS(group, fRHS, feqnos);
		} /*loop over facets in sideset*/
		
	} /*loop over sidesets*/

}


