/* $Id: UpLagAdaptiveT.cpp,v 1.12 2006-05-30 20:13:22 tdnguye Exp $ */
#include "UpLagAdaptiveT.h"

/* requires cohesive surface elements */
#ifdef COHESIVE_SURFACE_ELEMENT

#include "ofstreamT.h"
#include "CSEAnisoT.h"
#include "AutoFill2DT.h"
#include "TiedNodesT.h"
#include "FieldT.h"
#include "ShapeFunctionT.h"
#include "SolidMaterialT.h"

using namespace Tahoe;

/* constructor */
UpLagAdaptiveT::UpLagAdaptiveT(const ElementSupportT& support):
	UpdatedLagrangianT(support),
	fCSE(NULL),
	fTied(NULL),
	fReleaseThreshold(-1.0)
{
	SetName("updated_lagrangian_adaptive_insertion");
}

/* describe the parameters needed by the interface */
void UpLagAdaptiveT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	UpdatedLagrangianT::DefineParameters(list);

	/* cohesive element group */
	list.AddParameter(ParameterT::Integer, "cohesive_element_group");

	/* release threshold */
	ParameterT release(ParameterT::Double, "release_threshold");
	release.AddLimit(0, LimitT::Lower);
	list.AddParameter(release);
}

/* accept parameter list */
void UpLagAdaptiveT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "UpLagAdaptiveT::TakeParameterList";

	/* check */
	if (ElementSupport().Size() != 1) ExceptionT::GeneralFail(caller, "serial only");

	/* inherited */
	UpdatedLagrangianT::TakeParameterList(list);
	
	/* resolve CSE class */
	int cse_group = list.GetParameter("cohesive_element_group");
	cse_group--;
	//ElementBaseT& element_base = ElementSupport().ElementGroup(1);
	ElementBaseT& element_base = ElementSupport().ElementGroup(cse_group);
	fCSE = TB_DYNAMIC_CAST(CSEAnisoT*, &element_base);
	if (!fCSE) ExceptionT::BadInputValue(caller, "could not resolve CSE group at %d", cse_group+1);

	/* release threshold (t > 0) */
	fReleaseThreshold = list.GetParameter("release_threshold");

	/* get nodes used by the CSE element group */
	fCSE->NodesUsed(fCSENodesUsed);
	fGlobalToLocal.SetMap(fCSENodesUsed);
	fGlobalToLocal.SetOutOfRange(InverseMapT::MinusOne);

	/* collect CSE connectivites */
	AutoArrayT<const iArray2DT*> connects;
	fCSE->ConnectsX(connects);
	fConnectivitiesCSELocal.Dimension(fCSE->NumElements(), fCSE->NumElementNodes());
	int count = 0;
	for (int i = 0; i < connects.Length(); i++) {
//		cout << "\nconnects: "<<*(connects[i]);
		fConnectivitiesCSELocal.BlockRowCopyAt(*(connects[i]), count);
		count += connects[i]->MajorDim();
	}
	/* renumber in locally */
	int *pconn = fConnectivitiesCSELocal.Pointer();
	for (int i = 0; i < fConnectivitiesCSELocal.Length(); i++) {
		*pconn = fGlobalToLocal.Map(*pconn);
		pconn++;
	}

//	cout << "\nLocal Connects: "<<fConnectivitiesCSELocal;
	
	fCSEActive.Dimension(fConnectivitiesCSELocal.MajorDim());
	fCSEActive = ElementCardT::kOFF;

	/* generated inverse connectivities */
	int num_chunks = fCSENodesUsed.Length()/1000;
	num_chunks = (num_chunks < 1) ? 1 : num_chunks;
	int max_minordim = (NumSD() > 2) ? 20 : 8;
	AutoFill2DT<int> inv_conn(fCSENodesUsed.Length(), num_chunks, 25, max_minordim);
	int nen_cse = fConnectivitiesCSELocal.MinorDim();
	for (int i = 0; i < fConnectivitiesCSELocal.MajorDim(); i++)
	{
		int* pelem = fConnectivitiesCSELocal(i);
		for (int j = 0; j < nen_cse; j++)
			inv_conn.Append(i, *pelem++);
	}
	fInverseConnectivitiesCSELocal.CopyCompressed(inv_conn);
	inv_conn.Free();
	
/*
	cout << "\nfInverseConnectivities: ";
	for (int i = 0; i < fInverseConnectivitiesCSELocal.MajorDim(); i++)
	{
		cout << "\n";
		for (int j = 0; j < fInverseConnectivitiesCSELocal.MinorDim(i); j++)
			cout <<"\t"<< fInverseConnectivitiesCSELocal(i,j);
	}
*/
	/* get field (cast away const-ness) */
	FieldT* field = (FieldT*) &(Field());

	/* create tied node constraint */
	fTied = new TiedNodesT(ElementSupport(), *field);
	const ArrayT<int>* ex_nodes = ElementSupport().ExternalNodes();
	if (ex_nodes) fTied->SetExternalNodes(*ex_nodes);
	field->AddKBCController(fTied);

	/* set tied node constraints */
	SetNetwork(fCSEActive);

	/* space needed to compute nodal values (during form RHS) */
	fAvgCount.Dimension(fCSENodesUsed.Length());
	fAvgCount = 0;
	fNodalValues.Dimension(fCSENodesUsed.Length(), dSymMatrixT::NumValues(NumSD()));
	cout << "\nNodal Vals Major: "<<fNodalValues.MajorDim();
	cout << "\nNodal Vals Minor: "<<fNodalValues.MinorDim();
	
	fNodalValues = 0.0;
	fNodalExtrapolation.Dimension(NumElementNodes(), dSymMatrixT::NumValues(NumSD()));
}

/* initialize current time increment */
void UpLagAdaptiveT::RHSDriver(void)
{
	/* reset averaging workspace */
	fAvgCount = 0;
	fNodalValues = 0.0;

	/* inherited */
	UpdatedLagrangianT::RHSDriver();
}

/* calculate the internal force contribution ("-k*d") */
void UpLagAdaptiveT::FormKd(double constK)
{
	const double* Det    = fCurrShapes->IPDets();
	const double* Weight = fCurrShapes->IPWeights();

	/* collect incremental heat */
	bool need_heat = fElementHeat.Length() == fShapes->NumIP();

	fNodalExtrapolation = 0.0;
	fCurrShapes->TopIP();
	while ( fCurrShapes->NextIP() )
	{
		/* strain displacement matrix */
		Set_B(fCurrShapes->Derivatives_U(), fB);

		/* B^T * Cauchy stress */
		const dSymMatrixT& Cauchy = fCurrMaterial->s_ij();
		fB.MultTx(Cauchy, fNEEvec);

		/* accumulate */
		fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);

		/* incremental heat generation */
		if (need_heat) 
			fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();

		/* extrapolate the nodal stresses */
		fShapes->Extrapolate(Cauchy, fNodalExtrapolation);
		int elem = CurrElementNumber();
		int ip = CurrIP();

/*		if ((elem==0 && ip ==0) || (elem == 1577 && ip ==1) || (elem==11068 && ip==2) || (elem==5099 && ip==3))
			cout <<"\nelem: "<<elem<<"\tip: "<<ip<< "\ts22: "<<Cauchy[1];
*/
	}
	
	const iArrayT& element_nodes = CurrentElement().NodesU();
	for (int i = 0; i < element_nodes.Length(); i++)
	{
		int loc_node = fGlobalToLocal.Map(element_nodes[i]);
		if (loc_node > -1) {
			fAvgCount[loc_node]++;
			fNodalValues.AddToRowScaled(loc_node, 1.0, fNodalExtrapolation(i));
		}
/*		if(element_nodes[i]==8079)
		{
			cout << "\nelem: "<<CurrElementNumber();
			cout << "\nNodalExtrapolation: "<<fNodalExtrapolation(i,1);
			cout << "\nNodalValues: "<<fNodalValues(loc_node,1);
		}
*/
	}

}

/* element level reconfiguration for the current time increment */
GlobalT::RelaxCodeT UpLagAdaptiveT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax_code = UpdatedLagrangianT::RelaxSystem();

	/* compute average */
	for (int i = 0; i < fNodalValues.MajorDim(); i++) {
		int count = fAvgCount[i];
		if (count > 0) fNodalValues.ScaleRow(i, 1.0/count);
	}


//TEMP - only 2D for now
if (NumSD() != 2) ExceptionT::GeneralFail("UpLagAdaptiveT::RelaxSystem", "2D only");

	/* get state variable storage array */
	RaggedArray2DT<double>& state_variables = fCSE->StateVariables();

	int release_count = 0;
	double t_mag2 = 0.0;
	const dArray2DT& current_coords = ElementSupport().CurrentCoordinates();
	dSymMatrixT Cauchy(NumSD());
	dArrayT traction(NumSD()), tangent(NumSD()), normal(NumSD());
	ostream& out = ElementSupport().Output();
	for (int i = 0; i < fCSEActive.Length(); i++)
	  {
		if (fCSEActive[i] == ElementCardT::kOFF) /* only test rigid surfaces */
		{
			int* pface = fConnectivitiesCSELocal(i);
			int n1 = fCSENodesUsed[pface[0]];
			int n2 = fCSENodesUsed[pface[1]];
			
			/* face tangent */
			tangent[0] = current_coords(n2,0) - current_coords(n1,0);
			tangent[1] = current_coords(n2,1) - current_coords(n1,1);
			
			double mag = sqrt(tangent[0]*tangent[0]+tangent[1]*tangent[1]);
			tangent[0] /= mag;
			tangent[1] /= mag;
			
			/* face normal */
			normal[0] =-tangent[1];
			normal[1] = tangent[0];
			
			double tmax = 0.0;
			int jmax = 0;
			for (int j = 0; j < fConnectivitiesCSELocal.MinorDim(); j++)
			{
				double* s = fNodalValues(pface[0]);
				traction[0] = s[0]*normal[0] + s[2]*normal[1];
				traction[1] = s[2]*normal[0] + s[1]*normal[1];
				double t_mag2 = traction[0]*traction[0] + traction[1]*traction[1];

				if(tmax < sqrt(t_mag2)) {
					tmax =  sqrt(t_mag2);
					jmax = j;
				}
			}
			
			/* compute average of stresses */
/*
			double* s0 = fNodalValues(pface[0]);
			double* s1 = fNodalValues(pface[1]);
			double* s2 = fNodalValues(pface[2]);
			double* s3 = fNodalValues(pface[3]);

			Cauchy[0] = (s0[0] + s1[0] + s2[0] + s3[0])*0.25;
			Cauchy[1] = (s0[1] + s1[1] + s2[1] + s3[1])*0.25;
			Cauchy[2] = (s0[2] + s1[2] + s2[2] + s3[2])*0.25;
			
			/* traction */
/*			traction[0] = Cauchy(0,0)*normal[0] + Cauchy(0,1)*normal[1];
			traction[1] = Cauchy(1,0)*normal[0] + Cauchy(1,1)*normal[1];
			double t_mag2 = traction[0]*traction[0] + traction[1]*traction[1];
*/		
			/* "sense" */
/*			double sense = traction[0]*normal[0] + traction[1]*normal[1];

			if(tmax < sqrt(t_mag2))
				tmax =  sqrt(t_mag2);
*/
			double* s = fNodalValues(pface[jmax]);
			traction[0] = s[0]*normal[0] + s[2]*normal[1];
			traction[1] = s[2]*normal[0] + s[1]*normal[1];
			t_mag2 = traction[0]*traction[0] + traction[1]*traction[1];
			double sense = traction[0]*normal[0] + traction[1]*normal[1];
			
			/* tensile release */
			if (t_mag2 > fReleaseThreshold*fReleaseThreshold && sense > 0.0) {
				fCSEActive[i] = ElementCardT::kMarkON;
				release_count++;

				/* write traction into state variables */
//				cout << "\nlength: "<<state_variables.Length();

				if (state_variables.Length() > 0)
				{
					state_variables(i,0) = traction[0];
					state_variables(i,1) = traction[1];
				}
			}
		}
	  }
				out << "\nrelease_count: "<<release_count;
				out << "\ntmax: "<<sqrt(t_mag2);
				cout << "\nrelease_count: "<<release_count;
				cout << "\ntmax: "<<sqrt(t_mag2);

		
	/* relaxation code for CSE release */
	GlobalT::RelaxCodeT cse_relax_code = (release_count > 0) ? 
		GlobalT::kReEQ : GlobalT::kNoRelax;

	/* reset the network */
	if (release_count > 0) SetNetwork(fCSEActive);

	return GlobalT::MaxPrecedence(relax_code, cse_relax_code);
}

/* determine the tied nodes and reset the constraints based on the list of
 * active elements */
void UpLagAdaptiveT::SetNetwork(const ArrayT<ElementCardT::StatusT>& active_elements)
{
	/* determine duplicate nodes */

	FindLeaders(fConnectivitiesCSELocal, active_elements, fSameAs);

	/* collect constrainted pairs (global node numbers) */
	int num_constraints = fSameAs.Length() - fSameAs.Count(-1);
	int pair_num = 0;
	iArrayT follower(num_constraints);
	iArrayT leader(num_constraints);
	for (int i = 0; i < fSameAs.Length(); i++)
	{
		if (fSameAs[i] != -1) /* followers in 1st column */
		{
			follower[pair_num] = fCSENodesUsed[i];
			leader[pair_num] = fCSENodesUsed[fSameAs[i]];
			pair_num++;
		}
	}
//	cout << "\nfollower: "<<follower;
//	cout << "\nleader: "<<leader;
	
	/* reset tied nodes */
	fTied->SetTiedPairs(follower, leader);
	
	/* set the element flags in the CSE group */
	fCSE->SetStatus(active_elements);
}

/* generate the list of leaders for all nodes based on the active element list */
void UpLagAdaptiveT::FindLeaders(const iArray2DT& connects, const ArrayT<ElementCardT::StatusT>& active, iArrayT& same_as) const
{
	const char caller[] = "UpLagAdaptiveT::FindLeaders";
	cout << "\nUpLagAdaptiveT::FindLeaders";
//TEMP assume 4-noded (2D) and 8-noded (3D) CSE only
	if (NumSD() == 2 && fConnectivitiesCSELocal.MinorDim() != 4)
		ExceptionT::GeneralFail(caller, "expecting 4-noded cse in 2D");
	else if (NumSD() == 3 && fConnectivitiesCSELocal.MinorDim() != 8)
		ExceptionT::GeneralFail(caller, "expecting 8-nodes cse in 3D");
	else if (NumSD() == 1)
		ExceptionT::GeneralFail(caller, "1D not implemented");
	
	/* pair node to the nodes in the 1st facet */
//	int pair_2D[] = {2,3};
	int pair_2D[] = {3,2};
	int pair_3D[] = {4,5,6,7};
	int* pair_node = (NumSD() == 2) ? pair_2D : pair_3D;

	/* dimension and initialize */
	int max_node = connects.Max();
	same_as.Dimension(max_node+1);
	same_as = -1;

	/* initialize leader list */
	int nfn = connects.MinorDim()/2;
//	cout << "\nactive length: "<<active.Length();
	for (int i = 0; i < active.Length(); i++)
	if (active[i] == ElementCardT::kOFF)
	{
		const int* pelem = connects(i);
		for (int j = 0; j < nfn; j++)
		{
			int n1 = pelem[j];
			int n2 = pelem[pair_node[j]];
//			cout << "\nn1: "<<n1<<"\nn2: "<<n2;
			if (n1 > n2)
				same_as[n1] = n2;
			else if (n1 < n2)
				same_as[n2] = n1;
		}
	}

	/* set leaders to node number */
	int changed_count = 0;
	bool changed = true;
	while (changed)
	{
		changed = false;
		changed_count++;
		for (int i = 0; i < same_as.Length(); i++)
		{
			int& leader = same_as[i];
			if (leader > -1) {
				int leader_leader = same_as[leader];
				if (leader_leader > -1 && leader_leader < leader)
				{
					leader = leader_leader;
					changed = true;
				}
			}
		}
	}
}

#endif /* COHESIVE_SURFACE_ELEMENT */
