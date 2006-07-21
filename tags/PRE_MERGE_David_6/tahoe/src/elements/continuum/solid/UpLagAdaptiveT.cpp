/* $Id: UpLagAdaptiveT.cpp,v 1.14 2006-06-18 01:08:34 tdnguye Exp $ */
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
	ftractions_elem(LocalArrayT::kUnspecified),
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

GlobalT::RelaxCodeT UpLagAdaptiveT::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = UpdatedLagrangianT::ResetStep();
	GlobalT::RelaxCodeT cse_relax = relax;
	cout << "\nreseting step: ";
	
	/* get state variable storage array */
	RaggedArray2DT<double>& state_variables = fCSE->StateVariables();
	RaggedArray2DT<double>& state_variables_last = fCSE->StateVariables_Last();

	/*set state to last converged step*/
	state_variables = state_variables_last;
	
	ostream& out = ElementSupport().Output();
	int nip = fCSE->NumIP();

	if (state_variables.Length() > 0)
	{
		cout << state_variables.Length();
		for (int i = 0; i < fCSEActive.Length() ; i++)
		{			
			if (fCSEActive[i] == ElementCardT::kON || fCSEActive[i]==ElementCardT::kMarkON) 
			{
				int num_state = state_variables.MinorDim(i)/nip;
				for (int j = 0; j < nip; j++)
				{
					double T1 = state_variables(i,j*num_state);
					double T2 = state_variables(i,j*num_state);
					double sigma_max = sqrt(T1*T1+T2*T2);
					if (sigma_max < kSmall) {
						fCSEActive[i] = ElementCardT::kOFF;
//						cout << "\nRe-tie element: "<< i;
						cse_relax = GlobalT::kReEQ;
					}
				}
			}
		}
		SetNetwork(fCSEActive);
	}
	return GlobalT::MaxPrecedence(relax, cse_relax);
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
	ftractions_elem.Dimension(fCSE->NumElementNodes(), NumSD());
	ftractions_elem = 0.0;
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
	dArrayT tangent(NumSD()), normal(NumSD());
	dArrayT ip_tracts;
	ostream& out = ElementSupport().Output();
	bool over = false;
	
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
				double* s = fNodalValues(pface[j]);

				ftractions_elem(j,0) = s[0]*normal[0] + s[2]*normal[1];
				ftractions_elem(j,1) = s[2]*normal[0] + s[1]*normal[1];
				double t_mag2 = ftractions_elem(j,0)*ftractions_elem(j,0) + ftractions_elem(j,1)*ftractions_elem(j,1);

				if(tmax < sqrt(t_mag2)) {
					tmax =  sqrt(t_mag2);
					jmax = j;
				}
			}
			
			t_mag2 = ftractions_elem(jmax,0)*ftractions_elem(jmax,0) + ftractions_elem(jmax,1)*ftractions_elem(jmax,1);
			double sense = ftractions_elem(jmax,0)*normal[0] + ftractions_elem(jmax,1)*normal[1];
		    
				
			int nip = fCSE->NumIP();
			/* tensile release */
			double ratio = sqrt(t_mag2/(fReleaseThreshold*fReleaseThreshold));
			
			if (ratio - 1.0 > 0.01){
				over = true;
				cout << "\nNodal stress exceeds release threshold by more than 1%";
			}

			if (fabs(ratio-1.0) < 0.01 && sense > kSmall) {
				fCSEActive[i] = ElementCardT::kMarkON;
				release_count++;

				/*assume that tractions are stored at the top of the state variable array and
				that every ip has the same number of state variables*/
				int num_state = state_variables.MinorDim(i)/nip;
				/* write traction into state variables */
				cout << "\nelement: "<<i;		
				if (state_variables.Length() > 0)
				{
					for (int j = 0; j < nip; j++)
					{
						fCSE->Interpolate(ip_tracts, ftractions_elem, j);
//						cout << "\niptracts: "<<ip_tracts;
						state_variables(i,j*num_state) = ip_tracts[0];
						state_variables(i,j*num_state+1) = ip_tracts[1];
					}
//					for (int j = 0; j < num_state*nip; j++)
//						cout << "\nstate: "<<state_variables(i,j);
				}
			}
		}
	  }
				out << "\nrelease_count: "<<release_count;
				out << "\ntmax: "<<sqrt(t_mag2);
				cout << "\nrelease_count: "<<release_count;
				cout << "\ntmax: "<<sqrt(t_mag2);

	/* relaxation code for CSE release */
/*	GlobalT::RelaxCodeT cse_relax_code = (release_count > 0) ? 
		GlobalT::kReEQ : GlobalT::kNoRelax;
	if (release_count > 0) SetNetwork (fCSEActive);
*/
		
	GlobalT::RelaxCodeT cse_relax_code;

	if (release_count > 0)
		if (relax_code == GlobalT::kRelax || relax_code == GlobalT::kReEQRelax)
			cse_relax_code = GlobalT::kFailReset;
		else {
			cse_relax_code = GlobalT::kReEQ;
			SetNetwork(fCSEActive);
		}
	else
		cse_relax_code = GlobalT::kNoRelax;
	
	if (over) cse_relax_code = GlobalT::kFailReset;

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
	
	/* reset tied nodes */
	fTied->SetTiedPairs(follower, leader);
	
	/* set the element flags in the CSE group */
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
