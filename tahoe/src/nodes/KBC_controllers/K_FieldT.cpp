/* $Id: K_FieldT.cpp,v 1.20 2004-07-15 08:31:21 paklein Exp $ */
/* created: paklein (09/05/2000) */
#include "K_FieldT.h"

#include "NodeManagerT.h"
#include "FEManagerT.h"


#include "IsotropicT.h"
#include "SolidMaterialT.h"
#include "ElementsConfig.h"
#ifdef CONTINUUM_ELEMENT
#include "MaterialListT.h"
#include "ContinuumMaterialT.h"
#include "ContinuumElementT.h"
#include "IsotropicT.h"
#else
#include "ElementBaseT.h"
#endif
#include "ParameterContainerT.h"
#include "ParameterUtils.h"

using namespace Tahoe;

/* parameters */
const double Pi = acos(-1.0);
//const double TipNoise = 1.0e-06;

/* constructor */
K_FieldT::K_FieldT(const BasicSupportT& support):
	KBC_ControllerT(support),
	fLTf1(NULL),
	fLTf2(NULL),
	fK1(0.0),
	fK2(0.0),
	fDummySchedule(1.0),
	fNearTipGroupNum(-1),
	fNearTipOutputCode(-1),
	fTipColumnNum(-1),
	fTrackingCode(kMaximum),
	fMaxGrowthDistance(-1),
	fMaxGrowthSteps(-1),
	fmu(-1.0), fnu(-1.0), fkappa(-1.0),	
	fGroupNumber(-1),
	fMaterialNumber(-1)
{
	SetName("K-field");
}

void K_FieldT::InitialCondition(void)
{
	/* set initial crack tip position */
	fTipCoords = fInitTipCoords;

	/* set displacement factors */
	ComputeDisplacementFactors(fTipCoords);
}

/* restart operations */
void K_FieldT::ReadRestart(istream& in)
{
	/* inherited */
	KBC_ControllerT::ReadRestart(in);

	/* read tip coordinates */
	in >> fTipCoords;
	fLastTipCoords = fTipCoords;
	
	/* reset field factors */
	ComputeDisplacementFactors(fTipCoords);
}

void K_FieldT::WriteRestart(ostream& out) const
{
	/* inherited */
	KBC_ControllerT::WriteRestart(out);

	/* write tip coordinates */
	out << fTipCoords << '\n';
}

/* initialize/finalize/reset step */
void K_FieldT::InitStep(void)
{
	/* inherited */
	KBC_ControllerT::InitStep();

	/* reset extension count */
	fGrowthCount = 0;

	/* update BC cards */
	SetBCCards();
}

void K_FieldT::CloseStep(void)
{
	/* inherited */
	KBC_ControllerT::CloseStep();

	/* has moving tip */
	if (fNearTipGroupNum > -1)
		fLastTipCoords = fTipCoords;
}

void K_FieldT::Reset(void)
{
	/* inherited */
	KBC_ControllerT::Reset();

	/* has moving tip */
	if (fNearTipGroupNum > -1)
	{
		/* move the tip back */
		fTipCoords = fLastTipCoords;
	
		/* reset field factors */
		ComputeDisplacementFactors(fTipCoords);
	}
}

/* returns true if the internal force has been changed since
* the last time step */
GlobalT::RelaxCodeT K_FieldT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = KBC_ControllerT::RelaxSystem();

	/* no fracture path group */
	if (fNearTipGroupNum == -1) return relax;
		
	/* new tip coordinates */
	dArrayT tip_coords(fTipCoords.Length());
	GetNewTipCoordinates(tip_coords);
	
	/* extension vector/distance */
	dArrayT extension(fTipCoords.Length());
	extension.DiffOf(tip_coords, fTipCoords);
	double advance = dArrayT::Dot(fGrowthDirection, extension);

	/* sizable move */
	if (fabs(advance) > 10.0*kSmall)
	{
		/* reverse tip direction */
		if (advance < 0.0)
		{
			/* tip moving backwards */
			cout << " K_FieldT::RelaxSystem: tip moving backwards: IGNORED: " << advance << '\n';
			cout << " current position: " << fTipCoords[0] << '\n';
			return relax;
		}
		else
		{
			/* too much unzipping */
			if (++fGrowthCount == fMaxGrowthSteps || advance > fMaxGrowthDistance)
			{
				cout << "\n K_FieldT::RelaxSystem: exceeded max growth per increment\n"
					 <<   "   count: " << fGrowthCount << '\n'
					 <<   "    dist: " << advance << endl;
				throw ExceptionT::kBadJacobianDet; // to trigger step cut
			}

			/* move the crack tip coords */
			fTipCoords.AddScaled(advance, fGrowthDirection);
			cout << "\n Crack extension increment: " << setw(kDoubleWidth) << advance << '\n';
			cout <<   " New crack tip coordinates: \n" << fTipCoords << '\n';
			
			/* compute field factors */
			ComputeDisplacementFactors(fTipCoords);
			
			/* update BC cards */
			SetBCCards();
		
			return GlobalT::MaxPrecedence(relax, GlobalT::kRelax);
		}
	}
	else
		return relax;
}

/* output current configuration */
void K_FieldT::WriteOutput(ostream& out) const
{
	/* inherited */
	KBC_ControllerT::WriteOutput(out);
	
	/* K-field information */
	out << "\n K - f i e l d   D a t a :\n\n";
	out << " Crack tip coordinates: \n" << fTipCoords << '\n';
	out << " K I . . . . . . . . . . . . . . . . . . . . . . = " << (fLTf1 ? fK1*fLTf1->Value() : 0.0) << '\n';
	out << " K II. . . . . . . . . . . . . . . . . . . . . . = " << (fLTf2 ? fK2*fLTf2->Value() : 0.0) << '\n';
}

/* information about subordinate parameter lists */
void K_FieldT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	KBC_ControllerT::DefineSubs(sub_list);

	/* applied stress intensity factors */
	sub_list.AddSub("K_I", ParameterListT::ZeroOrOnce);
	sub_list.AddSub("K_II", ParameterListT::ZeroOrOnce);

	/* initial tip coordinates */
	sub_list.AddSub("initial_tip_coordinates");
	
	/* crack extension direction */
	sub_list.AddSub("crack_extension_direction");

	/* far field elastic properties */
	sub_list.AddSub("elastic_properties_choice", ParameterListT::Once, true);

	/* list of affected nodes */
	sub_list.AddSub("node_ID_list");

	/* tip tracking parameters */
	sub_list.AddSub("tip_tracking", ParameterListT::ZeroOrOnce);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* K_FieldT::NewSub(const StringT& name) const
{
	if (name == "initial_tip_coordinates" || name == "crack_extension_direction")
		return new VectorParameterT(name, 'x', 2);
	else if (name == "K_I" || name == "K_II")
	{
		ParameterContainerT* K_spec = new ParameterContainerT(name);
		K_spec->AddParameter(ParameterT::Double, "K");
		K_spec->AddParameter(ParameterT::Integer, "schedule");
		return K_spec;
	}
	else if (name == "tip_tracking")
	{
		ParameterContainerT* tracking = new ParameterContainerT(name);
		tracking->SetSubSource(this);
	
		/* define tracking data */
		tracking->AddParameter(fNearTipGroupNum, "near_tip_group");
		tracking->AddParameter(fNearTipOutputCode, "near_tip_output_code");
		tracking->AddParameter(fTipColumnNum, "tip_output_column");

		/* growth limits */
		tracking->AddParameter(fMaxGrowthDistance, "max_growth_distance");
		tracking->AddParameter(fMaxGrowthSteps, "max_growth_steps");
	
		/* choice of tracking methods */
		tracking->AddSub("tip_tracking_method", ParameterListT::Once, true);

		return tracking;
	}
	else if (name == "tip_tracking_method")
	{
		ParameterContainerT* method = new ParameterContainerT(name);
		method->SetListOrder(ParameterListT::Choice);
	
		/* maximum value */
		ParameterContainerT max("location_of_maximum");
		max.AddParameter(ParameterT::Double, "noise_level");
		method->AddSub(max);
		
		/* theshold value */
		ParameterContainerT theshold("farthest_above_threshold");
		theshold.AddParameter(ParameterT::Double, "theshold");
		method->AddSub(theshold);	
	
		return method;
	}
	else if (name == "elastic_properties_choice")
	{
		ParameterContainerT* props = new ParameterContainerT(name);
		props->SetListOrder(ParameterListT::Choice);
		props->SetSubSource(this);

#ifdef CONTINUUM_ELEMENT
		/* define from material in far-field element group */
		ParameterContainerT group("far_field_element_group");
		group.AddParameter(ParameterT::Integer, "group_number");
		group.AddParameter(ParameterT::Integer, "material_number");
		props->AddSub(group);
#endif

		/* define moduli directly */
		ParameterContainerT moduli("far_field_elastic_properties");
		moduli.SetSubSource(this);
		moduli.AddSub("isotropic");
		ParameterT constraint(ParameterT::Enumeration, "constraint_2D");
		constraint.AddEnumeration("plane_stress", SolidMaterialT::kPlaneStress);
		constraint.AddEnumeration("plane_strain", SolidMaterialT::kPlaneStrain);
		constraint.SetDefault(SolidMaterialT::kPlaneStrain);
		moduli.AddParameter(constraint);
		props->AddSub(moduli);
	
		return props;
	}
	else if (name == "isotropic")
		return new IsotropicT;
	else /* inherited */
		return KBC_ControllerT::NewSub(name);
}

/* accept parameter list */
void K_FieldT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "K_FieldT::TakeParameterList";

	/* inherited */
	KBC_ControllerT::TakeParameterList(list);

	/* only 2D for now */
	int nsd = fSupport.NumSD();
	if (nsd != 2) ExceptionT::GeneralFail(caller, "must be 2D: %d", nsd);

	/* K_I */
	const ParameterListT* K1 = list.List("K_I");
	if (K1) {
		fK1 = K1->GetParameter("K");
		int num_LTf1 = K1->GetParameter("schedule");
		num_LTf1--;
		fLTf1 = fSupport.Schedule(num_LTf1);
		if (!fLTf1) ExceptionT::BadInputValue(caller, "could not resolve schedule %d", num_LTf1+1);	
	}

	/* K_II */
	const ParameterListT* K2 = list.List("K_II");
	if (K2) {
		fK2 = K2->GetParameter("K");
		int num_LTf2 = K2->GetParameter("schedule");
		num_LTf2--;
		fLTf2 = fSupport.Schedule(num_LTf2);
		if (!fLTf2) ExceptionT::BadInputValue(caller, "could not resolve schedule %d", num_LTf2+1);	
	}

	/* check */
	if (!fLTf1 && !fLTf2)
		ExceptionT::GeneralFail(caller, "neither \"K_I\" or \"K_II\" K-field defined");

	/* initial tip coordinates */
	VectorParameterT vec('x', 2);
	vec.SetName("initial_tip_coordinates");
	vec.TakeParameterList(list.GetList(vec.Name()));
	fInitTipCoords = vec;
	fLastTipCoords = fTipCoords = fInitTipCoords;

	/* extension direction */
	vec.SetName("crack_extension_direction");
	vec.TakeParameterList(list.GetList(vec.Name()));
	fGrowthDirection = vec; 
	fGrowthDirection.UnitVector();

	/* resolve elastic properties */
	ResolveElasticProperties(list, fGroupNumber, fMaterialNumber, fmu, fnu, fkappa);

	/* tip tracking */
	const ParameterListT* tracking = list.List("tip_tracking");
	if (tracking) {
	
		/* define tracking data */
		fNearTipGroupNum = tracking->GetParameter("near_tip_group"); fNearTipGroupNum--;
		fNearTipOutputCode = tracking->GetParameter("near_tip_output_code"); fNearTipOutputCode--;
		fTipColumnNum = tracking->GetParameter("tip_output_column"); fTipColumnNum--;

		/* growth limits */
		fMaxGrowthDistance = tracking->GetParameter("max_growth_distance");
		fMaxGrowthSteps = tracking->GetParameter("max_growth_steps");

		/* resolve tracking method */
		ParameterInterfaceT* tracking_info = NewSub(tracking->Name());
		if (!tracking_info) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", tracking->Name().Pointer());
		const ParameterListT& tracking_method = tracking->GetListChoice(*tracking_info, "tip_tracking_method");
		delete tracking_info;
		if (tracking_method.Name() == "location_of_maximum") {
			fTrackingCode = kMaximum;
			fTrackingParameters.Dimension(1);
			fTrackingParameters[0] = tracking_method.GetParameter("noise_level");				
		}
		else if (tracking_method.Name() == "farthest_above_threshold") {
			fTrackingCode = kThreshold;
			fTrackingParameters.Dimension(1);
			fTrackingParameters[0] = tracking_method.GetParameter("threshold");
		}
		else
			ExceptionT::GeneralFail(caller, "unrecognized tracking method \"%s\"",
				tracking_method.Name().Pointer());
	}

	/* nodes */
	const ParameterListT* nodes = list.List("node_ID_list");
	if (nodes) {
		StringListT::Extract(*nodes,  fID_List);	
		GetNodes(fID_List, fNodes);
	}

	/* generate BC cards */
	fKBC_Cards.Dimension(fNodes.Length()*nsd);
	KBC_CardT* pcard = fKBC_Cards.Pointer();
	for (int i = 0; i < fNodes.Length(); i++)
		for (int j = 0; j < nsd; j++)
		{
			/* set values */
			pcard->SetValues(fNodes[i], j, KBC_CardT::kDsp, &fDummySchedule, 0.0);
			pcard++;
		}	

	/* allocate displacement field factors */
	fK1Disp.Dimension(fNodes.Length(), nsd);
	fK2Disp.Dimension(fNodes.Length(), nsd);
	
//TEMP - tip tracking not supporting for parallel execution
	if (fNearTipGroupNum != -1 && fSupport.Size() > 1) 
		ExceptionT::BadInputValue(caller, "tip tracking not implemented in parallel");

#if 0
	if (fID_List.Length() > 0)
	{
		out << " Number of group node sets . . . . . . . . . . . = " << fID_List.Length() << '\n';
		int wrap = 0;
		for (int i = 0; i < fID_List.Length(); i++)
		{
			if (++wrap == 4) {
				out << '\n';
				wrap = 0;
			}
			out << setw(12) << fID_List[i];
		}	
		out << '\n';
	}
	out << " Number of group nodes . . . . . . . . . . . . . = " << fNodes.Length() << '\n';	
	iArrayT tmp;
	tmp.Alias(fNodes);
	tmp++;
	out << tmp.wrap(6) << '\n';
	tmp--;
#endif
}

/* extract elastic constants */
void K_FieldT::ResolveElasticProperties(const ParameterListT& list,
	int& group_number, int& material_number, double& mu, double& nu, double& kappa) const
{
	const char caller[] = "K_FieldT::ResolveElasticProperties";
	
	const ParameterListT* elastic = list.ListChoice(*this, "elastic_properties_choice");
	if (elastic)
	{
		if (elastic->Name() == "far_field_element_group") 
		{
#ifndef CONTINUUM_ELEMENT
			ExceptionT::BadInputValue(caller, "\"%s\" requires CONTINUUM_ELEMENT", elastic->Name().Pointer());
#endif
			/* extract element group information - the group won't be available until later */
			group_number = elastic->GetParameter("group_number"); group_number--;
			material_number = elastic->GetParameter("material_number"); material_number--;
		}
		else if (elastic->Name() == "far_field_elastic_properties")
		{
			IsotropicT iso;
			iso.TakeParameterList(elastic->GetList("isotropic"));
			mu = iso.Mu();
			nu = iso.Poisson();	
			kappa = 3.0 - 4.0*nu;
			int constraint = elastic->GetParameter("constraint_2D");
			if (constraint == SolidMaterialT::kPlaneStress)
				kappa = (3.0 - nu)/(1.0 + nu);
		}
		else
			ExceptionT::GeneralFail(caller, "unrecognized properties choice \"%s\"",
				elastic->Name().Pointer());
	}
	else {
		group_number = -1;
		material_number = -1;
		mu = -1.0;
		nu = -1.0;
		kappa = -1.0;
	}
}

/**********************************************************************
 * Protected
 **********************************************************************/

/* determine the new tip coordinates */
void K_FieldT::GetNewTipCoordinates(dArrayT& tip_coords)
{
	const char caller[] = "K_FieldT::GetNewTipCoordinates";

	/* near tip element group */
	const FEManagerT& fe_man = fSupport.FEManager();
	ElementBaseT* neartip_group = fe_man.ElementGroup(fNearTipGroupNum);
	if (!neartip_group)
		ExceptionT::GeneralFail(caller, "could not resolve near tip element group number %d", fNearTipGroupNum+1);

	/* signal to accumulate nodal values */
	neartip_group->SendOutput(fNearTipOutputCode);

	/* the nodes */
	NodeManagerT& node_manager = fSupport.NodeManager();

	/* find new tip coordinates */
	tip_coords = fTipCoords;
	switch (fTrackingCode)
	{
		case kMaximum:
		{
			/* find the node with maximum value */
			int maxrow;
			double maxval;
			node_manager.MaxInColumn(fTipColumnNum, maxrow, maxval);
			if (maxrow == -1) ExceptionT::GeneralFail(caller);

			/* get new tip coordinates */
			double tip_noise = fTrackingParameters[0];
			if (maxval > tip_noise)
				fSupport.InitialCoordinates().RowCopy(maxrow, tip_coords);

			break;
		}
		case kThreshold:
		{
			/* nodal coordinates */
			const dArray2DT& initial_coordinates = fSupport.InitialCoordinates();
		
			/* get all nodal values */
			const dArray2DT& nodal_values = node_manager.OutputAverage(); 
		
			/* test threshold */
			double threshold = fTrackingParameters[0];
			double max_growth = 0.0;
			dArrayT advance(fGrowthDirection.Length());
			dArrayT x_node;
			for (int i = 0; i < nodal_values.MajorDim(); i++)
				if (nodal_values(i,fTipColumnNum) > threshold) {
			
					/* crack extension increment */
					initial_coordinates.RowAlias(i, x_node);
					advance.DiffOf(x_node, fTipCoords);
					double growth = dArrayT::Dot(fGrowthDirection, advance);
			
					/* maximum extension */
					if (growth > max_growth) {
						max_growth = growth;
						tip_coords = x_node;
					}
				}

			break;
		}
		default:
			ExceptionT::GeneralFail(caller, "unrecognized tracking method %d", fTrackingCode);
	}
}

/* resolve element info to isotropic material */
void K_FieldT::ResolveMaterialReference(int element_group,
	int material_num, const IsotropicT** iso, const SolidMaterialT** mat) const
{
	const char caller[] = "K_FieldT::ResolveMaterialReference";

#ifdef CONTINUUM_ELEMENT
	/* resolve element group */
	const FEManagerT& fe_man = fSupport.FEManager();
	const ElementBaseT* element = fe_man.ElementGroup(element_group);
	if (!element) ExceptionT::GeneralFail(caller, "could not resolve element group");

	const ContinuumElementT* cont_element = TB_DYNAMIC_CAST(const ContinuumElementT*, element);
	if (!cont_element)
		ExceptionT::GeneralFail(caller, "could not cast element group %d to ContinuumElementT", element_group+1);

	/* resolve material reference */
	const MaterialListT& material_list = cont_element->MaterialsList();
	ContinuumMaterialT* cont_mat = material_list[material_num];
	if (!cont_mat) ExceptionT::GeneralFail(caller, "could not resolve continuum material");
	*iso = TB_DYNAMIC_CAST(IsotropicT*, cont_mat);
	if (!(*iso))
		ExceptionT::GeneralFail(caller, "could not cast material %d \"%s\" to IsotropicT", 
			material_num+1, cont_mat->Name().Pointer());

#ifdef __NO_RTTI__
	ExceptionT::GeneralFail("K_FieldT::ResolveMaterialReference", "requires RTTI");
#endif

	if (fSupport.NumSD() == 2)
	{
		*mat = TB_DYNAMIC_CAST(SolidMaterialT*, cont_mat);
		if (!(*mat))
			ExceptionT::GeneralFail(caller, "could not cast material %d \"%s\" to Material2DT",
				material_num+1, cont_mat->Name().Pointer());		
	}
#else /* CONTINUUM_ELEMENT */
#pragma unused(element_group)
#pragma unused(material_num)
#pragma unused(iso)
#pragma unused(mat)
ExceptionT::GeneralFail(caller);
#endif /* CONTINUUM_ELEMENT */
}

/* compute K-field displacement factors */
void K_FieldT::ComputeDisplacementFactors(const dArrayT& tip_coords)
{
	/* (initial) nodal coordinates */
	int nsd = fSupport.NumSD();
	const dArray2DT& init_coords = fSupport.InitialCoordinates();

	/* resolve elastic constants */
	if (fmu < 0.0 && fGroupNumber > -1) 
	{
		/* resolve material and isotropy information */
		const IsotropicT* iso = NULL;
		const SolidMaterialT* mat = NULL;
		ResolveMaterialReference(fGroupNumber, fMaterialNumber, &iso, &mat);
			
		/* compute elastic constants */
		fmu = iso->Mu();
		fnu = iso->Poisson();	
		fkappa = 3.0 - 4.0*fnu;
		if (fSupport.NumSD() == 2 && mat->Constraint() == SolidMaterialT::kPlaneStress)
			fkappa = (3.0 - fnu)/(1.0 + fnu);
	}
	else
		ExceptionT::GeneralFail("K_FieldT::ComputeDisplacementFactors", "elastic constants not resolved");	

	/* compute K-field displacement factors (Andersen Table 2.2): */
	dArrayT coords;
	dArrayT	rvec(nsd);
	dArrayT ey(nsd);
	ey[0] =-fGrowthDirection[1];
	ey[1] = fGrowthDirection[0];
	for (int i = 0; i < fNodes.Length(); i++)
	{
		/* fetch coords */
		init_coords.RowAlias(fNodes[i], coords);
		
		/* vector from the tip */	
		rvec.DiffOf(coords, tip_coords);
		
		/* polar coords (factors) */
		double rx = dArrayT::Dot(rvec, fGrowthDirection);
		double ry = dArrayT::Dot(rvec, ey);

		double r = sqrt(rvec.Magnitude()/(2.0*Pi));
		double t = atan2(ry, rx)/2.0;

		/* K I components */
		fK1Disp(i,0) = (0.5/fmu)*r*cos(t)*(fkappa - 1.0 + 2.0*pow(sin(t), 2.0));
		fK1Disp(i,1) = (0.5/fmu)*r*sin(t)*(fkappa + 1.0 - 2.0*pow(cos(t), 2.0));

		/* K II components */
		fK2Disp(i,0) = (0.5/fmu)*r*sin(t)*(fkappa + 1.0 + 2.0*pow(cos(t), 2.0));
		fK2Disp(i,1) =-(0.5/fmu)*r*cos(t)*(fkappa - 1.0 - 2.0*pow(sin(t), 2.0));
	}
}

/* set BC cards with current displacement field */
void K_FieldT::SetBCCards(void)
{
	/* field intensities */
	double K1 = fLTf1 ? fK1*fLTf1->Value() : 0.0;
	double K2 = fLTf2 ? fK2*fLTf2->Value() : 0.0;

	/* apply K-field displacement */
	dArrayT disp;
	dArrayT K1disp;
	dArrayT K2disp;
	int dex = 0;
	for (int i = 0; i < fNodes.Length(); i++)
	{
		/* K-field node */
		int node = fNodes[i];

		/* shallow copies */
		fK1Disp.RowAlias(i, K1disp);
		fK2Disp.RowAlias(i, K2disp);

		/* displacement */
		double d1 = K1*K1disp[0] + K2*K2disp[0];
		double d2 = K1*K1disp[1] + K2*K2disp[1];
	
		/* set cards */
		fKBC_Cards[dex++].SetValues(node, 0, KBC_CardT::kDsp, NULL, d1);
		fKBC_Cards[dex++].SetValues(node, 1, KBC_CardT::kDsp, NULL, d2);
	}
}
