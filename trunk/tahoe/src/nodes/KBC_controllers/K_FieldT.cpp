/* $Id: K_FieldT.cpp,v 1.19 2004-06-17 07:41:57 paklein Exp $ */
/* created: paklein (09/05/2000) */
#include "K_FieldT.h"

#include "NodeManagerT.h"
#include "FEManagerT.h"
#include "ifstreamT.h"

#include "ElementsConfig.h"
#ifdef CONTINUUM_ELEMENT
#include "MaterialListT.h"
#include "ContinuumMaterialT.h"
#include "ContinuumElementT.h"
#include "IsotropicT.h"
#include "Material2DT.h"
#else
#include "ElementBaseT.h"
#endif

using namespace Tahoe;

/* parameters */
const double Pi = acos(-1.0);
const double TipNoise = 1.0e-06;

/* constructor */
K_FieldT::K_FieldT(NodeManagerT& node_manager):
	KBC_ControllerT(node_manager),
	fLTf1(NULL),
	fLTf2(NULL),
	fIsotropic(NULL),
	fMaterial2D(NULL),
	fDummySchedule(1.0),
	fNearTipGroupNum(-1),
	fNearTipOutputCode(-1),
	fTipColumnNum(-1),
	fTrackingCode(kMaximum),
	fMaxGrowthDistance(-1),
	fMaxGrowthSteps(-1)
{
	SetName("K_field");

#ifndef CONTINUUM_ELEMENT
	ExceptionT::BadInputValue("K_FieldT::K_FieldT", "CONTINUUM_ELEMENT not enabled");
#endif
}

/* initialize data - called immediately after construction */
void K_FieldT::Initialize(ifstreamT& in)
{
	const char caller[] = "K_FieldT::Initialize";

	/* only 2D for now */
	int nsd = fNodeManager.NumSD();
	if (nsd != 2) ExceptionT::GeneralFail(caller, "must be 2D: %d", nsd);

	/* K1 */
	in >> fnumLTf1 >> fK1; fnumLTf1--;
	fLTf1 = fNodeManager.Schedule(fnumLTf1);	
	if (!fLTf1) ExceptionT::BadInputValue(caller);

	/* K2 */
	in >> fnumLTf2 >> fK2; fnumLTf2--;
	fLTf2 = fNodeManager.Schedule(fnumLTf2);	
	if (!fLTf2) ExceptionT::BadInputValue(caller);

	/* coordinates of the crack tip */
	fInitTipCoords.Dimension(nsd);
	in >> fInitTipCoords;
	fLastTipCoords = fTipCoords = fInitTipCoords;

	/* crack extension parameters */
	fGrowthDirection.Dimension(nsd);
	in >> fGrowthDirection; fGrowthDirection.UnitVector();

	/* near tip group */
	in >> fNearTipGroupNum;   // -1: no nearfield group
	in >> fNearTipOutputCode; // variable to locate crack tip
	in >> fTipColumnNum;      // column of output variable to locate tip

	/* tracking type */
	int tracking_type = -99;
	in >> tracking_type;
	switch (tracking_type)
	{
		case kMaximum:
		{
			fTrackingCode = kMaximum;
			break;
		}
		case kThreshold:
		{
			fTrackingCode = kThreshold;
		
			/* threshold value */
			fTrackingParameters.Dimension(1);
			in >> fTrackingParameters;
			break;
		}
		default:
			ExceptionT::BadInputValue(caller, "unrecognized tip tracking method: %d", tracking_type);
	}

	/* growth limiting */
	in >> fMaxGrowthDistance; if (fMaxGrowthDistance < 0.0) ExceptionT::BadInputValue(caller);
	in >> fMaxGrowthSteps; if (fMaxGrowthSteps < 1) ExceptionT::BadInputValue(caller);

	/* offsets and checks */
	fNearTipOutputCode--;
	if (fNearTipGroupNum != -1) fNearTipGroupNum--;
	fTipColumnNum--;
	if (fNearTipGroupNum <  -1) ExceptionT::BadInputValue(caller);

	/* nodes */
	in >> fFarFieldGroupNum;
	in >> fFarFieldMaterialNum;
	ReadNodes(in, fID_List, fNodes);

	/* offsets and checks */
	fFarFieldGroupNum--;
	fFarFieldMaterialNum--;
	if (fFarFieldGroupNum < 0) ExceptionT::BadInputValue(caller);
	if (fFarFieldMaterialNum < 0) ExceptionT::BadInputValue(caller);

	/* generate BC cards */
	fKBC_Cards.Dimension(fNodes.Length()*nsd);
	KBC_CardT* pcard = fKBC_Cards.Pointer();
	for (int i = 0; i < fNodes.Length(); i++)
		for (int j = 0; j < nsd; j++)
		{
			/* set values */
			pcard->SetValues(fNodes[i], j, KBC_CardT::kDsp, 0, 0.0);
	
			/* dummy schedule */
			pcard->SetSchedule(&fDummySchedule);
			pcard++;
		}	

	/* allocate displacement field factors */
	fK1Disp.Dimension(fNodes.Length(), nsd);
	fK2Disp.Dimension(fNodes.Length(), nsd);
	
//TEMP - tip tracking not supporting for parallel execution
	if (fNearTipGroupNum != -1 && fNodeManager.Size() > 1) 
		ExceptionT::BadInputValue(caller, "tip tracking not implemented in parallel");
}

void K_FieldT::WriteParameters(ostream& out) const
{
	/* inherited */
	KBC_ControllerT::WriteParameters(out);

	out << "\n K - f i e l d   p a r a m e t e r s :\n\n";
	out << " K I LTf . . . . . . . . . . . . . . . . . . . . = " << fnumLTf1 + 1 << '\n';
	out << " K I . . . . . . . . . . . . . . . . . . . . . . = " << fK1         << '\n';
	out << " K II LTf. . . . . . . . . . . . . . . . . . . . = " << fnumLTf2 + 1 << '\n';
	out << " K II. . . . . . . . . . . . . . . . . . . . . . = " << fK2         << '\n';
	out << " Crack tip coordinates (initial):\n" << fTipCoords << '\n';
	out << " Crack extension direction:\n" << fGrowthDirection << '\n';
	out << " Fracture path element group . . . . . . . . . . = ";
	if (fNearTipGroupNum > -1) out << fNearTipGroupNum + 1 << '\n';
	else out << "<none>\n";
	out << " Fracture path output code . . . . . . . . . . . = " << fNearTipOutputCode + 1 << '\n';
	out << " Fracture path test value column number. . . . . = " << fTipColumnNum + 1  << '\n';
	out << " Tip tracking method . . . . . . . . . . . . . . = " << fTrackingCode << '\n';
	if (fTrackingParameters.Length() > 0) out << " parameters: " << fTrackingParameters.no_wrap() << '\n';
	out << " Maximum extension distance per load step. . . . = " << fMaxGrowthDistance << '\n';
	out << " Maximum number of extensions per load step. . . = " << fMaxGrowthSteps    << '\n';
	out << " Far field element group number. . . . . . . . . = " << fFarFieldGroupNum + 1 << '\n';
	out << " Far field group material number . . . . . . . . = " << fFarFieldMaterialNum + 1 << '\n';
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
	out << " K I . . . . . . . . . . . . . . . . . . . . . . = " << fK1*fLTf1->Value() << '\n';
	out << " K II. . . . . . . . . . . . . . . . . . . . . . = " << fK2*fLTf2->Value() << '\n';
}

/* describe the parameters needed by the interface */
void K_FieldT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	KBC_ControllerT::DefineParameters(list);

	list.AddParameter(fNearTipGroupNum, "near_tip_group");
	list.AddParameter(fNearTipOutputCode, "near_tip_output_code");
	list.AddParameter(fTipColumnNum, "tip_output_column");
	list.AddParameter(fMaxGrowthDistance, "max_growth_distance");
	list.AddParameter(fMaxGrowthSteps, "max_growth_steps");
}

/**********************************************************************
* Protected
**********************************************************************/

/* determine the new tip coordinates */
void K_FieldT::GetNewTipCoordinates(dArrayT& tip_coords)
{
	const char caller[] = "K_FieldT::GetNewTipCoordinates";

	/* near tip element group */
	const FEManagerT& fe_man = fNodeManager.FEManager();
	ElementBaseT* neartip_group = fe_man.ElementGroup(fNearTipGroupNum);
	if (!neartip_group)
		ExceptionT::GeneralFail(caller, "could not resolve near tip element group number %d", fNearTipGroupNum+1);

	/* signal to accumulate nodal values */
	neartip_group->SendOutput(fNearTipOutputCode);

	/* find new tip coordinates */
	tip_coords = fTipCoords;
	switch (fTrackingCode)
	{
		case kMaximum:
		{
			/* find the node with maximum value */
			int maxrow;
			double maxval;
			fNodeManager.MaxInColumn(fTipColumnNum, maxrow, maxval);
			if (maxrow == -1) ExceptionT::GeneralFail(caller);

			/* get new tip coordinates */
			if (maxval > TipNoise)
				fNodeManager.InitialCoordinates().RowCopy(maxrow, tip_coords);

			break;
		}
		case kThreshold:
		{
			/* nodal coordinates */
			const dArray2DT& initial_coordinates = fNodeManager.InitialCoordinates();
		
			/* get all nodal values */
			const dArray2DT& nodal_values = fNodeManager.OutputAverage(); 
		
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
	int material_num, const IsotropicT** iso, const Material2DT** mat) const
{
	const char caller[] = "K_FieldT::ResolveMaterialReference";

#ifdef CONTINUUM_ELEMENT
	/* resolve element group */
	const FEManagerT& fe_man = fNodeManager.FEManager();
	const ElementBaseT* element = fe_man.ElementGroup(element_group);
	if (!element) throw ExceptionT::kGeneralFail;

	const ContinuumElementT* cont_element = TB_DYNAMIC_CAST(const ContinuumElementT*, element);
	if (!cont_element)
		ExceptionT::GeneralFail(caller, "could not cast element group %d to ContinuumElementT", element_group+1);

	/* resolve material reference */
	const MaterialListT& material_list = cont_element->MaterialsList();
	ContinuumMaterialT* cont_mat = material_list[material_num];
	if (!cont_mat) throw ExceptionT::kGeneralFail;
	*iso = TB_DYNAMIC_CAST(IsotropicT*, cont_mat);
	if (!(*iso))
	{
		cout << "\n K_FieldT::ResolveReference: could not cast material "
		     << material_num+1<< " to\n" <<   "     IsotropicT:\n";
		cont_mat->PrintName(cout);
		cout.flush();
		throw ExceptionT::kGeneralFail;
	}

#ifdef __NO_RTTI__
	ExceptionT::GeneralFail("K_FieldT::ResolveMaterialReference", "requires RTTI");
#endif

	if (fNodeManager.NumSD() == 2)
	{
		*mat = TB_DYNAMIC_CAST(Material2DT*, cont_mat);
		if (!(*mat))
		{
			cout << "\n K_FieldT::ResolveReference: could not cast material "
			     << material_num+1<< " to\n" <<   "     Material2DT:\n";
			cont_mat->PrintName(cout);
			cout.flush();
			throw ExceptionT::kGeneralFail;
		}
	}
#else /* CONTINUUM_ELEMENT */
#pragma unused(element_group)
#pragma unused(material_num)
#pragma unused(iso)
#pragma unused(mat)
ExceptionT::GeneralFail();
#endif /* CONTINUUM_ELEMENT */
}

/* compute K-field displacement factors */
void K_FieldT::ComputeDisplacementFactors(const dArrayT& tip_coords)
{
	/* (initial) nodal coordinates */
	int nsd = fNodeManager.NumSD();
	const dArray2DT& init_coords = fNodeManager.InitialCoordinates();

	/* resolve near tip and material reference */
	if (!fIsotropic)
		ResolveMaterialReference(fFarFieldGroupNum, fFarFieldMaterialNum,
			&fIsotropic, &fMaterial2D);

#ifdef CONTINUUM_ELEMENT
	/* moduli */
	double mu = fIsotropic->Mu();
	double nu = fIsotropic->Poisson();	
	double kappa = 3.0 - 4.0*nu;
	if (fNodeManager.NumSD() == 2)
	{
		if (!fMaterial2D) throw ExceptionT::kGeneralFail;
		if (fMaterial2D->ConstraintOption() == Material2DT::kPlaneStress)
			kappa = (3.0 - nu)/(1.0 + nu);
	}
#else
	double mu = 0.0;
	double nu = 0.0;	
	double kappa = 0.0;
#endif

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
		fK1Disp(i,0) = (0.5/mu)*r*cos(t)*(kappa - 1.0 + 2.0*pow(sin(t), 2.0));
		fK1Disp(i,1) = (0.5/mu)*r*sin(t)*(kappa + 1.0 - 2.0*pow(cos(t), 2.0));

		/* K II components */
		fK2Disp(i,0) = (0.5/mu)*r*sin(t)*(kappa + 1.0 + 2.0*pow(cos(t), 2.0));
		fK2Disp(i,1) =-(0.5/mu)*r*cos(t)*(kappa - 1.0 - 2.0*pow(sin(t), 2.0));
	}
}

/* set BC cards with current displacement field */
void K_FieldT::SetBCCards(void)
{
	/* field intensities */
	double K1 = fK1*fLTf1->Value();
	double K2 = fK2*fLTf2->Value();

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
		fKBC_Cards[dex++].SetValues(node, 0, KBC_CardT::kDsp, 0, d1);
		fKBC_Cards[dex++].SetValues(node, 1, KBC_CardT::kDsp, 0, d2);
	}
}
