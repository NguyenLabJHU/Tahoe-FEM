/* $Id: BimaterialK_FieldT.cpp,v 1.1.1.1 2001-01-29 08:20:40 paklein Exp $ */
/* created: paklein (09/05/2000)                                          */

#include "BimaterialK_FieldT.h"

#include "NodeManagerT.h"
#include "fstreamT.h"
#include "IsotropicT.h"
#include "Material2DT.h"

/* parameters */
const double Pi = acos(-1.0);

/* constructor */
BimaterialK_FieldT::BimaterialK_FieldT(NodeManagerT& node_manager):
	K_FieldT(node_manager),
	fIsotropic_2(NULL),
	fMaterial2D_2(NULL)
{

}

/* initialize data - called immediately after construction */
void BimaterialK_FieldT::Initialize(ifstreamT& in)
{
	/* only 2D for now */
	int nsd = fNodeManager.NumSD();
	if (nsd != 2)
	{
		cout << "\n BimaterialK_FieldT::Initialize: must be 2D: " << nsd << endl;
		throw eGeneralFail;
	}

	/* K1 */
	in >> fnumLTf1 >> fK1; fnumLTf1--;
	fLTf1 = fNodeManager.GetLTfPtr(fnumLTf1);	
	if (!fLTf1) throw eBadInputValue;

	/* K2 */
	in >> fnumLTf2 >> fK2; fnumLTf2--;
	fLTf2 = fNodeManager.GetLTfPtr(fnumLTf2);	
	if (!fLTf2) throw eBadInputValue;

	/* coordinates of the crack tip */
	fInitTipCoords.Allocate(nsd);
	in >> fInitTipCoords;
	fLastTipCoords = fTipCoords = fInitTipCoords;

	/* crack extension parameters */
	fGrowthDirection.Allocate(nsd);
	in >> fGrowthDirection; fGrowthDirection.UnitVector();

	/* near tip group */
	in >> fNearTipGroupNum;   // -1: no nearfield group
	in >> fNearTipOutputCode; // variable to locate crack tip
	in >> fTipColumnNum;      // column of output variable to locate tip
	in >> fMaxGrowthDistance; if (fMaxGrowthDistance < 0.0) throw eBadInputValue;
	in >> fMaxGrowthSteps; if (fMaxGrowthSteps < 1) throw eBadInputValue;

	/* offsets and checks */
	if (fNearTipGroupNum != -1) fNearTipGroupNum--;
	fTipColumnNum--;
	if (fNearTipGroupNum <  -1) throw eBadInputValue;

	/* nodes */
	in >> fFarFieldGroupNum;
	if (fFarFieldGroupNum != -1)
	{
		in >> fFarFieldMaterialNum;
		ReadNodes(in, fID_List_1, fNodes_1);

		/* checks and offsets */
		fFarFieldGroupNum--;
		fFarFieldMaterialNum--;
		if (fFarFieldGroupNum < 0) throw eBadInputValue;
		if (fFarFieldMaterialNum < 0) throw eBadInputValue;
	}
	else
		fFarFieldMaterialNum = -1;

	in >> fFarFieldGroupNum_2;
	if (fFarFieldGroupNum_2 != -1)
	{
		in >> fFarFieldMaterialNum_2;
		ReadNodes(in, fID_List_2, fNodes_2);

		/* offsets and checks */
		fFarFieldGroupNum_2--;
		fFarFieldMaterialNum_2--;
		if (fFarFieldGroupNum_2 < 0) throw eBadInputValue;
		if (fFarFieldMaterialNum_2 < 0) throw eBadInputValue;
	}
	else
		fFarFieldMaterialNum_2 = -1;
		
	/* checks */
	if (fFarFieldGroupNum == -1 && fFarFieldGroupNum_2 == -1)
	{
		cout << "\n BimaterialK_FieldT::Initialize: at least one half plane must be specified" << endl;
		throw eBadInputValue;
	}

	/* create overall ID lists */
	if (fID_List_1.Length() == 0)
		fID_List.Alias(fID_List_2);
	else if (fID_List_2.Length() == 0)
		fID_List.Alias(fID_List_1);
	else
	{
		fID_List.Allocate(fID_List_1.Length() + fID_List_2.Length());
		fID_List.CopyIn(0, fID_List_1);
		fID_List.CopyIn(fID_List_1.Length(), fID_List_2);
		fID_List_1.Set(fID_List_1.Length(), fID_List.Pointer());
		fID_List_2.Set(fID_List_2.Length(), fID_List.Pointer(fID_List_1.Length()));
	}

	/* create overall node lists */
	if (fNodes_1.Length() == 0)
		fNodes.Alias(fNodes_2);
	else if (fNodes_2.Length() == 0)
		fNodes.Alias(fNodes_1);
	else
	{
		fNodes.Allocate(fNodes_1.Length() + fNodes_2.Length());
		fNodes.CopyIn(0, fNodes_1);
		fNodes.CopyIn(fNodes_1.Length(), fNodes_2);
		fNodes_1.Set(fNodes_1.Length(), fNodes.Pointer());
		fNodes_2.Set(fNodes_2.Length(), fNodes.Pointer(fNodes_1.Length()));
	}

	/* generate BC cards */
	fKBC_Cards.Allocate(fNodes.Length()*nsd);
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
	fK1Disp.Allocate(fNodes.Length(), nsd);
	fK2Disp.Allocate(fNodes.Length(), nsd);
	if (fNodes_1.Length() == 0)
	{
		fK1Disp_2.Alias(fK1Disp);
		fK2Disp_2.Alias(fK2Disp);
	}
	else if (fNodes_2.Length() == 0)
	{
		fK1Disp_1.Alias(fK1Disp);
		fK2Disp_1.Alias(fK2Disp);
	}
	else
	{
		fK1Disp_1.Set(fNodes_1.Length(), nsd, fK1Disp(0));
		fK1Disp_2.Set(fNodes_2.Length(), nsd, fK1Disp(fNodes_1.Length()));

		fK2Disp_1.Set(fNodes_1.Length(), nsd, fK2Disp(0));
		fK2Disp_2.Set(fNodes_2.Length(), nsd, fK2Disp(fNodes_1.Length()));
	}
	
	/* resolve groups in UHP/LHP */
	fUHP = UpperHalfPlane();
}

void BimaterialK_FieldT::WriteParameters(ostream& out) const
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
	out << " Fracture path output code . . . . . . . . . . . = " << fNearTipOutputCode << '\n';
	out << " Fracture path test value column number. . . . . = " << fTipColumnNum + 1  << '\n';
	out << " Maximum extension distance per load step. . . . = " << fMaxGrowthDistance << '\n';
	out << " Maximum number of extensions per load step. . . = " << fMaxGrowthSteps    << '\n';

	iArrayT tmp;
	/* set 1 */
	if (fFarFieldGroupNum != -1)
	{
		out << " Far field element group 1 number. . . . . . . . = " << fFarFieldGroupNum + 1 << '\n';
		out << " Far field group 1 material number . . . . . . . = " << fFarFieldMaterialNum + 1 << '\n';
		if (fID_List_1.Length() > 0)
		{
			out << " Number of group 1 node sets . . . . . . . . . . = " << fID_List.Length() << '\n';
			out << fID_List_1.wrap(6) << '\n';
		}
		out << " Number of group 1 nodes . . . . . . . . . . . . = " << fNodes.Length() << '\n';	
		tmp.Alias(fNodes_1);
		tmp++;
		out << tmp.wrap(6) << '\n';
		tmp--;
	}
	else
		out << " Far field element group 1 number. . . . . . . . = RIGID\n";
	
	/* set 2 */
	if (fFarFieldGroupNum_2 != -1)
	{
		out << " Far field element group 2 number. . . . . . . . = " << fFarFieldGroupNum_2 + 1 << '\n';
		out << " Far field group 2 material number . . . . . . . = " << fFarFieldMaterialNum_2 + 1 << '\n';
		if (fID_List_2.Length() > 0)
		{
			out << " Number of group 2 node sets . . . . . . . . . . = " << fID_List_2.Length() << '\n';
			out << fID_List_2.wrap(6) << '\n';
		}
		out << " Number of group 2 nodes . . . . . . . . . . . . = " << fNodes_2.Length() << '\n';	
		tmp.Alias(fNodes_2);
		tmp++;
		out << tmp.wrap(6) << '\n';
		tmp--;
	}
	else
		out << " Far field element group 2 number. . . . . . . . = RIGID\n";

	out << " Group in the upper half plane . . . . . . . . . = " << fUHP << '\n';
}

/***********************************************************************
* Protected
***********************************************************************/

/* compute K-field displacement factors */
void BimaterialK_FieldT::ComputeDisplacementFactors(const dArrayT& tip_coords)
{
	/* resolve near tip and material reference */
	if (fFarFieldGroupNum != -1 && !fIsotropic)
		ResolveMaterialReference(fFarFieldGroupNum, fFarFieldMaterialNum,
			&fIsotropic, &fMaterial2D);
	else
	{
		fIsotropic = NULL;
		fMaterial2D = NULL;
	}

	if (fFarFieldGroupNum_2 != -1 && !fIsotropic_2)		
		ResolveMaterialReference(fFarFieldGroupNum_2, fFarFieldMaterialNum_2,
			&fIsotropic_2, &fMaterial2D_2);
	else
	{
		fIsotropic_2 = NULL;
		fMaterial2D_2 = NULL;
	}

	/* moduli */
	double G_1, nu_1, mu_1;
	if (fIsotropic)
	{
		G_1 = fIsotropic->Mu();
		nu_1 = fIsotropic->Poisson();	
		mu_1 = 3.0 - 4.0*nu_1;	
	}
	double G_2, nu_2, mu_2;
	if (fIsotropic_2)
	{
		G_2 = fIsotropic_2->Mu();
		nu_2 = fIsotropic_2->Poisson();	
		mu_2 = 3.0 - 4.0*nu_2;
	}
	
	if (fNodeManager.NumSD() == 2)
	{
		if (fMaterial2D && fMaterial2D->ConstraintOption() == Material2DT::kPlaneStress)
			mu_1 = (3.0 - nu_1)/(1.0 + nu_1);
		if (fMaterial2D_2 && fMaterial2D_2->ConstraintOption() == Material2DT::kPlaneStress)
			mu_2 = (3.0 - nu_2)/(1.0 + nu_2);
	}

	if (fUHP == 1)
	{
		double eps;
		if (!fIsotropic)
			eps =-log(mu_2)/(2.0*Pi);
		else if (!fIsotropic_2)
			eps = log(mu_1)/(2.0*Pi);		
		else
			eps = log((mu_1/G_1 + 1.0/G_2)/
	                  (mu_2/G_2 + 1.0/G_1))/(2.0*Pi);
	
		/* UHP */
		if (fIsotropic)
			SetFieldFactors( 1, eps, mu_1, G_1, tip_coords, fNodes_1, fK1Disp_1, fK2Disp_1);

		/* LHP */
		if (fIsotropic_2)
			SetFieldFactors(-1, eps, mu_2, G_2, tip_coords, fNodes_2, fK1Disp_2, fK2Disp_2);
	}
	else
	{
		double eps;
		if (!fIsotropic)
			eps = log(mu_2)/(2.0*Pi);
		else if (!fIsotropic_2)
			eps =-log(mu_1)/(2.0*Pi);
		else
			eps = log((mu_2/G_2 + 1.0/G_1)/
	                  (mu_1/G_1 + 1.0/G_2))/(2.0*Pi);

		/* UHP */
		if (fIsotropic_2)
			SetFieldFactors( 1, eps, mu_2, G_2, tip_coords, fNodes_2, fK1Disp_2, fK2Disp_2);

		/* LHP */
		if (fIsotropic)
			SetFieldFactors(-1, eps, mu_1, G_1, tip_coords, fNodes_1, fK1Disp_1, fK2Disp_1);
	}
}

/***********************************************************************
* Private
***********************************************************************/

/* bimaterial displacement field factors */
void BimaterialK_FieldT::SetFieldFactors(int side, double eps, double mu,
	double G, const dArrayT& tip_coords, const iArrayT& nodes,
	dArray2DT& K1_disp, dArray2DT& K2_disp)
{
	if (side != 1 && side != -1) throw eGeneralFail;

	/* (initial) nodal coordinates */
	int nsd = fNodeManager.NumSD();
	const dArray2DT& init_coords = fNodeManager.InitialCoordinates();

	/* coefficient */
	double a = exp(Pi*eps)/(1.0 + exp(2.0*Pi*eps))/(2.0*G)/sqrt(2.0*Pi);

	/* compute K-field displacement factors */
	dArrayT coords;
	dArrayT	rvec(nsd);
	dArrayT ey(nsd);
	ey[0] =-fGrowthDirection[1];
	ey[1] = fGrowthDirection[0];	
	for (int i = 0; i < nodes.Length(); i++)
	{
		/* fetch coords */
		init_coords.RowAlias(nodes[i], coords);
		
		/* vector from the tip */	
		rvec.DiffOf(coords, tip_coords);
		
		/* (local) polar coords */
		double rx = dArrayT::Dot(rvec, fGrowthDirection);
		double ry = dArrayT::Dot(rvec, ey);

		double r = rvec.Magnitude();
		double t = atan2(ry, rx);
		
		double delta = (side == 1) ?
			exp((t - Pi)*eps) :
			exp((Pi + t)*eps);
		
		/* factors */
		double  gamma = mu*delta - 1.0/delta;
		double gamma1 = mu*delta + 1.0/delta;
		
		double eps_logr = eps*log(r); // log_e or log_10??
		double   beta = (0.5*cos(eps_logr) + eps*sin(eps_logr))/(0.25 + eps*eps);
		double  beta1 = (0.5*sin(eps_logr) - eps*cos(eps_logr))/(0.25 + eps*eps);

		double cos_tby2 = cos(0.5*t);
		double sin_tby2 = sin(0.5*t);
		double D = beta*gamma*cos_tby2 + beta1*gamma1*sin_tby2;
		double C = beta1*gamma*cos_tby2 - beta*gamma1*sin_tby2;

		double psi = eps_logr + 0.5*t;
		double sin_psi = sin(psi);
		double cos_psi = cos(psi);
		double a_sqrtr = a*sqrt(r);
		double d2sin_t = 2.0*delta*sin(t);

		/* K I factors */
		double f1I = D + d2sin_t*sin_psi; // (A1)
		double f2I =-C - d2sin_t*cos_psi; // (A2)

		K1_disp(i,0) = a_sqrtr*f1I; // (9.0)
		K1_disp(i,1) = a_sqrtr*f2I; // (9.1)

		/* K II factors */
		double f1II =-C + d2sin_t*cos_psi; // (A3)
		double f2II =-D + d2sin_t*sin_psi; // (A4)
		
		K2_disp(i,0) = a_sqrtr*f1II; // (12.0)
		K2_disp(i,1) = a_sqrtr*f2II; // (12.1)
	}
}

/* group in the "upper half plane" */
int BimaterialK_FieldT::UpperHalfPlane(void) const
{
	/* no nodes */
	if (fNodes_1.Length() == 0 && fNodes_2.Length() == 0) return 1;

	/* undeformed coordinates */
	const dArray2DT& init_coords = fNodeManager.InitialCoordinates();
	int nsd = init_coords.MinorDim();

	/* compute polar angles */
	double t_1 = 0.0;
	if (fNodes_1.Length() > 0)
	{
		dArrayT x_1(nsd), ey(nsd), tmp;
		x_1 = 0.0;
		ey[0] =-fGrowthDirection[1];
		ey[1] = fGrowthDirection[0];
		for (int i = 0; i < fNodes_1.Length(); i++)
		{
			init_coords.RowAlias(fNodes_1[i], tmp);
			x_1 += tmp;
		}
		x_1 /= fNodes_1.Length();

		/* vector from the tip */	
		x_1 -= fInitTipCoords;
		
		/* (local) polar coords */
		double rx = dArrayT::Dot(x_1, fGrowthDirection);
		double ry = dArrayT::Dot(x_1, ey);
		t_1 = atan2(ry, rx);
	}

	double t_2 = 0.0;
	if (fNodes_2.Length() > 0)
	{
		dArrayT x_2(nsd), ey(nsd), tmp;
		x_2 = 0.0;
		ey[0] =-fGrowthDirection[1];
		ey[1] = fGrowthDirection[0];
		for (int i = 0; i < fNodes_2.Length(); i++)
		{
			init_coords.RowAlias(fNodes_2[i], tmp);
			x_2 += tmp;
		}
		x_2 /= fNodes_2.Length();

		/* vector from the tip */	
		x_2 -= fInitTipCoords;
		
		/* (local) polar coords */
		double rx = dArrayT::Dot(x_2, fGrowthDirection);
		double ry = dArrayT::Dot(x_2, ey);
		t_2 = atan2(ry, rx);
	}
	
	/* resolve group in UHP */
	if (fNodes_1.Length() > 0 && fNodes_2.Length() > 0)
	{
		/* check */
		if (t_1*t_2 >= 0.0)
		{
			cout << "\n BimaterialK_FieldT::UpperHalfPlane: could not determine group in the\n"
			     <<   "     upper half plane" << endl;		
			throw eGeneralFail;
		}
	
		if (t_1 >= 0.0)
			return 1;
		else
			return 2;
	}
	else if (fNodes_1.Length() > 0)
	{
		if (t_1 >= 0.0)
			return 1;
		else
			return 2;
	}
	else
	{
		if (t_2 >= 0.0)
			return 2;
		else
			return 1;
	}
}
