/* $Id: BimaterialK_FieldT.cpp,v 1.10.2.3 2004-07-13 16:42:45 paklein Exp $ */
/* created: paklein (09/05/2000) */
#include "BimaterialK_FieldT.h"

#include "NodeManagerT.h"
#include "ParameterUtils.h"
#include "ParameterContainerT.h"

/* build options */
#include "ElementsConfig.h"
#ifdef CONTINUUM_ELEMENT
#include "IsotropicT.h"
#endif

/* parameters */
using namespace Tahoe;

/* parameters */
const double Pi = acos(-1.0);

/* constructor */
BimaterialK_FieldT::BimaterialK_FieldT(const BasicSupportT& support):
	K_FieldT(support),
	fmu_1(-1.0), fnu_1(-1.0), fkappa_1(-1.0),	
	fGroupNumber_1(-1), fMaterialNumber_1(-1),
	fmu_2(-1.0), fnu_2(-1.0), fkappa_2(-1.0),	
	fGroupNumber_2(-1), fMaterialNumber_2(-1)
{
	SetName("bi-material_K-field");
}

/* information about subordinate parameter lists */
void BimaterialK_FieldT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	K_FieldT::DefineSubs(sub_list);

	/* remove previous definition of elastic properties */
	sub_list.RemoveSub("elastic_properties_choice");
	sub_list.RemoveSub("node_ID_list");

	/* boundary regions */
	sub_list.AddSub("boundary_region_1", ParameterListT::ZeroOrOnce);
	sub_list.AddSub("boundary_region_2", ParameterListT::ZeroOrOnce);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* BimaterialK_FieldT::NewSub(const StringT& name) const
{
	if (name == "boundary_region_1" || name == "boundary_region_2")
	{
		ParameterContainerT* region = new ParameterContainerT(name);
		region->SetSubSource(this);
		region->AddSub("elastic_properties_choice", ParameterListT::Once, true);
		region->AddSub("node_ID_list");
		return region;
	}
	else /* inherited */
		return K_FieldT::NewSub(name);
}

/* accept parameter list */
void BimaterialK_FieldT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "BimaterialK_FieldT::TakeParameterList";

	/* inherited */
	K_FieldT::TakeParameterList(list);

	/* only 2D for now */
	int nsd = fSupport.NumSD();
	if (nsd != 2) ExceptionT::GeneralFail(caller, "must be 2D: %d", nsd);

	/* resolve elastic properties */
	const ParameterListT* region_1 = list.List("boundary_region_1");
	if (region_1) {
	
		/* moduli */
		ResolveElasticProperties(*region_1, fGroupNumber_1, fMaterialNumber_1,  fmu_1, fnu_1, fkappa_1);
	
		/* nodes */
		StringListT::Extract(region_1->GetList("node_ID_list"),  fID_List_1);	
		GetNodes(fID_List_1, fNodes_1);
	}

	const ParameterListT* region_2 = list.List("boundary_region_2");
	if (region_2) {
	
		/* moduli */
		ResolveElasticProperties(*region_2, fGroupNumber_2, fMaterialNumber_2,  fmu_2, fnu_2, fkappa_2);
	
		/* nodes */
		StringListT::Extract(region_2->GetList("node_ID_list"),  fID_List_2);	
		GetNodes(fID_List_2, fNodes_2);
	}
	
	/* check */
	if (!region_1 && !region_2)
		ExceptionT::GeneralFail(caller, "expecting \"boundary_region_1\" and/or \"boundary_region_2\"");	

	/* create overall ID lists */
	if (fID_List_1.Length() == 0)
		fID_List.Alias(fID_List_2);
	else if (fID_List_2.Length() == 0)
		fID_List.Alias(fID_List_1);
	else
	{
		fID_List.Dimension(fID_List_1.Length() + fID_List_2.Length());
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
		fNodes.Dimension(fNodes_1.Length() + fNodes_2.Length());
		fNodes.CopyIn(0, fNodes_1);
		fNodes.CopyIn(fNodes_1.Length(), fNodes_2);
		fNodes_1.Set(fNodes_1.Length(), fNodes.Pointer());
		fNodes_2.Set(fNodes_2.Length(), fNodes.Pointer(fNodes_1.Length()));
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
	
//TEMP - tip tracking not supporting for parallel execution
	if (fNearTipGroupNum != -1 && fSupport.Size() > 1) 
		ExceptionT::BadInputValue(caller, "tip tracking not implemented in parallel");
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* compute K-field displacement factors */
void BimaterialK_FieldT::ComputeDisplacementFactors(const dArrayT& tip_coords)
{
	const char caller[] = "BimaterialK_FieldT::ComputeDisplacementFactors";

	/* resolve elastic constants */
	if (fmu_1 < 0.0 && fGroupNumber_1 > -1) 
	{
		/* resolve material and isotropy information */
		const IsotropicT* iso = NULL;
		const SolidMaterialT* mat = NULL;
		ResolveMaterialReference(fGroupNumber_1, fMaterialNumber_1, &iso, &mat);
			
		/* compute elastic constants */
		fmu_1 = iso->Mu();
		fnu_1 = iso->Poisson();	
		fkappa_1 = 3.0 - 4.0*fnu_1;
		if (fSupport.NumSD() == 2 && mat->Constraint() == SolidMaterialT::kPlaneStress)
			fkappa = (3.0 - fnu_1)/(1.0 + fnu_1);
	}

	/* resolve elastic constants */
	if (fmu_2 < 0.0 && fGroupNumber_2 > -1) 
	{
		/* resolve material and isotropy information */
		const IsotropicT* iso = NULL;
		const SolidMaterialT* mat = NULL;
		ResolveMaterialReference(fGroupNumber_2, fMaterialNumber_2, &iso, &mat);
			
		/* compute elastic constants */
		fmu_2 = iso->Mu();
		fnu_2 = iso->Poisson();	
		fkappa_2 = 3.0 - 4.0*fnu_2;
		if (fSupport.NumSD() == 2 && mat->Constraint() == SolidMaterialT::kPlaneStress)
			fkappa = (3.0 - fnu_2)/(1.0 + fnu_2);
	}

	/* rename variables */
	double  G_1 = fmu_1;
	double nu_1 = fnu_1;
	double mu_1 = fkappa_1;
	double  G_2 = fmu_2;
	double nu_2 = fnu_2;
	double mu_2 = fkappa_2;

	if (fUHP == 1)
	{
		double eps;
		if (G_2 > 0.0 && G_1 > 0.0)
			eps = log((mu_1/G_1 + 1.0/G_2)/
	                  (mu_2/G_2 + 1.0/G_1))/(2.0*Pi);
		else if (G_2 > 0.0)
			eps =-log(mu_2)/(2.0*Pi);
		else if (G_1 > 0.0)
			eps = log(mu_1)/(2.0*Pi);		
	
		/* UHP */
		if (G_1 > 0.0)
			SetFieldFactors( 1, eps, mu_1, G_1, tip_coords, fNodes_1, fK1Disp_1, fK2Disp_1);

		/* LHP */
		if (G_2 > 0.0)
			SetFieldFactors(-1, eps, mu_2, G_2, tip_coords, fNodes_2, fK1Disp_2, fK2Disp_2);
	}
	else
	{
		double eps;
		if (G_1 > 0.0 && G_2 > 0.0)
			eps = log((mu_2/G_2 + 1.0/G_1)/
	                  (mu_1/G_1 + 1.0/G_2))/(2.0*Pi);
		else if (G_2 > 0.0)
			eps = log(mu_2)/(2.0*Pi);
		else if (G_1 > 0.0)
			eps =-log(mu_1)/(2.0*Pi);

		/* UHP */
		if (G_2 > 0.0)
			SetFieldFactors( 1, eps, mu_2, G_2, tip_coords, fNodes_2, fK1Disp_2, fK2Disp_2);

		/* LHP */
		if (G_1 > 0.0)
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
	if (side != 1 && side != -1) ExceptionT::GeneralFail("BimaterialK_FieldT::SetFieldFactors");

	/* (initial) nodal coordinates */
	int nsd = fSupport.NumSD();
	const dArray2DT& init_coords = fSupport.InitialCoordinates();

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
	const dArray2DT& init_coords = fSupport.InitialCoordinates();
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
			ExceptionT::GeneralFail("BimaterialK_FieldT::UpperHalfPlane",
				" could not determine group in the upper half plane");
	
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
