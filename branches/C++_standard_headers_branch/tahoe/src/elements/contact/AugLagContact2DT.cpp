/* $Id: AugLagContact2DT.cpp,v 1.21 2005-07-20 16:18:23 paklein Exp $ */
/* created: paklein (05/31/1998) */
#include "AugLagContact2DT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ifstreamT.h"
#include "eIntegratorT.h"
#include "ElementSupportT.h"
#include "XDOF_ManagerT.h"

using namespace Tahoe;

/* parameters */
const int kNumAugLagDOF  = 1;

/* constructor */
AugLagContact2DT::AugLagContact2DT(const ElementSupportT& support):
	Contact2DT(support),
	fr(0.0)
{
	SetName("contact_2D_multiplier");
}

/* append element equations numbers to the list */
void AugLagContact2DT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_2)

	/* collect using method allowing mixed node/tag numbers */
	ElementSupport().XDOF_Manager().XDOF_SetLocalEqnos(Group(), fXDOFConnectivities, fXDOFEqnos);

	/* add to list */
	eq_1.Append(&fXDOFEqnos);
}

/* returns the array for the DOF tags needed for the current config */
void AugLagContact2DT::SetDOFTags(void)
{
	/* DOF space about to be reset. store history */
	dArrayT constraints;
	constraints.Alias(ElementSupport().XDOF_Manager().XDOF(this, 0));
	fLastDOF = constraints;

	/* resize DOF tags array */
	fContactDOFtags.Dimension(fActiveStrikers.Length());
}

iArrayT& AugLagContact2DT::DOFTags(int tag_set)
{
#pragma unused(tag_set)
#if __option(extended_errorcheck)
	/* check */
	if (tag_set != 0)
		ExceptionT::OutOfRange("ugLagContact2DT::DOFTags",
			"group only has 1 tag set: %d", tag_set);
#endif
	
	return fContactDOFtags;
}

/* generate element data (based on current striker/body data) */
void AugLagContact2DT::GenerateElementData(void)
{
	/* inherited - set nodal connectivities */
	Contact2DT::SetConnectivities();

	/* dimension */
	int num_active = fConnectivities[0]->MajorDim();

	/* resize work space */
	fXDOFConnectivities_man.SetMajorDimension(num_active, false);
	fXDOFEqnos_man.SetMajorDimension(num_active, false);
	const int *pelem = fConnectivities[0]->Pointer();
	int rowlength = fConnectivities[0]->MinorDim();
	for (int i = 0; i < num_active; i++, pelem += rowlength)
	{	
		int* pxelem = fXDOFConnectivities(i);

		/* XDOF element tags */
		pxelem[0] = pelem[0]; // 1st facet node
		pxelem[1] = pelem[1]; // 2nd facet node
		pxelem[2] = pelem[2]; // striker node
		pxelem[3] = fContactDOFtags[i]; // contact DOF tag
	}
}

/* return the contact elements */
const iArray2DT& AugLagContact2DT::DOFConnects(int tag_set) const
{
#pragma unused(tag_set)
#if __option(extended_errorcheck)
	/* check */
	if (tag_set != 0)
		ExceptionT::OutOfRange("ugLagContact2DT::DOFConnects",
			"group only has 1 tag set: %d", tag_set);
#endif

	return fXDOFConnectivities;
}

/* restore the DOF values to the last converged solution */
void AugLagContact2DT::ResetDOF(dArray2DT& DOF, int tag_set) const
{
#pragma unused(tag_set)
#if __option(extended_errorcheck)
	/* check */
	if (tag_set != 0)
		ExceptionT::OutOfRange("ugLagContact2DT::ResetDOF",
			"group only has 1 tag set: %d", tag_set);
#endif

	/* alias */
	dArrayT constraints;
	constraints.Alias(DOF);
	constraints = 0.0;
	for (int i = 0; i < fLastActiveMap.Length(); i++)
	{
		int old_map = fLastActiveMap[i];
		int new_map = fActiveMap[i];
		if (old_map > -1 && new_map > -1)
			constraints[new_map] = fLastDOF[old_map];
	}
}

/* returns 1 if group needs to reconfigure DOF's, else 0 */
int AugLagContact2DT::Reconfigure(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = Contact2DT::RelaxSystem();
	if (relax != GlobalT::kNoRelax)
		return 1;
	else
		return 0;
}

/* element level reconfiguration for the current solution */
GlobalT::RelaxCodeT AugLagContact2DT::RelaxSystem(void)
{
	/* override all inherited - relaxation handled through
	 * DOFElementT interface */
	return GlobalT::kNoRelax;
}

/* appends group connectivities to the array */
void AugLagContact2DT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
	/* inherited */
	Contact2DT::ConnectsU(connects_1, connects_2);

	/* replace contact connects */
	bool found = false;
	for (int i = connects_1.Length() - 1; i > -1 && !found; i--)
		if (connects_1[i] == fConnectivities[0])
		{
			connects_1[i] = &fXDOFConnectivities;
			found = true;
		}

	/* check */
	if (!found) connects_1.AppendUnique(&fXDOFConnectivities);
}

/* restart functions */
void AugLagContact2DT::ReadRestart(istream& in)
{
#pragma unused(in)
	ExceptionT::GeneralFail("AugLagContact2DT::ReadRestart", "not tested");
}

void AugLagContact2DT::WriteRestart(ostream& out) const
{
#pragma unused(out)
	ExceptionT::GeneralFail("AugLagContact2DT::WriteRestart", "not tested");
}
//TEMP - restarts have not been tested. these functions
//       throw ExceptionT::xceptions

/* describe the parameters needed by the interface */
void AugLagContact2DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	Contact2DT::DefineParameters(list);

	/* regularization */
	ParameterT regularization(ParameterT::Double, "regularization");
	regularization.SetDefault(1.0);
	regularization.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(regularization);
}

/* accept parameter list */
void AugLagContact2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	Contact2DT::TakeParameterList(list);

	/* regularization */
	fr = list.GetParameter("regularization");

	/* reset base class parameters */
	int neq = NumElementNodes()*NumDOF() + 1; // 1 additional dof

	/* re-size element results */
	fLHS.Dimension(neq); // or make new variables?
	fRHS.Dimension(neq);

	/* dynamic work space managers for element arrays */
	fXDOFConnectivities_man.SetWard(0, fXDOFConnectivities, NumElementNodes() + 1);		
	fXDOFEqnos_man.SetWard(0, fXDOFEqnos, neq);

	/* only 1 tag set for the group */
	iArrayT numDOF(1);
	numDOF = kNumAugLagDOF;

	/* register with node manager - sets initial fContactDOFtags */
	ElementSupport().XDOF_Manager().XDOF_Register(this, numDOF);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* step in setting contact configuration. */
bool AugLagContact2DT::SetActiveInteractions(void)
{
	/* striker map about to be reset. store history */
	fLastActiveMap = fActiveMap;

	/* inherited */
	return Contact2DT::SetActiveInteractions();
}

/* called by FormRHS and FormLHS */
void AugLagContact2DT::LHSDriver(GlobalT::SystemTypeT)
{
	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	/* get reference to global coordinates and constrain force vector */
	const dArray2DT& coords = ElementSupport().CurrentCoordinates();
	const dArray2DT& constr = ElementSupport().XDOF_Manager().XDOF(this, 0);
	const dArrayT force(constr.MajorDim(),constr.Pointer());

	/* loop over active elements */
	int neq = NumElementNodes()*NumDOF() + 1;
	dArrayT tangent(NumSD());
	iArrayT eqnos;
	dMatrixT uLHS(neq - kNumAugLagDOF);
	dArrayT  uRHS(neq - kNumAugLagDOF, fRHS.Pointer());
	for (int i = 0; i < fXDOFConnectivities.MajorDim(); i++)
	{
		int* pelem = fXDOFConnectivities(i);
	
		/* get facet and striker coords */
		coords.RowAlias(pelem[0], fx1);
		coords.RowAlias(pelem[1], fx2);
		coords.RowAlias(pelem[2], fStriker);

		/* penetration vectors */
		fv1.DiffOf(fStriker, fx1);
		fv2.DiffOf(fStriker, fx2);

		/* tangent vector */
		tangent.DiffOf(fx2, fx1);

		/* distance to facet (could store some of this) */
		double magtan = tangent.Magnitude();				
		double      h = (fv2[0]*fv1[1] - fv1[0]*fv2[1])/magtan;

		/* augmented Lagragian multiplier */
		double g = force[i] + fr*h;

		/* tributary area */
		int striker_index = fStrikerTags_map.Map(pelem[2]);
		double area = fStrikerArea[striker_index];

		/* contact */
		if (g < 0.0)
		{
			/* initialize */
			fRHS = 0.0; // to hold d h/d d_i
			fLHS = 0.0;
			uLHS = 0.0;
					
			/* d_tan_j = tan_i d tan_i/d d_j */
			fdtanT.Multx(tangent,fNEEvec);
						
			/* compute  d h/d d_i*/
			uRHS.AddScaled(-h*area/(magtan*magtan),fNEEvec);

			fColtemp1.Set(neq - 1, fdv1T(0));
			fColtemp2.Set(neq - 1, fdv2T(1));
			uRHS.AddCombination(-fv2[1]*area/magtan,fColtemp1,
				                -fv1[0]*area/magtan,fColtemp2);
			
			fColtemp1.Set(neq - 1, fdv1T(1));
			fColtemp2.Set(neq - 1, fdv2T(0));
			uRHS.AddCombination(fv2[0]*area/magtan, fColtemp1,
				                fv1[1]*area/magtan, fColtemp2);

			/* dd_ e_3jk v2_j v1_k */
			double Kh_by_t = g*area/magtan;
			uLHS(0,3) = uLHS(3,0) =-Kh_by_t;
			uLHS(0,5) = uLHS(5,0) = Kh_by_t;
			uLHS(1,2) = uLHS(2,1) = Kh_by_t;
			uLHS(1,4) = uLHS(4,1) =-Kh_by_t;
			uLHS(2,5) = uLHS(5,2) =-Kh_by_t;
			uLHS(3,4) = uLHS(4,3) = Kh_by_t;

			/* (d_tan_ (x) d h/d d_)^s */
			fNEEmat.Outer(uRHS, fNEEvec);
			fNEEmat.Symmetrize();
			uLHS.AddScaled(-2.0*g*area/(magtan*magtan), fNEEmat);

			/* d_tan_ (x) d_tan_ */
			fNEEmat.Outer(fNEEvec, fNEEvec);
			uLHS.AddScaled(g*h*area/pow(magtan,4), fNEEmat);

			/* tan_k/d d_i tan_k/d d_j */
			fNEEmat.MultABT(fdtanT, fdtanT);
			uLHS.AddScaled(-g*h*area/(magtan*magtan), fNEEmat);

			/* d h/d d_ (x) d h/d d_ */
			fNEEmat.Outer(uRHS, uRHS);
			uLHS.AddScaled(fr*area, fNEEmat);

			/* assemble sub-block */
			fLHS.AddBlock(0, 0, uLHS);
			
			/* augmented Lagrangian DOF */
			int dex = neq - 1;
			fLHS.SetRow(dex, fRHS);
			fLHS.SetCol(dex, fRHS);
		}
		else /* gap */
		{
			/* initialize */
			fLHS = 0.0;
		
			/* augmented Lagrangian DOF */
			int dex = neq - 1;
			fLHS(dex,dex) = -area/fr;							
		}

		/* get equation numbers */
		fXDOFEqnos.RowAlias(i, eqnos);
			
		/* assemble */
		ElementSupport().AssembleLHS(Group(), fLHS, eqnos);
	}
}

void AugLagContact2DT::RHSDriver(void)
{
	/* time-stepping parameters */
	double constKd = 0.0;
	int     formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* get reference to global coordinates and constrain force vector */
	const dArray2DT& coords = ElementSupport().CurrentCoordinates(); //EFFECTIVE_DVA
	const dArray2DT& constr = ElementSupport().XDOF_Manager().XDOF(this, 0);
	const dArrayT force(constr.MajorDim(), constr.Pointer()); // general for all
	                                                          // value of kNumAugLagDOF
	/* reset tracking data */
	int num_contact = 0;
	double h_max = 0.0;

	/* loop over active elements */
	int neq = NumElementNodes()*NumDOF() + 1;
	dArrayT tangent(NumSD());
	iArrayT eqnos;
	dArrayT uRHS(neq - kNumAugLagDOF, fRHS.Pointer());
	for (int i = 0; i < fXDOFConnectivities.MajorDim(); i++)
	{
		int* pelem = fXDOFConnectivities(i);
	
		/* get facet and striker coords */
		coords.RowAlias(pelem[0], fx1);
		coords.RowAlias(pelem[1], fx2);
		coords.RowAlias(pelem[2], fStriker);

		/* penetration vectors */
		fv1.DiffOf(fStriker, fx1);
		fv2.DiffOf(fStriker, fx2);

		/* tangent vector */
		tangent.DiffOf(fx2, fx1);

		/* distance to facet (could store some of this) */
		double magtan = tangent.Magnitude();
		double      h = (fv2[0]*fv1[1] - fv1[0]*fv2[1])/magtan;
//double  max_d =-magtan/10; //max penetration
	  //double      b = Dot(tangent,fv1)/magtan;
	  //don't worry about "on facet" or max penetration for now

		/* augmented Lagrangian multiplier */
		double g = force[i] + fr*h;

		/* tributary area */
		int striker_index = fStrikerTags_map.Map(pelem[2]);
		double area = fStrikerArea[striker_index];

		/* store force vector output */
		double f_by_t = -force[i]*area/magtan;
		fStrikerForce2D(striker_index,0) = f_by_t*tangent[1];
		fStrikerForce2D(striker_index,1) =-f_by_t*tangent[0];

		/* contact */
		if (g < 0.0)
		{
			/* tracking data */
			num_contact++;
			h_max = (h < h_max) ? h : h_max;

			/* initialize */
			fRHS = 0.0;
			
			/* grad_disp contribution */
					
			/* d_tan contribution */
			fdtanT.Multx(tangent,fNEEvec);
			uRHS.AddScaled(area*g*h/(magtan*magtan),fNEEvec);
						
			/* d_area */
			fColtemp1.Set(neq - 1, fdv1T(0));
			fColtemp2.Set(neq - 1, fdv2T(1));
			uRHS.AddCombination(area*g*fv2[1]/magtan, fColtemp1,
				                area*g*fv1[0]/magtan, fColtemp2);
			
			fColtemp1.Set(neq - 1, fdv1T(1));
			fColtemp2.Set(neq - 1, fdv2T(0));
			uRHS.AddCombination(-area*g*fv2[0]/magtan, fColtemp1,
				                -area*g*fv1[1]/magtan, fColtemp2);
				
			/* augmented Lagrangian DOF */				
			fRHS[neq - 1] = -h*area;
		}
		/* gap */
		else
		{
			/* grad_disp contribution */
			fRHS = 0.0;

			/* augmented Lagrangian DOF */				
			fRHS[neq - 1] = area*force[i]/fr;							
		}

		/* get equation numbers */
		fXDOFEqnos.RowAlias(i, eqnos);

		/* assemble */
		ElementSupport().AssembleRHS(Group(), fRHS, eqnos);
	}

	/* set tracking */
	SetTrackingData(num_contact, h_max);
}
