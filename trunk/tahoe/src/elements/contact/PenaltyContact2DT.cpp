/* $Id: PenaltyContact2DT.cpp,v 1.15 2005-03-12 08:38:09 paklein Exp $ */
/* created: paklein (12/11/1997) */
#include "PenaltyContact2DT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ifstreamT.h"
#include "eIntegratorT.h"

using namespace Tahoe;

/* constructor */
PenaltyContact2DT::PenaltyContact2DT(const ElementSupportT& support):
	Contact2DT(support),
	fK(0.0)
{
	SetName("contact_2D_penalty");
}

/* describe the parameters needed by the interface */
void PenaltyContact2DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	Contact2DT::DefineParameters(list);

	/* penalty stiffness */
	ParameterT stiffness(ParameterT::Double, "penalty_stiffness");
	stiffness.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(stiffness);
}

/* accept parameter list */
void PenaltyContact2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	Contact2DT::TakeParameterList(list);

	/* contact stiffness */
	fK = list.GetParameter("penalty_stiffness");

	/* dimension work space */
	fElCoord.Dimension(fNumFacetNodes + 1, NumSD());
	fElDisp.Dimension(fNumFacetNodes + 1, NumDOF());	
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* called by FormRHS and FormLHS */
void PenaltyContact2DT::LHSDriver(GlobalT::SystemTypeT)
{
	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	/* get reference to global coordinates */
	const dArray2DT& coords = ElementSupport().CurrentCoordinates(); //EFFECTIVE_DVA

	/* loop over active elements */
	dArrayT tangent(NumSD());
	iArrayT eqnos;
	const int* pelem = fConnectivities[0]->Pointer();
	int rowlength = fConnectivities[0]->MinorDim();
	for (int i = 0; i < fConnectivities[0]->MajorDim(); i++, pelem += rowlength)
	{
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
//		double  max_d =-magtan/10; //max penetration

		/* contact */
		if (h < 0.0)
		{
			/* initialize */
			fRHS = 0.0; // to hold d h/d d_i
			fLHS = 0.0;
					
			/* d_tan_j = tan_i d tan_i/d d_j */
			fdtanT.Multx(tangent, fNEEvec);
						
			/* compute  d h/d d_i*/
			fRHS.AddScaled(-h/(magtan*magtan), fNEEvec);

			fColtemp1.Set(fdv1T.Rows(), fdv1T(0));
			fColtemp2.Set(fdv2T.Rows(), fdv2T(1));
			fRHS.AddCombination(-fv2[1]/magtan,fColtemp1,
				                -fv1[0]/magtan,fColtemp2);
			
			fColtemp1.Set(fdv1T.Rows(), fdv1T(1));
			fColtemp2.Set(fdv2T.Rows(), fdv2T(0));
			fRHS.AddCombination(fv2[0]/magtan,fColtemp1,
				                fv1[1]/magtan,fColtemp2);

			/* dd_ e_3jk v2_j v1_k */
			double Kh_by_t = fK*h/magtan;
			fLHS(0,3) = fLHS(3,0) =-Kh_by_t;
			fLHS(0,5) = fLHS(5,0) = Kh_by_t;
			fLHS(1,2) = fLHS(2,1) = Kh_by_t;
			fLHS(1,4) = fLHS(4,1) =-Kh_by_t;
			fLHS(2,5) = fLHS(5,2) =-Kh_by_t;
			fLHS(3,4) = fLHS(4,3) = Kh_by_t;

			/* d h/d d_ (x) d h/d d_ */
			fNEEmat.Outer(fRHS, fRHS);
			fLHS.AddScaled(fK, fNEEmat);
			
			/* (d_tan_ (x) d h/d d_)^s */
			fNEEmat.Outer(fRHS, fNEEvec);
			fNEEmat.Symmetrize();
			fLHS.AddScaled(-2.0*fK*h/(magtan*magtan), fNEEmat);

			/* d_tan_ (x) d_tan_ */
			fNEEmat.Outer(fNEEvec, fNEEvec);
			fLHS.AddScaled(fK*h*h/pow(magtan,4), fNEEmat);

			/* tan_k/d d_i tan_k/d d_j */
			fNEEmat.MultABT(fdtanT, fdtanT);
			fLHS.AddScaled(-fK*h*h/(magtan*magtan), fNEEmat);

			/* get equation numbers */
			fEqnos[0].RowAlias(i, eqnos);
			
			/* time integration factor */
			fLHS *= constK;
			
			/* assemble */
			ElementSupport().AssembleLHS(Group(), fLHS, eqnos);
		}
	}
}

void PenaltyContact2DT::RHSDriver(void)
{
	/* time integration parameters */
	double constKd = 0.0;
	int     formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* references to global nodal data */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& disp = Field()[0]; /* displacements */

	/* reset tracking data */
	int num_contact = 0;
	double h_max = 0.0;

	/* clear force */
	fStrikerForce2D = 0.0;

	/* loop over active elements */
	dArrayT tangent(NumSD());
	iArrayT eqnos;
	const int* pelem = fConnectivities[0]->Pointer();
	int rowlength = fConnectivities[0]->MinorDim();
	for (int i = 0; i < fConnectivities[0]->MajorDim(); i++, pelem += rowlength)
	{
		/* collect element configuration */
		fElCoord.RowCollect(pelem, init_coords);
		fElDisp.RowCollect(pelem, disp);

		/* current configuration using effective displacement */
		fElCoord.AddScaled(constKd, fElDisp); //EFFECTIVE_DVA
	
		/* get facet and striker coords */
		fElCoord.RowAlias(0, fx1);
		fElCoord.RowAlias(1, fx2);
		fElCoord.RowAlias(2, fStriker);

		/* penetration vectors */
		fv1.DiffOf(fStriker, fx1);
		fv2.DiffOf(fStriker, fx2);

		/* tangent vector */
		tangent.DiffOf(fx2, fx1);

		/* distance to facet (could store some of this) */
		double magtan = tangent.Magnitude();				
		double      h = (fv2[0]*fv1[1] - fv1[0]*fv2[1])/magtan;
//		double  max_d =-magtan/10; //max penetration

		/* contact */
		if (h < 0.0)
		{
			/* tracking data */
			num_contact++;
			h_max = (h < h_max) ? h : h_max;

			/* penetration force */
			double dphi =-fK*h;
			
			/* initialize */
			fRHS = 0.0;
					
			/* d_tan contribution */
			fdtanT.Multx(tangent, fNEEvec);
			fRHS.AddScaled(-dphi*h/(magtan*magtan), fNEEvec);
						
			/* d_area */
			fColtemp1.Set(fdv1T.Rows(), fdv1T(0));
			fColtemp2.Set(fdv2T.Rows(), fdv2T(1));
			fRHS.AddCombination(-dphi*fv2[1]/magtan, fColtemp1,
				                -dphi*fv1[0]/magtan, fColtemp2);
			
			fColtemp1.Set(fdv1T.Rows(), fdv1T(1));
			fColtemp2.Set(fdv2T.Rows(), fdv2T(0));
			fRHS.AddCombination(dphi*fv2[0]/magtan, fColtemp1,
				                dphi*fv1[1]/magtan, fColtemp2);
					
			/* get equation numbers */
			fEqnos[0].RowAlias(i, eqnos);

			/* assemble */
			ElementSupport().AssembleRHS(Group(), fRHS, eqnos);

			/* store force vector output */
			int index = fStrikerTags_map.Map(pelem[2]);
			fStrikerForce2D(index,0) = dphi*tangent[1];
			fStrikerForce2D(index,1) =-dphi*tangent[0];
		}
	}

	/* set tracking */
	SetTrackingData(num_contact, h_max);
}
