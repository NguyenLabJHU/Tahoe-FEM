/* $Id: PenaltyContact2DT.cpp,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (12/11/1997)                                          */

#include "PenaltyContact2DT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
#include "FEManagerT.h"
#include "eControllerT.h"
#include "NodeManagerT.h"

/* parameters (duplicated from Contact2DT) */
const int kNumFacetNodes = 2;

/* constructor */
PenaltyContact2DT::PenaltyContact2DT(FEManagerT& fe_manager):
	Contact2DT(fe_manager),
	fElCoord(fNumFacetNodes + 1, fNumSD),
	fElDisp(fNumFacetNodes + 1, fNumDOF)
	
{
	fFEManager.Input() >> fK;
	if (fK < 0.0)
	{
		cout << "\n PenaltyContact2DT::PenaltyContact2DT: reguralization must be > 0: "
		     << fK << endl;
		throw eBadInputValue;
	}
}

/* print/compute element output quantities */
void PenaltyContact2DT::WriteOutput(IOBaseT::OutputModeT mode)
{
	/* inherited */
	Contact2DT::WriteOutput(mode);

	/* contact statistics */
	ostream& out = fFEManager.Output();
	out << " Number of contact interactions = " << fnum_contact << '\n';
	out << " Maximum penetration depth      = " << fh_max << '\n';
}

/***********************************************************************
* Protected
***********************************************************************/

/* print element group data */
void PenaltyContact2DT::PrintControlData(ostream& out) const
{
	/* inherited */
	Contact2DT::PrintControlData(out);

	/* regularization */
	out << " Regularization parameter. . . . . . . . . . . . = " << fK << '\n';	
}

/* called by FormRHS and FormLHS */
void PenaltyContact2DT::LHSDriver(void)
{
	double constK = 0.0;
	int formK = fController->FormK(constK);
	if (!formK) return;

	/* get reference to global coordinates */
	const dArray2DT& coords = fNodes->CurrentCoordinates(); //EFFECTIVE_DVA

	/* loop over active elements */
	dArrayT tangent(fNumSD);
	iArrayT eqnos;
	for (int i = 0; i < fConnectivities.MajorDim(); i++)
	{
		int* pelem = fConnectivities(i);
	
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

			fColtemp1.Set(fNumElemEqnos, fdv1T(0));
			fColtemp2.Set(fNumElemEqnos, fdv2T(1));
			fRHS.AddCombination(-fv2[1]/magtan,fColtemp1,
				                -fv1[0]/magtan,fColtemp2);
			
			fColtemp1.Set(fNumElemEqnos, fdv1T(1));
			fColtemp2.Set(fNumElemEqnos, fdv2T(0));
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
			fEqnos.RowAlias(i, eqnos);
			
			/* time integration factor */
			fLHS *= constK;
			
			/* assemble */
			fFEManager.AssembleLHS(fLHS, eqnos);
		}
	}
}

void PenaltyContact2DT::RHSDriver(void)
{
	/* time integration parameters */
	double constKd = 0.0;
	int     formKd = fController->FormKd(constKd);
	if (!formKd) return;

	/* references to global nodal data */
	const dArray2DT& init_coords = fNodes->InitialCoordinates();
	const dArray2DT& disp = fNodes->Displacements();

	/* reset tracking data */
	fnum_contact = 0;
	fh_max = 0.0;

	/* loop over active elements */
	dArrayT tangent(fNumSD);
	iArrayT eqnos;
	for (int i = 0; i < fConnectivities.MajorDim(); i++)
	{
		int* pelem = fConnectivities(i);

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
			fnum_contact++;
			fh_max = (h < fh_max) ? h : fh_max;

			/* penetration force */
			double dphi =-fK*h;
			
			/* initialize */
			fRHS = 0.0;
					
			/* d_tan contribution */
			fdtanT.Multx(tangent, fNEEvec);
			fRHS.AddScaled(-dphi*h/(magtan*magtan), fNEEvec);
						
			/* d_area */
			fColtemp1.Set(fNumElemEqnos,fdv1T(0));
			fColtemp2.Set(fNumElemEqnos,fdv2T(1));
			fRHS.AddCombination(-dphi*fv2[1]/magtan, fColtemp1,
				                -dphi*fv1[0]/magtan, fColtemp2);
			
			fColtemp1.Set(fNumElemEqnos,fdv1T(1));
			fColtemp2.Set(fNumElemEqnos,fdv2T(0));
			fRHS.AddCombination(dphi*fv2[0]/magtan, fColtemp1,
				                dphi*fv1[1]/magtan, fColtemp2);
					
			/* get equation numbers */
			fEqnos.RowAlias(i, eqnos);
			
			/* assemble */
			fFEManager.AssembleRHS(fRHS, eqnos);
		}
	}
}

/***********************************************************************
* Private
***********************************************************************/
