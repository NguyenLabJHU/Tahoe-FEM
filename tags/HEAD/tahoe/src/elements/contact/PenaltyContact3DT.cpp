/* $Id: PenaltyContact3DT.cpp,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (02/09/2000)                                          */

#include "PenaltyContact3DT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
#include "FEManagerT.h"
#include "eControllerT.h"
#include "NodeManagerT.h"

/* vector functions */
inline static void CrossProduct(const double* A, const double* B, double* AxB)
{   AxB[0] = A[1]*B[2] - A[2]*B[1];
	AxB[1] = A[2]*B[0] - A[0]*B[2];
	AxB[2] = A[0]*B[1] - A[1]*B[0];
};

inline static double Dot(const double* A, const double* B)
{ return A[0]*B[0] + A[1]*B[1] + A[2]*B[2]; };

inline static void Vector(const double* start, const double* end, double* v)
{
	v[0] = end[0] - start[0];
	v[1] = end[1] - start[1];
	v[2] = end[2] - start[2];
};

/* constructor */
PenaltyContact3DT::PenaltyContact3DT(FEManagerT& fe_manager):
	Contact3DT(fe_manager),
	fElCoord(fNumFacetNodes + 1, fNumSD),
	fElDisp(fNumFacetNodes + 1, fNumDOF),
	fdc_du(fNumSD, fElDisp.Length()),
	fdn_du(fNumSD, fElDisp.Length()),
	fM1(fNumSD),
	fM2(fNumSD, fElDisp.Length()),
	fV1(fElDisp.Length()),
	fnum_contact(0),
	fh_max(0.0)
{
	fFEManager.Input() >> fK;
	if (fK < 0.0)
	{
		cout << "\n PenaltyContact3DT::PenaltyContact3DT: reguralization must be > 0: "
		     << fK << endl;
		throw eBadInputValue;
	}

	double third = 1.0/3.0;
	double* p = fdc_du.Pointer();
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ = 1;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ = 1;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p   = 1;
	
	/* set console access */
	iAddVariable("penalty_parameter", fK);
}

/* print/compute element output quantities */
void PenaltyContact3DT::WriteOutput(IOBaseT::OutputModeT mode)
{
	/* inherited */
	Contact3DT::WriteOutput(mode);

	/* contact statistics */
	ostream& out = fFEManager.Output();
	out << " Number of contact interactions = " << fnum_contact << '\n';
	out << " Maximum penetration depth      = " << fh_max       << '\n';
}

/***********************************************************************
* Protected
***********************************************************************/

/* print element group data */
void PenaltyContact3DT::PrintControlData(ostream& out) const
{
	/* inherited */
	Contact3DT::PrintControlData(out);

	/* regularization */
	out << " Regularization parameter. . . . . . . . . . . . = " << fK << '\n';	
}

/* called by FormRHS and FormLHS */
void PenaltyContact3DT::LHSDriver(void)
{
	double constK = 0.0;
	int formK = fController->FormK(constK);
	if (!formK) return;

//TEMP - consistent tangent not implemented
	fLHS.Identity(fK);
	
	/* loop over active elements */
	iArrayT eqnos;
	for (int i = 0; i < fConnectivities.MajorDim(); i++)
	{
		int* pelem = fConnectivities(i);

		/* contact */
		if (fDists[i] < 0.0)
		{
			/* get equation numbers */
			fEqnos.RowAlias(i, eqnos);
			
			/* assemble */
			fFEManager.AssembleLHS(fLHS, eqnos);
		}
	}
}

void PenaltyContact3DT::RHSDriver(void)
{
	/* time integration parameters */
	double constKd = 0.0;
	int     formKd = fController->FormKd(constKd);
	if (!formKd) return;

	/* references to global nodal data */
	const dArray2DT& init_coords = fNodes->InitialCoordinates();
	const dArray2DT& disp = fNodes->Displacements();

	/* loop over active elements */
	dArrayT c(3), n(3);
	iArrayT eqnos;
	double a[3], b[3];

	/* reset tracking data */
	fnum_contact = 0;
	fh_max = 0.0;

	fDists.Allocate(fConnectivities.MajorDim());
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
		fElCoord.RowAlias(2, fx3);
		fElCoord.RowAlias(3, fStriker);

		/* facet normal (direction) = a x b */
		Vector(fx1.Pointer(), fx2.Pointer(), a);
		Vector(fx1.Pointer(), fx3.Pointer(), b);
		CrossProduct(a, b, n.Pointer());
		double mag = sqrt(Dot(n.Pointer(), n.Pointer()));
		n[0] /= mag;
		n[1] /= mag;
		n[2] /= mag;

		/* (store) distance to facet */
		c[0] = fStriker[0] - (fx1[0] + fx2[0] + fx3[0])/3.0;
		c[1] = fStriker[1] - (fx1[1] + fx2[1] + fx3[1])/3.0;
		c[2] = fStriker[2] - (fx1[2] + fx2[2] + fx3[2])/3.0;
		double h = fDists[i] = Dot(n.Pointer(), c.Pointer());

		/* contact */
		if (h < 0.0)
		{
			/* tracking data */
			fnum_contact++;
			fh_max = (h < fh_max) ? h : fh_max;
		
			/* penetration force */
			double dphi =-fK*h;

			/* d_c */
			fdc_du.MultTx(n, fV1);
			fRHS.SetToScaled(dphi, fV1);

			/* d_normal */
			Set_dn_du(fElCoord, fdn_du);					
			fM1.Outer(n, n);
			fM1.PlusIdentity(-1.0);
			fM2.MultATB(fM1, fdn_du);
			fM2.MultTx(c, fV1);
			fRHS.AddScaled(-dphi/mag, fV1);
								
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

/* set surface normal derivative matrix */
void PenaltyContact3DT::Set_dn_du(const dArray2DT& curr_coords,
	dMatrixT& dn_du) const
{
	double* p = dn_du.Pointer();
	double* x1 = curr_coords(0);
	double* x2 = curr_coords(1);
	double* x3 = curr_coords(2);

	*p++ = 0;
	*p++ =-x2[2] + x3[2];
	*p++ = x2[1] - x3[1];
	*p++ = x2[2] - x3[2];
	*p++ = 0;
	*p++ =-x2[0] + x3[0];
	*p++ =-x2[1] + x3[1];
	*p++ = x2[0] - x3[0];
	*p++ = 0;
	*p++ = 0;
	*p++ = x1[2] - x3[2];
	*p++ =-x1[1] + x3[1];
	*p++ =-x1[2] + x3[2];
	*p++ = 0;
	*p++ = x1[0] - x3[0];
	*p++ = x1[1] - x3[1];
	*p++ =-x1[0] + x3[0];
	*p++ = 0;
	*p++ = 0;
	*p++ =-x1[2] + x2[2];
	*p++ = x1[1] - x2[1];
	*p++ = x1[2] - x2[2];
	*p++ = 0;
	*p++ =-x1[0] + x2[0];
	*p++ =-x1[1] + x2[1];
	*p++ = x1[0] - x2[0];
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p   = 0;
}
