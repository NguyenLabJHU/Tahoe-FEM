/* $Id: AugLagContact3DT.cpp,v 1.4 2004-06-17 07:13:38 paklein Exp $ */
#include "AugLagContact3DT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ifstreamT.h"
#include "eIntegratorT.h"
#include "ElementSupportT.h"
#include "XDOF_ManagerT.h"

using namespace Tahoe;

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

/* parameters */
const int kNumAugLagDOF = 1;

/* constructor */
AugLagContact3DT::AugLagContact3DT(const ElementSupportT& support, const FieldT& field):
	Contact3DT(support, field),
	fElCoord(fNumFacetNodes + 1, NumSD()),
	fElRefCoord(fNumFacetNodes + 1, NumSD()),
	fElDisp(fNumFacetNodes + 1, NumDOF()),
	fdc_du(NumSD(), fElDisp.Length()),
	fdn_du(NumSD(), fElDisp.Length()),
	fM1(NumSD()),
	fM2(NumSD(), fElDisp.Length()),
	fV1(fElDisp.Length())
{
	/* regularization parameter */
	const char caller[] = "AugLagContact3DT::AugLagContact3DT";
	ElementSupport().Input() >> fr;
	if (fr < 0.0)
		ExceptionT::BadInputValue(caller, "regularization must be > 0: %g", fr);

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
}

/* allocates space and reads connectivity data */
void AugLagContact3DT::Initialize(void)
{
	/* inherited */
	Contact3DT::Initialize();

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

/* append element equations numbers to the list */
void AugLagContact3DT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_2)

	/* collect using method allowing mixed node/tag numbers */
	ElementSupport().XDOF_Manager().XDOF_SetLocalEqnos(Group(), fXDOFConnectivities, fXDOFEqnos);

	/* add to list */
	eq_1.Append(&fXDOFEqnos);
}

/* returns the array for the DOF tags needed for the current config */
void AugLagContact3DT::SetDOFTags(void)
{
	/* DOF space about to be reset. store history */
	dArrayT constraints;
	constraints.Alias(ElementSupport().XDOF_Manager().XDOF(this, 0));
	fLastDOF = constraints;

	/* resize DOF tags array */
	fContactDOFtags.Dimension(fActiveStrikers.Length());
}

iArrayT& AugLagContact3DT::DOFTags(int tag_set)
{
#pragma unused(tag_set)
#if __option(extended_errorcheck)
	/* check */
	if (tag_set != 0)
		ExceptionT::OutOfRange("AugLagContact3DT::DOFTags",
			"group only has 1 tag set: %d", tag_set);
#endif
	
	return fContactDOFtags;
}

/* generate element data (based on current striker/body data) */
void AugLagContact3DT::GenerateElementData(void)
{
	/* inherited - set nodal connectivities */
	Contact3DT::SetConnectivities();

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
		pxelem[2] = pelem[2]; // 3rd facet node
		pxelem[3] = pelem[3]; // striker node
		pxelem[4] = fContactDOFtags[i]; // contact DOF tag
	}
}

/* return the contact elements */
const iArray2DT& AugLagContact3DT::DOFConnects(int tag_set) const
{
#pragma unused(tag_set)
#if __option(extended_errorcheck)
	/* check */
	if (tag_set != 0)
		ExceptionT::OutOfRange("AugLagContact3DT::DOFConnects", 
			"group only has 1 tag set: %d", tag_set);
#endif

	return fXDOFConnectivities;
}

/* restore the DOF values to the last converged solution */
void AugLagContact3DT::ResetDOF(dArray2DT& DOF, int tag_set) const
{
#pragma unused(tag_set)
#if __option(extended_errorcheck)
	/* check */
	if (tag_set != 0)
		ExceptionT::OutOfRange("AugLagContact3DT::ResetDOF",
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
int AugLagContact3DT::Reconfigure(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = Contact3DT::RelaxSystem();
	if (relax != GlobalT::kNoRelax)
		return 1;
	else
		return 0;
}

/* element level reconfiguration for the current solution */
GlobalT::RelaxCodeT AugLagContact3DT::RelaxSystem(void)
{
	/* override all inherited - relaxation handled through
	 * DOFElementT interface */
	return GlobalT::kNoRelax;
}

/* appends group connectivities to the array */
void AugLagContact3DT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
	/* inherited */
	Contact3DT::ConnectsU(connects_1, connects_2);

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
void AugLagContact3DT::ReadRestart(istream& in)
{
#pragma unused(in)
	ExceptionT::GeneralFail("AugLagContact3DT::ReadRestart", "has not been tested");
}

void AugLagContact3DT::WriteRestart(ostream& out) const
{
#pragma unused(out)
	ExceptionT::GeneralFail("AugLagContact3DT::WriteRestart", "has not been tested");
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* step in setting contact configuration. */
bool AugLagContact3DT::SetActiveInteractions(void)
{
	/* striker map about to be reset. store history */
	fLastActiveMap = fActiveMap;

	/* inherited */
	return Contact3DT::SetActiveInteractions();
}

/* print element group data */
void AugLagContact3DT::PrintControlData(ostream& out) const
{
	/* inherited */
	Contact3DT::PrintControlData(out);

	/* regularization */
	out << " Regularization parameter. . . . . . . . . . . . = " << fr << '\n';	
}

void AugLagContact3DT::RHSDriver(void)
{
	/* time integration parameters */
	double constKd = 0.0;
	int     formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* references to global nodal data */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& disp = Field()[0]; /* displacements */
	const dArray2DT& constr = ElementSupport().XDOF_Manager().XDOF(this, 0);
	const dArrayT force(constr.MajorDim(), constr.Pointer()); // general for all
	                                                          // value of kNumAugLagDOF
	/* work space */
	dArrayT c(3), n(3);
	double a[3], b[3];
	iArrayT eqnos;

	/* reset tracking data */
	int num_contact = 0;
	double h_max = 0.0;

	int neq = NumElementNodes()*NumDOF() + 1;
	dArrayT uRHS(neq - kNumAugLagDOF, fRHS.Pointer());
	for (int i = 0; i < fXDOFConnectivities.MajorDim(); i++)
	{
		int* pelem = fXDOFConnectivities(i);

		/* collect element configuration */
		fElRefCoord.RowCollect(pelem, init_coords);
		fElDisp.RowCollect(pelem, disp);

		/* current configuration */
		fElCoord.SumOf(fElRefCoord, fElDisp);
	
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
		double h = Dot(n.Pointer(), c.Pointer());

		/* augmented Lagrangian */
		double g = force[i] + fr*h;

		/* store for output */
		fActiveStrikersForce[i] = force[i];
		
		/* contact */
		if (g < 0.0)
		{
			/* initialize */
			fRHS = 0.0;

			/* tracking data */
			num_contact++;
			h_max = (h < h_max) ? h : h_max;

			/* d_c */
			fdc_du.MultTx(n, fV1);
			uRHS.SetToScaled(-g, fV1);

			/* d_normal */
			Set_dn_du(fElCoord, fdn_du);					
			fM1.Outer(n, n);
			fM1.PlusIdentity(-1.0);
			fM2.MultATB(fM1, fdn_du);
			fM2.MultTx(c, fV1);
			uRHS.AddScaled(g/mag, fV1);
								
			/* augmented Lagrangian DOF */				
			fRHS[neq - 1] = -h;					
		}
		else /* gap */
		{
			/* grad_disp contribution */
			fRHS = 0.0;

			/* augmented Lagrangian DOF */				
			fRHS[neq - 1] = force[i]/fr;							
		}

		/* get equation numbers */
		fXDOFEqnos.RowAlias(i, eqnos);

		/* assemble */
		ElementSupport().AssembleRHS(Group(), fRHS, eqnos);
	}

	/* set tracking */
	SetTrackingData(num_contact, h_max);
}

/* called by FormRHS and FormLHS */
void AugLagContact3DT::LHSDriver(GlobalT::SystemTypeT)
{
	/* time integrator parameters */
	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	/* references to global nodal data */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& disp = Field()[0]; /* displacements */
	const dArray2DT& constr = ElementSupport().XDOF_Manager().XDOF(this, 0);
	const dArrayT force(constr.MajorDim(),constr.Pointer());

	/* work space */
	dArrayT c(3), n(3);
	double a[3], b[3];
	iArrayT eqnos;
	
	/* loop over active elements */
	int neq = NumElementNodes()*NumDOF() + 1;
	dMatrixT uLHS(neq - kNumAugLagDOF);
	dArrayT  uRHS(neq - kNumAugLagDOF, fRHS.Pointer());
	for (int i = 0; i < fXDOFConnectivities.MajorDim(); i++)
	{
		int* pelem = fXDOFConnectivities(i);

		/* collect element configuration */
		fElRefCoord.RowCollect(pelem, init_coords);
		fElDisp.RowCollect(pelem, disp);

		/* current configuration using effective displacement */
		fElCoord.SumOf(fElRefCoord, fElDisp);
	
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
		double h = Dot(n.Pointer(), c.Pointer());

		/* augmented Lagrangian */
		double g = force[i] + fr*h;

		/* contact */
		if (g < 0.0)
		{
			/* initialize */
			fRHS = 0.0;
			fLHS = 0.0;
			uLHS = 0.0;
	
			/* second variation of gap */
			DDg_tri_facet(
				fElRefCoord(0), fElRefCoord(1), fElRefCoord(2), fElRefCoord(3),
				fElDisp(0), fElDisp(1), fElDisp(2), fElDisp(3),
				uLHS);
			uLHS *= g*constK; 

			/* d_c */
			fdc_du.MultTx(n, fV1);
			uRHS.SetToScaled(1.0, fV1);

			/* d_normal */
			Set_dn_du(fElCoord, fdn_du);					
			fM1.Outer(n, n);
			fM1.PlusIdentity(-1.0);
			fM2.MultATB(fM1, fdn_du);
			fM2.MultTx(c, fV1);
			uRHS.AddScaled(-1.0/mag, fV1);

			/* add term g^T g */
			uLHS.Outer(uRHS, uRHS, fr*constK, dMatrixT::kAccumulate);

			/* assemble sub-block */
			fLHS.AddBlock(0, 0, uLHS);
			
			/* augmented Lagrangian DOF */
			int dex = neq - 1;
			fRHS *= constK;
			fLHS.SetRow(dex, fRHS);
			fLHS.SetCol(dex, fRHS);
		}
		else /* gap */
		{
			/* initialize */
			fLHS = 0.0;
		
			/* augmented Lagrangian DOF */
			int dex = neq - 1;
			fLHS(dex,dex) = -1.0/fr;							
		}		

		/* get equation numbers */
		fXDOFEqnos.RowAlias(i, eqnos);

		/* assemble */
		ElementSupport().AssembleLHS(Group(), fLHS, eqnos);
	}
}
