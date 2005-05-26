/* $Id: PenaltyContactDrag3DT.cpp,v 1.6 2005-01-13 01:50:47 paklein Exp $ */
/* created: paklein (12/11/1997) */
#include "PenaltyContactDrag3DT.h"

#include "eIntegratorT.h"
#include "ModelManagerT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

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

/* constructor */
PenaltyContactDrag3DT::PenaltyContactDrag3DT(const ElementSupportT& support):
	PenaltyContact3DT(support),
	fDrag(0),
	fGapTolerance(0),
	fSlipTolerance(0)
{
	SetName("contact_drag_3D_penalty");
}

/* describe the parameters needed by the interface */
void PenaltyContactDrag3DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	PenaltyContact3DT::DefineParameters(list);

	list.AddParameter(fDrag, "drag_traction");
	list.AddParameter(fGapTolerance, "gap_tolerance");
	list.AddParameter(fSlipTolerance, "slip_tolerance");
}

/* accept parameter list */
void PenaltyContactDrag3DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	PenaltyContact3DT::TakeParameterList(list);

	fDrag = list.GetParameter("drag_traction");
	fGapTolerance = list.GetParameter("gap_tolerance");
	fSlipTolerance = list.GetParameter("slip_tolerance");

	/* collect volume element block ID's containing the strikers */
	ModelManagerT& model = ElementSupport().ModelManager();
	ArrayT<StringT> element_id_all;
	model.ElementGroupIDsWithNodes(fStrikerTags, element_id_all);
	iArrayT volume_element(element_id_all.Length());
	for (int i = 0; i < element_id_all.Length(); i++) {
		GeometryT::CodeT geom = model.ElementGroupGeometry(element_id_all[i]);
		volume_element[i] = (GeometryT::GeometryToNumSD(geom) == 3) ? 1 : 0;
	}
	int count = 0;
	ArrayT<StringT> element_id(volume_element.Count(1));
	for (int i = 0; i < element_id_all.Length(); i++)
		if (volume_element[i])
			element_id[count++] = element_id_all[i];

	
	/* compute associated nodal area */
	ComputeNodalArea(element_id, fNodalArea, fStrikerLocNumber);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

void PenaltyContactDrag3DT::RHSDriver(void)
{
	/* time integration parameters */
	double constKd = 0.0;
	int     formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* references to global nodal data */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& disp = Field()(0,0); /* displacements */
	const dArray2DT& disp_last = Field()(-1,0); /* displacements from last step */

	/* work space */
	dArrayT c(3), n(3);
	iArrayT eqnos;
	double a[3], b[3];

	/* reset tracking data */
	int num_contact = 0;
	double h_max = 0.0;

	const int* pelem = fConnectivities[0]->Pointer();
	int rowlength = fConnectivities[0]->MinorDim();
	for (int i = 0; i < fConnectivities[0]->MajorDim(); i++, pelem += rowlength)
	{
		/* collect element configuration */
		fElCoord.RowCollect(pelem, init_coords);
		fElDisp.RowCollect(pelem, disp);

		/* current configuration using effective displacement */
		fElCoord.AddScaled(constKd, fElDisp);
	
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

		/* contact */
		bool has_contact = false;
		if (h < 0.0)
		{
			has_contact = true;

			/* tracking data */
			num_contact++;
			h_max = (h < h_max) ? h : h_max;
		
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
		}

		/* drag */
		bool has_drag = false;
		if (h < fGapTolerance)
		{		
			/* reference facet normal (direction) = a x b */
			double n_ref[3];
			Vector(init_coords(pelem[0]), init_coords(pelem[1]), a);
			Vector(init_coords(pelem[0]), init_coords(pelem[2]), b);
			CrossProduct(a, b, n_ref);
			double ref_norm_mag = sqrt(Dot(n_ref, n_ref));
			n_ref[0] /= ref_norm_mag;
			n_ref[1] /= ref_norm_mag;
			n_ref[2] /= ref_norm_mag;

			/* displacement from the last increment */
			double drag[3];
			int striker_node = pelem[3];
			Vector(disp_last(striker_node), disp(striker_node), drag);

			/* just the tangent part of the drag */
			double inc_norm = Dot(drag, n_ref);
			drag[0] -= inc_norm*n_ref[0];
			drag[1] -= inc_norm*n_ref[1];
			drag[2] -= inc_norm*n_ref[2];
			double mag_slip = sqrt(Dot(drag, drag));
			if (mag_slip > fSlipTolerance)
			{
				has_drag = true;
				if (!has_contact) fRHS = 0.0;

				/* drag force */
				int striker_index = fStrikerLocNumber.Map(striker_node);
				double drag_force = -fNodalArea[striker_index]*fDrag;
				double f_x = drag_force*drag[0]/mag_slip;
				double f_y = drag_force*drag[1]/mag_slip;
				double f_z = drag_force*drag[2]/mag_slip;

				double f_x_by_3 = -f_x/3.0;
				double f_y_by_3 = -f_y/3.0;
				double f_z_by_3 = -f_z/3.0;
			
				/* assemble - equal and opposite force on facet nodes */
				fRHS[0] += f_x_by_3;
				fRHS[1] += f_y_by_3;
				fRHS[2] += f_z_by_3;

				fRHS[3] += f_x_by_3;
				fRHS[4] += f_y_by_3;
				fRHS[5] += f_z_by_3;

				fRHS[6] += f_x_by_3;
				fRHS[7] += f_y_by_3;
				fRHS[8] += f_z_by_3;

				fRHS[9]  += f_x;
				fRHS[10] += f_y;
				fRHS[11] += f_z;
			}
		}

		/* assemble */
		if (has_contact || has_drag)
		{
			/* get equation numbers */
			fEqnos[0].RowAlias(i, eqnos);

			/* assemble */
			ElementSupport().AssembleRHS(Group(), fRHS, eqnos);
		}
	}

	/* set tracking */
	SetTrackingData(num_contact, h_max);
}

void PenaltyContactDrag3DT::LHSDriver(GlobalT::SystemTypeT)
{
	/* time integrator parameters */
	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	/* references to global nodal data */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& disp = Field()(0,0); /* displacements */
	const dArray2DT& disp_last = Field()(-1,0); /* displacements from last step */

	/* work space */
	dArrayT c(3), n(3);
	double a[3], b[3];
	iArrayT eqnos;
	dMatrixT K_drag(3);
	
	/* loop over active elements */
	const iArray2DT& connects = *(fConnectivities[0]);
	for (int i = 0; i < connects.MajorDim(); i++)
	{
		const int* pelem = connects(i);

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

		/* contact */
		bool has_contact = false;
		if (h < 0.0)
		{
			has_contact = true;
		
			/* initialize */
			fRHS = 0.0;
			fLHS = 0.0;
	
			/* second variation of gap */
			DDg_tri_facet(
				fElRefCoord(0), fElRefCoord(1), fElRefCoord(2), fElRefCoord(3),
				fElDisp(0), fElDisp(1), fElDisp(2), fElDisp(3),
				fLHS);
			fLHS *= fK*h*constK; 

			/* d_c */
			fdc_du.MultTx(n, fV1);
			fRHS.SetToScaled(1.0, fV1);

			/* d_normal */
			Set_dn_du(fElCoord, fdn_du);					
			fM1.Outer(n, n);
			fM1.PlusIdentity(-1.0);
			fM2.MultATB(fM1, fdn_du);
			fM2.MultTx(c, fV1);
			fRHS.AddScaled(-1.0/mag, fV1);

			/* add term g^T g */
			fLHS.Outer(fRHS, fRHS, fK*constK, dMatrixT::kAccumulate);
		}
		
		/* drag */
		bool has_drag = false;
		if (h < fGapTolerance)
		{
			/* reference facet normal (direction) = a x b */
			double n_ref[3];
			Vector(init_coords(pelem[0]), init_coords(pelem[1]), a);
			Vector(init_coords(pelem[0]), init_coords(pelem[2]), b);
			CrossProduct(a, b, n_ref);
			double ref_norm_mag = sqrt(Dot(n_ref, n_ref));
			n_ref[0] /= ref_norm_mag;
			n_ref[1] /= ref_norm_mag;
			n_ref[2] /= ref_norm_mag;

			/* displacement from the last increment */
			double drag[3];
			int striker_node = pelem[3];
			Vector(disp_last(striker_node), disp(striker_node), drag);

			/* just the tangent part of the drag */
			double inc_norm = Dot(drag, n_ref);
			drag[0] -= inc_norm*n_ref[0];
			drag[1] -= inc_norm*n_ref[1];
			drag[2] -= inc_norm*n_ref[2];
			double mag_slip = sqrt(Dot(drag, drag));
			if (mag_slip > fSlipTolerance)
			{
				has_drag = true;
				if (!has_contact) fLHS = 0.0;

				/* drag force */
				int striker_index = fStrikerLocNumber.Map(striker_node);
				double drag_force = fNodalArea[striker_index]*fDrag;

				K_drag.Identity();
				K_drag.Outer(n_ref, n_ref, -1.0, dMatrixT::kAccumulate);
				K_drag.Outer(drag, drag, -1.0/(mag_slip*mag_slip), dMatrixT::kAccumulate);
				K_drag *= constK*drag_force/mag_slip;

				/* assemble */
				fLHS.AddBlock(9, 9, K_drag);
			}
		}
		
		/* assemble */
		if (has_contact || has_drag)
		{
			/* get equation numbers */
			fEqnos[0].RowAlias(i, eqnos);
			
			/* assemble */
			ElementSupport().AssembleLHS(Group(), fLHS, eqnos);
		}
	}
}
