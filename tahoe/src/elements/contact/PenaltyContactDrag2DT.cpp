/* $Id: PenaltyContactDrag2DT.cpp,v 1.1 2003-08-14 05:50:41 paklein Exp $ */
/* created: paklein (12/11/1997) */
#include "PenaltyContactDrag2DT.h"
#include "ContinuumElementT.h"
#include "ModelManagerT.h"
#include "ParentDomainT.h"
#include "InverseMapT.h"
#include "fstreamT.h"
#include "eIntegratorT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

using namespace Tahoe;

/* constructor */
PenaltyContactDrag2DT::PenaltyContactDrag2DT(const ElementSupportT& support, const FieldT& field):
	PenaltyContact2DT(support, field),
	fDrag(0),
	fGapTolerance(0),
	fSlipTolerance(0)
{

}

/* initialization after constructor */
void PenaltyContactDrag2DT::Initialize(void)
{
	/* inherited */
	PenaltyContact2DT::Initialize();

	ifstreamT& in = ElementSupport().Input();

	/* drag parameters */
	in >> fDrag
	   >> fGapTolerance
	   >> fSlipTolerance;

	/* compute associated nodal area (using all element blocks) */
	const ArrayT<StringT>& element_id = ElementSupport().Model().ElementGroupIDs();
	ComputeNodalArea(element_id, fNodalArea);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* print element group data */
void PenaltyContactDrag2DT::PrintControlData(ostream& out) const
{
	/* inherited */
	Contact2DT::PrintControlData(out);

	/* regularization */
	out << " Regularization parameter. . . . . . . . . . . . = " << fK << '\n';	
}

void PenaltyContactDrag2DT::RHSDriver(void)
{
	/* time integration parameters */
	double constKd = 0.0;
	int     formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* references to global nodal data */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& disp = Field()(0,0); /* displacements */
	const dArray2DT& disp_last = Field()(-1,0); /* displacements from last step */

	/* reset tracking data */
	int num_contact = 0;
	double h_max = 0.0;

	/* loop over active elements */
	dArrayT tangent(NumSD()), drag(NumDOF());
	iArrayT eqnos;
	int* pelem = fConnectivities[0]->Pointer();
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
		bool has_contact = false;
		if (h < 0.0)
		{
			has_contact = true;
		
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
		}
		
		/* drag */
		bool has_drag = false;
		if (h < fGapTolerance)
		{
			/* displacement from the last increment */
			int striker_node = pelem[2];
			drag.DiffOf(disp_last(striker_node), disp(striker_node));
			
			/* striker is "sliding" */
			double mag_slip = dArrayT::Dot(drag, tangent)/magtan;
			if (fabs(mag_slip) > fSlipTolerance)
			{
				has_drag = true;
			
				/* drag force */
				double f_x = -fDrag*tangent[0]/magtan;
				double f_y = -fDrag*tangent[1]/magtan;
			
				/* assemble - equal and opposite force on facet nodes */
				fRHS[0] += -0.5*f_x;
				fRHS[1] += -0.5*f_y;
				fRHS[2] += -0.5*f_x;
				fRHS[3] += -0.5*f_y;
				fRHS[4] += f_x;
				fRHS[5] += f_y;
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

/* compute the nodal area associated with each striker node */
void PenaltyContactDrag2DT::ComputeNodalArea(const ArrayT<StringT>& striker_blocks, 
	dArrayT& nodal_area)
{
	/* initialize nodal area */
	nodal_area.Dimension(fStrikerTags.Length());
	nodal_area = 0.0;

	/* get surface faces */
	GeometryT::CodeT geometry;
	ArrayT<iArray2DT> surfaces;
	iArrayT surface_nodes;
	ElementSupport().Model().SurfaceFacets(striker_blocks, geometry, surfaces, surface_nodes);

	/* no surfaces */
	if (surfaces.Length() == 0) return;

	/* map to local id of striker nodes */
	InverseMapT inverse_map;
	inverse_map.SetOutOfRange(InverseMapT::MinusOne);
	inverse_map.SetMap(fStrikerTags);

	/* shape functions over the faces */
	int nip = 1;
	int nfn = surfaces[0].MinorDim();
	ParentDomainT surf_shape(geometry, nip, nfn);
	surf_shape.Initialize();

	/* coordinates over the face */
	int nsd = NumSD();
	LocalArrayT ref_coords(LocalArrayT::kInitCoords, nfn, nsd);
	ElementSupport().RegisterCoordinates(ref_coords);
	dMatrixT jacobian(nsd, nsd-1);

	/* loop over surfaces */
	const double* Na = surf_shape.Shape(0);
	const double* w  = surf_shape.Weight();
	iArrayT facet_nodes;
	for (int i = 0; i < surfaces.Length(); i++)
	{
		const iArray2DT& surface = surfaces[i];

		/* loop over faces */
		for (int j = 0; j < surface.MajorDim(); j++)
		{
			/* face nodes */
			surface.RowAlias(j, facet_nodes);
		
			/* gather coordinates */
			ref_coords.SetLocal(facet_nodes);
		
			/* coordinate mapping */
			surf_shape.DomainJacobian(ref_coords, 0, jacobian);
			double detj = surf_shape.SurfaceJacobian(jacobian);	
		
			/* loop over face nodes */
			for (int k = 0; k < facet_nodes.Length(); k++)
			{
				/* striker node index */
				int index = inverse_map.Map(facet_nodes[k]);
				if (index != -1)
					nodal_area[index] += w[0]*detj*Na[k];
			}
		}
	}
}
