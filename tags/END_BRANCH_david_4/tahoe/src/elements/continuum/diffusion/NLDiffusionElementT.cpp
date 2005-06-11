/* $Id: NLDiffusionElementT.cpp,v 1.8 2005-01-24 06:52:28 paklein Exp $ */
#include "NLDiffusionElementT.h"

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include "ElementCardT.h"
#include "ShapeFunctionT.h"
#include "eIntegratorT.h"
#include "iAutoArrayT.h"
#include "ModelManagerT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"

/* materials */
#include "NLDiffusionMaterialT.h"
#include "DiffusionMatSupportT.h"

using namespace Tahoe;

/* constructor */
NLDiffusionElementT::NLDiffusionElementT(const ElementSupportT& support):
	DiffusionElementT(support),
	feps(0.0), 
	fT0(0.0), 
	falpha(0.0)
{
	SetName("nonlinear_diffusion");
}

/* collecting element group equation numbers */
void NLDiffusionElementT::Equations(AutoArrayT<const iArray2DT*>& eq_1, AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	/* inherited */
	DiffusionElementT::Equations(eq_1, eq_2);
	
	/* collect equation numbers for nonlinear boundary conditions */
	fBCEqnos.Dimension(fBCFaces);
	Field().SetLocalEqnos(fBCFaces, fBCEqnos);
}

/* information about subordinate parameter lists */
void NLDiffusionElementT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	DiffusionElementT::DefineSubs(sub_list);

	/* mixed boundary condition */
	sub_list.AddSub("mixed_bc", ParameterListT::ZeroOrOnce);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* NLDiffusionElementT::NewSub(const StringT& name) const
{
	if (name == "mixed_bc") {
		ParameterContainerT* mixed_bc = new ParameterContainerT(name);
		mixed_bc->SetDescription("outward flux: h = epsilon*(T - T_0)^alpha");
		mixed_bc->SetSubSource(this);
	
		/* side sets */
		mixed_bc->AddSub("side_set_ID_list", ParameterListT::Once);
	
		/* flux parameters */
		mixed_bc->AddParameter(  feps, "epsilon");
		mixed_bc->AddParameter(   fT0, "T0");
		mixed_bc->AddParameter(falpha, "alpha");
		return mixed_bc;
	}
	else /* inherited */
		return DiffusionElementT::NewSub(name);
}

/* accept parameter list */
void NLDiffusionElementT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	DiffusionElementT::TakeParameterList(list);

	/* dimension work space */
	fField_list.Dimension(NumIP());

	/* extract information about mixed bc's */
	TakeTractionBC(list);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct the effective mass matrix */
void NLDiffusionElementT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
	/* inherited - skip DiffusionElementT implementation */
	ContinuumElementT::LHSDriver(sys_type);

	/* set components and weights */
	double constC = 0.0;
	double constK = 0.0;
	
	int formC = fIntegrator->FormC(constC);
	int formK = fIntegrator->FormK(constK);

	/* time integration order */
	int order = fIntegrator->Order();

	/* dimensions */
	int nen = NumElementNodes();	

	/* loop over elements */
	dArrayT Na;
	Top();
	while (NextElement())
	{
		/* initialize */
		fLHS = 0.0;
		
		/* set shape function derivatives */
		SetGlobalShape();

		/* collect local field values */
		SetLocalU(fLocDisp);
		if (order == 1) SetLocalU(fLocVel);

		/* element stiffness - calculates/stored ip temperatures */
		if (formK) FormStiffness(constK);

		/* element mass */
		if (formC)
		{
			/* integration parameters */
			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();
		
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				double scale = constC*(*Det++)*(*Weight++)*fCurrMaterial->Capacity();
			
				/* shape function array */
				Na.Alias(nen, fShapes->IPShapeU());
		
				/* accumulate */
				fLHS.Outer(Na, Na, scale, dMatrixT::kAccumulate);
			}
		}
	
		/* add to global equations */
		AssembleLHS();		
	}

	/* compute contribution to LHS from mixed BC's */
	TractionBC_LHS();
}

/* form the residual force vector */
void NLDiffusionElementT::RHSDriver(void)
{
	/* inherited - skip DiffusionElementT implementation */
	ContinuumElementT::RHSDriver();

	/* set components and weights */
	double constCv = 0.0;
	double constKd = 0.0;
	
	/* components dicated by the algorithm */
	int formCv = fIntegrator->FormCv(constCv);
	int formKd = fIntegrator->FormKd(constKd);

	/* block info - needed for source terms */
	int block_dex = 0;
	int block_count = 0;
	dArray2DT ip_source;
	const ElementBlockDataT* block_data = fBlockData.Pointer(block_dex);
	const dArray2DT* block_source = Field().Source(block_data->ID());

	/* body forces */
	int formBody = 0;
	if ((fBodySchedule && fBody.Magnitude() > kSmall) || block_source) {	
		formBody = 1;
		if (!formCv) constCv = 1.0; // correct value ??
	}

	int nen = NumElementNodes();
	double dt = ElementSupport().TimeStep();
	double by_dt = (fabs(dt) > kSmall) ? 1.0/dt : 0.0; /* for dt -> 0 */
	Top();
	while (NextElement())
	{
		/* reset block info (skip empty) */
		while (block_count == block_data->Dimension()) {
			block_data = fBlockData.Pointer(++block_dex);
			block_source = Field().Source(block_data->ID());
			block_count = 0;
		}
		
		/* convert heat increment/volume to rate */
		if (block_source) ip_source.Alias(NumIP(), 1, (*block_source)(block_count));
		block_count++;
		
		/* initialize */
		fRHS = 0.0;

		/* global shape function values */
		SetGlobalShape();

		/* conduction term */
		if (formKd) 
		{
			SetLocalU(fLocDisp);
			FormKd(-constKd);
		}

		/* capacity term */
		if (formCv || formBody)
		{
			if (formCv) SetLocalU(fLocVel);
			else fLocVel = 0.0;
			if (formBody) AddBodyForce(fLocVel);

			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();

			fShapes->TopIP();
			while (fShapes->NextIP())
			{					
				/* capacity */
				double pc = fCurrMaterial->Capacity();

				/* interpolate nodal values to ip */
				fShapes->InterpolateU(fLocVel, fDOFvec);
					
				/* ip sources */
				if (block_source) fDOFvec.AddScaled(-by_dt/pc, ip_source(fShapes->CurrIP()));

				/* accumulate in element residual force vector */				
				double*	res = fRHS.Pointer();
				const double* Na = fShapes->IPShapeU();
				
				double temp = -constCv*(*Weight++)*(*Det++)*pc;
				for (int lnd = 0; lnd < nen; lnd++)
					*res++ += temp*(*Na++)*fDOFvec[0];
			}
		}
				
		/* assemble */
		AssembleRHS();
	}

	/* compute contribution to RHS from mixed BC's */
	TractionBC_RHS();
}

/* calculate the internal force contribution ("-k*d") */
void NLDiffusionElementT::FormKd(double constK)
{
	/* integration parameters */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	
	int nsd = NumSD();
	dMatrixT grad;
	dArrayT field;
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		/* set field gradient */
		grad.Set(1, nsd, fGradient_list[CurrIP()].Pointer());
		IP_ComputeGradient(fLocDisp, grad);
		
		/* interpolate field */
		field.Set(1, fField_list.Pointer(CurrIP()));
		IP_Interpolate(fLocDisp, field);

		/* get strain-displacement matrix */
		B(fShapes->CurrIP(), fB);

		/* compute heat flow */
		fB.MultTx(fCurrMaterial->q_i(), fNEEvec);

		/* accumulate */
		fRHS.AddScaled(-constK*(*Weight++)*(*Det++), fNEEvec);
	}	
}

/* form the element stiffness matrix */
void NLDiffusionElementT::FormStiffness(double constK)
{
	/* must be nonsymmetric */
	if (fLHS.Format() != ElementMatrixT::kNonSymmetric)
		ExceptionT::GeneralFail("NLDiffusionElementT::FormStiffness",
			"LHS matrix must be nonsymmetric");

	/* integration parameters */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/* dimensions */
	int nsd = NumSD();
	int nen = NumElementNodes();	

	/* time integration order */
	int order = fIntegrator->Order();
	
	/* integrate element stiffness */
	dMatrixT grad;
	dArrayT field;
	dArrayT Na;
	fShapes->TopIP();
	dArrayT dfield(1);
	while (fShapes->NextIP())
	{
		double scale = constK*(*Det++)*(*Weight++);

		/* set field gradient */
		grad.Set(1, nsd, fGradient_list[CurrIP()].Pointer());
		IP_ComputeGradient(fLocDisp, grad);
		
		/* interpolate field */
		field.Set(1, fField_list.Pointer(CurrIP()));
		IP_Interpolate(fLocDisp, field);

		/* shape function array */
		Na.Set(nen, (double*) fShapes->IPShapeU());

		/* strain displacement matrix */
		B(fShapes->CurrIP(), fB);
	
		/* nonlinear term */
		fB.MultTx(fCurrMaterial->dq_i_dT(), fNEEvec);
		fLHS.Outer(fNEEvec, Na, -scale, dMatrixT::kAccumulate);

		/* get D matrix */
		fD.SetToScaled(scale, fCurrMaterial->k_ij());
							
		/* multiply b(transpose) * db, taking account of symmetry, */
		/* and accumulate in elstif */
		fLHS.MultQTBQ(fB, fD, dMatrixT::kWhole, dMatrixT::kAccumulate);	

		/* variation in capacity */
		if (order == 1)
		{
			/* field rate */
			IP_Interpolate(fLocVel, dfield);
		
			/* accumulate */
			fLHS.Outer(Na, Na, scale*fCurrMaterial->dCapacity_dT()*dfield[0], dMatrixT::kAccumulate);
		}
	}
}

/* construct a new material support and return a pointer */
MaterialSupportT* NLDiffusionElementT::NewMaterialSupport(MaterialSupportT* p) const
{
	/* allocate */
	if (!p) p = new DiffusionMatSupportT(NumDOF(), NumIP());

	/* inherited initializations */
	DiffusionElementT::NewMaterialSupport(p);
	
	/* set DiffusionMatSupportT fields */
	DiffusionMatSupportT* ps = TB_DYNAMIC_CAST(DiffusionMatSupportT*, p);
	if (ps) {
		ps->SetField(&fField_list);
	}

	return p;
}

/***********************************************************************
 * Private
 ***********************************************************************/

void NLDiffusionElementT::TakeTractionBC(const ParameterListT& list)
{
	const char caller[] = "NLDiffusionElementT::TakeTractionBC";

	/* quick exit */
	const ParameterListT* mixed_bc = list.List("mixed_bc");
	if (!mixed_bc)
		return;

	/* extract BC parameters */
	feps   = mixed_bc->GetParameter("epsilon");
	fT0    = mixed_bc->GetParameter("T0");
	falpha = mixed_bc->GetParameter("alpha");
	
	const ParameterListT& ss_ID_list = mixed_bc->GetList("side_set_ID_list");
	int num_sides = ss_ID_list.NumLists("String");
	if (num_sides > 0)
	{
		/* model manager */
		ModelManagerT& model = ElementSupport().ModelManager();

		/* total number of faces */
		int num_faces = 0;
		ArrayT<StringT> side_ID(num_sides);
		for (int i = 0; i < side_ID.Length(); i++) {
			side_ID[i] = ss_ID_list.GetList("String", i).GetParameter("value");
			num_faces += model.SideSetLength(side_ID[i]);
		}

		/* element topology */
		iArrayT nodes_on_faces(fShapes->NumFacets());
		fShapes->NumNodesOnFacets(nodes_on_faces);
		int min, max;
		nodes_on_faces.MinMax(min, max);
		if (min != max) ExceptionT::GeneralFail(caller, "all faces must have same shape");
		
		/* collect nodes on faces */
		int face_num = 0;
		fBCFaces.Dimension(num_faces, nodes_on_faces[0]);
		for (int i = 0; i < side_ID.Length(); i++)
		{
			int num_sides = model.SideSetLength(side_ID[i]);
			iArray2DT faces(num_sides, fBCFaces.MinorDim(), fBCFaces(face_num));
		
			/* read side set */
			ArrayT<GeometryT::CodeT> facet_geom;
			iArrayT facet_nodes;
			model.SideSet(side_ID[i], facet_geom, facet_nodes, faces);		
		
			/* next set */
			face_num += num_sides;
		}		
	}
}

/* compute contribution to RHS from mixed BC's */
void NLDiffusionElementT::TractionBC_RHS(void)
{
	/* quick exit */
	if (fBCFaces.MajorDim() == 0) return;

	/* dimensions */
	int nsd = NumSD();
	int ndof = NumDOF();
	int nfn = fBCFaces.MinorDim();
	
	/* force vector */
	dArrayT rhs(nfn*ndof);
		
	/* local coordinates */
	LocalArrayT coords(LocalArrayT::kInitCoords, nfn, nsd);
	ElementSupport().RegisterCoordinates(coords);
		
	/* nodal field values */
	LocalArrayT field(LocalArrayT::kDisp, nfn, ndof);
	Field().RegisterLocal(field);

	/* boundary shape functions - using face 0 */
	const ParentDomainT& surf_shape = ShapeFunction().FacetShapeFunction(0);
	int nip = surf_shape.NumIP();
	const dArray2DT& Na_all = surf_shape.Na();
	dArray2DT ip_field(nip, ndof);

	/* Jacobian of the surface mapping */
	dMatrixT jacobian(nsd, nsd-1);

	/* loop over faces */
	iArrayT nodes, eqnos;
	dArrayT Na;
	for (int i = 0; i < fBCFaces.MajorDim(); i++)
	{
		/* face info */
		fBCFaces.RowAlias(i, nodes);
		fBCEqnos.RowAlias(i, eqnos);
			
		/* local values */
		coords.SetLocal(nodes);
		field.SetLocal(nodes);
		
		/* all ip field values: (nip x ndof) */
		surf_shape.Interpolate(field, ip_field);

		/* integrate */			
		rhs = 0.0;
		const double* w = surf_shape.Weight();
		for (int j = 0; j < nip; j++)
		{
			/* coordinate mapping */
			surf_shape.DomainJacobian(coords, j, jacobian);
			double detj = surf_shape.SurfaceJacobian(jacobian);
	
			/* ip weight */
			double jw = detj*w[j];
					
			/* flux */
			double qn = -feps*pow(ip_field[j] - fT0, falpha);
			
			/* shape functions */
			Na_all.RowAlias(j, Na);

			/* accumulate */
			rhs.AddScaled(jw*qn, Na);			
		}

		/* assemble */
		ElementSupport().AssembleRHS(Group(), rhs, eqnos);
	}
}

/* compute contribution to LHS from BC's */
void NLDiffusionElementT::TractionBC_LHS(void)
{
	/* quick exit */
	if (fBCFaces.MajorDim() == 0) return;

	/* dimensions */
	int nsd = NumSD();
	int ndof = NumDOF();
	int nfn = fBCFaces.MinorDim();
	
	/* stiffness matrix */
	ElementMatrixT lhs(nfn*ndof, ElementMatrixT::kNonSymmetric);
		
	/* local coordinates */
	LocalArrayT coords(LocalArrayT::kInitCoords, nfn, nsd);
	ElementSupport().RegisterCoordinates(coords);
		
	/* nodal field values */
	LocalArrayT field(LocalArrayT::kDisp, nfn, ndof);
	Field().RegisterLocal(field);

	/* boundary shape functions - using face 0 */
	const ParentDomainT& surf_shape = ShapeFunction().FacetShapeFunction(0);
	int nip = surf_shape.NumIP();
	const dArray2DT& Na_all = surf_shape.Na();
	dArray2DT ip_field(nip, ndof);

	/* Jacobian of the surface mapping */
	dMatrixT jacobian(nsd, nsd-1);

	/* loop over faces */
	iArrayT nodes, eqnos;
	dArrayT Na;
	for (int i = 0; i < fBCFaces.MajorDim(); i++)
	{
		/* face info */
		fBCFaces.RowAlias(i, nodes);
		fBCEqnos.RowAlias(i, eqnos);
			
		/* local values */
		coords.SetLocal(nodes);
		field.SetLocal(nodes);
		
		/* all ip field values: (nip x ndof) */
		surf_shape.Interpolate(field, ip_field);

		/* integrate */			
		lhs = 0.0;
		const double* w = surf_shape.Weight();
		for (int j = 0; j < nip; j++)
		{
			/* coordinate mapping */
			surf_shape.DomainJacobian(coords, j, jacobian);
			double detj = surf_shape.SurfaceJacobian(jacobian);
	
			/* ip weight */
			double jw = detj*w[j];
					
			/* d_flux */
			double d_qn = feps*falpha*pow(ip_field[j] - fT0, falpha-1);
			
			/* shape functions */
			Na_all.RowAlias(j, Na);

			/* accumulate */
			lhs.Outer(Na, Na, jw*d_qn, dMatrixT::kAccumulate);			
		}

		/* assemble */
		ElementSupport().AssembleLHS(Group(), lhs, eqnos);
	}
}
