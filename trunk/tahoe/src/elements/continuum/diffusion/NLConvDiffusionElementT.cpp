/* $Id: NLConvDiffusionElementT.cpp,v 1.2 2016-11-05 15:46:19 tdnguye Exp $ */
#include "NLConvDiffusionElementT.h"

#include <iostream>
#include <iomanip>
#include <cmath>

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
#include "ScheduleT.h"

using namespace Tahoe;
const int kRadialDirection = 0; /* the x direction is radial */
/*TDN: Why didn't you multiply integration by 2Pi?*/
const double Pi2 = 2.0*acos(-1.0);
/* constructor */
NLConvDiffusionElementT::NLConvDiffusionElementT(const ElementSupportT& support):
	DiffusionElementT(support),
    fSchedule(NULL),
    fLocCurrCoords(LocalArrayT::kCurrCoords),
    fLocDisplacement(NULL),
   fLocDisplacement_last(NULL),
	feps(0.0),
	falpha(0.0)
{
	//cout <<" \n  NLConvDiffusionElementT::NLConvDiffusionElementT: "<<endl;
    SetName("nonlinear_diffusion_general_convection");
}

/* collecting element group equation numbers */
void NLConvDiffusionElementT::Equations(AutoArrayT<const iArray2DT*>& eq_1, AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	//cout <<" \n  NLConvDiffusionElementT::Equations "<<endl;
    /* inherited */
	DiffusionElementT::Equations(eq_1, eq_2);
	
	/* collect equation numbers for nonlinear boundary conditions */
	fBCEqnos.Dimension(fBCFaces);
	Field().SetLocalEqnos(fBCFaces, fBCEqnos);
}

/* information about subordinate parameter lists */
void NLConvDiffusionElementT::DefineSubs(SubListT& sub_list) const
{
	//cout <<" \n  NLConvDiffusionElementT::DefineSubs: "<<endl;
    /* inherited */
	DiffusionElementT::DefineSubs(sub_list);

	/* mixed boundary condition */
	sub_list.AddSub("convection_bc", ParameterListT::Any);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* NLConvDiffusionElementT::NewSub(const StringT& name) const
{
	//cout <<" \n  NLConvDiffusionElementT::NewSub: "<<endl;
    if (name == "convection_bc") {
		ParameterContainerT* convection_bc = new ParameterContainerT(name);
		convection_bc->SetDescription("outward flux: h = epsilon*(T - T_wall)^alpha");
		convection_bc->SetSubSource(this);
	
		/* side sets */
		convection_bc->AddSub("side_set_ID_list", ParameterListT::Once);
		convection_bc->AddParameter(ParameterT::Integer, "schedule");
	
		/* flux parameters */
        ParameterT eps(ParameterT::Double, "epsilon");
        eps.AddLimit(0.0, LimitT::Lower);
		convection_bc->AddParameter(eps);
        
        ParameterT alpha(ParameterT::Double, "alpha");
        alpha.AddLimit(0.0, LimitT::Lower);
		convection_bc->AddParameter(alpha);
        
		return convection_bc;
	}
	else /* inherited */
		return DiffusionElementT::NewSub(name);
}

/* accept parameter list */
void NLConvDiffusionElementT::TakeParameterList(const ParameterListT& list)
{
	//cout <<" \n  NLConvDiffusionElementT::TakeParameterList: "<<endl;
    /* inherited */
	DiffusionElementT::TakeParameterList(list);
    /* dimensions */
    int nip = NumIP();
    
    /* integration point radii over the current element */
    fRadius_X.Dimension(nip);
    fRadius_x.Dimension(nip);
	/* dimension work space */
	fField_list.Dimension(NumIP());

	/* extract information about mixed bc's */
	TakeTractionBC(list);
}

/***********************************************************************
 * Protected
 ***********************************************************************/
//R Xiao added SetLocalArrays
void NLConvDiffusionElementT::SetLocalArrays(void)
{
    /* inherited */
    DiffusionElementT::SetLocalArrays();
    fLocCurrCoords.Dimension(NumElementNodes(), NumSD());
    ElementSupport().RegisterCoordinates(fLocCurrCoords);
    
}

/* construct the effective mass matrix */
void NLConvDiffusionElementT::LHSDriver(GlobalT::SystemTypeT sys_type)
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
    /*R Xiao added */
  // bool axisymmetric = Axisymmetric();
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
               /* R Xiao added */
                int ip = fShapes->CurrIP();
                double r = fRadius_x[ip];
                double scale = constC*(*Det++)*(*Weight++)*fCurrMaterial->Capacity()*r*Pi2;
			
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
void NLConvDiffusionElementT::RHSDriver(void)
{
 //  cout <<" \n  I reached NLConvDiffusionElementT::RHSDriver: "<<endl;
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
//   cout <<" \n  I reached NLConvDiffusionElementT::RHSDriver2: "<<endl;
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
                /* R Xiao added */
                int ip = fShapes->CurrIP();
                double r = fRadius_x[ip];
                /* capacity */
				double pc = fCurrMaterial->Capacity();

				/* interpolate nodal values to ip */
				fShapes->InterpolateU(fLocVel, fDOFvec);
					
				/* ip sources */
				if (block_source) fDOFvec.AddScaled(-by_dt/pc, ip_source(fShapes->CurrIP()));

				/* accumulate in element residual force vector */				
				double*	res = fRHS.Pointer();
				const double* Na = fShapes->IPShapeU();
				
				double temp = -constCv*(*Weight++)*(*Det++)*pc*r*Pi2;
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
void NLConvDiffusionElementT::FormKd(double constK)
{
//	cout <<" \n  NLConvDiffusionElementT::FormKd: "<<endl;
    /* integration parameters */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	
	int nsd = NumSD();
	dMatrixT grad;
	dArrayT field;
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
        /* R Xiao added */
        int ip = fShapes->CurrIP();
        double r = fRadius_x[ip];
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
		fRHS.AddScaled(-constK*(*Weight++)*(*Det++)*r*Pi2, fNEEvec);
	}	
}

/* form the element stiffness matrix */
void NLConvDiffusionElementT::FormStiffness(double constK)
{
//	cout <<" \n  NLConvDiffusionElementT::FormStiffness: "<<endl;
    /* must be nonsymmetric */
	if (fLHS.Format() != ElementMatrixT::kNonSymmetric)
		ExceptionT::GeneralFail("NLConvDiffusionElementT::FormStiffness",
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
        /*R Xiao added*/
        int ip = fShapes->CurrIP();
        double r = fRadius_x[ip];
        double scale = constK*(*Det++)*(*Weight++)*r*Pi2;

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
MaterialSupportT* NLConvDiffusionElementT::NewMaterialSupport(MaterialSupportT* p) const
{
	//cout <<" \n  NLConvDiffusionElementT::NewMaterialSupport: "<<endl;
    /* allocate */
	if (!p) p = new DiffusionMatSupportT(NumDOF(), NumIP());

	/* inherited initializations */
	DiffusionElementT::NewMaterialSupport(p);
	
	/* set DiffusionMatSupportT fields */
	DiffusionMatSupportT* ps = TB_DYNAMIC_CAST(DiffusionMatSupportT*, p);
	if (ps) {
     /////   ps->SetLocalArray(fLocLastDisp);
		ps->SetField(&fField_list);
	}

	return p;
}

/*void NLConvDiffusionElementT::SetLocalArrays(void)
{
    /* inherited */
 //   ContinuumElementT::SetLocalArrays();
    
    /* allocate */
//    int nen = NumElementNodes();
//    fLocLastDisp.Dimension(nen, NumDOF());
    
    /* register */
//    Field().RegisterLocal(fLocLastDisp);
//}


/***********************************************************************
 * Private
 ***********************************************************************/

void NLConvDiffusionElementT::TakeTractionBC(const ParameterListT& list)
{
//	cout <<" \n  NLConvDiffusionElementT::TakeTractionBC: "<<endl;
    const char caller[] = "NLConvDiffusionElementT::TakeTractionBC";

	/* quick exit */
    int num_conv_bc = list.NumLists("convection_bc");
//    cout << "\nnum_conv_bc: "<<num_conv_bc<<endl;
	if (num_conv_bc == 0)
		return;
/*count the total number of faces*/
    int tot_num_faces = 0;
    int num_face_nodes;
    for (int j = 0; j<num_conv_bc; j++)
    {
        const ParameterListT& convection_bc = list.GetList("convection_bc",j);
        const ParameterListT& ss_ID_list = convection_bc.GetList("side_set_ID_list");
        int num_sides = ss_ID_list.NumLists("String");
        if (num_sides > 0)
        {
            /* model manager */
            ModelManagerT& model = ElementSupport().ModelManager();
            
            /* total number of faces */
            int num_faces = 0;
            ArrayT<StringT> side_ID(num_sides);
            for (int i = 0; i < side_ID.Length(); i++)
            {
                side_ID[i] = ss_ID_list.GetList("String", i).GetParameter("value");
                num_faces += model.SideSetLength(side_ID[i]);
            }
            tot_num_faces += num_faces;
            /* element topology */
            iArrayT nodes_on_faces(fShapes->NumFacets());
            fShapes->NumNodesOnFacets(nodes_on_faces);
            int min, max;
            nodes_on_faces.MinMax(min, max);
            if (min != max) ExceptionT::GeneralFail(caller, "all faces must have same shape");
            if (j==0) num_face_nodes = min;
            if (min != num_face_nodes) ExceptionT::GeneralFail(caller, "all surfaces must have same facet shape");
        }
    }
    feps.Dimension(tot_num_faces);
    falpha.Dimension(tot_num_faces);
    fschedulenum.Dimension(tot_num_faces);
    fBCFaces.Dimension(tot_num_faces, num_face_nodes);

    /*read in parameters for each convection surface bc*/
    int face_num = 0;
    for (int j = 0; j<num_conv_bc; j++)
    {
        const ParameterListT& convection_bc = list.GetList("convection_bc",j);
        const ParameterListT& ss_ID_list = convection_bc.GetList("side_set_ID_list");
        /* extract BC parameters */
        double eps   = convection_bc.GetParameter("epsilon");
        double alpha = convection_bc.GetParameter("alpha");
        double schedule_num = convection_bc.GetParameter("schedule");
        
        int num_sides = ss_ID_list.NumLists("String");
        if (num_sides > 0)
        {
            /* model manager */
            ModelManagerT& model = ElementSupport().ModelManager();

            /* collect nodes on faces */
            ArrayT<StringT> side_ID(num_sides);
            for (int i = 0; i < side_ID.Length(); i++)
            {
                side_ID[i] = ss_ID_list.GetList("String", i).GetParameter("value");
                int num_sides = model.SideSetLength(side_ID[i]);
                for (int k = 0; k<num_sides; k++)
                {
                    int index = face_num+k;
                    feps[index] = eps;
                    falpha[index] = alpha;
                    fschedulenum[index] = schedule_num;
                }
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
}

/* compute contribution to RHS from mixed BC's */
void NLConvDiffusionElementT::TractionBC_RHS(void)
{
//	cout <<" \n  NLConvDiffusionElementT::TractionBC_RHS: "<<endl;
    /* quick exit */
	if (fBCFaces.MajorDim() == 0) return;

	/* dimensions */
	int nsd = NumSD();
	int ndof = NumDOF();
	int nfn = fBCFaces.MinorDim();
	
	/* force vector */
	dArrayT rhs(nfn*ndof);
		
/*TDN: Rui, I don't understand what you're doing here.  I think that you're getting the inital coordinates and then trying to get the displacements.  But the tag kDisp here does not refer to the dispalcement, rather the temperature.  The Field() here refers to the temperature field*/
    
	/* local coordinates */
	LocalArrayT coords(LocalArrayT::kInitCoords, nfn, nsd);
	ElementSupport().RegisterCoordinates(coords);
		
	/* nodal field values */
	LocalArrayT field(LocalArrayT::kDisp, nfn, ndof);
	Field().RegisterLocal(field);
    
    /*local nodal displacements */
/*
    LocalArrayT locdisp(LocalArrayT::kDisp, nfn, nsd);
    Field().RegisterLocal(locdisp);
*/
    /*TDN: Replacing your code above*/
    const FieldT* displacement = ElementSupport().Field("displacement");
    LocalArrayT locdisp(LocalArrayT::kDisp, nfn, nsd);
    if (displacement)
        displacement->RegisterLocal(locdisp);
    
	/* boundary shape functions - using face 0 */
	const ParentDomainT& surf_shape = ShapeFunction().FacetShapeFunction(0);
	int nip = surf_shape.NumIP();
	const dArray2DT& Na_all = surf_shape.Na();
	dArray2DT ip_field(nip, ndof);
      dArrayT ip_coords(2);
      dArrayT ip_displ(2);
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

        int schedule_num = fschedulenum[i];
        fSchedule = ElementSupport().Schedule(--schedule_num);
        
        double eps = feps[i];
        double alpha = falpha[i];
        
		/* local values */
		coords.SetLocal(nodes);
		field.SetLocal(nodes);
//        cout << "\nfield: "<<field;

        if(displacement)
            locdisp.SetLocal(nodes);
        else locdisp=0.0;
//        cout << "\nlocdisp: "<<locdisp;
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
            /* R Xiao added */
            surf_shape.Interpolate(coords, ip_coords,j);
            surf_shape.Interpolate(locdisp, ip_displ,j);
            double r = ip_coords[kRadialDirection]+ip_displ[kRadialDirection];
			/* ip weight */
			double jw = detj*w[j]*r*Pi2;
					
			/* flux */
			double qn = -eps*pow(ip_field[j] - fSchedule->Value(), alpha);
//            cout <<"\nval: "<<fSchedule->Value()<< "\tqn: "<<qn;
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
void NLConvDiffusionElementT::TractionBC_LHS(void)
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
    /*local nodal displacements */
/*    LocalArrayT locdisp(LocalArrayT::kDisp, nfn, nsd);
    Field().RegisterLocal(locdisp);
*/
	/* nodal field values */
	LocalArrayT field(LocalArrayT::kDisp, nfn, ndof);
	Field().RegisterLocal(field);

    /*TDN: Replacing your code above*/
    const FieldT* displacement = ElementSupport().Field("displacement");
    LocalArrayT locdisp(LocalArrayT::kDisp, nfn, nsd);
    if (displacement)
        displacement->RegisterLocal(locdisp);


	/* boundary shape functions - using face 0 */
	const ParentDomainT& surf_shape = ShapeFunction().FacetShapeFunction(0);
	int nip = surf_shape.NumIP();
	const dArray2DT& Na_all = surf_shape.Na();
	dArray2DT ip_field(nip, ndof);
     dArrayT ip_coords(2);
    dArrayT ip_displ(2);

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
			
        int schedule_num = fschedulenum[i];
        fSchedule = ElementSupport().Schedule(--schedule_num);
        
        double eps = feps[i];
        double alpha = falpha[i];

		/* local values */
		coords.SetLocal(nodes);
		field.SetLocal(nodes);
        if(displacement)
            locdisp.SetLocal(nodes);
        else locdisp=0.0;
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
            /* R Xiao added */
            surf_shape.Interpolate(coords, ip_coords,j);
            surf_shape.Interpolate(locdisp, ip_displ,j);
            double r = ip_coords[0]+ip_displ[0];
			/* ip weight */
			double jw = detj*w[j]*r*Pi2;
					
			/* d_flux */
			double d_qn = eps*alpha*pow(ip_field[j] - fSchedule->Value(), alpha-1.);
			
			/* shape functions */
			Na_all.RowAlias(j, Na);

			/* accumulate */
			lhs.Outer(Na, Na, jw*d_qn, dMatrixT::kAccumulate);			
		}

		/* assemble */
		ElementSupport().AssembleLHS(Group(), lhs, eqnos);
	}
}

/* form shape functions and derivatives */
void NLConvDiffusionElementT::SetGlobalShape()
{
//    cout << "\n NLConvDiffusionElementT::SetGlobalShape"<<endl;
    ContinuumElementT::SetGlobalShape();
    /* what needs to get computed */
    
    /* get current element coordinates and temperature */
    SetLocalX(fLocCurrCoords);
    SetLocalU(fLocDisp);
    
    const LocalArrayT& Temp = *fLocDisplacement;

    /* get nodal temperatures if available */
	if (fLocDisplacement)SetLocalU(*fLocDisplacement);
	if (fLocDisplacement_last) SetLocalU(*fLocDisplacement_last);

    
 ////   SetLocalU(fLocLastDisp);
    int nen = fLocCurrCoords.NumberOfNodes();
    /* loop over integration points */
    for (int i = 0; i < NumIP(); i++)
    {
        /* compute radii */
        const double* NaX = fShapes->IPShapeX(i);
//        cout << "\n Loc Init: "<<fLocInitCoords<<endl;
//        cout << "\nr: "<<kRadialDirection<<endl;
        
        const double* X_r = fLocInitCoords(kRadialDirection);
        
        double R=0.0;
        double u=0.0;
        double u_last =0.0;
        /*TDN: Add this so it will run if there are no displacements*/
        if (fLocDisplacement)
        {
            const double* NaU = fShapes->IPShapeU(i);
            int nun = Temp.NumberOfNodes();
            const double* u_r = Temp(kRadialDirection);
//            cout << "n: Xf: "<<*X_r<<endl;
     //   const double* u_r_last = (needs_F_last) ? fLocLastDisp(kRadialDirection) : u_r; /* fLocLastDisp not used */
 ////       const double* u_r_last=fLocLastDisp(kRadialDirection);
            if (nen == nun)
            {
                for (int a = 0; a < nen; a++) {
                    R += (*NaX)*(*X_r++);
                    u += (*NaU)*(*u_r++);
                    /////       u_last += (*NaU)*(*u_r_last++);
                    NaX++;
                    NaU++;
                }
            }
            else /* separate loops for field and geometry */
            {
                for (int a = 0; a < nen; a++) {
                    R += (*NaX)*(*X_r++);
                    NaX++;
                }
                for (int a = 0; a < nun; a++) {
                    u += (*NaU)*(*u_r++);
         ////       u_last += (*NaU)*(*u_r_last++);
                    NaU++;
                }
            }
        }
        else
        {
            for (int a = 0; a < nen; a++) {
                R += (*NaX)*(*X_r++);
                NaX++;
            }
        }
        double r = R + u;
//          cout <<" \n R =  "<<R<<endl;
//        cout <<" \n u =  "<<u<<endl;
        fRadius_X[i] = R;
        fRadius_x[i] = r;
    }
}

