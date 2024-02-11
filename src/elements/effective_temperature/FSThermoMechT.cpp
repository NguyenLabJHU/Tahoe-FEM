#include "FSThermoMechT.h"
#include "FSThermoMechSupportT.h"
#include "FSThermoMechMatT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "OutputSetT.h"
#include "ShapeFunctionT.h"
#include "eIntegratorT.h"
#include "FSThermoMechMatListT.h"
#include "ModelManagerT.h"


/* ISSUES:
1.  AddNodalForce - calls SetGlobalShape, which is modified in the Q1P0 formulation.
	Could be a problem for reaction force output accuracy
2.  ComputeOutput - all based on TL (reference shape functions), not ULagrangian
3.  Integration factor for K terms multiplied by beta*dt^2 - correct!
*/

// materials lists (3D only)
//#include "FSSolidMatList2DT.h"
//#include "FSSolidMatList3DT.h"

namespace Tahoe {

  FSThermoMechT::FSThermoMechT(const ElementSupportT& support):
    UpdatedLagrangianT(support),
	fFSThermoMechSupport(NULL),
    //  fLocTe(LocalArrayT::kUnspecified),
    fLocTe_last(LocalArrayT::kUnspecified),
    fTemperatureField(NULL),
	fLocTemperatureRate(NULL),
		fCurrMaterial(0)
  {
    SetName("thermo_mech_coupled");
  }


/* Destructor */
  FSThermoMechT::~FSThermoMechT()
  {
	  delete fFSThermoMechSupport;
      delete fTemperatureField;
      delete fLocTemperatureRate;
  }

    
  // specify parameters needed by the interface  
  void FSThermoMechT::DefineParameters(ParameterListT& list) const
  {
    // inherited
	  // additional fields
      FiniteStrainT::DefineParameters(list);

	  ParameterT temperature_field_name(ParameterT::Word, "temperature_field_name");
	  temperature_field_name.SetDefault("temperature");
	  list.AddParameter(temperature_field_name); 

  }

    /* information about subordinate parameter lists */
    void FSThermoMechT::DefineSubs(SubListT& sub_list) const
    {
        
        /* inherited */
        SolidElementT::DefineSubs(sub_list);
        
        /* element block/material specification */
        sub_list.AddSub("fsthermomech_element_block", ParameterListT::OnePlus);
        sub_list.AddSub("mixed_bc", ParameterListT::ZeroOrOnce);
    }
    
    /* a pointer to the ParameterInterfaceT of the given subordinate */
    ParameterInterfaceT* FSThermoMechT::NewSub(const StringT& name) const
    {
        /* inherited */
        if (name == "fsthermomech_element_block")
        {
            ParameterContainerT* block = new ParameterContainerT(name);
            
            /* list of element block ID's (defined by ElementBaseT) */
            block->AddSub("block_ID_list", ParameterListT::Once);
            
            /* choice of materials lists */
            block->AddSub("fsthermomech_material", ParameterListT::Once);
            
            /* set this as source of subs */
            block->SetSubSource(this);
            
            return block;
        }
        else if (name == "mixed_bc")
        {
            ParameterContainerT* mixedbc = new ParameterContainerT(name);
            mixedbc->SetDescription("outward flux: h = epsilon*(T - T_0)^alpha");
            mixedbc->SetSubSource(this);
            
            /* side sets */
            mixedbc->AddSub("side_set_ID_list", ParameterListT::Once);
            
            /* flux parameters */
            mixedbc->AddParameter(  feps, "epsilon");
            mixedbc->AddParameter(   fT0, "T0");
            mixedbc->AddParameter(falpha, "alpha");
            return mixedbc;
        }
        else /* inherited */
            return SolidElementT::NewSub(name);        
    }
    
    
    
  // accept parameter list
  void FSThermoMechT::TakeParameterList(const ParameterListT& list)
  {
      const char caller[] = "FSThermoMechT::TakeParameterList";

      /*inherited*/
      UpdatedLagrangianT::TakeParameterList(list);
 
      fTemperatureField = ElementSupport().Field("temperature");
 //     cout << "\ngot the temperture field"<<endl;
	if (!fTemperatureField) {
        ExceptionT::GeneralFail("FSThermoMechT::TakeParameterList",
                             	"could not resolve temperature field ");
	}
      /* extract information about mixed bc's */
      TakeTractionBC(list);
      
 
      
    /* Define matrix sizes */
    ModelManagerT& model = ElementSupport().ModelManager();
    int nnd = model.NumNodes();
    int nen = NumElementNodes();
    int nsd = NumSD();
    int nt = nen;	// # temperature DOFs per element
    int nme = nen * nsd;	// # of mechanical DOFs per element
    int dof = nsd + 1;	// total # of DOFs per node (mech + temp)
    int neq = nen * dof;	// total # of DOFs per element (mech + temp)

      /*TDN: define and dimension workspaces to to extrapolate Te at ip to nodes*/
      /*assume all materials within elemblock have the same internal dissipation variables*/
      ContinuumMaterialT* pmat0 = (*fMaterialList)[ElementCard(0).MaterialNumber()];
      fCurrMaterial = dynamic_cast<FSThermoMechMatT*>(pmat0);
      if (!fCurrMaterial) throw ExceptionT::kGeneralFail;
      
      fNumTe = fCurrMaterial->NumStructuralRelaxation();
      
      fGradEffectiveTemp_last.Dimension(fNumTe,NumSD());
      fGlobalMass.Dimension(nnd);
      fGlobalVal_last.Dimension(nnd, fNumTe);
      felem_mass.Dimension(nen);
      felem_val_last.Dimension(nen, fNumTe);
      fGlobalVal_last = 0.0;
      fGlobalMass = 0.0;

      fLocTe_last.Dimension(nen,fNumTe);

      /*      fGradEffectiveTemp.Dimension(fNumTe,NumSD());
      fGlobalVal.Dimension(nnd, fNumTe);
      felem_val.Dimension(nen, fNumTe);
      fGlobalVal = 0.0;
       fLocTe.Dimension(nen,fNumTe);
*/
 
    /* Dimension temperature gradient arrays */
	const int nip = NumIP();
	fTempGrad_all.Dimension(nip*nsd);
	fTempGrad_all = 0.0;
	fTempGrad_List.Dimension(nip);
	
    for (int i = 0; i < nip; ++i) {
      fTempGrad_List[i].Alias(nsd, fTempGrad_all.Pointer(i * nsd));
    }	

      fGradTe_last_List.Dimension(nip);
      //      fGradTe_List.Dimension(nip);
      for (int i = 0; i < nip; ++i) {
          fGradTe_last_List[i].Dimension(fNumTe,NumSD());
          fGradTe_last_List[i] = 0;
          //    fGradTe_List[i].Dimension(fNumTe,NumSD());
          //    fGradTe_List[i] = 0;
      }

	/* Tangent moduli for LHS */
    fAmm_mat.Dimension(nme, nme);
    fAmm_geo.Dimension(nen, nen);	// dimensions changed for Q1P0!
    fAme.Dimension(nme, nt);
    fAem.Dimension(nt, nme);
    fAee.Dimension(nt, nt);

	/* Initialize thermal mass matrix */
	fMassMatrix.Dimension(nt, nt);
      /* Initialize the matrix for thermal mechanical coupled */
      fdij.Dimension(nsd);
      fHij.Dimension(nsd);
    //  ip_temperaturerate.Dimension(1);
	/*TDN: If 2 fields are defined (temperature and displacement), then each has an integration order associated with it.  
	 *fIntegrator->Order() refers to the integration order specified for displacement while fTempIntegrator->Order() refers to 
	 integration order for temperature.  For right now, require that they be the same order.*/
	fTempIntegrator = &(fTemperatureField->Integrator().eIntegrator());

      
	int order_disp = fIntegrator->Order();
	int order_temp = fTempIntegrator->Order();
	
	if (order_disp != order_temp) 
		ExceptionT::GeneralFail(caller,
								"both temperature and displacement fields must have the same integration order");
	if (order_disp == 1)
		fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
	else
		ExceptionT::GeneralFail(caller,
								"must specify first order integrator for the temperature and displacement fields");
	
    fLHS.Dimension(neq);
    fRHS.Dimension(neq);
  //     cout << "\n got FSThermoMechT::TakeParmeterList()"<<endl;
	
  }
    
/*TDN: initialize/finalize step. Te from last converged step is extrapolated to nodes*/
void FSThermoMechT::InitStep(void)
{
    /* inherited */
    UpdatedLagrangianT::InitStep();
        
    /*extrapolate ip effective temperatues to nodes*/
    Extrapolate();
}

    
/* extract the list of material parameters */
void FSThermoMechT::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{
    const char caller[] = "FSThermoMechT::CollectMaterialInfo";
        
        /* initialize */
        mat_params.Clear();
        
        /* set materials list name */
        mat_params.SetName("fsthermomech_material");
        
        /* collected material parameters */
        int num_blocks = all_params.NumLists("fsthermomech_element_block");
        for (int i = 0; i < num_blocks; i++) {
            
            /* block information */
            const ParameterListT& block = all_params.GetList("fsthermomech_element_block", i);
            
            /* collect material parameters */
            const ParameterListT& mat_list = block.GetList(mat_params.Name());
            const ArrayT<ParameterListT>& mat = mat_list.Lists();
            mat_params.AddList(mat[0]);
    }
}
    

	  

/* form of tangent matrix */
GlobalT::SystemTypeT FSThermoMechT::TangentType(void) const
{
	/* Define LHS type based upon analysis type, i.e. static vs. dynamic */
	int order_disp = fIntegrator->Order();
 //   cout << "\n good with order_disp"<<endl;
/*temperoy disabled the following restrictions for debugging reasons "*/
//	int order_temp = fTempIntegrator->Order();
//    int order_temp = &(fTemperatureField->Integrator().eIntegrator())->Order();

//	cout << "\n good with order_temp"<<endl;
//	if (order_disp != order_temp)
//		ExceptionT::GeneralFail("FSThermoMechT::TangentType",
//								"both temperature and displacement fields must have the same integration order");
	if (order_disp == 1)
		return GlobalT::kNonSymmetric;
	else
		ExceptionT::GeneralFail("FSThermoMechT::TangentType",
								"must specify first order integrator for the temperature and displacement fields");	
}

/* finalize current step - step is solved */
void FSThermoMechT::CloseStep(void)
{
	/* inherited */
	UpdatedLagrangianT::CloseStep();
	
}
	
/* restore last converged state */
GlobalT::RelaxCodeT FSThermoMechT::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = UpdatedLagrangianT::ResetStep();
	
	return relax;
}

/* read restart information from stream */
void FSThermoMechT::ReadRestart(istream& in)
{
	/* inherited */
	UpdatedLagrangianT::ReadRestart(in);
	
}

/* write restart information from stream */
void FSThermoMechT::WriteRestart(ostream& out) const
{
	/* inherited */
	UpdatedLagrangianT::WriteRestart(out);
}

  // PROTECTED
  // construct a new material support and return a pointer
  MaterialSupportT* FSThermoMechT::NewMaterialSupport(MaterialSupportT* p) const
  {
    // allocate
    if (!p) p = new FSThermoMechSupportT(1, NumIP());
    
    // inherited initializations
    UpdatedLagrangianT::NewMaterialSupport(p);
    FSThermoMechSupportT* ps = dynamic_cast<FSThermoMechSupportT*> (p);

    if (ps != 0) {
      ps->SetTemperatureGradient(&fTempGrad_List);
      ps->SetTeGradient_last(&fGradTe_last_List);
    }
	 
    return p;
  }

  // construct materials manager and read data
  MaterialListT* FSThermoMechT::NewMaterialList(const StringT& name, int size)
  {
 	
  	// resolve number of spatial dimensions 
 /*	int nsd = -1;
	if (name == "large_strain_material_2D")
  		nsd = 2;
  	else if (name == "large_strain_material_3D")
  		nsd = 3;
  	
  	// no match 
 	if (nsd == -1) return NULL; */
      
      if (name != "fsthermomech_material")
          return NULL;
  	
    if (size > 0) {
   	  	// material support 
      	if (!fFSThermoMechSupport) {
        	fFSThermoMechSupport = dynamic_cast<FSThermoMechSupportT*> (NewMaterialSupport());
            if (!fFSThermoMechSupport) ExceptionT::GeneralFail("FSThermoMechT::NewMaterialList");
      }
        return new FSThermoMechMatListT(size, *fFSThermoMechSupport);
    }
	else
	 		return new FSThermoMechMatListT;
  }


  // form shape functions and derivatives; for ComputeOutput
  void FSThermoMechT::SetGlobalShape()
  {
    // inherited
    UpdatedLagrangianT::SetGlobalShape();

	/*temperatures and displacements set in SolidElementT*/
	/* get nodal temperature rate*/
      SetLocalU(*fLocTemperatureRate);

    /*TDN: Calculate gradient of the effective temperature, last converged time step*/
    /*get local effective temperatures*/
    /*get local element nodes*/
    iArrayT& nodes = CurrentElement().NodesU();
    /*initialize workspace*/
    //  felem_val = 0;
    felem_val_last = 0;
      
    /*Collect values from element nodes*/
    for (int i = 0; i<nodes.Length(); i++)
    {
        int node = nodes[i];
        double* pelem_val_last = felem_val_last(i);
        double* pglobal_val_last = fGlobalVal_last(node);
        //  double* pelem_val = felem_val(i);
        //  double* pglobal_val = fGlobalVal(node);
        for (int j = 0; j<fNumTe; j++)
            *pelem_val_last++ = *pglobal_val_last++;
            //  *pelem_val++ = *pglobal_val++;
    }
    /*store them in localarray for gradient calculation*/
    fLocTe_last.SetGlobal(felem_val_last);
      //    fLocTe.SetGlobal(felem_val);

    /*calculate gradients*/
	for (int i = 0; i < NumIP(); i++)
    {
		/*  temperature*/
		dArrayT& TempGrad = fTempGrad_List[i];
		fgrad.Dimension(1, NumSD());
		fShapes->GradU(*fLocTemp, fgrad, i);
        
        
        double* ptemp_grad = fgrad.Pointer();
        double* ptemp_grad_list = TempGrad.Pointer();
        for (int j = 0; j < NumSD(); j++)
            *ptemp_grad_list++ = *ptemp_grad++;

        /*TDN: effective temperature from last time step*/
        dArray2DT& TeGrad_last = fGradTe_last_List[i];
        fShapes->GradU(fLocTe_last, fGradEffectiveTemp_last);
        //  dArray2DT& TeGrad = fGradTe_List[i];
        //  fShapes->GradU(fLocTe,fGradEffectiveTemp);
        
        double* pte_grad_last = fGradEffectiveTemp_last.Pointer();
        double* pte_grad_last_list = TeGrad_last.Pointer();
        //  double* pte_grad = fGradEffectiveTemp.Pointer();
        //  double* pte_grad_list = TeGrad.Pointer();
		for (int j = 0; j < fNumTe; j++)
            for (int k = 0; k < NumSD(); k++)
                *pte_grad_last_list++ = *pte_grad_last++;
                //  *pte_grad_list++ = *pte_grad++;
	}
      
  }

  // write all current element information to the stream
  void FSThermoMechT::CurrElementInfo(ostream& out) const
  {
    // inherited
    UpdatedLagrangianT::CurrElementInfo(out);

    // write deformation gradients
    out << std::endl;
    out << "Temperature gradient field at IP:";
    out << std::endl;

    for (int i = 0; i < fTempGrad_List.Length(); ++i) {
      out << " ip: " << i + 1 << std::endl << fTempGrad_List[i] << std::endl;
    }

    out << std::endl;
  }

//  increment current element - for ComputeOutput
  bool FSThermoMechT::NextElement()
  {
    bool isThereNext = UpdatedLagrangianT::NextElement();

    if (isThereNext == true) {

      const int index = CurrentElement().MaterialNumber();

      ContinuumMaterialT* pMaterial = (*fMaterialList)[index];

     fCurrMaterial = dynamic_cast<FSThermoMechMatT*> (pMaterial);
    }

    return isThereNext;
  }

  // Initialize local arrays
  void FSThermoMechT::SetLocalArrays()
  {
  //   cout << "\n I have reached FSThermoMechT::SetLocalArrays()"<<endl;
	/*Check that temperature field is defined*/
    fTemperatureField = ElementSupport().Field("temperature");

    if (0 == fTemperatureField) {

      std::cout << std::endl;
      std::cout << "FSThermoMechT::SetLocalArrays: ";
      std::cout << "temprature field not found.";
      std::cout << std::endl;

      throw ExceptionT::kGeneralFail;
    }
 
	  int nen = NumElementNodes();

      /*fLocTemp and fLocTemp_last registered with temperature field in SolidElementT*/
      if (fTemperatureField) {
	  fLocTemperatureRate = new LocalArrayT(LocalArrayT::kVel, nen, fTemperatureField->NumDOF());
	  fTemperatureField->RegisterLocal(*fLocTemperatureRate);

      }

	  /* Inherited*/
    UpdatedLagrangianT::SetLocalArrays();
}

  void FSThermoMechT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
      AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
  {
    //  cout<<"\n fEqnos.Major is "<<fEqnos[0].MajorDim();
    //  cout<<"\n fEqnos.Minor is "<<fEqnos[0].MinorDim();
    //   cout<<"\n fEqnos[i] is "<<fEqnos[0];
      for (int i = 0; i < fEqnos.Length(); ++i) {

      const int ndf = NumDOF();	// displacement DOFs e.g. ndf = nsd
      const int nen = fConnectivities[i]->MinorDim();

	  /*offset in equation number for multiple fields. DOFs are arranged {u_for_all_nodes, T_for_all_nodes}*/
      const int offset = ndf * nen;  
     // cout<<"\n fConnectivities[i]->MinorDim() is :"<<fConnectivities[i]->MinorDim();
      //  cout<<"\n fEqnos.Length() is "<<fEqnos.Length();
      // cout<<"\n offset is :"<<offset;
 
    /*assign equation numbers for temperature*/
      fTemperatureField->SetLocalEqnos(*fConnectivities[i],
          fEqnos[i], offset);
     //  cout<<"\n fEqnos[i] is "<<fEqnos[0];
    }
     ElementBaseT::Equations(eq_1, eq_2);
    /*the following is added to have the mixed BC */
     fBCEqnos.Dimension(fBCFaces);
     // Field().SetLocalEqnos(fBCFaces, fBCEqnos);
     fTemperatureField->SetLocalEqnos(fBCFaces, fBCEqnos);
    //  cout<<"\n fBCFaces is:" <<fBCFaces;
   // cout<<"\n fBCEqnos is:" <<fBCEqnos;
    
  }

  //
  int FSThermoMechT::TotalNumDOF() const
  {
 	int mechdof = 3;
 	int temperaturedof = 1;
    return (mechdof+temperaturedof);
  }

  const dArrayT&
  FSThermoMechT::TemperatureGradient() const
  {
    return fTempGrad_List[CurrIP()];
  }

  //
  const dArrayT&
  FSThermoMechT::TemperatureGradient(int ip) const
  {
    return fTempGrad_List[ip];
  }

/* accumulate the residual force on the specified node */
void FSThermoMechT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
	/*implement body force term of equilibrium equations,	*
	 *viscoplastic heating term, thermoelastic heat terms,	*
	 *\dot{Tf} term.
     /* quick exit */
	bool hasnode = false;
	for (int i=0; i < fBlockData.Length() && !hasnode; i++)
		if (fConnectivities[i]->HasValue(node)) hasnode = true;
	if (!hasnode) return;
    // cout<<"\n  fBlockData.Length()  is:" << fBlockData.Length() ;
	/* set components and weights */
	double constMa = 0.0;
	double constKd = 0.0;
    
	/* components dicated by the algorithm */
	int formMa = fIntegrator->FormMa(constMa);
	int formKd = fIntegrator->FormKd(constKd);
    
	/* body forces */
	int formBody = 0;
	if (fMassType != kNoMass &&
        (fBodySchedule && fBody.Magnitude() > kSmall))
	{
		formBody = 1;
		if (!formMa) constMa = 1.0; /* override */
	}
    /* override controller */
	if (fMassType == kNoMass) formMa = 0;
    
	/* temp for nodal force */
	dArrayT nodalforce;
    
    bool axisymmetric = Axisymmetric();
	Top();
    while (NextElement())
	{
        int nodeposition;
		const iArrayT& nodes_u = CurrentElement().NodesU();
		if (nodes_u.HasValue(node, nodeposition))
		{
			/* initialize */
			fRHS = 0.0;
            /* global shape function values */
			SetGlobalShape();
            
			/* internal force contribution */
			if (formKd) FormKd(constKd);
            /* inertia forces */
			/*if (formMa)
             {
             SetLocalU(fLocAcc);
             FormMa(fMassType, constMa*fCurrMaterial->Density(), axisymmetric, &fLocAcc, NULL, NULL);
             }*/
            /* mechanical and thermal reaction forces */
			double mr1, mr2, mr3, er1;
			dArrayT react(3);
			
			/* loop over nodes (double-noding OK) */
			int dex = 0;
			int dex2 = 0;
			int whichdof;	// new
			for (int i = 0; i < nodes_u.Length(); i++)
			{
				if (nodes_u[i] == node)
				{
					/* not my field - thermal */
					if (&field != &(Field()))
					{
						er1 = fRHS[dex2+3*NumElementNodes()];
						react[0] = er1;
						whichdof = 1;
					}
					else	// otherwise do mechanical
					{
						mr1 = fRHS[dex];
						mr2 = fRHS[dex+1];
						mr3 = fRHS[dex+2];
						react[0] = mr1;
						react[1] = mr2;
						react[2] = mr3;
						whichdof = NumDOF();
					}
					
					/* components for node - mechanical + electrical DOFs */
                    //					nodalforce.Set(TotalNumDOF(), react.Pointer(0));
					nodalforce.Set(whichdof, react.Pointer(0));
                    
					/* accumulate */
					force += nodalforce;
				}
                
				dex += NumDOF();
				dex2 += 1;
            }
            
        }
    }
        
    }
    
/* construct the effective mass matrix */
    void FSThermoMechT::LHSDriver(GlobalT::SystemTypeT sys_type)
    {
        //SolidElementT::LHSDriver(sys_type);
        /* Basically there is a bug in the diffusion problem if we have the initial diffusion condition */
        ContinuumElementT::LHSDriver(sys_type);
        double constM = 0.0;
        double constK = 0.0;
        
        int formM = fIntegrator->FormM(constM);
        int formK = fIntegrator->FormK(constK);
        bool axisymmetric = Axisymmetric();
        Top();
        while (NextElement())
            if (CurrentElement().Flag() != ElementCardT::kOFF)
            {
                double constKe = constK;
                double constMe = constM;
                /* initialize */
                fLHS = 0.0;
                
                /* set shape function derivatives */
                SetGlobalShape();
                if (fabs(constMe) > kSmall)
                {
                    if (!fCurrMaterial->HasChangingDensity())
                        FormMass(fMassType, constMe*(fCurrMaterial->Density()), axisymmetric, NULL);
                    else
                    {
                        /* collect densities */
                        fShapes->TopIP();
                        while (fShapes->NextIP())
                            fDensity[fShapes->CurrIP()] = fCurrMaterial->Density();
                        
                        FormMass(fMassType, constMe, axisymmetric, fDensity.Pointer());
                    }
                    //           cout<<"\n add FormMass in SolidElementT is"<<endl;
                }
                /* element stiffness */
               // if (fabs(constKe) > kSmall)
                    FormStiffness(constKe);
                /* add to global equations */
                AssembleLHS();
                //      cout<<" \n I have reached the final SolidElementT::ElementLHSDriver"<<endl;
            }

          TractionBC_LHS();
      

    }
    
    void FSThermoMechT::RHSDriver(void)
    {
        /* inherited */
        SolidElementT::RHSDriver();
    	
        TractionBC_RHS();
   //   cout<<"\n fRHS2:"<<fRHS;
        
    }
/* calculate the LHS of residual, or element stiffness matrix */
  void FSThermoMechT::FormStiffness(double constK)
    {
	  /*implement the element K matrix for coupled problem*/
      /* Time integrator info for dynamic problems */
      int order = fIntegrator->Order();
      
      /* Matrix format - depends upon time integration order */
      dMatrixT::SymmetryFlagT format = (fLHS.Format()
                                        == ElementMatrixT::kNonSymmetric)
      ? dMatrixT::kWhole
      : dMatrixT::kUpperOnly;
      
      int nen = NumElementNodes();
      int nsd = NumSD();
      int nt = nen;	// # temperature DOFs per element
      int nme = nen * nsd;	// # of mechanical DOFs per element
      int dof = nsd + 1;	// total # of DOFs per node (mech + temp)
      int neq = nen * dof;	// total # of DOFs per element (mech + temp)
            
      fAmm_mat = 0.0;
      fAmm_geo = 0.0;
      fAme = 0.0;
      fAem = 0.0;
      fAee = 0.0;
      
      /* integration */
      const double* Det    = fCurrShapes->IPDets();
      const double* Weight = fCurrShapes->IPWeights();
      
      
      fShapes->TopIP();
      while (fShapes->NextIP() )
      {
          double scale2 = (*Det++)*(*Weight++);
          double scale=scale2*constK;
         double scale1=scale*fCurrMaterial->b1();
          scale2 *=fCurrMaterial->Capacity();

          
          const double* Na = fShapes->IPShapeU();
          fCurrShapes->GradNa(fGradNa);
          Set_B(fCurrShapes->Derivatives_U(), fB);
          
          /* S T R E S S   S T I F F N E S S */
          /* compute Cauchy stress */
          const dSymMatrixT& cauchy = fCurrMaterial->s_ij();
       //   cout<<"\n s_ij in FSThermoMechT is "<<cauchy<<endl;
          cauchy.ToMatrix(fCauchyStress);
          fCauchyStress *= scale;
          /* using the stress symmetry */
          fAmm_geo.MultQTBQ(fGradNa, fCauchyStress, format, dMatrixT::kAccumulate);
          
          /* M A T E R I A L   S T I F F N E S S */
          /* get D matrix */
          fD.SetToScaled(scale, fCurrMaterial->c_ijkl());
          /* accumulate */
          fAmm_mat.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);
          
          int a,b,i;
          double temp3;
          dMatrixT temp4(nen,nsd);
          dArrayT temp5(nen);
        //  dArrayT temp6(nen);
       //   cout<<"\n b1 is:"<<fCurrMaterial->b1()<<endl;
          const dMatrixT kij=fCurrMaterial->Conductivity();
          const dSymMatrixT& dij = fCurrMaterial->d_ij();
     //    cout<<"\n d_ij in FSThermoMechT is "<<dij<<endl;
          dij.ToMatrix(fdij);
          const dSymMatrixT& hij = fCurrMaterial->h_ij();
     //     cout<<"\n h_ij in FSThermoMechT is "<<hij<<endl;
          hij.ToMatrix(fHij);
          const dArrayT& GradT=TemperatureGradient();
        //  scale2=scale/constK*fCurrMaterial->Capacity();
     //   cout<<"\n capacity is "<<fCurrMaterial->Capacity()<<endl;
        // cout<<"\n scale2 is "<<scale2<<endl;
          // cout<<"\n scale is "<<scale<<endl;
          // cout<<"\n scale1 is "<<scale1<<endl;
           /*fAee*/
          for (a=0;a<nen;a++)
          {
             for(b=0;b<nen;b++)         
             {
                 fAee(a,b)+=Na[a]*Na[b]*scale1;
           /* if (a==0&b==0)
                 {
               //  cout<<"\n Na[0] is:"<<Na[0]<<endl;
                 cout<<"\n scale1 is:"<<scale1<<endl;
                 cout<<"\n Na[0]*Na[0]*scale1 is:"<<Na[a]*Na[b]*scale1<<endl;
                 cout<<"the first fAee[0,0] is:"<<fAee(0,0)<<endl;
                 } */
                 temp3=fGradNa(0,a)*fGradNa(0,b)*kij(0,0)+fGradNa(1,a)*fGradNa(1,b)*kij(1,1)+fGradNa(2,a)*fGradNa(2,b)*kij(2,2);
                 fAee(a,b)+=temp3*scale;
              /*  if (a==0&b==0)
                  {
                      cout<<"\n temp3  is "<<temp3<<endl;
                      cout<<"\n scale is:"<<scale<<endl;
                cout<<"\n temp3*scale is:"<<temp3*scale<<endl;
                 cout<<"the second fAee[0,0] is:"<<fAee(0,0)<<endl;
                  }  */
                 /*this is adding FormMass*/
                 fAee(a,b)+= Na[a]*Na[b]*scale2;
              /*  if (a==0&b==0)
                 {
                     cout<<"\n scale2  is "<<scale2<<endl;
                     cout<<"\n Na[0] is:"<<Na[0]<<endl;
                     cout<<"\n Na[0]*Na[0]*scale2 is:"<<Na[a]*Na[b]*scale2<<endl;
                 cout<<"the third fAee[0,0] is:"<<fAee(0,0)<<endl;
                 } */
             }
            }
    

      //  cout<<"\n the first fAee is "<<fAee<<endl;
          /*fAme*/
           for (a=0;a<nen;a++)
           {
                for(b=0;b<nen;b++)
               {
               for(i=0;i<nsd;i++)
               {
                temp4(a,i) =fGradNa(0,a)*fdij(i,0)+fGradNa(1,a)*fdij(i,1)+fGradNa(2,a)*fdij(i,2);
                fAme(3*a+i,b)+=temp4(a,i)*Na[b]*scale;
               }
               }
           }

           /*fAem*/
          for (a=0;a<nen;a++)
          {
              for(b=0;b<nen;b++)
              {
                  temp5[b]=GradT[0]*fGradNa(0,b)*kij(0,0)+GradT[1]*fGradNa(1,b)*kij(1,1)+GradT[2]*fGradNa(2,b)*kij(2,2);
                  temp3 =fGradNa(0,a)*fGradNa(0,b)*kij(0,0)+fGradNa(1,a)*fGradNa(1,b)*kij(1,1)+fGradNa(2,a)*fGradNa(2,b)*kij(2,2);
                  for(i=0;i<nsd;i++)
                  {   temp4(b,i) =fGradNa(0,b)*fHij(i,0)+fGradNa(1,b)*fHij(i,1)+fGradNa(2,b)*fHij(i,2);
                      fAem(a,3*b+i)+=temp4(b,i)*Na[a]*scale;
                      fAem(a,3*b+i)+=-fGradNa(i,a)*temp5[b]*scale;
                      fAem(a,3*b+i)+=-GradT[i]*temp3*scale;
                  }
                }
            }
      //    cout<<"\n GradT is "<<GradT<<endl;
      //    cout<<"\n kij is :"<<kij<<endl;
      //    cout<<"\n temp5 is :"<<temp5<<endl;
       //   cout<<"\n temp3 is :"<<temp3<<endl;
       //   cout<<"\n temp4 is :"<<temp4<<endl;

      }
        
      
      /* stress stiffness into fLHS (i.e. fAmm_mat) */
      fAmm_mat.Expand(fAmm_geo, NumDOF(), dMatrixT::kAccumulate);
      
      /* Assemble into fLHS, or element stiffness matrix */
      fLHS.AddBlock(0, 0, fAmm_mat);
    //  fAme = 0.0;
    //   fAem = 0.0;
      fLHS.AddBlock(0, fAmm_mat.Cols(), fAme);
      fLHS.AddBlock(fAmm_mat.Rows(),0,fAem);
        //fAee *=0.5;
      fLHS.AddBlock(fAmm_mat.Rows(), fAmm_mat.Cols(), fAee);
      //  cout<<"\n left hand side is :"<<fLHS<<endl;
  //   cout<<"\n fAee is "<<fAee<<endl;
 //    cout<<"\n fAem is "<<fAem<<endl;
     //  cout<<"\n fAme is "<<fAme<<endl;

  }
/* Compute RHS, or residual of element equations */
  void FSThermoMechT::FormKd(double constK)
  {  	
	  /*implement residual vector R for coupled problem*/
      /* element preliminaries */
      const int nsd = NumSD();
      const int nen = NumElementNodes();
      
      /* Define mechanical and electrical residuals */
      dArrayT Rtotal((nsd+1)*nen);
      Rtotal = 0.0;
      dArrayT Rmech(nen*nsd);
      Rmech = 0.0;
      dArrayT Relec(nen);
      Relec = 0.0; 
      
      const double* Det    = fCurrShapes->IPDets();
      const double* Weight = fCurrShapes->IPWeights();
      
      /* current element number */
      int elem = CurrElementNumber();
      

    // cout<< "\n  fLocTemperatureRate is at all the nodal: "<<*fLocTemperatureRate<<endl;
    //   cout<< "\n  fLocTemperature is at all the nodal: "<<*fLocTemp<<endl;
      fCurrShapes->TopIP();
      dArrayT tempqi;
      dArrayT nodal_temperaturerate;
      nodal_temperaturerate.Dimension(nen);
      nodal_temperaturerate=*fLocTemperatureRate;
      while (fCurrShapes->NextIP() )
      {
          const double* Na = fShapes->IPShapeU();
          fCurrShapes->GradNa(fGradNa);
          Set_B(fCurrShapes->Derivatives_U(), fB);
          
          /* B^T * Cauchy stress */
          const dSymMatrixT& cauchy = fCurrMaterial->s_ij();
          fB.MultTx(cauchy, fNEEvec);
          
          /* double scale factor */
          double scale = constK*(*Det++)*(*Weight++);
          double scale1=scale*fCurrMaterial->Capacity();
        //  cout<<"\n constK in FromKd is :"<<constK<<endl; constK is -1 in the above function, this may be the reason to have the error.
          /* accumulate - use Rmech instead of fRHS */
          Rmech.AddScaled(scale, fNEEvec);
          
     //     dArrayT& GradT=TemperatureGradient();
         const dArrayT& qi=fCurrMaterial->q_i();
     /*     if(CurrIP()==0)
     cout<<"\n the heat flux in the first integration point is :"<<qi<<endl;*/
    //      if(CurrIP()==1)
     //         cout<<"\n the heat flux in the second integration point is :"<<qi<<endl;
          int a,b;
          double ip_temperaturerate;
      // cout<<"\n heat residule is :"<<fCurrMaterial->heatres()<<endl;
          
          for (a=0;a<nen;a++)
          {
              Relec[a]+= scale*Na[a]*fCurrMaterial->heatres();
              Relec[a]+=(fGradNa(0,a)*qi[0]+fGradNa(1,a)*qi[1]+fGradNa(2,a)*qi[2])*scale;
              ip_temperaturerate=0.0;
              for (b=0;b<nen;b++)
              {
                  ip_temperaturerate +=Na[b]*(nodal_temperaturerate[b]);
              }
             // fFSThermoMechSupport->Interpolate(*fLocTemperatureRate, ip_temperaturerate);
            //  cout<<"\n flocTempraturerate at all the nodal is :"<<*fLocTemperatureRate<<endl;
           //   cout<<"\n flocTempraturerate at all the nodal by nodaltemperaturerate is :"<<nodal_temperaturerate<<endl;
          

              Relec[a]+=Na[a]*ip_temperaturerate*scale1;
            /*  if (a==0)
              {
                  cout<<"\n flocTempraturerate at the integration point is :"<<ip_temperaturerate<<endl;
                  cout<<"\n Na[0]*ip_temperaturerate*scale at the current intergration point is:"<<Na[a]*ip_temperaturerate*scale1;
                  cout<<"\n the accumlative residual is:"<<Relec[0];
              } */
          }
      //    cout<<"\n flocTempraturerate at all the nodal is :"<<*fLocTemperatureRate<<endl;
          /* temperature resudial in current configuration */
       
      }
     
      Relec *= 1.0;
   //   cout<<"\n Relec is:"<<Relec<<endl;
      Rtotal.CopyIn(0, Rmech);
      Rtotal.CopyIn(Rmech.Length(), Relec);
      fRHS += Rtotal; 	
  //   cout<<"\n thermal residual is:"<<Relec<<endl;
   //   cout<<"\n check the mechanical residual:"<<Rmech<<endl;
      /* volume averaged */
     // p_bar /= fElementVolume[CurrElementNumber()];
  }

/* Dummy mass matrix for dynamic calculations */
void FSThermoMechT::FormMass(MassTypeT mass_type, double constM, bool axisymmetric, const double* ip_weight)
{
  //  int order = fIntegrator->Order();
    
    /*The FormMass is not inclued in the LHSDriver in SolidElements if the integrator is first order */
    /*So I will add the formmass directly from the FormKd */
    /*
    int nen = NumElementNodes();
    int nt = nen;	// # temperature DOFs per element
    fMassMatrix=0.0;
    
    const double* Det    = fCurrShapes->IPDets();
    const double* Weight = fCurrShapes->IPWeights();
    
    
    fShapes->TopIP();
    double temp1;
    while (fShapes->NextIP() )
    {
        double scale = (*Det++)*(*Weight++)*fCurrMaterial->Capacity();
       const double* Na = fShapes->IPShapeU();
        int a,b;
        for (a=0;a<nt;a++)
        {
            for(b=0;b<nt;b++)
                temp1 =Na[a]*Na[b]*scale;
            fMassMatrix(a,b)+=temp1;
            }
        }
        
        
   
    fLHS.AddBlock(fAmm_mat.Rows(), fAmm_mat.Cols(), fMassMatrix);*/
    fLHS +=0.0;
}
    

/* Calculate inertial force for dynamic calculations */
void FSThermoMechT::FormMa(MassTypeT mass_type, double constM, bool axisymmetric, 
	const LocalArrayT* nodal_values, const dArray2DT* ip_values, const double* ip_weight)
{
    /*we may not be able to add deltac*dotT here because the intergrotor order is 1 and FormMa may not add the  right hand side, so I suggest to add this part at FormKd */
  /*  const char caller[] = "FSThermoMechT::FormMa";
    
    const int nsd = NumSD();
    const int nen = NumElementNodes();
	dArrayT Relec(nen);
	dArrayT Rtotal((nsd+1)*nen);
    dArrayT Rmech(nen*nsd);
    Relec = 0.0;
	Rmech = 0.0;
	Rtotal = 0.0;
    const double* Det    = fCurrShapes->IPDets();
    const double* Weight = fCurrShapes->IPWeights();
    dArrayT ip_temperaturerate;
     while (fCurrShapes->NextIP() )
     {
          const double* Na = fShapes->IPShapeU();
         double scale = constM*(*Det++)*(*Weight++);
         int a;
         fFSThermoMechSupport->Interpolate(*fLocTemperatureRate, ip_temperaturerate);
         for (a=0;a<nen;a++)
         {
           Relec[a]=scale*Na[a]*ip_temperaturerate[0]*fCurrMaterial->Capacity();
         }
     }
    Rtotal.CopyIn(Rmech.Length(), Relec);
	fRHS += Rtotal; */
    fRHS +=0.0;
    
}


/***********************************************************************
 * Private
 ***********************************************************************/

    void FSThermoMechT::TakeTractionBC(const ParameterListT& list)
    {
        const char caller[] = "FSThermoMechT::TakeTractionBC";
        
        /* quick exit */
        const ParameterListT* mixedbc = list.List("mixed_bc");
        if (!mixedbc)
            return;
        
        /* extract BC parameters */
        feps   = mixedbc->GetParameter("epsilon");
        fT0    = mixedbc->GetParameter("T0");
        falpha = mixedbc->GetParameter("alpha");
        
        const ParameterListT& ss_ID_list = mixedbc->GetList("side_set_ID_list");
        int num_sides = ss_ID_list.NumLists("String");
      //   cout<<"\n num_sides is:"<<num_sides;
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
        //      cout<<"\n num_faces is:"<<num_faces;
            /* element topology */
            iArrayT nodes_on_faces(fShapes->NumFacets());
            fShapes->NumNodesOnFacets(nodes_on_faces);
            int min, max;
            nodes_on_faces.MinMax(min, max);
            if (min != max) ExceptionT::GeneralFail(caller, "all faces must have same shape");
            
            /* collect nodes on faces */
            int face_num = 0;
            fBCFaces.Dimension(num_faces, nodes_on_faces[0]);
           //  cout<<"\n nodes_on_faces is:"<<nodes_on_faces[0];
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
       //         cout<<"\n fBCFaces at the second place is:"<<fBCFaces;
            }		
        }
    }
    
    void FSThermoMechT::TractionBC_RHS(void)
    {
        /* quick exit */
        if (fBCFaces.MajorDim() == 0) return;
        
        /* dimensions */
        int nsd = NumSD();
        //choose mechanical dof +temperature dof or just temperature dof?
       // int ndof = NumDOF();
        int ndof =1;
        int nfn = fBCFaces.MinorDim();
        
        /* force vector */
        dArrayT rhs(nfn*ndof);
		
        /* local coordinates */
        LocalArrayT coords(LocalArrayT::kInitCoords, nfn, nsd);
        ElementSupport().RegisterCoordinates(coords);
		
        /* nodal field values */
        //RX modify this by assuming field is *fLocalTemp, need to registerlocal?
         LocalArrayT field(LocalArrayT::kDisp, nfn, ndof);
        //Field().RegisterLocal(field);
         fTemperatureField->RegisterLocal(field);
        
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
		//	cout<<"\n nodes is:"<<nodes;
         //   cout<<"\n eqnos is:"<<eqnos;
            /* local values */
            coords.SetLocal(nodes);
            field.SetLocal(nodes);
            //RX modify this
           // field.SetLocal(nodes);
           
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
       //     cout<<"\n rhs is "<<rhs;
            /* assemble */
            ElementSupport().AssembleRHS(Group(), rhs, eqnos);
        }
    }
    
    /* compute contribution to LHS from BC's */
    void FSThermoMechT::TractionBC_LHS(void)
    {
        /* quick exit */
        if (fBCFaces.MajorDim() == 0) return;
        
        /* dimensions */
        int nsd = NumSD();
        //RX modify this 
      //  int ndof = NumDOF();
        int ndof=1;
        int nfn = fBCFaces.MinorDim();
        
        /* stiffness matrix */
        ElementMatrixT lhs(nfn*ndof, ElementMatrixT::kNonSymmetric);
		
        /* local coordinates */
        LocalArrayT coords(LocalArrayT::kInitCoords, nfn, nsd);
        ElementSupport().RegisterCoordinates(coords);
		
        /* nodal field values */
        LocalArrayT field(LocalArrayT::kDisp, nfn, ndof);
        // Field().RegisterLocal(field);
        fTemperatureField->RegisterLocal(field);

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
        double constK = 0.0;
        int formK = fIntegrator->FormK(constK);
        for (int i = 0; i < fBCFaces.MajorDim(); i++)
        {
            /* face info */
            fBCFaces.RowAlias(i, nodes);
            fBCEqnos.RowAlias(i, eqnos);
		//	cout<<"\n nodes is:"<<nodes;
        //    cout<<"\n eqnos is:"<<eqnos;
            /* local values */
            coords.SetLocal(nodes);
           field.SetLocal(nodes);
            
            /* all ip field values: (nip x ndof) */
            surf_shape.Interpolate(field, ip_field);
      //      cout<<"\n ip_field is:"<<ip_field;
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
                lhs.Outer(Na, Na, jw*d_qn*constK, dMatrixT::kAccumulate);			
            }
        //   cout<<"\n lhs is "<<lhs;
      //      cout<<"\n I reached this point"<<endl;

            /* assemble */
            ElementSupport().AssembleLHS(Group(), lhs, eqnos);
        }
    }

    /*TDN: Extrapolate effective temepratures from last converged time step at IP to nodes*/
    void FSThermoMechT::Extrapolate(void)
    {
        const char caller[] = "FSThermoMechT::Extrapolate";
        
        int nen = NumElementNodes();
        
        Top();
        while (NextElement())
        {
            /*initialize element matrices*/
            felem_mass = 0.0;
            
            if (CurrentElement().Flag() != ElementCardT::kOFF)
            {
                ContinuumMaterialT* pmat0 = (*fMaterialList)[ElementCard(0).MaterialNumber()];
                fCurrMaterial = dynamic_cast<FSThermoMechMatT*>(pmat0);
                if (!fCurrMaterial) throw ExceptionT::kGeneralFail;
                UpdatedLagrangianT::SetGlobalShape();
                
                //  felem_val = 0.0;
                felem_val_last = 0.0;
                
                const double* jac = fShapes->IPDets();;
                const double* weight = fShapes->IPWeights();
                fShapes->TopIP();
                while(fShapes->NextIP())
                {
                    const dArrayT& Te_last = fCurrMaterial->Effective_Temperature_last();
                    //  const dArrayT& Te = fCurrMaterial->Effective_Temperature();
                    const double* pQbU = fShapes->IPShapeU();
                    
                    for (int i=0; i<nen; i++)
                    {
                        const double* pQaU = fShapes->IPShapeU();
                        
                        /*lumped element mass matrix*/
                        for (int j = 0; j<NumElementNodes(); j++)
                            felem_mass[i] += (*pQaU++)*(*pQbU)*(*jac)*(*weight);
                        
                        /*element forcing function*/
                        for (int cnt = 0; cnt < fNumTe; cnt++)
                            felem_val_last(i,cnt) += (*pQbU)*(*jac)*(*weight)*Te_last[cnt];
                        //  felem_val(i,cnt) += (*pQbU)*(*jac)*(*weight)*Te[cnt];

                        pQbU++;
                    }
                    weight++;
                    jac++;
                } 
                AssembleArray(felem_mass, fGlobalMass, CurrentElement().NodesX());
                AssembleArray2D(felem_val_last, fGlobalVal_last, CurrentElement().NodesX());
//                AssembleArray2D(felem_val, fGlobalVal, CurrentElement().NodesX());
            }
        }
        
        for (int i = 0; i< fGlobalMass.Length(); i++)
            for (int j = 0; j< fNumTe; j++)
                fGlobalVal_last(i,j) /= fGlobalMass[i];
                //  fGlobalVal(i,j) /= fGlobalMass[i];
    }

    /*TDN: utility function to assemble lumped mass values*/
    void FSThermoMechT::AssembleArray(const dArrayT& elem_val, dArrayT& global_val, const iArrayT& elem_nodes)
    {
        int nen = elem_nodes.Length();
        int nsd = NumSD();
        int index;
        int numval = elem_val.Length();
        
        if (numval != nen*nsd && numval != nen) throw ExceptionT::kGeneralFail;
        
        for (int i = 0; i< nen; i++)
        {
            if (numval == nen*nsd)
            {
                index=elem_nodes[i]*nsd;
                for (int j = 0; j<nsd; j++)
                    global_val[index+j] += elem_val[i*nsd+j];
            }
            else
            {
                index = elem_nodes[i];
                global_val[index] += elem_val[i];
            }
        }
    }
    
    /*TDN: utility function to assemble ip values*/
    void FSThermoMechT::AssembleArray2D(const dArray2DT& elem_val, dArray2DT& global_val, const iArrayT& elem_nodes)
    {
        int nen = elem_nodes.Length();
        int numval = elem_val.MinorDim();
        
        if (global_val.MinorDim() != numval) throw ExceptionT::kGeneralFail;
        if (elem_val.MajorDim() != nen) throw ExceptionT::kGeneralFail;
        
        int index;
        for(int i = 0; i < nen; i++)
        {
            index=elem_nodes[i];
            for (int j = 0; j < numval; j++)
                global_val(index,j) += elem_val(i,j);
        }
    }
    


// extrapolate from integration points and compute output nodal/element values
/*void FSThermoMechT::ComputeOutput(const iArrayT& n_codes,
      dArray2DT& n_values, const iArrayT& e_codes, dArray2DT& e_values)
{
    SolidElementT::ComputeOutput(n_codes,n_values,  e_codes, e_values);
   }
*/
    
} // namespace Tahoe
