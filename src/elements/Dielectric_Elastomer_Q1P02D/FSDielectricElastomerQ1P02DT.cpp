#include "FSDielectricElastomerQ1P02DT.h"
#include "FSDEMatSupportQ1P02DT.h"
#include "FSDEMatQ1P02DT.h"
#include "ParameterContainerT.h"
#include "OutputSetT.h"
#include "ShapeFunctionT.h"
#include "eIntegratorT.h"
#include "SolidMaterialT.h"
#include "SolidMatListT.h"
#include "FSSolidMatT.h"
#include "FSMatSupportT.h"

#include <iostream>
#include <fstream>

/* ISSUES:
1.  AddNodalForce - calls SetGlobalShape, which is modified in the Q1P0 formulation.
	Could be a problem for reaction force output accuracy
2.  ComputeOutput - all based on TL (reference shape functions), not ULagrangian
3.  Integration factor for K terms multiplied by beta*dt^2 - correct!
*/

// materials lists (3D only)
#include "FSSolidMatList2DT.h"
using namespace std;

namespace Tahoe {

  FSDielectricElastomerQ1P02DT::FSDielectricElastomerQ1P02DT(
      const ElementSupportT& support) :
    FiniteStrainT(support), fFSDEMatSupportQ1P02D(0), fCurrMaterial(0),
    fLocScalarPotential(LocalArrayT::kESP), fElectricScalarPotentialField(0),
    fLocCurrCoords(LocalArrayT::kCurrCoords)
  {
    SetName("dielectric_elastomer_Q1P02D");
  }


/* Destructor */
  FSDielectricElastomerQ1P02DT::~FSDielectricElastomerQ1P02DT()
  {
  	  delete fCurrShapes;
  	  fCurrShapes = NULL;
	  if (0 != fFSDEMatSupportQ1P02D) delete fFSDEMatSupportQ1P02D;
  }

  // specify parameters needed by the interface  
  void FSDielectricElastomerQ1P02DT::DefineParameters(ParameterListT& list) const
  {
    // inherited
    FiniteStrainT::DefineParameters(list);

    // additional fields
    list.AddParameter(ParameterT::Word, "electric_field_name");  
    
	/* remove option to store shape functions - for Q1P0 - WHAT IS WRONG?  */
//	list.RemoveParameter("store_shapefunctions");   
  }

  // accept parameter list
  void FSDielectricElastomerQ1P02DT::TakeParameterList(const ParameterListT& list)
  {
    // inherited
    FiniteStrainT::TakeParameterList(list);

    // get electric scalar potential field
    // for now use same integration and interpolation schemes as primary field
    const StringT& electric_field_name = list.GetParameter(
        "electric_field_name");

    fElectricScalarPotentialField = ElementSupport().Field(electric_field_name);
    if (!fElectricScalarPotentialField) {
      ExceptionT::GeneralFail("FSDielectricElastomerQ1P02DT::TakeParameterList",
          "could not resolve \"%s\" field", electric_field_name.Pointer());
    }

    /* Define matrix sizes */
    int nen = NumElementNodes();
    int nsd = NumSD();
    int nel = nen;	// # electrical DOFs per element
    int nme = nen * nsd;	// # of mechanical DOFs per element
    int dof = nsd + 1;	// total # of DOFs per node (mech + elec)
    int neq = nen * dof;	// total # of DOFs per element (mech + elec)

	/* Dimension electric field arrays */
	const int nip = NumIP();
	fE_all.Dimension(nip*nsd);
	fE_all = 0.0;	// testing HSP
	fE_List.Dimension(nip);
	
	/* what does this do? */
    for (int i = 0; i < nip; ++i) {
      fE_List[i].Alias(nsd, fE_all.Pointer(i * nsd));
    }	

	/* Tangent moduli for LHS */
    fAmm_mat.Dimension(nme, nme);
    fAmm_geo.Dimension(nen, nen);	// dimensions changed for Q1P0!
    fAme.Dimension(nme, nel);
    fAem.Dimension(nel, nme);
    fAee.Dimension(nel, nel);

	/* Initialize mass matrix */
	fMassMatrix.Dimension(nme, nme);

	/* allocate workspace - from UpdatedLagrangianT.cpp */
	fCauchyStress.Dimension(nsd);
	fGradNa.Dimension(nsd, nen);
	fStressStiff.Dimension(nen);

	/* Define LHS type based upon analysis type, i.e. static vs. dynamic */
	int order = fIntegrator->Order();
	if (order == 2)
		fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
	else
		fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);
	
    fLHS.Dimension(neq);
    fRHS.Dimension(neq);
    
	/* Q1P0 STUFF */
	const char caller[] = "FSDielectricElastomerQ1P02DT::TakeParameterList";	
	
	/* check geometry code and number of element nodes -> Q1 */
	if (GeometryCode() == GeometryT::kQuadrilateral) {
		if (NumElementNodes() != 4) 
			ExceptionT::BadInputValue(caller, "expecting 4 node quad: %d", NumElementNodes());
	}
	else if (GeometryCode() == GeometryT::kHexahedron) {
		if (NumElementNodes() != 8) 
			ExceptionT::BadInputValue(caller, "expecting 8 node hex: %d", NumElementNodes());
	}
	else
		ExceptionT::BadInputValue(caller, "expecting hex or quad geometry: %d", GeometryCode());
	
	/* need to store last deformed element volume */
	fElementVolume.Dimension(NumElements());	
	fElementVolume = 0.0;
	fElementVolume_last.Dimension(NumElements());
	fElementVolume_last = 0.0;
	
	/* element pressure */
	fPressure.Dimension(NumElements());
	fPressure = 0.0;
	
	/* determinant of the deformation gradient */
	fJacobian.Dimension(NumIP());
	fJacobian = 1.0;
	
	/* dimension work space */
	fMeanGradient.Dimension(NumSD(), NumElementNodes());
	fNEEmat.Dimension(nme);
	fdiff_b.Dimension(NumSD(), NumElementNodes());
	fb_bar.Dimension(NumSD(), NumElementNodes());
	fb_sig.Dimension(NumSD(), NumElementNodes());
	fF_tmp.Dimension(NumSD());

	/* need to initialize previous volume */
	Top();

	while (NextElement())
	{
		/* inherited - computes gradients and standard 
		 * deformation gradients */
		FiniteStrainT::SetGlobalShape();

		/* compute mean of shape function gradients */
		double H; /* reference volume */
		double& v = fElementVolume_last[CurrElementNumber()];

		SetMeanGradient(fMeanGradient, H, v);
	}	
  }

/* form of tangent matrix */
GlobalT::SystemTypeT FSDielectricElastomerQ1P02DT::TangentType(void) const
{
	/* Define LHS type based upon analysis type, i.e. static vs. dynamic */
 	int order = fIntegrator->Order();
 	if (order == 2)
 		return GlobalT::kNonSymmetric;
 	else
 		return GlobalT::kSymmetric;

}

/* finalize current step - step is solved */
void FSDielectricElastomerQ1P02DT::CloseStep(void)
{
	/* inherited */
	FiniteStrainT::CloseStep();
	
	/* store converged solution */
	fElementVolume_last = fElementVolume;
}
	
/* restore last converged state */
GlobalT::RelaxCodeT FSDielectricElastomerQ1P02DT::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = FiniteStrainT::ResetStep();
	
	/* store converged solution */
	fElementVolume = fElementVolume_last;

	return relax;
}

/* read restart information from stream */
void FSDielectricElastomerQ1P02DT::ReadRestart(istream& in)
{
	/* inherited */
	FiniteStrainT::ReadRestart(in);
	
	/* read restart data */
	in >> fElementVolume;
	
	/* reset last state */
	fElementVolume_last = fElementVolume;
}

/* write restart information from stream */
void FSDielectricElastomerQ1P02DT::WriteRestart(ostream& out) const
{
	/* inherited */
	FiniteStrainT::WriteRestart(out);
	
	/* read restart data */
	out << fElementVolume << '\n';
}

  // PROTECTED
  // construct a new material support and return a pointer
  MaterialSupportT*
  FSDielectricElastomerQ1P02DT::NewMaterialSupport(MaterialSupportT* p) const
  {
    // allocate
    if (!p) p = new FSDEMatSupportQ1P02DT(1, NumIP());
    
    // inherited initializations
    FiniteStrainT::NewMaterialSupport(p);

    // set parent class fields
    FSDEMatSupportQ1P02DT* ps = dynamic_cast<FSDEMatSupportQ1P02DT*> (p);

    if (ps != 0) {
      ps->SetElectricField(&fE_List);
    }

    return p;
  }

  // construct materials manager and read data
  MaterialListT*
  FSDielectricElastomerQ1P02DT::NewMaterialList(const StringT& name, int size)
  {
    if (name != "large_strain_material_2D") {
      return 0;
    }
    MaterialListT* mlp = 0;

    if (size > 0) {

   	  	/* material support */
      	if (0 == fFSDEMatSupportQ1P02D) {
        	fFSDEMatSupportQ1P02D = dynamic_cast<FSDEMatSupportQ1P02DT*> (NewMaterialSupport());

        if (0 == fFSDEMatSupportQ1P02D) {
          ExceptionT::GeneralFail("FSDielectricElastomerQ1P02DT::NewMaterialList");
        }
      }
      
      mlp = new FSSolidMatList2DT(size, *fFSDEMatSupportQ1P02D);

    }
	else {
		mlp = new FSSolidMatList2DT;
    } 
    return mlp;

  }

/* initialization functions */
void FSDielectricElastomerQ1P02DT::SetShape(void)
{
	/* inherited */
	FiniteStrainT::SetShape();

	/* linked shape functions */
	fCurrShapes = new ShapeFunctionT(*fShapes, fLocCurrCoords);
	if (!fCurrShapes) throw ExceptionT::kOutOfMemory ;

	fCurrShapes->Initialize();
}

  // form shape functions and derivatives; for ComputeOutput
  void FSDielectricElastomerQ1P02DT::SetGlobalShape()
  {
    // inherited
    FiniteStrainT::SetGlobalShape();

    // what needs to be computed
    SetLocalU(fLocScalarPotential);

    for (int i = 0; i < NumIP(); i++) {
    
      // electric field
        dArrayT& E = fE_List[i];
		dMatrixT E1(1, NumSD());

		fShapes->GradU(fLocScalarPotential, E1, i);
		E1 *= -1.0;
		for (int i = 0; i < NumSD(); i++)
			E[i] = E1(0,i);
      }
      
	/* shape function wrt current config - Q1P0 */
	SetLocalX(fLocCurrCoords);
	fCurrShapes->SetDerivatives();

	/* MORE Q1P0 STUFF */
	/* current element number */
	int elem = CurrElementNumber();

	/* compute mean of shape function gradients */
	double H; /* reference volume */
	double& v = fElementVolume[elem];
	SetMeanGradient(fMeanGradient, H, v);
	
	/* last deformed volume */
	double& v_last = fElementVolume_last[elem];

	/* what needs to get computed */
	int material_number = CurrentElement().MaterialNumber();
	bool needs_F = Needs_F(material_number);
	bool needs_F_last = Needs_F_last(material_number);

	/* Getting ready for calculating F_0 */
	fNa_0.Dimension(ElementSupport().NumNodes());
	fDNa_0.Dimension(NumSD(), ElementSupport().NumNodes());
	fJ.Dimension(NumSD());
	fJ_0_inv.Dimension(NumSD());
	fGrad_U.Dimension(2, NumSD());


	double px[2] = {0.0, 0.0};
	dArrayT coords_0(NumSD(), px);
	fShapes->GradU(fLocDisp, fGrad_U, coords_0, fNa_0, fDNa_0);
	fGrad_U.PlusIdentity(); // Computing F_0 Need to be confirmed!!!
	double J_0 = fGrad_U.Det();

	/* loop over integration points */
	for (int i = 0; i < NumIP(); i++)
	{
		/* deformation gradient */
		if (needs_F)
		{
			/* "replace" dilatation */
			dMatrixT& F = fF_List[i];
			double J = F.Det();
			F *= pow((v*J_0)/(H*J), 1.0/2.0);
			
			/* store Jacobian */
			fJacobian[i] = J;
		}

		/* "last" deformation gradient */
		if (needs_F_last)
		{
			/* "replace" dilatation */
			dMatrixT& F = fF_last_List[i];

			double J = F.Det();
			F *= pow((v_last*J_0)/(H*J), 1.0/2.0);
		}
	}	
  }

  // write all current element information to the stream
  void FSDielectricElastomerQ1P02DT::CurrElementInfo(ostream& out) const
  {
    // inherited
    FiniteStrainT::CurrElementInfo(out);

    // write deformation gradients
    out << std::endl;
    out << "electric field at IP:";
    out << std::endl;

    for (int i = 0; i < fE_List.Length(); ++i) {
      out << " ip: " << i + 1 << std::endl << fE_List[i] << std::endl;
    }

    out << std::endl;
  }

//  increment current element - for ComputeOutput
  bool FSDielectricElastomerQ1P02DT::NextElement()
  {
    bool isThereNext = FiniteStrainT::NextElement();

    if (isThereNext == true) {

      const int index = CurrentElement().MaterialNumber();

      ContinuumMaterialT* pMaterial = (*fMaterialList)[index];

      fCurrMaterial = dynamic_cast<FSDEMatQ1P02DT*> (pMaterial);
    }

    return isThereNext;
  }

  // Initialize local arrays
  void FSDielectricElastomerQ1P02DT::SetLocalArrays()
  {
    // look for an electric scalar potential field
    const FieldT* esp = 0;

    if (0 == fElectricScalarPotentialField) {
      esp = ElementSupport().Field("electric_scalar_potential");
      fElectricScalarPotentialField = esp;
    } else {
      esp = fElectricScalarPotentialField;
    }

    if (0 == esp) {

      std::cout << std::endl;
      std::cout << "FSDielectricElastomerQ1P02DT::SetLocalArrays: ";
      std::cout << "Electric scalar potential field not found.";
      std::cout << std::endl;

      throw ExceptionT::kGeneralFail;
    }

    // Inherited
    FiniteStrainT::SetLocalArrays();

    // Allocate storage space:  1 DOF/node for scalar potential
    const int nen = NumElementNodes();

    fLocScalarPotential.Dimension(nen, 1);

    // Register fields
    esp->RegisterLocal(fLocScalarPotential);
    
	/* allocate and set source - for Q1P0 */
	fLocCurrCoords.Dimension(NumElementNodes(), NumSD());
	ElementSupport().RegisterCoordinates(fLocCurrCoords);      
  }

  //
  void FSDielectricElastomerQ1P02DT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
      AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
  {
    for (int i = 0; i < fEqnos.Length(); ++i) {

      const int ndf = NumDOF();	// scalar potential
      const int nen = fConnectivities[i]->MinorDim();
      const int offset = ndf * nen;

      fElectricScalarPotentialField->SetLocalEqnos(*fConnectivities[i],
          fEqnos[i], offset);

    }
	
    ElementBaseT::Equations(eq_1, eq_2);
  }

  //
  int FSDielectricElastomerQ1P02DT::TotalNumDOF() const
  {
 	int mechdof = 2;
 	int elecdof = 1;
    return (mechdof+elecdof);
  }

  const dArrayT&
  FSDielectricElastomerQ1P02DT::ElectricField() const
  {
    return fE_List[CurrIP()];
  }

  //
  const dArrayT&
  FSDielectricElastomerQ1P02DT::ElectricField(int ip) const
  {
    return fE_List[ip];
  }

/* accumulate the residual force on the specified node */
void FSDielectricElastomerQ1P02DT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
	/* not my field */
//	if (&field != &(Field())) return;
	
	/* quick exit */
	bool hasnode = false;
	for (int i=0; i < fBlockData.Length() && !hasnode; i++)
		if (fConnectivities[i]->HasValue(node)) hasnode = true;
	if (!hasnode) return;

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
		cout << "\nWarning: Body forces not yet implemented in DielectricElastomerT";
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
			if (formMa)
			{
				SetLocalU(fLocAcc);
				FormMa(fMassType, constMa*fCurrMaterial->Density(), axisymmetric, &fLocAcc, NULL, NULL);
			}
	
			/* mechanical and electrical reaction forces */
			double mr1, mr2, er1;
			dArrayT react(2);
			
			/* loop over nodes (double-noding OK) */
			int dex = 0;
			int dex2 = 0;
			int whichdof;
			for (int i = 0; i < nodes_u.Length(); i++)
			{
				if (nodes_u[i] == node)
				{
					/* not my field - electrical */
					if (&field != &(Field()))
					{
						er1 = fRHS[dex2+2*NumElementNodes()];
						react[0] = er1;
						whichdof = 1;
					}
					else	// otherwise do mechanical
					{
						mr1 = fRHS[dex];
						mr2 = fRHS[dex+1];	
						react[0] = mr1;
						react[1] = mr2;
						whichdof = NumDOF();
					}
					
					/* components for node - mechanical + electrical DOFs */
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

/* calculate the LHS of residual, or element stiffness matrix */
  void FSDielectricElastomerQ1P02DT::FormStiffness(double constK)
  {
	/* Time integrator info for dynamic problems */
 	int order = fIntegrator->Order();  
  
	/* Matrix format - depends upon time integration order */
    dMatrixT::SymmetryFlagT format = (fLHS.Format()
        == ElementMatrixT::kNonSymmetric)
        ? dMatrixT::kWhole
        : dMatrixT::kUpperOnly;

	/* current element info */
	int el = CurrElementNumber();
	double v = fElementVolume[el];
	double p_bar = fPressure[el];

    fAmm_mat = 0.0;
    fAmm_geo = 0.0;
    fAme = 0.0;
    fAem = 0.0;
    fAee = 0.0;

	/* integration */
	const double* Det    = fCurrShapes->IPDets();
	const double* Weight = fCurrShapes->IPWeights();

	/* initialize */
//	fStressStiff = 0.0;
	fCurrShapes->GradNa(fMeanGradient, fb_bar);	
	
    fShapes->TopIP();
    while (fShapes->NextIP() ) 
    {
		/* double scale factor (for dynamic problems) */
		/* NOTE:  constK = beta * dt^2 */
		double scale = constK*(*Det++)*(*Weight++);
		/* scale factor for K matrix terms without beta * dt^2 factor */
 		double scale1 = scale/constK;	
	
	/* S T R E S S   S T I F F N E S S */			
		/* compute Cauchy stress */
		const dSymMatrixT& cauchy = fCurrMaterial->s_ij();
		const dMatrixT& a = fCurrMaterial->a_ijkl();
		cauchy.ToMatrix(fCauchyStress);
		/* determinant of modified deformation gradient */
		double J_bar = DeformationGradient().Det();
		
		/* detF correction */
		double J_correction = J_bar/fJacobian[CurrIP()];
		double p = J_correction*cauchy.Trace()/2.0;

		/* get shape function gradients matrix */
		fCurrShapes->GradNa(fGradNa);
		fb_sig.MultAB(fCauchyStress, fGradNa);

		/* integration constants */		
		fCauchyStress *= scale*J_correction;

	
		/* using the stress symmetry */
		fAmm_geo.MultQTBQ(fGradNa, a, format, dMatrixT::kAccumulate);
	/* M A T E R I A L   S T I F F N E S S */									
		/* strain displacement matrix */
		Set_B_bar(fCurrShapes->Derivatives_U(), fMeanGradient, fB);

		/* get D matrix */
		fD.SetToScaled(scale*J_correction, fCurrMaterial->c_ijkl());
						
		/* accumulate */
		fAmm_mat.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);
		
	/* Electromechanical Coupling Stiffnesses in current configuration */
	/* May need to modify integration constants (scale) for BIJ and EIJK as compared to CIJKL */
	/* J_correction for eijk terms? */
		dMatrixT bij = fCurrMaterial->b_ij();
		dMatrixT eijk = fCurrMaterial->e_ijk();
		eijk *= scale*J_correction;
		bij *= scale1*J_correction;	// integration constant
		
		/* mechanical-electrical stiffness (24 x 8 matrix for 8-node 3D element) */
		/* Need similar for EIJK1 though with different integration constant */
       	fAme.MultATBC(fB, eijk, fGradNa, dMatrixT::kWhole, dMatrixT::kAccumulate); 
		
		/* mechanical-electrical stiffness (24 x 8 matrix for 8-node 3D element) */
		/* Need similar for EIJK1 though with different integration constant */
//		fAem.MultATBC(fB, eijk1, fGradNa, dMatrixT::kWhole, dMatrixT::kAccumulate); 		
		
		/* electrical-electrical stiffness (8 x 8 matrix for 8-node 3D element) */
  		fAee.MultQTBQ(fGradNa, bij, format, dMatrixT::kAccumulate);
	}
	
	/* stress stiffness into fLHS (i.e. fAmm_mat) */
	fAmm_mat.Expand(fAmm_geo, NumDOF(), dMatrixT::kAccumulate);
	fAem.Transpose();
	
	/* Add mass matrix and non-symmetric electromechanical tangent if dynamic problem */
	if (order == 2)
	{
		/* Calculate mass matrix for dynamic problems */
		MassMatrix();		
		fLHS.AddBlock(0, 0, fMassMatrix);

		/* Need non-symmetric EM coupling tangent for dynamic problems */
		fLHS.AddBlock(fAmm_mat.Rows(), 0, fAem);
	}

	/* Assemble into fLHS, or element stiffness matrix */
	fLHS.AddBlock(0, 0, fAmm_mat);
	fLHS.AddBlock(fAmm_mat.Rows(), fAmm_mat.Cols(), fAee);
	fLHS.AddBlock(0, fAmm_mat.Cols(), fAme);
			// Saving the fLHS matrix
		/* ofstream myLHS;
		myLHS.open("fLHS_bulk.txt");
		for (int i = 0; i < fLHS.Rows(); i++)
		{
			for (int j = 0; j < fLHS.Cols(); j++)
			{
				// myLHS << "fLHS(" << i << "," << j << ")= " << fLHS(i, j); // List the values of fLHS(i,j)
				if (fLHS(i, j) == 0)
					myLHS << "0.00000" << " ";
				else
					myLHS << fLHS(i, j) << " "; // See the matrix form of fLHS for a quick look
			}
			myLHS << endl;
		}
		myLHS.close(); */
  }

/* Compute RHS, or residual of element equations */
  void FSDielectricElastomerQ1P02DT::FormKd(double constK)
  {  	
	//  cout << "FormKd" << endl;
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

	/* constant pressure */
	double& p_bar = fPressure[elem];
	p_bar = 0.0;

    fCurrShapes->TopIP();
    while (fCurrShapes->NextIP() )
    {
		/* strain displacement matrix */
		Set_B_bar(fCurrShapes->Derivatives_U(), fMeanGradient, fB);

		/* B^T * Cauchy stress */
		const dSymMatrixT& cauchy = fCurrMaterial->s_ij();
		fB.MultTx(cauchy, fNEEvec);
		
		/* determinant of modified deformation gradient */
		double J_bar = DeformationGradient().Det();
		
		/* detF correction */
		double J_correction = J_bar/fJacobian[CurrIP()];
		
		/* integrate pressure */
		p_bar += (*Weight)*(*Det)*J_correction*cauchy.Trace()/2.0;

		/* double scale factor */
		double scale = constK*(*Det++)*(*Weight++);   
		
		/* accumulate - use Rmech instead of fRHS */
		Rmech.AddScaled(scale*J_correction, fNEEvec);
    
 	  	/* electrical stress in current configuration */
 	  	dArrayT di = fCurrMaterial->d_i();
 	  	//for (int i = 0; i < di.Length(); i++)
 	  	// 	  		cout << di[i] << endl;
 	  	//cout << scale << endl;
 	  	//cout << J_correction << endl;
 	  	di *= scale*J_correction;
 	  	//for (int i = 0; i < di.Length(); i++)
 	  	//  		cout << di[i] << endl;
 	  	/* get shape function gradients matrix */
		fCurrShapes->GradNa(fGradNa);
		fGradNa.MultTx(di, Relec, 1.0, dMatrixT::kAccumulate);  
	}

  	Relec *= -1.0;	
 	Rtotal.CopyIn(0, Rmech);
 	Rtotal.CopyIn(Rmech.Length(), Relec);
 	fRHS += Rtotal; 
	
            /* Saving RHS to a file */
           // ofstream myRHS;
           // myRHS.open("fRHS_bulk.txt");
           // for (int i = 0; i < fRHS.Length(); i++)
            //{
            //      myRHS << "fRHS(" << i << ") = " << fRHS[i] << endl;
           // }
           // myRHS.close();
 
	/* volume averaged */
	p_bar /= fElementVolume[CurrElementNumber()]; 
  }

/* Dummy mass matrix for dynamic calculations */
void FSDielectricElastomerQ1P02DT::FormMass(MassTypeT mass_type, double constM, bool axisymmetric, const double* ip_weight)
{
	/* SHOULD THIS BE INHERITED? */
	
	/* Do nothing but add 0 to fLHS - implement mass in FormStiffness */
	fLHS += 0.0;	
}

/* Calculate inertial force for dynamic calculations */
void FSDielectricElastomerQ1P02DT::FormMa(MassTypeT mass_type, double constM, bool axisymmetric, 
	const LocalArrayT* nodal_values, const dArray2DT* ip_values, const double* ip_weight)
{
	const char caller[] = "FSDielectricElastomerQ1P02DT::FormMa";
    
    /* Define mechanical contribution to inertial force only */
    /* Currently assuming a lumped mass matrix - good for implicit dynamics convergence */
	dArrayT Rmech(NumSD()*NumElementNodes());
	dArrayT Rtotal((NumSD()+1)*NumElementNodes());
	Rmech = 0.0;
	Rtotal = 0.0;

	/* Compute lumped mass matrix */
	MassMatrix();	
	/* Comply with negative factor in ContinuumElementT::FormMa */
	fMassMatrix *= -1.0;	
	
	/* init nodal values */
	if (nodal_values)
		nodal_values->ReturnTranspose(fNEEvec);
	else {
		ExceptionT::GeneralFail(caller, "expecting nodal values for lumped mass");
	}
	
	double* pAcc = fNEEvec.Pointer();
	double* pRes = Rmech.Pointer();
	int     massdex = 0;

	int nee = nodal_values->Length();

	for (int i = 0; i < nee; i++)
	{
		*pRes++ += (*pAcc++)*fMassMatrix(massdex,massdex);
		massdex++;
	}
	
 	Rtotal.CopyIn(0, Rmech);
	fRHS += Rtotal;
}

void FSDielectricElastomerQ1P02DT::MassMatrix()
{
	/* Calculate mass matrix for mechanical portion of LHS for FormStiffness */
	/* Implement lumped mass matrix only - better convergence for implicit dynamics 
		as compared to consistent mass - see Hughes FEM book */
	int nen = NumElementNodes();
	int ndof = NumDOF();
	int nsd = NumSD();
    int nme = nen * nsd;	// # of mechanical DOFs per element
	dArrayT NEEvec(nme);
	NEEvec = 0.0;
	double dsum = 0.0;
	double totmas = 0.0;
	fMassMatrix = 0.0;

	fShapes->TopIP();	
	while (fShapes->NextIP() != 0) {
		
		/* integration factor - ignoring constM factor */
		double temp1 = fShapes->IPDet() * fShapes->IPWeight();
//		if (ip_weight) temp1 *= *ip_weight++;

		const double* Na = fShapes->IPShapeU();
		totmas += temp1;
		for (int lnd = 0; lnd < nen; lnd++) {
			double temp2 = temp1*Na[lnd]*Na[lnd];
			dsum += temp2;
			NEEvec[lnd] += temp2;
		}
	}
	
	/* scale diagonal to conserve total mass */
	double diagmass = totmas/dsum;
	
	/* lump mass onto diagonal */
	double* pmass = fMassMatrix.Pointer();
	int inc = fMassMatrix.Rows() + 1;
	for (int lnd = 0; lnd < nen; lnd++)
	{
		double temp = diagmass*NEEvec[lnd];
		for (int ed = 0; ed < ndof; ed++)
		{
			*pmass += temp;
			pmass += inc;
		}
	}

}

/***********************************************************************
 * Private
 ***********************************************************************/

/* compute mean shape function gradient, Hughes (4.5.23) */
void FSDielectricElastomerQ1P02DT::SetMeanGradient(dArray2DT& mean_gradient, double& H, double& v) const
{
	/* assume same integration rule defined for current and references
	 * shape functions */
	int nip = NumIP();
	const double*   det = fCurrShapes->IPDets();
	const double* det_0 = fShapes->IPDets();
	const double*     w = fShapes->IPWeights();

	/* H and current volume */
	H = 0.0;
	v = 0.0;

	for (int i = 0; i < nip; i++)
	{
		H += w[i]*det_0[i];
		v += w[i]*det[i];
	}

	/* initialize */
	mean_gradient = 0.0;			

	/* integrate */
	for (int i = 0; i < nip; i++)
		mean_gradient.AddScaled(w[i]*det[i]/v, fCurrShapes->Derivatives_U(i));
}

void FSDielectricElastomerQ1P02DT::bSp_bRq_to_KSqRp(const dMatrixT& b, dMatrixT& K) const
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (b.Length() != K.Rows() ||
	    K.Rows() != K.Cols()) ExceptionT::SizeMismatch("SimoQ1P0::bSp_bRq_to_KSqRp");
#endif

	int dim = K.Rows();
	int sub_dim = b.Rows();
	int S = 0;
	int p = 0;
	for (int i = 0; i < dim; i++)
	{
		int R = 0;
		int q = 0;
		for (int j = 0; j < dim; j++)
		{
			K(i,j) = b(q,S)*b(p,R);
		
			q++;
			if (q == sub_dim) {
				R++;
				q = 0;
			}
		}
		p++;
		if (p == sub_dim) {
			S++;
			p = 0;
		}
	}	
}

  // extrapolate from integration points and compute output nodal/element values
  void FSDielectricElastomerQ1P02DT::ComputeOutput(const iArrayT& n_codes,
      dArray2DT& n_values, const iArrayT& e_codes, dArray2DT& e_values)
  {
    //
    // number of output values
    //
    int n_out = n_codes.Sum();
    int e_out = e_codes.Sum();

    // nothing to output
    if (n_out == 0 && e_out == 0) return;

    // dimensions
    int nsd = NumSD();
    int ndof = NumDOF();
    int nen = NumElementNodes();
    int nnd = ElementSupport().NumNodes();

    // reset averaging work space
    ElementSupport().ResetAverage(n_out);

    // allocate element results space
    e_values.Dimension(NumElements(), e_out);

    // nodal work arrays
    dArray2DT nodal_space(nen, n_out);
    dArray2DT nodal_all(nen, n_out);
    dArray2DT coords, disp;
    dArray2DT nodalstress, princstress, matdat;
    dArray2DT energy, speed;
    dArray2DT ndElectricScalarPotential;
//	dArrayT ndElectricScalarPotential;
    dArray2DT ndElectricDisplacement;
    dArray2DT ndElectricField;

    // ip values
    dArrayT ipmat(n_codes[iMaterialData]), ipenergy(1);
    dArrayT ipspeed(nsd), ipprincipal(nsd);
    dMatrixT ippvector(nsd);
    dArrayT Eall(nsd+1);

    // set shallow copies
    double* pall = nodal_space.Pointer();
    coords.Alias(nen, n_codes[iNodalCoord], pall);
    pall += coords.Length();
    disp.Alias(nen, n_codes[iNodalDisp], pall);
    pall += disp.Length();

    nodalstress.Alias(nen, n_codes[iNodalStress], pall);
    pall += nodalstress.Length();
    princstress.Alias(nen, n_codes[iPrincipal], pall);
    pall += princstress.Length();
    energy.Alias(nen, n_codes[iEnergyDensity], pall);
    pall += energy.Length();
    speed.Alias(nen, n_codes[iWaveSpeeds], pall);
    pall += speed.Length();
    matdat.Alias(nen, n_codes[iMaterialData], pall);
    pall += matdat.Length();

    ndElectricDisplacement.Alias(nen, n_codes[ND_ELEC_DISP], pall);
    pall += ndElectricDisplacement.Length();

    ndElectricField.Alias(nen, n_codes[ND_ELEC_FLD], pall);
    pall += ndElectricField.Length();
// 	ndElectricField.Alias(nen, NumSD(), pall);
// 	pall += ndElectricField.Length();
    
    ndElectricScalarPotential.Alias(nen, n_codes[ND_ELEC_POT_SCALAR], pall);
    pall += ndElectricScalarPotential.Length();

    // element work arrays
    dArrayT element_values(e_values.MinorDim());
    pall = element_values.Pointer();
    dArrayT centroid, ip_centroid, ip_mass;
    dArrayT ip_coords(nsd);
    if (e_codes[iCentroid]) {
      centroid.Alias(nsd, pall);
      pall += nsd;
      ip_centroid.Dimension(nsd);
    }
    if (e_codes[iMass]) {
      ip_mass.Alias(NumIP(), pall);
      pall += NumIP();
    }
    double w_tmp, ke_tmp;
    double mass;
    double& strain_energy = (e_codes[iStrainEnergy])
        ? *pall++
        : w_tmp;
    double& kinetic_energy = (e_codes[iKineticEnergy])
        ? *pall++
        : ke_tmp;
    dArrayT linear_momentum, ip_velocity;

    if (e_codes[iLinearMomentum]) {
      linear_momentum.Alias(ndof, pall);
      pall += ndof;
      ip_velocity.Dimension(ndof);
    } else if (e_codes[iKineticEnergy]) ip_velocity.Dimension(ndof);

    dArray2DT ip_stress;
    if (e_codes[iIPStress]) {
      ip_stress.Alias(NumIP(), e_codes[iIPStress] / NumIP(), pall);
      pall += ip_stress.Length();
    }
    dArray2DT ip_material_data;
    if (e_codes[iIPMaterialData]) {
      ip_material_data.Alias(NumIP(), e_codes[iIPMaterialData] / NumIP(), pall);
      pall += ip_material_data.Length();
      ipmat.Dimension(ip_material_data.MinorDim());
    }

    dArray2DT ipElectricDisplacement;
    if (e_codes[IP_ELEC_DISP]) {
      ipElectricDisplacement.Alias(NumIP(), NumSD(), pall);
      pall += NumIP() * NumSD();
    }

    dArray2DT ipElectricField;
    if (e_codes[IP_ELEC_FLD]) {
      ipElectricField.Alias(NumIP(), NumSD(), pall);
      pall += NumIP() * NumSD();
    }

    // check that degrees are displacements
    int interpolant_DOF = InterpolantDOFs();

    Top();
    while (NextElement()) {

      if (CurrentElement().Flag() == ElementCardT::kOFF) continue;

      // initialize
      nodal_space = 0.0;

      // global shape function values
      SetGlobalShape();

      // collect nodal values
      if (e_codes[iKineticEnergy] || e_codes[iLinearMomentum]) {
        if (fLocVel.IsRegistered())
          SetLocalU(fLocVel);
        else
          fLocVel = 0.0;
      }

      // coordinates and displacements all at once
      if (n_codes[iNodalCoord]) fLocInitCoords.ReturnTranspose(coords);
      if (n_codes[iNodalDisp]) {
        if (interpolant_DOF)
          fLocDisp.ReturnTranspose(disp);
        else
          NodalDOFs(CurrentElement().NodesX(), disp);
      }

      if (n_codes[ND_ELEC_POT_SCALAR]) {
        if (interpolant_DOF) {
          fLocScalarPotential.ReturnTranspose(ndElectricScalarPotential);
        } else {
          NodalDOFs(CurrentElement().NodesX(), ndElectricScalarPotential);
        }
      }

      // initialize element values
      mass = strain_energy = kinetic_energy = 0;
      if (e_codes[iCentroid]) centroid = 0.0;
      if (e_codes[iLinearMomentum]) linear_momentum = 0.0;
      const double* j = fShapes->IPDets();
      const double* w = fShapes->IPWeights();

      // integrate
      dArray2DT Na_X_ip_w;
      fShapes->TopIP();
      while (fShapes->NextIP() != 0) {

        // density may change with integration point
        double density = fCurrMaterial->Density();

        // element integration weight
        double ip_w = (*j++) * (*w++);

        if (qNoExtrap) {
          Na_X_ip_w.Dimension(nen, 1);
          for (int k = 0; k < nen; k++) {
            Na_X_ip_w(k, 0) = 1.;
          }
        }

        // get Cauchy stress
        const dSymMatrixT& stress = fCurrMaterial->s_ij();
        dSymMatrixT strain;

        // stresses
        if (n_codes[iNodalStress]) {
          if (qNoExtrap) {
            for (int k = 0; k < nen; k++) {
              nodalstress.AddToRowScaled(k, Na_X_ip_w(k, 0), stress);
            }
          } else {
            fShapes->Extrapolate(stress, nodalstress);
          }
        }

        if (e_codes[iIPStress]) {
          double* row = ip_stress(fShapes->CurrIP());
          strain.Set(nsd, row);
          strain = stress;
          row += stress.Length();
          strain.Set(nsd, row);
          fCurrMaterial->Strain(strain);
        }

        // wave speeds
        if (n_codes[iWaveSpeeds]) {
          // acoustic wave speeds
          fCurrMaterial->WaveSpeeds(fNormal, ipspeed);
          if (qNoExtrap) {
            for (int k = 0; k < nen; k++) {
              speed.AddToRowScaled(k, Na_X_ip_w(k, 0), ipspeed);
            }
          } else {
            fShapes->Extrapolate(ipspeed, speed);
          }
        }

        // principal values - compute principal before smoothing
        if (n_codes[iPrincipal]) {
          // compute eigenvalues
          stress.PrincipalValues(ipprincipal);

          if (qNoExtrap) {
            for (int k = 0; k < nen; k++) {
              princstress.AddToRowScaled(k, Na_X_ip_w(k, 0), ipprincipal);
            }
          } else {
            fShapes->Extrapolate(ipprincipal, princstress);
          }
        }

        // strain energy density
        if (n_codes[iEnergyDensity] || e_codes[iStrainEnergy]) {
          double ip_strain_energy = fCurrMaterial->StrainEnergyDensity();

          // nodal average
          if (n_codes[iEnergyDensity]) {
            ipenergy[0] = ip_strain_energy;
            if (qNoExtrap) {
              for (int k = 0; k < nen; k++) {
                energy.AddToRowScaled(k, Na_X_ip_w(k, 0), ipenergy);
              }
            } else {
              fShapes->Extrapolate(ipenergy, energy);
            }
          }

          // integrate over element
          if (e_codes[iStrainEnergy]) {
            strain_energy += ip_w * ip_strain_energy;
          }

        }

        // material stuff
        if (n_codes[iMaterialData] || e_codes[iIPMaterialData]) {
          // compute material output
          fCurrMaterial->ComputeOutput(ipmat);

          // store nodal data
          if (n_codes[iMaterialData]) {
            if (qNoExtrap) {
              for (int k = 0; k < nen; k++) {
                matdat.AddToRowScaled(k, Na_X_ip_w(k, 0), ipmat);
              }
            } else {
              fShapes->Extrapolate(ipmat, matdat);
            }
          }

          // store element data
          if (e_codes[iIPMaterialData]) {
            ip_material_data.SetRow(fShapes->CurrIP(), ipmat);
          }
        }

        // mass averaged centroid
        if (e_codes[iCentroid] || e_codes[iMass]) {
          // mass
          mass += ip_w * density;

          // integration point mass
          if (e_codes[iMass]) ip_mass[fShapes->CurrIP()] = ip_w * density;

          // moment
          if (e_codes[iCentroid]) {
            fShapes->IPCoords(ip_centroid);
            centroid.AddScaled(ip_w * density, ip_centroid);
          }
        }

        // kinetic energy/linear momentum
        if (e_codes[iKineticEnergy] || e_codes[iLinearMomentum]) {
          // velocity at integration point
          fShapes->InterpolateU(fLocVel, ip_velocity);
          double ke_density = 0.5 * density * dArrayT::Dot(ip_velocity,
              ip_velocity);

          // kinetic energy
          if (e_codes[iKineticEnergy]) {
            kinetic_energy += ip_w * ke_density;
          }

          // linear momentum
          if (e_codes[iLinearMomentum]) {
            linear_momentum.AddScaled(ip_w * density, ip_velocity);
          }

        }

        // electric displacements
        const dArrayT& D = fCurrMaterial->D_I();
        if (n_codes[ND_ELEC_DISP]) {
          if (qNoExtrap) {
            for (int k = 0; k < nen; k++) {
              ndElectricDisplacement.AddToRowScaled(k, Na_X_ip_w(k, 0), D);
            }
          } else {
            fShapes->Extrapolate(D, ndElectricDisplacement);
          }
        }

        // electric field
        const dArrayT& E = fCurrMaterial->E_I();
        double emag = E.Magnitude();
        Eall[0] = E[0];
        Eall[1] = E[1];
        Eall[2] = emag;
        //Eall = E;
        if (n_codes[ND_ELEC_FLD]) {
          if (qNoExtrap) {
            for (int k = 0; k < nen; k++) {
              ndElectricField.AddToRowScaled(k, Na_X_ip_w(k, 0), Eall);
            }
          } else {
            fShapes->Extrapolate(Eall, ndElectricField);
          }
        }

      }

      // copy in the cols
      int colcount = 0;
      nodal_all.BlockColumnCopyAt(disp, colcount);
      colcount += disp.MinorDim();

      nodal_all.BlockColumnCopyAt(coords, colcount);
      colcount += coords.MinorDim();

      if (qNoExtrap) {
        double nip(fShapes->NumIP());
        nodalstress /= nip;
        princstress /= nip;
        energy /= nip;
        speed /= nip;
        matdat /= nip;
        ndElectricDisplacement /= nip;
        ndElectricField /= nip;
        ndElectricScalarPotential /= nip;
      }
      nodal_all.BlockColumnCopyAt(nodalstress, colcount);
      colcount += nodalstress.MinorDim();

      nodal_all.BlockColumnCopyAt(princstress, colcount);
      colcount += princstress.MinorDim();

      nodal_all.BlockColumnCopyAt(energy, colcount);
      colcount += energy.MinorDim();

      nodal_all.BlockColumnCopyAt(speed, colcount);
      colcount += speed.MinorDim();

      nodal_all.BlockColumnCopyAt(matdat, colcount);
      colcount += matdat.MinorDim();

      nodal_all.BlockColumnCopyAt(ndElectricDisplacement, colcount);
      colcount += ndElectricDisplacement.MinorDim();

      nodal_all.BlockColumnCopyAt(ndElectricField, colcount);
      colcount += ndElectricField.MinorDim();

      nodal_all.BlockColumnCopyAt(ndElectricScalarPotential, colcount);
      colcount += ndElectricScalarPotential.MinorDim();

      // accumulate - extrapolation done from ip's to corners => X nodes
      ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_all);

      // element values
      if (e_codes[iCentroid]) centroid /= mass;

      // store results
      e_values.SetRow(CurrElementNumber(), element_values);

    }

    // get nodally averaged values
    const OutputSetT& output_set = ElementSupport().OutputSet(fOutputID);
    const iArrayT& nodes_used = output_set.NodesUsed();
    dArray2DT extrap_values(nodes_used.Length(), n_out);
    extrap_values.RowCollect(nodes_used, ElementSupport().OutputAverage());

    int tmpDim = extrap_values.MajorDim();
    n_values.Dimension(tmpDim, n_out);
    n_values.BlockColumnCopyAt(extrap_values, 0);
   }

} // namespace Tahoe
