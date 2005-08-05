/* $Id: MFGPAssemblyT.cpp,v 1.3 2005-08-05 22:24:04 kyonten Exp $ */
#include "MFGPAssemblyT.h"
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include "ifstreamT.h"
#include "ofstreamT.h"
#include "ModelManagerT.h"
#include "D3MeshfreeShapeFunctionT.h"
#include "eIntegratorT.h"
#include "Traction_CardT.h"
#include "iAutoArrayT.h"
#include "OutputSetT.h"
#include "ScheduleT.h"
#include "ParameterContainerT.h"
#include "CommunicatorT.h"
#include "ParameterUtils.h"
#include "OutputSetT.h"

//TEMP: all this for general traction BC implementation?
#include "VariArrayT.h"
#include "nVariArray2DT.h"
#include "VariLocalArrayT.h"

/* materials lists */
#include "MFGPMatSupportT.h"
#include "MFGPMatListT.h"
#include "MFGPSolidMatListT.h"
 
const double Pi = acos(-1.0);

using namespace Tahoe;

/* constructor */
MFGPAssemblyT::MFGPAssemblyT(const ElementSupportT& support):
	ElementBaseT(support), //pass the displacement field to the base class
	u(LocalArrayT::kDisp),
	u_n(LocalArrayT::kLastDisp),
	DDu(LocalArrayT::kAcc),
	lambda(LocalArrayT::kDisp),
	lambda_n(LocalArrayT::kLastDisp),
	fInitCoords_displ(LocalArrayT::kInitCoords),
	fCurrCoords_displ(LocalArrayT::kCurrCoords),
	fInitCoords_plast(LocalArrayT::kInitCoords),
	fCurrCoords_plast(LocalArrayT::kCurrCoords),
	fBodySchedule(NULL),
	fTractionBCSet(0),
	fDispl(NULL),
	fPlast(NULL),
	fKuu(ElementMatrixT::kNonSymmetric),
	fKulambda(ElementMatrixT::kNonSymmetric),
	fKlambdau(ElementMatrixT::kNonSymmetric),
	fKlambdalambda(ElementMatrixT::kNonSymmetric),
	fKuu_temp(ElementMatrixT::kNonSymmetric),
	fKulambda_temp(ElementMatrixT::kNonSymmetric),
	fKlambdau_temp(ElementMatrixT::kNonSymmetric),
	fKlambdalambda_temp(ElementMatrixT::kNonSymmetric),
	fGroupCommunicator(NULL),
	fShapes_displ(NULL),
	fShapes_plast(NULL),
	fNumIP_displ(0),
	fNumIP_plast(0),
	fOutputID(),
	fGeometryCode_displ(GeometryT::kNone),
	fGeometryCode_plast(GeometryT::kNone),
	fMFGPMatList(NULL)
{
	SetName("mfgp_assembly");
}

/* destructor */
MFGPAssemblyT::~MFGPAssemblyT(void) 
{  
	delete fGroupCommunicator;
	delete fBodySchedule;
	delete fShapes_displ;
	delete fShapes_plast;  
}

/* form LHS */
void MFGPAssemblyT::LHSDriver(GlobalT::SystemTypeT)
{
	/* Everything done in RHSDriver for efficiency */
}

/* choose between monolithic or staggered coupled solution scheme */
void MFGPAssemblyT::RHSDriver(void)	
{
	int curr_group = ElementSupport().CurrentGroup();

	/* traction boundary conditions acting on displacement equations */
	if (curr_group == fDispl->Group()) 
		ApplyTractionBC();

	/* choose solution method */
	if (fDispl->Group() == fPlast->Group())
	  RHSDriver_monolithic();
	else
	  RHSDriver_staggered();
}

/* relate local and global equation numbers */
void MFGPAssemblyT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
  /* inherited */
	ElementBaseT::Equations(eq_1, eq_2);

	/* mark traction BC data as old */
	fTractionBCSet = 0;
}

/* form of tangent matrix */
GlobalT::SystemTypeT MFGPAssemblyT::TangentType(void) const
{
	return GlobalT::kNonSymmetric; 
}

/* return true if the element contributes to the solution of a given group. */
bool MFGPAssemblyT::InGroup(int group) const
{
	return group == fDispl->Group() ||
	       group == fPlast->Group();
}

/* accessors */
const int& MFGPAssemblyT::CurrIP(void) const
{
	return ShapeFunction().CurrIP();
}

/* the coordinates of the current integration points */
void MFGPAssemblyT::IP_Coords(dArrayT& ip_coords) const
{
	/* computed by shape functions */
	ShapeFunction().IPCoords(ip_coords);
}

/* interpolate the nodal field values to the current integration point */
void MFGPAssemblyT::IP_Interpolate(const LocalArrayT& nodal_u, dArrayT& ip_u) const
{
    /* computed by shape functions */
    ShapeFunction().InterpolateU(nodal_u, ip_u);
}

void MFGPAssemblyT::IP_Interpolate(const LocalArrayT& nodal_u, dArrayT& ip_u, int ip) const
{
    /* computed by shape functions */
    ShapeFunction().InterpolateU(nodal_u, ip_u, ip);
}

void MFGPAssemblyT::IP_ComputeLaplacian(const LocalArrayT& field, dArrayT& laplacian) const
{
#if __option(extended_errorcheck)
	if (laplacian.Length() != field.MinorDim())
    	throw ExceptionT::kSizeMismatch;
#endif
    /* computed by shape functions */
    //ShapeFunction().LaplaceU(field, laplacian);
    dMatrixT gradgrad_u;
    ShapeFunction().GradGradU(field, gradgrad_u);
    laplacian = 0.0;
    
    for (int j = 0; j < field.MinorDim(); j++) 
    	for (int i = 0; i < NumSD(); i++) 
    		laplacian[j] += gradgrad_u(j,i);
}

void MFGPAssemblyT::IP_ComputeLaplacian(const LocalArrayT& field, dArrayT& laplacian, int ip) const
{
#if __option(extended_errorcheck)
	if (laplacian.Length() != field.MinorDim())
    	throw ExceptionT::kSizeMismatch;
#endif
    /* computed by shape functions */
    //ShapeFunction().LaplaceU(field, laplacian);
    dMatrixT gradgrad_u;
    ShapeFunction().GradGradU(field, gradgrad_u, ip);
    laplacian = 0.0;
    // cout << endl << "gradgrad_u" << gradgrad_u << endl;
    for (int j = 0; j < field.MinorDim(); j++) 
    	for (int i = 0; i < NumSD(); i++) 
    		laplacian[j] += gradgrad_u(j,i);
}

/* field gradients */
void MFGPAssemblyT::IP_ComputeGradient(const LocalArrayT& field, 
	dMatrixT& gradient) const
{
	/* computed by shape functions */
	ShapeFunction().GradU(field, gradient);
}

void MFGPAssemblyT::IP_ComputeGradient(const LocalArrayT& field, 
	dMatrixT& gradient, int ip) const
{
	/* computed by shape functions */
	ShapeFunction().GradU(field, gradient, ip);
}

/* extrapolate all integration point values to the nodes */
void MFGPAssemblyT::IP_ExtrapolateAll(const dArrayT& ip_values, 
	dArrayT& nodal_values) const
{
	/* computed by shape functions */
	ShapeFunction().ExtrapolateAll(ip_values, nodal_values);
}

/* extrapolate the integration point stresses and internal variables, 
   check the yield condition on the nodes and pass the flag whether 
   the nodes are elastically or plastically loaded
*/ 
void MFGPAssemblyT::CheckNodalYield()
{
	// do nothing
}

/* initialize/finalize step */
void MFGPAssemblyT::InitStep(void)
{
	/* inherited */
	ElementBaseT::InitStep();
	
	/* set material variables */
	if (fMFGPMatList)  fMFGPMatList->InitStep();
}

/* initialize/finalize step */
void MFGPAssemblyT::CloseStep(void)
{
	/* inherited */
	ElementBaseT::CloseStep();

	if (fMFGPMatList) 
	{
		/* set material variables */
		fMFGPMatList->CloseStep();

		/* update element level internal variables */
		if (fMFGPMatList->HasHistoryMaterials())
		{
			Top();
			while (NextElement())
			{
				const ElementCardT& element = CurrentElement();
				if (element.IsAllocated())
				{
					MFGPMaterialT* pmat = (*fMFGPMatList)[element.MaterialNumber()];

				/* material update function */
					pmat->UpdateHistory();
				}
			}
		}
	} 
}

/* write element group parameters to out */
void MFGPAssemblyT::PrintControlData(ostream& out) const
{
	/* inherited */
	//ElementBaseT::PrintControlData(out);

	out << " Displacement field. . . . . . . . . . . . . . . = \"" << fDispl->Name() << "\"\n";
	out << " Plastic gradient field. . . . . . . . . . . . . = \"" << fPlast->Name() << "\"\n";
	out << " Element geometry code . . . . . . . . . . . . . = " << fGeometryCode_displ << '\n';
	out << "    eq." << GeometryT::kPoint         << ", point\n";
	out << "    eq." << GeometryT::kLine          << ", line\n";
	out << "    eq." << GeometryT::kQuadrilateral << ", quadrilateral\n";
	out << "    eq." << GeometryT::kTriangle	  << ", triangle\n";
	out << "    eq." << GeometryT::kHexahedron	  << ", hexahedron\n";
	out << "    eq." << GeometryT::kTetrahedron   << ", tetrahedron\n";
	out << " Number of integration points. . . . . . . . . . = " << fNumIP_plast    << '\n';
}


/* resets to the last converged solution */
GlobalT::RelaxCodeT MFGPAssemblyT::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::ResetStep();

	/* update material internal variables */
	if (fMFGPMatList && fMFGPMatList->HasHistoryMaterials())
	{
		Top();
		while (NextElement())
		{
			const ElementCardT& element = CurrentElement();		
			if (element.IsAllocated())
			{
				MFGPMaterialT* pmat = (*fMFGPMatList)[element.MaterialNumber()];

				/* material reset function */
				pmat->ResetHistory();
			}
		}
	}
	return relax;
}

/* write restart data to the output stream */
void MFGPAssemblyT::WriteRestart(ostream& out) const
{
	/* inherited */
	ElementBaseT::WriteRestart(out);

	/* update element level internal variables */
	if (fMFGPMatList && fMFGPMatList->HasHistoryMaterials())
	{
		for (int i = 0; i < fElementCards.Length(); i++)
		{
			const ElementCardT& element = fElementCards[i];
			out << element.IsAllocated() << '\n';
			if (element.IsAllocated()) element.WriteRestart(out);
		}
	}
}

/* element level reconfiguration for the current time increment */
GlobalT::RelaxCodeT MFGPAssemblyT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	/* loop over materials */
	for (int i = 0; i < fMFGPMatList->Length(); i++)
		relax = GlobalT::MaxPrecedence((*fMFGPMatList)[i]->RelaxCode(), relax);

	return relax;
}

/* read restart data to the output stream */
void MFGPAssemblyT::ReadRestart(istream& in)
{
	/* inherited */
	ElementBaseT::ReadRestart(in);

	/* update element level internal variables */
	if (fMFGPMatList && fMFGPMatList->HasHistoryMaterials())
	{
		for (int i = 0; i < fElementCards.Length(); i++)
		{
			int isallocated;
			in >> isallocated;
			if (isallocated) fElementCards[i].ReadRestart(in);
		}
	}
}

/* register output and output labels */
void MFGPAssemblyT::RegisterOutput(void)
{
	/* nodal output */
	iArrayT n_counts;
	SetNodalOutputCodes(IOBaseT::kAtInc, fNodalOutputCodes, n_counts);

	/* element output */
	iArrayT e_counts;
	SetElementOutputCodes(IOBaseT::kAtInc, fElementOutputCodes, e_counts);

	/* inherited */
	if (n_counts.Length() == 0 && e_counts.Length() == 0)
		ElementBaseT::RegisterOutput();
	else
	{
		/* collect variable labels */
		ArrayT<StringT> n_labels(n_counts.Sum());
		ArrayT<StringT> e_labels(e_counts.Sum());
		GenerateOutputLabels(n_counts, n_labels, e_counts, e_labels);

		/* block ID's used by the group */
		ArrayT<StringT> block_ID(fBlockData.Length());
		for (int i = 0; i < block_ID.Length(); i++)
			block_ID[i] = fBlockData[i].ID();

		/* set output specifier */
		OutputSetT output_set(fGeometryCode_displ, block_ID, fConnectivities, n_labels, e_labels, false);
		
		/* register and get output ID */
		fOutputID = ElementSupport().RegisterOutput(output_set);
	}
}

void MFGPAssemblyT::WriteOutput(void)
{
	/* regular output */
	IOBaseT::OutputModeT mode = IOBaseT::kAtInc;

	/* map output flags to count of values */
	iArrayT n_counts;
	SetNodalOutputCodes(mode, fNodalOutputCodes, n_counts);
	iArrayT e_counts;
	SetElementOutputCodes(mode, fElementOutputCodes, e_counts);

	/* inherited */
	if (n_counts.Length() == 0 && e_counts.Length() == 0)
		ElementBaseT::WriteOutput();
	else
	{
		/* calculate output values */
		dArray2DT n_values;
		dArray2DT e_values;
		ComputeOutput(n_counts, n_values, e_counts, e_values);

		/* send to output */
		ElementSupport().WriteOutput(fOutputID, n_values, e_values);
	}
}	

/* return geometry and number of nodes on each facet */
void MFGPAssemblyT::FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geometry, 
	iArrayT& num_facet_nodes) const
{
	/* from integration domain */
	ShapeFunction().FacetGeometry(facet_geometry, num_facet_nodes);
}

void MFGPAssemblyT::SetStatus(const ArrayT<ElementCardT::StatusT>& status) 
{
  /* work space */
  dArrayT state;
  dArrayT t_in;

  /* loop over elements and initial state variables */
  int elem_num = 0;
  Top();
  while (NextElement())
    {
      /* current element */
      ElementCardT::StatusT& flag = CurrentElement().Flag();
      flag = status[elem_num++];
      /* material pointer */
      MFGPMaterialT* pmat = (*fMFGPMatList)[CurrentElement().MaterialNumber()];

      if (flag == ElementCardT::kMarkON) {
	  	if (pmat->NeedsPointInitialization()) {
	  /* global shape function values */
	  	SetGlobalShape();

	  	fShapes_displ->TopIP();
	  	while (fShapes_displ->NextIP())
	    	pmat->PointInitialize();
		}
		flag = ElementCardT::kON;
      }
      else if (flag == ElementCardT::kMarkOFF)
	  	flag = ElementCardT::kOFF;
    }
}

/* initial condition/restart functions (per time sequence) */
void MFGPAssemblyT::InitialCondition(void)
{
	/* inherited */
	ElementBaseT::InitialCondition();

	/* check for initialization materials */
	bool need_init = false;
	if (fMFGPMatList)
		for (int i = 0; i < fMFGPMatList->Length() && !need_init; i++)
			need_init = (*fMFGPMatList)[i]->NeedsPointInitialization();

	/* need to run through elements */
	if (fStoreShape || need_init)
	{
		/* initialize storage */
		if (fStoreShape) fShapes_displ->InitStore(NumElements(), &(fElementCards.Position()));

		/* loop over elements */
		Top();
		while (NextElement())
		{
			/* compute shape function derivarives */
			SetGlobalShape();

			/* store */
			if (fStoreShape) fShapes_displ->Store();
		
			/* initialize material */
			if (need_init) 
			{
				/* material pointer */
				MFGPMaterialT* pmat = (*fMFGPMatList)[CurrentElement().MaterialNumber()];
				if (pmat->NeedsPointInitialization())
				{
					/* loop over integration points */
					fShapes_displ->TopIP();
					while (fShapes_displ->NextIP())
						pmat->PointInitialize();
				}
			}
		}
		
		/* finalize storage */
		if (fStoreShape) fShapes_displ->CloseStore();
	}
}

MFGPAssemblyT::MassTypeT MFGPAssemblyT::int2MassTypeT(int i)
{
	if (i == kNoMass)
		return kNoMass;
	else if (i == kConsistentMass)
		return kConsistentMass;
	else if (i == kLumpedMass)
		return kLumpedMass;
	else
		ExceptionT::GeneralFail("MFGPAssemblyT::int2MassTypeT", 
			"could not translate %d", i);	
	return kNoMass;
}

/* extract the list of material parameters */
void MFGPAssemblyT::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{
#pragma unused(all_params)
	mat_params.Clear();
}

/***********************************************************************
 * Protected
 ***********************************************************************/

namespace Tahoe {

/* stream extraction operator */
istream& operator>>(istream& in, MFGPAssemblyT::MassTypeT& mtype)
{
	int i_type;
	in >> i_type;
	switch (i_type)
	{
		case MFGPAssemblyT::kNoMass:
			mtype = MFGPAssemblyT::kNoMass;
			break;
		case MFGPAssemblyT::kConsistentMass:
			mtype = MFGPAssemblyT::kConsistentMass;
			break;
		case MFGPAssemblyT::kLumpedMass:
			mtype = MFGPAssemblyT::kLumpedMass;
			break;
		default:
			cout << "\n MFGPAssemblyT::MassTypeT: unknown type: "
			<< i_type<< endl;
			throw ExceptionT::kBadInputValue;	
	}
	return in;
}

} // namespace Tahoe

void MFGPAssemblyT::SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
#pragma unused(mode)
#pragma unused(flags)
	counts.Dimension(0);
}

void MFGPAssemblyT::SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
#pragma unused(mode)
#pragma unused(flags)
	counts.Dimension(0);
}

void MFGPAssemblyT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	const iArrayT& e_codes, dArray2DT& e_values)
{
#pragma unused(n_codes)
#pragma unused(e_codes)
	n_values.Dimension(0, n_values.MinorDim());
	e_values.Dimension(0, e_values.MinorDim());
}

void MFGPAssemblyT::GenerateOutputLabels( const iArrayT& n_codes, 
	ArrayT<StringT>& n_labels, const iArrayT& e_codes, ArrayT<StringT>& e_labels) const
{
#pragma unused(n_codes)
#pragma unused(e_codes)
	n_labels.Dimension(0);
	e_labels.Dimension(0);
}

/* construct a new material support and return a pointer */
MFGPMatSupportT* MFGPAssemblyT::NewMFGPMatSupport(MFGPMatSupportT* p) const
{
	if (!p) p = new MFGPMatSupportT(fDispl->NumDOF(), fPlast->NumDOF(), NumIP(), NumIP());

	/* MFGPAssemblyT sources */
	p->SetMFGPAssembly(this);
	p->SetElementCards(const_cast<AutoArrayT<ElementCardT>* >(&fElementCards));
	p->SetCurrIP(fShapes_displ->CurrIP());
	//p->SetCurrIP(fShapes_plast->CurrIP());
	p->SetGroup(fDispl->Group());
	//p->SetGroup(fPlast->Group());
	
	/* set pointer to local array */
	p->SetLocalArray(u);
	p->SetLocalArray(u_n);
	//p->SetLocalArray(Du);
	p->SetLocalArray(DDu);
	p->SetLocalArray(lambda);
	p->SetLocalArray(lambda_n);
	return p;
}

/* return a pointer to a new material list */
MFGPMatListT* MFGPAssemblyT::NewMFGPMatList(const StringT& name, int size)
{
#pragma unused(name)
#pragma unused(size)
	return NULL;
}

/* form group contribution to the stiffness matrix and RHS */
void MFGPAssemblyT::RHSDriver_staggered(void)
{
	const char caller[] = "MFGPAssemblyT::RHSDriver_staggered";
	if (fDispl->Group() == fPlast->Group())
		ExceptionT::GeneralFail(caller, "displacement and plastic multiplier group must be different: %d == %d",
			fDispl->Group(), fPlast->Group());

	int curr_group = ElementSupport().CurrentGroup();

 	/* set components and weights */
	double constK = 0.0;
	double constKd = 0.0;
	int formK = fIntegrator->FormK(constK);
	int formKd = fIntegrator->FormKd(constKd);

	/* quick exit */
	if ((formK == 0) || fabs(constK) < kSmall) return;
 	/* has (displacement) body forces */
	int formBody = 0;
	if (fBodySchedule && fBody.Magnitude() > kSmall)
		formBody = 1;

	/* loop over "elements" */ 
	Top();
	while (NextElement())
	{
		double constKe = constK;
		
		/* initialize */
		fKuu=0.0;
		fKulambda=0.0;
		fKlambdau=0.0;
		fKlambdalambda=0.0;
		fFu_int=0.0;
		fFlambda=0.0;
		
		/* global shape function derivatives, local arrays */ 
	    SetGlobalShape();

		/* which field */
	  	//SolverGroup 1 (gets field 1) <-- u (obtained by a rearranged Equation_u)
		if ( curr_group == fDispl->Group() )	
		{
			/* element stiffness */
			if (fabs(constKe) > kSmall)
				FormStiffness(constKe);
				
			/* internal force contribution */	
			if (formKd) FormKd(-constKd);
				
			/** Set displacement LHS */
			fLHS = fKuu;

			/** Compute displacement RHS */
			fKulambda.Multx (del_lambda_vec, fRHS );
			fRHS += fFu_int; 
			fRHS *= -1.0; 

			/** Compute Traction B.C. and Body Forces */
			//Get_Fu_ext ( fFu_ext );
			//fRHS += fFu_ext;
				
			/* add body forces */
			if (formBody) 
			{
				//double density = fBalLinMomMaterial->Retrieve(Iso_MatlT::kDensity);
				double density = 1.0;
				DDu = 0.0;
				AddBodyForce(DDu);
				//FormMa(kConsistentMass, -density, &DDu, NULL);				
			}
			
			/* add to global equations */
			ElementSupport().AssembleLHS ( fDispl->Group(), fLHS, CurrentElement().Equations() );
			ElementSupport().AssembleRHS ( fDispl->Group(), fRHS, CurrentElement().Equations() );
		} //if ( curr_group == fDispl->Group() )

		// SolverGroup 2 (gets field 2) <-- lambda (obtained by a rearranged Equation_eps)
		else if (curr_group == fPlast->Group() )	
		{
			/* element stiffness */
			if (fabs(constKe) > kSmall)
				FormStiffness(constKe);
				
			/* internal force contribution */	
			if (formKd) FormKd(-constKd);
				
			/** Set LHS */
			fLHS = fKlambdalambda;	
		
			/** Compute plasticity RHS  */
			fKlambdau.Multx (del_u_vec, fRHS);
			fRHS += fFlambda; 
			fRHS *= -1.0; 
		
			int e = CurrElementNumber();
			/* plasticity equation numbers */
			const iArrayT& plast_eq = fElementCards_plast[e].Equations();

			/* add to global equations */
			ElementSupport().AssembleLHS (fPlast->Group(), fLHS, plast_eq );
			ElementSupport().AssembleRHS (fPlast->Group(), fRHS, plast_eq );
		}// else if (curr_group == fPlast->Group() )
		else ExceptionT::GeneralFail(caller);
	}// while (NextElement())
}

//---------------------------------------------------------------------
/* form group contribution to the stiffness matrix and RHS */
void MFGPAssemblyT::RHSDriver_monolithic(void)
{
	const char caller[] = "MFGPAssemblyT::RHSDriver_monolithic";
	if (fDispl->Group() != fPlast->Group())
		ExceptionT::GeneralFail(caller, "displacment and plastic multiplier group must be the same: %d != %d",
			fDispl->Group(), fPlast->Group());

	int curr_group = ElementSupport().CurrentGroup();
	
 	/* set components and weights */
	double constK = 0.0;
	double constKd = 0.0;
	int formK = fIntegrator->FormK(constK);
	int formKd = fIntegrator->FormKd(constKd);

	/* quick exit */
	if ((formK == 0) || fabs(constK) < kSmall) return;
	     
 	/* has (displacement) body forces */
	int formBody = 0;
	if (fBodySchedule && fBody.Magnitude() > kSmall)
		formBody = 1;
	
	/* check yield conditions on the nodes and set up yield flags */
	CheckNodalYield();
	
	/* dimension and initialize penalty flags array */
	iArrayT fPenaltyFlags(fNodalYieldFlags.Length()); 
	fPenaltyFlags = 0; 
	
	/* loop over nodes */
	Top();
	while (NextElement())
		if (CurrentElement().Flag() != ElementCardT::kOFF)
		{
			double constKe = constK;
			
			/* initialize */
			fKuu = 0.0;
			fKulambda = 0.0;
			fKlambdau = 0.0;
			fKlambdalambda = 0.0;
			fFu_int = 0.0;
			fFlambda = 0.0;
		
			/* global shape function derivatives, local arrays */ 
	    	SetGlobalShape();
	    	
	    	/* element stiffness */
			if (fabs(constKe) > kSmall)
				FormStiffness(constKe);
				
			/* internal force contribution */	
			if (formKd) FormKd(-constKd);
			
			/* add body force */
			if (formBody) 
			{
				double density = 1.0;
				DDu = 0.0;
				AddBodyForce(DDu);
				
				/* add body force to fRHS */
				fRHS = 0.0;
				//FormMa(kConsistentMass, -density, &DDu, NULL);
				//fFu_int += fRHS;
				fRHS = fFu_int;
			}

			int e = CurrElementNumber();
			/* equations numbers */
			const iArrayT& displ_eq = fElementCards_displ[e].Equations();
			const iArrayT& plast_eq = fElementCards_plast[e].Equations();
			
			/* plastic multiplier field node arrays */
			const iArrayT& nodes_plast = fElementCards_plast[e].NodesU();
			// node numbers are local or global??
			/*
			for (int i = 0; i < nodes_plast.Length(); i++)
				cout << "element # " << e << endl;
				cout << "contributing nodes: " << endl;
				cout << nodes_plast[i]+1 << endl;
				cout << "displacement equation number:" << endl;
				cout << displ_eq << endl << endl;
				cout << "plastic multiplier equation number:" << endl;
				cout << plast_eq << endl;
			*/
			PrintStiffness("before_penalty");
			
			/* impose boundary conditions via penalty method on plastic multiplier; 
			   penalty number added to diagonal terms of Klambdalambda of elastic nodes
			*/ 
			for (int i = 0; i < nodes_plast.Length(); i++) 
			{
					int j = nodes_plast[i]; // local and global number are same ??
					if(fNodalYieldFlags[j] == 0 && fPenaltyFlags[j] == 0) {
						fKlambdalambda(i,i) = fKlambdalambda(i,i) + penalty_num;
						fPenaltyFlags[j] = 1;
					}
			}
			
			PrintStiffness("after_penalty");
			
			/* assemble residuals */
			ElementSupport().AssembleRHS(curr_group, fFu_int, displ_eq);
			ElementSupport().AssembleRHS(curr_group, fFlambda, plast_eq);

			/* assemble components of the tangent */
			ElementSupport().AssembleLHS(curr_group, fKuu, displ_eq);
			ElementSupport().AssembleLHS(curr_group, fKlambdalambda, plast_eq);
			ElementSupport().AssembleLHS(curr_group, fKulambda, displ_eq, plast_eq);
			ElementSupport().AssembleLHS(curr_group, fKlambdau, plast_eq, displ_eq); 
		}// if (CurrentElement().Flag() != ElementCardT::kOFF)
}

/* initialize local arrays */
void MFGPAssemblyT::SetLocalArrays(void)
{
	int nen_displ = NumElementNodes();
	int nen_plast = nen_displ;
	int nsd = NumSD(); 
	int ndof_displ = fDispl->NumDOF();
	int ndof_plast = fPlast->NumDOF();
	
	/* set local arrays for geometry */
	fInitCoords_displ.Dimension(nen_displ, nsd);
	fInitCoords_plast.Dimension(nen_plast, nsd);
	//fCurrCoords_displ.Dimension(n_en_displ, n_sd);
	//fCurrCoords_plast.Dimension(n_en_plast, n_sd);
	
	/* register local arrays */
	ElementSupport().RegisterCoordinates(fInitCoords_displ);
	ElementSupport().RegisterCoordinates(fInitCoords_plast);
	
	/* set local arrays for displacement field */
	u.Dimension (nen_displ, ndof_displ);
	u_n.Dimension (nen_displ, ndof_displ);
	DDu.Dimension (nen_displ, ndof_displ);
	del_u.Dimension (nen_displ, ndof_displ);

	/* register local arrays */
	fDispl->RegisterLocal(u);
	fDispl->RegisterLocal(u_n);	

	/* set local arrays for plastic multiplier field */
	lambda.Dimension (nen_plast, ndof_plast);
	lambda_n.Dimension (nen_plast, ndof_plast);
	del_lambda.Dimension (nen_plast, ndof_plast);
	del_lambda_vec.Dimension (nen_plast);
	
	/* register local arrays */
	fPlast->RegisterLocal(lambda);
	fPlast->RegisterLocal(lambda_n);
}

/* form global shape function derivatives */
void MFGPAssemblyT::SetGlobalShape(void)
{
	/* fetch (initial) coordinates */
	SetLocalX(fInitCoords_displ);
	SetLocalX(fInitCoords_plast);
	
	/* compute shape function derivatives */
	fShapes_displ->SetDerivatives();
	fShapes_plast->SetDerivatives();
	
	/* field node arrays */
	int e = CurrElementNumber();
	const iArrayT& nodes_displ = fElementCards_displ[e].NodesU();
	const iArrayT& nodes_plast = fElementCards_plast[e].NodesU();
	
	u.SetLocal(nodes_displ);
	u_n.SetLocal(nodes_displ);
	lambda.SetLocal(nodes_plast);
	lambda_n.SetLocal(nodes_plast);
	del_u.DiffOf(u, u_n);
	del_lambda.DiffOf(lambda, lambda_n);
	
	/* material dependent local arrays */
	/*
	SetLocalU(u);	
	SetLocalU(u_n);	
	SetLocalU(lambda);	
	SetLocalU(lambda_n);
	*/
	//cout << endl << "displacement" << endl << u << endl; //check displ passed
	
	/* if lambda is negative set lambda to its previous value */
	// should negative lambdas be set to zeros or previous values??
	for (int i = 0; i < lambda.Length(); i++)
	{
		if (lambda[i] < 0.) {
		lambda[i] = lambda_n[i];
		}	
	}
}

/* describe the parameters needed by the interface */
void MFGPAssemblyT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ElementBaseT::DefineParameters(list);

	/* displacement field */
	//already done in ElementBaseT
	//list.AddParameter(ParameterT::Word, "displ_field_name");
	
	/* plastic gradient field */
	list.AddParameter(ParameterT::Word, "plastic_mult_field_name");
	
	ParameterT penalty_num(ParameterT::Double, "penalty_number"); //penalty number for lambda b.c.
	penalty_num.SetDefault(1.0e+20);  //set default value to 1e+20
	list.AddParameter(penalty_num); 
}

/* information about subordinate parameter lists */
void MFGPAssemblyT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ElementBaseT::DefineSubs(sub_list);
	
	/* geometry and integration rule (inline) */
	sub_list.AddSub("element_geometry", ParameterListT::Once, true);
	
	/* optional body force */
	sub_list.AddSub("body_force", ParameterListT::ZeroOrOnce);
	
	/* tractions */
	sub_list.AddSub("natural_bc", ParameterListT::Any);
}

/* return the description of the given inline subordinate parameter list */
void MFGPAssemblyT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	/* geometry and integration rule (inline) */
	if (name == "element_geometry")
	{
		/* choice */
		order = ParameterListT::Choice;
	
		/* element geometries */
		sub_lists.AddSub(GeometryT::ToString(GeometryT::kQuadrilateral));
		sub_lists.AddSub(GeometryT::ToString(GeometryT::kTriangle));
		sub_lists.AddSub(GeometryT::ToString(GeometryT::kHexahedron));
		sub_lists.AddSub(GeometryT::ToString(GeometryT::kTetrahedron));
		sub_lists.AddSub(GeometryT::ToString(GeometryT::kLine));
	}
	
	else
		ElementBaseT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* MFGPAssemblyT::NewSub(const StringT& name) const
{
		/* create non-const this */
		MFGPAssemblyT* non_const_this = const_cast<MFGPAssemblyT*>(this);
 		/* try material list */
	    MFGPMatListT* material_list = non_const_this->NewMFGPMatList(name, 0);
		if (material_list)
			return material_list;
			
		/* try geometry */
		ParameterInterfaceT* geometry = GeometryT::New(name);
		if (geometry)
			return geometry;

		/* body force */
		if (name == "body_force")
		{
			ParameterContainerT* body_force = new ParameterContainerT(name);
	
			/* schedule number */
			body_force->AddParameter(ParameterT::Integer, "schedule");
	
			/* body force vector */
			body_force->AddSub("Double", ParameterListT::OnePlus); 		
		
			return body_force;
		}
		else if (name == "natural_bc") /* traction bc */
		{
			ParameterContainerT* natural_bc = new ParameterContainerT(name);

			natural_bc->AddParameter(ParameterT::Word, "side_set_ID");
			natural_bc->AddParameter(ParameterT::Integer, "schedule");

			ParameterT coord_sys(ParameterT::Enumeration, "coordinate_system");
			coord_sys.AddEnumeration("global", Traction_CardT::kCartesian);
			coord_sys.AddEnumeration( "local", Traction_CardT::kLocal);
			coord_sys.SetDefault(Traction_CardT::kCartesian);
			natural_bc->AddParameter(coord_sys);

			natural_bc->AddSub("DoubleList", ParameterListT::OnePlus); 		
		
			return natural_bc;
		}
		else /* inherited */
			return ElementBaseT::NewSub(name);	
}

/* accept parameter list */
void MFGPAssemblyT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "MFGPAssemblyT::TakeParameterList";
	
	/* resolve geometry before calling inherited method - geometry code
	 * may be needed while reading connectivities */
	const ParameterListT& integration_domain = list.GetListChoice(*this, "element_geometry");
	fGeometryCode_displ = GeometryT::string2CodeT(integration_domain.Name());
	fGeometryCode_plast = fGeometryCode_displ;  // use same geometry for now
	
	/* inherited */
	ElementBaseT::TakeParameterList(list);
	
	/* get the fields */
	/* get displacement field */
	/*
	const StringT& displ_field_name = list.GetParameter("displ_field_name");
	fDispl = ElementSupport().Field(displ_field_name);
	if (!fDispl)
		ExceptionT::GeneralFail(caller, "could not resolve \"%s\" displ_field", 
		displ_field_name.Pointer());
	*/
	const StringT& displ_field_name = list.GetParameter("field_name");
	fDispl = ElementSupport().Field(displ_field_name);
	if (!fDispl)
		ExceptionT::GeneralFail(caller, "could not resolve \"%s\" displ_field", 
		displ_field_name.Pointer());	

	/* get plastic multiplier field */
	const StringT& plastic_mult_field_name = list.GetParameter("plastic_mult_field_name");
	fPlast = ElementSupport().Field(plastic_mult_field_name);
	if (!fPlast)
		ExceptionT::GeneralFail(caller, "could not resolve \"%s\" plastic_mult_field", 
		plastic_mult_field_name.Pointer());
	
	/* construct group communicator */
	const CommunicatorT& comm = ElementSupport().Communicator();
	int color = (NumElements() > 0) ? 1 : CommunicatorT::kNoColor;
	fGroupCommunicator = new CommunicatorT(comm, color, comm.Rank());
	
	/* initialize connectivities */
	fConnectivities_displ.Alias(fConnectivities);
	fConnectivities_plast.Alias(fConnectivities);
	
	/* initialize element cards before calling SetShape() */
	fElementCards_displ.Alias(fElementCards);
	fElementCards_plast = fElementCards_displ; // use same element card for now
	
	/* dimensions */
	fNEEvec.Dimension(NumElementNodes()*NumDOF());
	fDOFvec.Dimension(NumDOF());
	fNodalYieldFlags.Dimension(ElementSupport().NumNodes());
	
	
	/* initialize/set up local arrays */
	SetLocalArrays();
	
	/* construct shape functions */
	fNumIP_displ = integration_domain.GetParameter("num_ip");
	fNumIP_plast = fNumIP_displ;
	SetShape();
	
	/* construct material list */
	ParameterListT mat_params;
	CollectMaterialInfo(list, mat_params);
	if (mat_params.NumLists() > 0)
	{
		fMFGPMatList = NewMFGPMatList(mat_params.Name(), mat_params.NumLists());
		if (!fMFGPMatList) ExceptionT::GeneralFail(caller, "could not construct material list \"%s\"", 
			mat_params.Name().Pointer());
		fMFGPMatList->TakeParameterList(mat_params);
	}
	
	/* get form of tangent */
	GlobalT::SystemTypeT type = TangentType();
	
	/* set form of element stiffness matrix */
	if (type == GlobalT::kSymmetric)
		fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);
	else if (type == GlobalT::kNonSymmetric){
		fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
		fKuu.SetFormat(ElementMatrixT::kNonSymmetric);
		fKulambda.SetFormat(ElementMatrixT::kNonSymmetric);
		fKlambdau.SetFormat(ElementMatrixT::kNonSymmetric);
		fKlambdalambda.SetFormat(ElementMatrixT::kNonSymmetric);
	}
	else if (type == GlobalT::kDiagonal)
		fLHS.SetFormat(ElementMatrixT::kDiagonal);
	
	/* body force */
	const ParameterListT* body_force = list.List("body_force");
	if (body_force) 
	{
		int schedule = body_force->GetParameter("schedule");
		fBodySchedule = ElementSupport().Schedule(--schedule);

		/* body force vector */
		const ArrayT<ParameterListT>& body_force_vector = body_force->Lists();
		if (body_force_vector.Length() != NumDOF())
			ExceptionT::BadInputValue(caller, "body force is length %d not %d",
				body_force_vector.Length(), NumDOF());
		fBody.Dimension(NumDOF());
		for (int i = 0; i < fBody.Length(); i++)
			fBody[i] = body_force_vector[i].GetParameter("value");
	}
		
	/* extract natural boundary conditions */
	TakeNaturalBC(list);
	
	/* get penalty number */
	penalty_num = list.GetParameter("penalty_number");  
	
	// enable the model manager
	ModelManagerT& model = ElementSupport().ModelManager();
	
	/* dimension */
	int nen_displ = NumElementNodes();
	int nen_plast = nen_displ; 
	int ndof_displ = fDispl->NumDOF();
	int ndof_plast = fPlast->NumDOF();
	int nen_displ_x_ndof = nen_displ*ndof_displ;
	fKuu.Dimension 			( nen_displ_x_ndof, nen_displ_x_ndof );
	fKulambda.Dimension		( nen_displ_x_ndof, nen_plast );
	fKlambdau.Dimension		( nen_plast, nen_displ_x_ndof );
	fKlambdalambda.Dimension( nen_plast, nen_plast );
	fFu_int.Dimension 		( nen_displ_x_ndof);
	fFu_ext.Dimension 		( nen_displ_x_ndof );
	fFlambda.Dimension 	( nen_plast );
	
	fKuu_temp.Dimension 			( nen_displ_x_ndof, nen_displ_x_ndof );
	fKulambda_temp.Dimension		( nen_displ_x_ndof, nen_plast );
	fKlambdau_temp.Dimension		( nen_plast, nen_displ_x_ndof );
	fKlambdalambda_temp.Dimension	( nen_plast, nen_plast );
	fFu_int_temp.Dimension 		( nen_displ_x_ndof );
	fFlambda_temp.Dimension 	( nen_plast );
	
	/* streams */
	ofstreamT& out = ElementSupport().Output();

} //MFGPAssemblyT::TakeParameterList(const ParameterListT& list)

/* reset loop */
void MFGPAssemblyT::Top(void)
{
	/* inherited */
	ElementBaseT::Top();

	/* reset element card arrays */
	fElementCards_displ.Top();
	fElementCards_plast.Top();
}

/* current element operations */
bool MFGPAssemblyT::NextElement(void)
{
	/* inherited */
	bool result = ElementBaseT::NextElement();
	
	/* configure for current element */
	if (result)
	{
		/* reset element card arrays */
		fElementCards_displ.Next();
		fElementCards_plast.Next();
	}

	return result;
}

/* return the default number of element nodes */
int MFGPAssemblyT::DefaultNumElemNodes(void) const
{
	switch (fGeometryCode_displ)
	{
		case GeometryT::kLine:
			return 2;
		case GeometryT::kQuadrilateral:
			return 4;
		case GeometryT::kTriangle:
			return 3;
		case GeometryT::kHexahedron:
			return 8;
		case GeometryT::kTetrahedron:
			return 4;
		case GeometryT::kPentahedron:
			return 6;
		default:
			cout << "\n MFGPAssemblyT::DefaultNumElemNodes: unknown geometry code: "
			     << fGeometryCode_displ << endl;
			return 0;
	}
}
//NOTE: needed because ExodusII does not store ANY information about
//      empty element groups, which causes trouble for parallel execution
//      when a partition contains no element from a group.

/* update traction BC data */
void MFGPAssemblyT::SetTractionBC(void)
{
//NOTE: With the possibility of variable global node numbers and
//		and equations, we assume as little as possible here with
//      regard to the validity of the node/equation numbers, requiring
//      only that NodesX in the element cards has the correct global
//      node numbers.

	/* dimensions */
	int ndof = NumDOF();

	/* echo values */
	iArray2DT nd_tmp, eq_tmp;
	for (int i = 0; i < fTractionList.Length(); i++)
	{
		Traction_CardT& BC_card = fTractionList[i];
			
		/* traction element/facet */
		int elem, facet;
		BC_card.Destination(elem, facet);

		/* set global node numbers */
		const iArrayT& loc_nodes = BC_card.LocalNodeNumbers();
		int nnd = loc_nodes.Length();
		
		iArrayT& nodes = BC_card.Nodes();
		nodes.Dimension(nnd);
		nodes.Collect(loc_nodes, fElementCards[elem].NodesX());
		
		/* set global equation numbers */
		iArrayT& eqnos = BC_card.Eqnos();
		eqnos.Dimension(ndof*nnd);
		
		/* get from node manager */
		nd_tmp.Set(1, nnd, nodes.Pointer());
		eq_tmp.Set(1, ndof*nnd, eqnos.Pointer());
		fDispl->SetLocalEqnos(nd_tmp, eq_tmp);
	}

	/* set flag */
	fTractionBCSet = 1;
}

/* add contribution from the body force */
void MFGPAssemblyT::AddBodyForce(LocalArrayT& body_force) const
{
	if (fBodySchedule)
	{
		int ndof = NumDOF();
		int nen = body_force.NumberOfNodes();
		double loadfactor = fBodySchedule->Value();
		double* p = body_force.Pointer();

		for (int i = 0; i < ndof; i++)
		{
			double temp = -fBody[i]*loadfactor;
			for (int j = 0; j < nen; j++)
				*p++ = temp;
		}
	}
}

/* calculate the body force contribution */
void MFGPAssemblyT::FormMa(MassTypeT mass_type, double constM, bool axisymmetric,
	const LocalArrayT* nodal_values, const dArray2DT* ip_values)
{
	const char caller[] = "MFGPAssemblyT::FormMa";

	/* quick exit */
	if (!nodal_values && !ip_values) return;

#if __option(extended_errorcheck)
	/* dimension checks */
	if (nodal_values && 
		fRHS.Length() != nodal_values->Length()) 
			ExceptionT::SizeMismatch(caller);

	if (ip_values &&
		(ip_values->MajorDim() != fShapes_displ->NumIP() ||
		 ip_values->MinorDim() != NumDOF()))
			ExceptionT::SizeMismatch(caller);
#endif

	switch (mass_type)
	{
		case kConsistentMass:
		{
			int ndof = NumDOF();
			//int  nen = MeshFreeElementSupportT::NumElementNodes();
			int  nen = NumElementNodes();
			int  nun = nodal_values->NumberOfNodes();

			const double* Det    = fShapes_displ->IPDets();
			const double* Weight = fShapes_displ->IPWeights();

			if (axisymmetric)
			{
				const LocalArrayT& coords = fShapes_displ->Coordinates();
				fShapes_displ->TopIP();
				while (fShapes_displ->NextIP())
				{
					/* compute radius */
					const double* NaX = fShapes_displ->IPShapeX();
					const double* x_r = coords(0); /* r is x-coordinate */
					double r = 0.0;
					for (int a = 0; a < nen; a++)
						r += (*NaX++)*(*x_r++);
				
					/* interpolate nodal values to ip */
					if (nodal_values)
						fShapes_displ->InterpolateU(*nodal_values, fDOFvec);
					
					/* ip sources */
					if (ip_values)
						fDOFvec -= (*ip_values)(fShapes_displ->CurrIP());

					/* accumulate in element residual force vector */				
					double*	res      = fRHS.Pointer();
					const double* Na = fShapes_displ->IPShapeU();
				
					double temp = 2.0*Pi*r*constM*(*Weight++)*(*Det++);				
					for (int lnd = 0; lnd < nun; lnd++)
					{
						double temp2 = temp*(*Na++);
						double* pacc = fDOFvec.Pointer();
						for (int dof = 0; dof < ndof; dof++)			
							*res++ += temp2*(*pacc++);
					}
				}
			}
			else /* not axisymmetric */
			{
				fShapes_displ->TopIP();
				while (fShapes_displ->NextIP())
				{					
					/* interpolate nodal values to ip */
					if (nodal_values)
						fShapes_displ->InterpolateU(*nodal_values, fDOFvec);
					
					/* ip sources */
					if (ip_values)
						fDOFvec -= (*ip_values)(fShapes_displ->CurrIP());

					/* accumulate in element residual force vector */				
					double*	res      = fRHS.Pointer();
					const double* Na = fShapes_displ->IPShapeU();
				
					double temp = constM*(*Weight++)*(*Det++);				
					for (int lnd = 0; lnd < nun; lnd++)
					{
						double temp2 = temp*(*Na++);
						double* pacc = fDOFvec.Pointer();
						for (int dof = 0; dof < ndof; dof++)			
							*res++ += temp2*(*pacc++);
					}
				}
			}
			break;
		}	
		case kLumpedMass:
		{
			fLHS = 0.0; //hope there's nothing in there!
			FormMass(kLumpedMass, constM, axisymmetric);

			/* init nodal values */
			if (nodal_values)
				nodal_values->ReturnTranspose(fNEEvec);
			else 
			{
				ExceptionT::GeneralFail(caller, "expecting nodal values for lumped mass");
			}
				
//TEMP - what to do with ip values?
if (ip_values)
	ExceptionT::GeneralFail(caller, "lumped mass not implemented for ip sources");

			double* pAcc = fNEEvec.Pointer();
			double* pRes = fRHS.Pointer();
			int     massdex = 0;
			
			int nee = nodal_values->Length();
			for (int i = 0; i < nee; i++)
			{
				*pRes++ += (*pAcc++)*fLHS(massdex,massdex);
				massdex++;
			}
			
			break;
		}
		
	}//switch (mass_type)
}

/* form the element mass matrix */
void MFGPAssemblyT::FormMass(int mass_type, double constM, bool axisymmetric)
{
#if __option(extended_errorcheck)
	if (u.Length() != fLHS.Rows()) throw ExceptionT::kSizeMismatch;
#endif

	switch (mass_type)
	{
		case kNoMass:			/* no mass matrix */
		
			break;
		
		case kConsistentMass:	/* consistent mass	*/
		{
			// integration of the element mass is done
			// in the reference configuration since density
			// is mass/(undeformed volume)
			const double* Det    = fShapes_displ->IPDets();
			const double* Weight = fShapes_displ->IPWeights();
			
			//int nen = MeshFreeElementSupportT::NumElementNodes();
			int nen = NumElementNodes();
			int nun = u.NumberOfNodes();
			int ndof = NumDOF();
			
			/* matrix form */
			int a = 0, zero = 0;
			int& b_start = (fLHS.Format() == ElementMatrixT::kSymmetricUpper) ? a : zero;
			
			if (axisymmetric)
			{
				const LocalArrayT& coords = fShapes_displ->Coordinates();
				fShapes_displ->TopIP();	
				while ( fShapes_displ->NextIP() )
				{
					/* compute radius */
					const double* NaX = fShapes_displ->IPShapeX();
					const double* x_r = coords(0); /* r is x-coordinate */
					double r = 0.0;
					for (a = 0; a < nen; a++)
						r += (*NaX++)*(*x_r++);
				
					double temp = 2.0*Pi*r*constM*(*Weight++)*(*Det++);
					const double* Na = fShapes_displ->IPShapeU();		
					for (a = 0; a < nun; a++)
						for (int i = 0; i < ndof; i++)
						{
							int p = a*ndof + i;
							
							/* upper triangle only */
							for (int b = b_start; b < nun; b++) //TEMP - interpolate at the same time?
								for (int j = 0; j < ndof; j++)
									if(i == j) 
									{
										int q = b*ndof + j;
										fLHS(p,q) += temp*Na[a]*Na[b];
									}
						}
				}
			}			
			else /* not axisymmetric */
			{
				fShapes_displ->TopIP();	
				while ( fShapes_displ->NextIP() )
				{
					double temp = constM*(*Weight++)*(*Det++);
					const double* Na = fShapes_displ->IPShapeU();		
					for (a = 0; a < nun; a++)
						for (int i = 0; i < ndof; i++)
						{
							int p = a*ndof + i;
							
							/* upper triangle only */
							for (int b = b_start; b < nun; b++)
								for (int j = 0; j < ndof; j++)
									if(i == j) 
									{
										int q = b*ndof + j;
										fLHS(p,q) += temp*Na[a]*Na[b];
									}
						}
				}
			}
			break;
		}

		case kLumpedMass:	/* lumped mass */
		{
			//int nen = MeshFreeElementSupportT::NumElementNodes();
			int nen = NumElementNodes();
			int nun = u.NumberOfNodes();
			int ndof = NumDOF();

		    double dsum   = 0.0;
		    double totmas = 0.0;
		    fNEEvec = 0.0;

			const double* Det    = fShapes_displ->IPDets();
			const double* Weight = fShapes_displ->IPWeights();

			/* total mass and diagonal sum */
			if (axisymmetric)
			{
				const LocalArrayT& coords = fShapes_displ->Coordinates();
				fShapes_displ->TopIP();
				while (fShapes_displ->NextIP()) 
				{

					/* compute radius */
					const double* NaX = fShapes_displ->IPShapeX();
					const double* x_r = coords(0); /* r is x-coordinate */
					double r = 0.0;
					for (int a = 0; a < nen; a++)
						r += (*NaX++)*(*x_r++);

					double temp1     = 2.0*Pi*r*constM*(*Weight++)*(*Det++);
					const double* Na = fShapes_displ->IPShapeU();
					totmas += temp1;
					for (int lnd = 0; lnd < nun; lnd++) 
					{
						double temp2 = temp1*Na[lnd]*Na[lnd];
						dsum += temp2;
						fNEEvec[lnd] += temp2;
					}
				}
			}
			else /* not axisymmetric */
			{
				fShapes_displ->TopIP();
				while (fShapes_displ->NextIP()) 
				{
					double temp1     = constM*(*Weight++)*(*Det++);
					const double* Na = fShapes_displ->IPShapeU();
					totmas += temp1;
					for (int lnd = 0; lnd < nun; lnd++) 
					{
						double temp2 = temp1*Na[lnd]*Na[lnd];
						dsum += temp2;
						fNEEvec[lnd] += temp2;
					}
				}
			}	
				
			/* scale diagonal to conserve total mass */
			double diagmass = totmas/dsum;
			
			/* lump mass onto diagonal */
			double* pmass = fLHS.Pointer();
			int inc = fLHS.Rows() + 1;
			for (int lnd = 0; lnd < nun; lnd++)
			{
				double temp = diagmass*fNEEvec[lnd];
				for (int ed = 0; ed < ndof; ed++)
				{
					*pmass += temp;
					pmass += inc;	
				}
			}
			break;
		}			
		default:
			ExceptionT::BadInputValue("Elastic::FormMass", "unknown mass matrix code");
	}//switch (mass_type)
}

void MFGPAssemblyT::EchoBodyForce(ifstreamT& in, ostream& out)
{
	/* schedule number and body force vector */
	int n_sched;
	in >> n_sched >> fBody;		
	n_sched--;

	/* no LTf => no body force */
	if (n_sched < 0) 
		fBody = 0.0;
	else
	{
		fBodySchedule = ElementSupport().Schedule(n_sched);
		if (!fBodySchedule)
			ExceptionT::BadInputValue("MFGPAssemblyT::EchoBodyForce", 
				"could not resolve schedule %d", n_sched+1);
	}
	
	out << "\n Body force vector:\n";
	out << " Body force load-time function number. . . . . . = " << n_sched + 1<< '\n';
	out << " Body force vector components:\n";
	for (int j = 0 ; j < NumDOF(); j++)
	{
		out << "   x[" << j+1 << "] direction. . . . . . . . . . . . . . . . = ";
		out << fBody[j] << '\n';
	}
	out.flush();   	   	
}

/* extract natural boundary condition information */
void MFGPAssemblyT::TakeNaturalBC(const ParameterListT& list)
{
	const char caller[] = "MFGPAssemblyT::TakeTractionBC";

	int num_natural_bc = list.NumLists("natural_bc");
	if (num_natural_bc > 0)
	{
		/* model manager */
		ModelManagerT& model = ElementSupport().ModelManager();
	
		/* temp space */
		ArrayT<StringT> block_ID(num_natural_bc);
	    ArrayT<iArray2DT> localsides(num_natural_bc);
	    iArrayT LTf(num_natural_bc);
	    ArrayT<Traction_CardT::CoordSystemT> coord_sys(num_natural_bc);
	    ArrayT<dArray2DT> values(num_natural_bc);

	    /* nodes on element facets */
	    iArrayT num_facet_nodes;
	    fShapes_displ->NumNodesOnFacets(num_facet_nodes);
	    
	    /* loop over natural BC's */
	    int tot_num_sides = 0;
	    for (int i = 0; i < num_natural_bc; i++) 
	   	{
	    	const ParameterListT& natural_bc = list.GetList("natural_bc", i);
	    
	    	/* side set */
	    	const StringT& ss_ID = natural_bc.GetParameter("side_set_ID");
			localsides[i] = model.SideSet(ss_ID);
			int num_sides = localsides[i].MajorDim();
			tot_num_sides += num_sides;
			if (num_sides > 0)
			{
				block_ID[i] = model.SideSetGroupID(ss_ID);
				LTf[i] = natural_bc.GetParameter("schedule");
				coord_sys[i] = Traction_CardT::int2CoordSystemT(natural_bc.GetParameter("coordinate_system"));

				/* switch to elements numbering within the group */
				iArray2DT& side_set = localsides[i];
				iArrayT elems(num_sides);
				side_set.ColumnCopy(0, elems);
				BlockToGroupElementNumbers(elems, block_ID[i]);
				side_set.SetColumn(0, elems);

				/* all facets in set must have the same number of nodes */
				int num_nodes = num_facet_nodes[side_set(0,1)];
				for (int f = 0; f < num_sides; f++)
					if (num_facet_nodes[side_set(f,1)] != num_nodes)
						ExceptionT::BadInputValue(caller, "faces side set \"%s\" have different numbers of nodes",
							ss_ID.Pointer());

				/* read traction nodal values */
				dArray2DT& nodal_values = values[i];
				nodal_values.Dimension(num_nodes, NumDOF());
				int num_traction_vectors = natural_bc.NumLists("DoubleList");
				if (num_traction_vectors != 1 && num_traction_vectors != num_nodes)
					ExceptionT::GeneralFail(caller, "expecting 1 or %d vectors not %d",
						num_nodes, num_traction_vectors);
						
				/* constant over the face */
				if (num_traction_vectors == 1) 
				{
					const ParameterListT& traction_vector = natural_bc.GetList("DoubleList");
					int dim = traction_vector.NumLists("Double");
					if (dim != NumDOF())
						ExceptionT::GeneralFail(caller, "expecting traction vector length %d not %d",
							NumDOF(), dim);

					/* same for all face nodes */
					for (int f = 0; f < NumDOF(); f++) 
					{
						double t = traction_vector.GetList("Double", f).GetParameter("value");
						nodal_values.SetColumn(f, t);
					}
				}
				else
				{
					/* read separate vector for each face node */
					dArrayT t;
					for (int f = 0; f < num_nodes; f++) 
					{
						const ParameterListT& traction_vector = natural_bc.GetList("DoubleList", f);
					int dim = traction_vector.NumLists("Double");
						if (dim != NumDOF())
							ExceptionT::GeneralFail(caller, "expecting traction vector length %d not %d",
								NumDOF(), dim);

						nodal_values.RowAlias(f, t);
						for (int j = 0; j < NumDOF(); j++)
							t[j] = traction_vector.GetList("Double", j).GetParameter("value");
					}
				}// if (num_traction_vectors == 1)
			}// if (num_sides > 0)
	    }// for (int i = 0; i < num_natural_bc; i++)
#pragma message("OK with empty side sets?")

		/* allocate all traction BC cards */
	    fTractionList.Dimension(tot_num_sides);

	    /* correct numbering offset */
	    LTf--;

		/* define traction cards */
		if (tot_num_sides > 0)
		{
			iArrayT loc_node_nums;
			int dex = 0;
			for (int i = 0; i < num_natural_bc; i++)
			{
				/* set traction BC cards */
				iArray2DT& side_set = localsides[i];
				int num_sides = side_set.MajorDim();
				for (int j = 0; j < num_sides; j++)
				{					
					/* get facet local node numbers */
					fShapes_displ->NodesOnFacet(side_set(j, 1), loc_node_nums);
					
					/* set and echo */
					fTractionList[dex++].SetValues(ElementSupport(), side_set(j,0), side_set (j,1), LTf[i],
						 coord_sys[i], loc_node_nums, values[i]);
				}
			}
		}

		/* check coordinate system specifications */
		if (NumSD() != NumDOF())
			for (int i = 0; i < fTractionList.Length(); i++)
				if (fTractionList[i].CoordSystem() != Traction_CardT::kCartesian)
					ExceptionT::BadInputValue(caller, "coordinate system must be Cartesian if (nsd != ndof) for card %d", i+1);
	}// if (num_natural_bc > 0)
}

/* compute contribution to RHS from traction BC's */
void MFGPAssemblyT::ApplyTractionBC(void)
{
	if (fTractionList.Length() > 0)
	{		
		/* dimensions */
		int nsd = NumSD();
		int ndof = NumDOF();
		bool is_axi = Axisymmetric();
		if (is_axi && nsd != 2) ExceptionT::GeneralFail();
	
		/* update equation numbers */
		if (!fTractionBCSet) SetTractionBC();
	
		/* force vector */
		dArrayT rhs;
		VariArrayT<double> rhs_man(25, rhs);
		
		/* local coordinates */
		LocalArrayT coords(LocalArrayT::kInitCoords);
		VariLocalArrayT coord_man(25, coords, nsd);
		ElementSupport().RegisterCoordinates(coords);
		
		/* nodal tractions */
		LocalArrayT tract(LocalArrayT::kUnspecified);
		VariLocalArrayT tract_man(25, tract, ndof);

		/* integration point tractions */
		dArray2DT ip_tract;
		nVariArray2DT<double> ip_tract_man(25, ip_tract, ndof);
		dArrayT tract_loc, tract_glb(ndof);
		dMatrixT Q(ndof);
		dArrayT ip_coords(2);
		
		/* Jacobian of the surface mapping */
		dMatrixT jacobian(nsd, nsd-1);
		double Pi2 = 2.0*Pi;
		for (int i = 0; i < fTractionList.Length(); i++)
		{
			const Traction_CardT& BC_card = fTractionList[i];

			/* dimension */
			const iArrayT& nodes = BC_card.Nodes();
			int nnd = nodes.Length();
			rhs_man.SetLength(nnd*ndof, false);
			coord_man.SetNumberOfNodes(nnd);
			tract_man.SetNumberOfNodes(nnd);
			
			/* local coordinates */
			coords.SetLocal(nodes);

			/* nodal traction vectors: (ndof x nnd) */
			BC_card.CurrentValue(tract);
			
			/* BC destination */
			int elem, facet;
			BC_card.Destination(elem, facet);
			
			/* boundary shape functions */
			//const ParentDomainT& surf_shape = ShapeFunction().FacetShapeFunction(facet);
			const ParentDomainT& surf_shape = fShapes_displ->FacetShapeFunction(facet);
			int nip = surf_shape.NumIP();
			
			/* all ip tractions: (nip x ndof) */
			ip_tract_man.SetMajorDimension(nip, false);
			surf_shape.Interpolate(tract, ip_tract);

			/* traction vector coordinate system */
			if (BC_card.CoordSystem() == Traction_CardT::kCartesian)
			{
				/* integrate */			
				rhs = 0.0;
				const double* w = surf_shape.Weight();
				for (int j = 0; j < nip; j++)
				{
					/* coordinate mapping */
					surf_shape.DomainJacobian(coords, j, jacobian);
					double detj = surf_shape.SurfaceJacobian(jacobian);
	
					/* ip weight */
					double jwt = detj*w[j];
					if (is_axi) 
					{
						surf_shape.Interpolate(coords, ip_coords, j);
						jwt *= Pi2*ip_coords[0];
					}
					
					/* ip traction */
					const double* tj = ip_tract(j);
					
					/* accumulate */
					for (int l = 0; l < ndof; l++)
					{
						/* nodal shape function */
						const double* Na = surf_shape.Shape(j);
					
						double* prhs = rhs.Pointer(l);
						double  fact = jwt*(*tj++);
						for (int k = 0; k < nnd; k++)
						{
							*prhs += fact*(*Na++);
							prhs += ndof;
						}
					}				
				}
			}
			else if (BC_card.CoordSystem() == Traction_CardT::kLocal)
			{
				/* integrate */			
				rhs = 0.0;
				const double* w = surf_shape.Weight();
				for (int j = 0; j < nip; j++)
				{
					/* coordinate mapping */
					surf_shape.DomainJacobian(coords, j, jacobian);
					double detj = surf_shape.SurfaceJacobian(jacobian, Q);
	
					/* ip weight */
					double jwt = detj*w[j];
					if (is_axi) 
					{
						surf_shape.Interpolate(coords, ip_coords, j);
						jwt *= Pi2*ip_coords[0];
					}
					
					/* transform ip traction out of local frame */
					ip_tract.RowAlias(j, tract_loc);
					Q.Multx(tract_loc, tract_glb);

					/* ip traction */
					const double* tj = tract_glb.Pointer();
					
					/* accumulate */
					for (int l = 0; l < ndof; l++)
					{
						/* nodal shape function */
						const double* Na = surf_shape.Shape(j);
					
						double* prhs = rhs.Pointer(l);
						double  fact = jwt*(*tj++);
						for (int k = 0; k < nnd; k++)
						{
							*prhs += fact*(*Na++);
							prhs += ndof;
						}
					}				
				}
			}//if (BC_card.CoordSystem() == Traction_CardT::kCartesian)
			else
				throw ExceptionT::kGeneralFail;

			/* assemble */
			ElementSupport().AssembleRHS(fDispl->Group(), rhs, BC_card.Eqnos());
		}// for (int i = 0; i < fTractionList.Length(); i++)
	}// if (fTractionList.Length() > 0)
}

/* write all current element information to the stream */
void MFGPAssemblyT::CurrElementInfo(ostream& out) const
{
	/* inherited */
	ElementBaseT::CurrElementInfo(out);
	dArray2DT temp;
	temp.Dimension(fInitCoords_displ.NumberOfNodes(), fInitCoords_displ.MinorDim());
	
	out <<   " initial coords:\n";
	temp.Dimension(fInitCoords_displ.NumberOfNodes(), fInitCoords_displ.MinorDim());
	fInitCoords_displ.ReturnTranspose(temp);
	temp.WriteNumbered(out);

	out <<   " displacements:\n";
	temp.Dimension(u.NumberOfNodes(), u.MinorDim());
	u.ReturnTranspose(temp);
	temp.WriteNumbered(out);
	
	out <<   " plastic multiplier:\n";
	temp.Dimension(lambda.NumberOfNodes(), lambda.MinorDim());
	lambda.ReturnTranspose(temp);
	temp.WriteNumbered(out);
}

/* check material outputs - return true if OK */
bool MFGPAssemblyT::CheckMaterialOutput(void) const
{
	/* check compatibility of output */
	if (fMFGPMatList && fMFGPMatList->Length() > 1)
	{
		/* check compatibility of material outputs */
		bool OK = true;
		int i, j;
		for (i = 0; OK && i < fMFGPMatList->Length(); i++)
		{
			MFGPMaterialT* m_i = (*fMFGPMatList)[i];
			for (j = i+1; OK && j < fMFGPMatList->Length(); j++)
			{
				MFGPMaterialT* m_j = (*fMFGPMatList)[j];
				OK = MFGPMaterialT::CompatibleOutput(*m_i, *m_j);
			}
		}
		i--; j--;
			
		/* output not compatible */
		if (!OK)	
		{
#pragma message("report names")
			cout << "\n MFGPAssemblyT::CheckMaterialOutput: incompatible output\n"
			    <<    "     between materials " << i+1 << " and " << j+1 << ":\n";
//			(*fMFGPMatList)[i]->PrintName(cout);
			cout << '\n';
//			(*fMFGPMatList)[j]->PrintName(cout);
			cout << endl;
			return false;
		}
	}
	
	/* no problems */
	return true;
}

/***********************************************************************
 * Private
 ***********************************************************************/
 /* write stiffness matrices to the output files 
 before or after penalty number is added */
void MFGPAssemblyT::PrintStiffness(StringT before_after) const
{
	StringT file_name; int e = CurrElementNumber(); 
	if(before_after == "before_penalty") {			
		/* one output for each element */
		file_name = "C:/Documents and Settings/kyonten/My Documents/tahoe_xml/bp_stiffness.";
		//file_name = "C:/Documents and Settings/Administrator/My Documents/tahoe/bp_stiffness.";
	}
	else if (before_after == "after_penalty") {
		file_name = "C:/Documents and Settings/kyonten/My Documents/tahoe_xml/ap_stiffness.";
		//file_name = "C:/Documents and Settings/Administrator/My Documents/tahoe/ap_stiffness.";
	}
	file_name.Append(e); // append element number to output string
	file_name.Append(".txt");
	ofstream output(file_name);
	/* print displ_eq and plast_eq numbers */
	if (!output) {
		cout << "Error opening output file" << endl;
	}
	output << "element # " << e << endl;
			
		/* print element stiffness matrices */
		output << "accumulated stiffness matrices " << endl;
		for (int i = 0; i < fKuu.Rows(); i++)
		{
			for (int j = 0; j < fKuu.Cols(); j++)
				output << "KUU("<< i << ","<< j <<"): " << fKuu(i,j) << endl;
			output << endl;	
			for (int j = 0; j < fKlambdalambda.Rows(); j++) 
				output << "KULambda("<< i << ","<< j <<"): " << fKulambda(i,j) << endl;
			output << endl;
		}
				
		for (int i = 0; i < fKlambdalambda.Rows(); i++)
		{
			for (int j = 0; j < fKlambdalambda.Cols(); j++)
				output << "KLambdaLambda("<< i << ","<< j <<"): " << fKlambdalambda(i,j) << endl;
			cout << endl;
			for (int j = 0; j < fKlambdalambda.Cols(); j++)
				output << "KLambdaU("<< i << ","<< j <<"): " << fKlambdau(i,j) << endl;
			output << endl;
		}	
		output.close();
}
		