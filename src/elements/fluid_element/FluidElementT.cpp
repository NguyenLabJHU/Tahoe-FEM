/* $Header: /home/regueiro/tahoe_cloudforge_repo_snapshots/development/src/elements/fluid_element/FluidElementT.cpp,v 1.5 2006-07-13 17:55:11 a-kopacz Exp $ */
/* created: a-kopacz (07/04/2006) */
#include "FluidElementT.h"

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include "ifstreamT.h"
#include "eIntegratorT.h"
#include "ShapeFunctionT.h"
#include "ParameterContainerT.h"

/* materials */
#include "FluidMaterialT.h"
#include "FluidMatSupportT.h"
#include "FluidMatListT.h"

const double Pi = acos(-1.0);

using namespace Tahoe;

/* initialize static data */
const int FluidElementT::NumNodalOutputCodes = 4;
static const char* NodalOutputNames[] = {
  "coordinates",
  "velocities",
  "accelerations",
  "pressures"};

const int FluidElementT::NumElementOutputCodes = 0;
static const char* ElementOutputNames[] = {
  "NONE"};

const int FluidElementT::NumStabParamCodes = 1;
static const char* StabParamNames[] = {
  "tau_m_is_tau_c"}; /* T.E. Tezduyar, Stabilized Finite Element Formulations for Incompressible Flow Computations, Adv. Appl. Mech. 28(1991) 1-44. */
  
/* parameters */
const int FluidElementT::kPressureNDOF = 1;

/* constructor */
FluidElementT::FluidElementT(const ElementSupportT& support):
  ContinuumElementT(support),
  /*dofs: vels and press*/
  fLocLastDisp(LocalArrayT::kLastDisp),
  fLocVel(LocalArrayT::kVel), /* velocity */
  fLocCurVel(LocalArrayT::kUnspecified),  
  fLocOldVel(LocalArrayT::kUnspecified),
  /* accelaration */
  fLocCurAcc(LocalArrayT::kUnspecified),
  /* pressure */   
  fLocPrs(LocalArrayT::kUnspecified)
{
  SetName("incompressible_newtonian_fluid_element");
}

/* destructor */
FluidElementT::~FluidElementT(void)
{
	delete fFluidMatSupport;
}

/* compute nodal force */
void FluidElementT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
  WriteCallLocation("AddNodalForce: not implemented"); //DEBUG  
  //not implemented
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
}

/* returns the energy as defined by the derived class types */
double FluidElementT::InternalEnergy(void)
{
  WriteCallLocation("InternalEnergy: not implemented"); //DEBUG
  //not implemented
  double energy = 0.0;
  return energy;
}

/** compute specified output parameter and send for smoothing */
void FluidElementT::SendOutput(int kincode)
{
  WriteCallLocation("SendOutput"); //DEBUG 
  /* output flags */
  iArrayT flags(fNodalOutputCodes.Length());

  /* set flags to get desired output */
  flags = IOBaseT::kAtNever;
  
  switch (kincode)
  {
    case iNodalCrd:
      flags[iNodalCrd] = 1;
    break;
    case iNodalVel:
      flags[iNodalVel] = 1;
    break;
    case iNodalAcc:
      flags[iNodalAcc] = 1;
    break;
    case iNodalPrs:
      flags[iNodalPrs] = 1;
    break;
    default:
      cout << "\n FluidElementT::SendOutput: invalid output code: "
        << kincode << endl;
    }

  /* number of output values */
  iArrayT n_counts;
  SetNodalOutputCodes(IOBaseT::kAtInc, flags, n_counts);

  /* reset averaging workspace */
  ElementSupport().ResetAverage(n_counts.Sum());

  /* no element output */
  iArrayT e_counts(fElementOutputCodes.Length());
  e_counts = 0;

  /* generate output */
  dArray2DT n_values, e_values;
  ComputeOutput(n_counts, n_values, e_counts, e_values);
}
/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct a new material support and return a pointer */
MaterialSupportT* FluidElementT::NewMaterialSupport(MaterialSupportT* p) const
{
	/* allocate */
	if (!p) p = new FluidMatSupportT(NumDOF(), NumIP());

	/* inherited initializations */
	ContinuumElementT::NewMaterialSupport(p);
	
	/* set FluidMatSupportT fields */
	FluidMatSupportT* ps = TB_DYNAMIC_CAST(FluidMatSupportT*, p);
	if (ps) {
		ps->SetContinuumElement(this);
		ps->SetField(&fVel_list, &fPres_list);
		ps->SetGradient(&fGradVel_list, &fGradPres_list);
	}

	return p;
}

/* return a pointer to a new material list */
MaterialListT* FluidElementT::NewMaterialList(const StringT& name, int size)
{
	/* no match */
	if (name != "fluid_material") 
		return NULL;

	if (size > 0)
	{
		if (!fFluidMatSupport) {
			fFluidMatSupport = TB_DYNAMIC_CAST(FluidMatSupportT*, NewMaterialSupport());
			if (!fFluidMatSupport) ExceptionT::GeneralFail("FluidElementT::NewMaterialList");
		}

		return new FluidMatListT(size, *fFluidMatSupport);
	}
	else
		return new FluidMatListT;

}

/* driver for calculating output values */
void FluidElementT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
  const iArrayT& e_codes, dArray2DT& e_values)
{
  WriteCallLocation("ComputeOutput"); //DEBUG 
  /* number of output values */
  int n_out = n_codes.Sum();
  int e_out = e_codes.Sum();

  /* nothing to output */
  if (n_out == 0 && e_out == 0) return;

  // not implemented
  return;
}  

/* current element operations */
bool FluidElementT::NextElement(void)
{
	/* inherited */
	bool result = ContinuumElementT::NextElement();
	
	/* get material pointer */
	if (result)
	{
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[CurrentElement().MaterialNumber()];
	
		/* cast is safe since class contructs materials list */
		fCurrMaterial = (FluidMaterialT*) pcont_mat;
	}
	
	return result;
}
/** initialize local arrays */
void FluidElementT::SetLocalArrays(void)
{
  WriteCallLocation("SetLocalArrays"); //DEBUG

  /* inherited */
  ContinuumElementT::SetLocalArrays();

  /* allocate */
  int nen = NumElementNodes();
  int ndof  = NumDOF();

  fLocLastDisp.Dimension(nen, ndof);

  /* set source */
  Field().RegisterLocal(fLocLastDisp);

  if (fIntegrator->Order() > 0)
    Field().RegisterLocal(fLocVel);
  
  const double* p_cur_vel = fLocDisp(0);  
  fLocCurVel.Alias(nen, ndof - 1, p_cur_vel);

  const double* p_cur_pres = fLocDisp(ndof-1);
  fLocPrs.Alias(nen, 1, p_cur_pres);
  
  const double* p_old_vel = fLocLastDisp(0);  
  fLocOldVel.Alias(nen, ndof - 1, p_old_vel);

  const double* p_cur_acc = fLocVel(0);  
  fLocCurAcc.Alias(nen, ndof - 1, p_cur_acc);
}

/* form shape functions and derivatives */
void FluidElementT::SetGlobalShape(void)
{
	/* inherited */
	ContinuumElementT::SetGlobalShape();


	/* material dependent local arrays */
	SetLocalU(fLocDisp);	
	SetLocalU(fLocLastDisp);	

	/* have velocity */
	if (fLocVel.IsRegistered())
		SetLocalU(fLocVel);
	
	/* loop over integration points */
	for (int i = 0; i < NumIP(); i++)
	{
		/* velocity gradient */
		dArrayT& vel = fVel_list[i];
		dMatrixT& gradV = fGradVel_list[i];
		dMatrixT gradP;
		gradP.Set(1,NumSD(), fGradPres_list[i].Pointer());
    dArrayT pres;
    double* p = fPres_list.Pointer();
    pres.Set(1,p+i);
		fShapes->GradU(fLocCurVel, gradV, i);
		fShapes->GradU(fLocPrs, gradP, i);
		fShapes->InterpolateU(fLocPrs, pres, i);
		fShapes->InterpolateU(fLocCurVel, vel, i);
	}
}

void FluidElementT::Set_B(const dArray2DT& DNa, dMatrixT& B) const
{
#if __option(extended_errorcheck)
	if (B.Rows() != dSymMatrixT::NumValues(DNa.MajorDim()) ||
	    B.Cols() != DNa.Length())
	    throw ExceptionT::kSizeMismatch;
#endif

	int nnd = DNa.MinorDim();
	double* pB = B.Pointer();

	/* 1D */
	if (DNa.MajorDim() == 1)
	{
		const double* pNax = DNa(0);
		for (int i = 0; i < nnd; i++)
			*pB++ = *pNax++;
	}
	/* 2D */
	else if (DNa.MajorDim() == 2)
	{
		const double* pNax = DNa(0);
		const double* pNay = DNa(1);
		for (int i = 0; i < nnd; i++)
		{
			/* see Hughes (2.8.20) */
			*pB++ = *pNax;
			*pB++ = 0.0;
			*pB++ = *pNay;

			*pB++ = 0.0;
			*pB++ = *pNay++;
			*pB++ = *pNax++;
		}
	}
	/* 3D */
	else		
	{
		const double* pNax = DNa(0);
		const double* pNay = DNa(1);
		const double* pNaz = DNa(2);
		for (int i = 0; i < nnd; i++)
		{
			/* see Hughes (2.8.21) */
			*pB++ = *pNax;
			*pB++ = 0.0;
			*pB++ = 0.0;
			*pB++ = 0.0;
			*pB++ = *pNaz;
			*pB++ = *pNay;

			*pB++ = 0.0;
			*pB++ = *pNay;
			*pB++ = 0.0;
			*pB++ = *pNaz;
			*pB++ = 0.0;
			*pB++ = *pNax;

			*pB++ = 0.0;
			*pB++ = 0.0;
			*pB++ = *pNaz++;
			*pB++ = *pNay++;
			*pB++ = *pNax++;
			*pB++ = 0.0;
		}
	}
}

/* set the \e B matrix for axisymmetric deformations */
void FluidElementT::Set_B_axi(const dArrayT& Na, const dArray2DT& DNa, 
	double r, dMatrixT& B) const
{
#if __option(extended_errorcheck)
	if (B.Rows() != 4 || /* (number of stress 2D components) + 1 */
	    B.Cols() != DNa.Length() ||
    DNa.MajorDim() != 2 ||
    DNa.MinorDim() != Na.Length())
			ExceptionT::SizeMismatch("SolidElementT::Set_B_axi");
#endif

	int nnd = DNa.MinorDim();
	double* pB = B.Pointer();

	const double* pNax = DNa(0);
	const double* pNay = DNa(1);
	const double* pNa  = Na.Pointer();
	for (int i = 0; i < nnd; i++)
	{
		/* see Hughes (2.12.8) */
		*pB++ = *pNax;
		*pB++ = 0.0;
		*pB++ = *pNay;
		*pB++ = *pNa++/r; /* about y-axis: u_r = u_x */

		*pB++ = 0.0;
		*pB++ = *pNay++;
		*pB++ = *pNax++;
		*pB++ = 0.0;
	}
}

/* form the element mass matrix */
void FluidElementT::FormMass(MassTypeT mass_type, double constM, bool axisymmetric, const double* ip_weight)
{
#if __option(extended_errorcheck)
	if (fLocDisp.Length() != fLHS.Rows()) ExceptionT::SizeMismatch("FluidElementT::FormMass");
#endif

	switch (mass_type)
	{
		case kNoMass:			/* no mass matrix */
		
			break;
		
		/*Only consider consistent mass for now.  Do we need a case for lumped mass? */
		case kConsistentMass:	/* consistent mass	*/
		{
			// integration of the element mass is done
			// in the reference configuration since density
			// is mass/(undeformed volume)
			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();
			
			int nen = NumElementNodes();
			int nun = fLocDisp.NumberOfNodes();
			int ndof = NumDOF();
			
			/* matrix form */
			int a = 0, zero = 0;
			int& b_start = (fLHS.Format() == ElementMatrixT::kSymmetricUpper) ? a : zero;
			
			if (axisymmetric)
			{
				const LocalArrayT& coords = fShapes->Coordinates();
				fShapes->TopIP();	
				while ( fShapes->NextIP() )
				{
					/* compute radius */
					const double* NaX = fShapes->IPShapeX();
					const double* x_r = coords(0); /* r is x-coordinate */
					double r = 0.0;
					for (a = 0; a < nen; a++)
						r += (*NaX++)*(*x_r++);

					/* integration factor */				
					double temp = 2.0*Pi*r*constM*(*Weight++)*(*Det++);
					if (ip_weight) temp *= *ip_weight++;

					const double* Na = fShapes->IPShapeU();		
					for (a = 0; a < nun; a++)
						for (int i = 0; i < ndof; i++)
						{
							int p = a*ndof + i;
							
							/* upper triangle only */
							for (int b = b_start; b < nun; b++) //TEMP - interpolate at the same time?
								for (int j = 0; j < ndof; j++)
									if(i == j) {
										int q = b*ndof + j;
										if (j < fPresIndex)
											fLHS(p,q) += temp*Na[a]*Na[b];
										else fLHS(p,q) = 0.0;
									}
						}
				}
			}			
			else /* not axisymmetric */
			{
				fShapes->TopIP();	
				while ( fShapes->NextIP() )
				{
					/* integration factor */
					double temp = constM*(*Weight++)*(*Det++);
					if (ip_weight) temp *= *ip_weight++;

					const double* Na = fShapes->IPShapeU();		
					for (a = 0; a < nun; a++)
						for (int i = 0; i < ndof; i++)
						{
							int p = a*ndof + i;
							
							/* upper triangle only */
							for (int b = b_start; b < nun; b++)
								for (int j = 0; j < ndof; j++)
									if(i == j) {
										int q = b*ndof + j;
										if (j < fPresIndex)
											fLHS(p,q) += temp*Na[a]*Na[b];
										else fLHS(p,q) = 0.0;
									}
						}
				}
			}
			break;
		}
		default:
			ExceptionT::BadInputValue("FluidElementT::FormMass", "unknown mass matrix code");
	}
}

void FluidElementT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
  WriteCallLocation("LHSDriver"); //DEBUG
	/* inherited */
	ContinuumElementT::LHSDriver(sys_type);
	/* set components and weights */
	double constC = 0.0;
	double constK = 0.0;
	
	int formC = fIntegrator->FormC(constC);
	int formK = fIntegrator->FormK(constK);

	/* loop over elements */
	bool axisymmetric = Axisymmetric();
	Top();
	while (NextElement())
	{
		/* initialize */
		fLHS = 0.0;
		
		/* set shape function derivatives */
		SetGlobalShape();

		/* element mass */
		double density =  fCurrMaterial->Density();
		if (formC) FormMass(kConsistentMass, constC*density, axisymmetric, NULL);

		/* element stiffness */
		if (formK) FormStiffness(constK);
	
		/* add to global equations */
		AssembleLHS();		
	}
}

void FluidElementT::RHSDriver(void)
{
  WriteCallLocation("RHSDriver"); //DEBUG  
	/* inherited */
	ContinuumElementT::RHSDriver();

	/* set components and weights */
	double constCv = 0.0;
	double constKd = 0.0;
	
	/* components dicated by the algorithm */
	int formCv = fIntegrator->FormCv(constCv);
	int formKd = fIntegrator->FormKd(constKd);


	/* body forces */
	int formBody = 0;
	if ((fBodySchedule && fBody.Magnitude() > kSmall)) {	
		formBody = 1;
		if (!formCv) constCv = 1.0; // correct value ??
	}


	bool axisymmetric = Axisymmetric();
	double dt = ElementSupport().TimeStep();
	double by_dt = (fabs(dt) > kSmall) ? 1.0/dt: 0.0; /* for dt -> 0 */
	Top();
	while (NextElement())
	{
		const ElementCardT& element = CurrentElement();

		if (element.Flag() != ElementCardT::kOFF)
		{
			/* initialize */
			fRHS = 0.0;

			/* global shape function values */
			SetGlobalShape();

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

				/* add internal contribution */
				double density = fCurrMaterial->Density();
				FormMa(kConsistentMass, -constCv*density, axisymmetric,
					&fLocVel,NULL, NULL);			  		
			}
				
			/* assemble */
			AssembleRHS();
		}

	}
}

/* calculate the body force contribution */
void FluidElementT::FormMa(MassTypeT mass_type, double constM, bool axisymmetric,
	const LocalArrayT* nodal_values,
	const dArray2DT* ip_values,
	const double* ip_weight)
{
	const char caller[] = "FluidElementT::FormMa";

	/* quick exit */
	if (!nodal_values && !ip_values) return;

#if __option(extended_errorcheck)
	/* dimension checks */
	if (nodal_values && 
		fRHS.Length() != nodal_values->Length()) 
			ExceptionT::SizeMismatch(caller);

	if (ip_values &&
		(ip_values->MajorDim() != fShapes->NumIP() ||
		 ip_values->MinorDim() != NumDOF()))
			ExceptionT::SizeMismatch(caller);
#endif

	switch (mass_type)
	{
		/*Only consider consistent mass for now.  Do we need a case for lumped mass? */
		case kConsistentMass:
		{
			int ndof = NumDOF();
			int  nen = NumElementNodes();
			int  nun = nodal_values->NumberOfNodes();

			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();

			if (axisymmetric)
			{
				const LocalArrayT& coords = fShapes->Coordinates();
				fShapes->TopIP();
				while (fShapes->NextIP())
				{
					/* compute radius */
					const double* NaX = fShapes->IPShapeX();
					const double* x_r = coords(0); /* r is x-coordinate */
					double r = 0.0;
					for (int a = 0; a < nen; a++)
						r += (*NaX++)*(*x_r++);
				
					/* interpolate nodal values to ip */
					if (nodal_values)
						fShapes->InterpolateU(*nodal_values, fDOFvec);
					
					/* ip sources */
					if (ip_values)
						fDOFvec -= (*ip_values)(fShapes->CurrIP());

					/* accumulate in element residual force vector */				
					double*	res      = fRHS.Pointer();
					const double* Na = fShapes->IPShapeU();
				
					/* integration factor */
					double temp = 2.0*Pi*r*constM*(*Weight++)*(*Det++);
					if (ip_weight) temp *= *ip_weight++;

					for (int lnd = 0; lnd < nun; lnd++)
					{
						double temp2 = temp*(*Na++);
						double* pacc = fDOFvec.Pointer();
						for (int dof = 0; dof < ndof; dof++){
							if (dof < fPresIndex)
								*res++ += temp2*(*pacc++);
							else {
								*res++ = 0;
								pacc++;
							}
						}
					}
				}
			}
			else /* not axisymmetric */
			{
				fShapes->TopIP();
				while (fShapes->NextIP())
				{					
					/* interpolate nodal values to ip */
					if (nodal_values)
						fShapes->InterpolateU(*nodal_values, fDOFvec);
					
					/* ip sources */
					if (ip_values)
						fDOFvec -= (*ip_values)(fShapes->CurrIP());

					/* accumulate in element residual force vector */				
					double*	res      = fRHS.Pointer();
					const double* Na = fShapes->IPShapeU();
				
					/* integration factor */
					double temp = constM*(*Weight++)*(*Det++);
					if (ip_weight) temp *= *ip_weight++;
					
					for (int lnd = 0; lnd < nun; lnd++)
					{
						double temp2 = temp*(*Na++);
						double* pacc = fDOFvec.Pointer();
						for (int dof = 0; dof < ndof; dof++){
							if (dof < fPresIndex)
								*res++ += temp2*(*pacc++);
							else {
								*res++ = 0;
								pacc++;
							}
						}
					}
				}
			}
			break;
		}	
		default:
			ExceptionT::BadInputValue("FluidElementT::FormMass", "unknown mass matrix code");
	}
}


/* form the element stiffness matrix */
void FluidElementT::FormStiffness(double constK)
{
	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	/* integration parameters */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	
	/* integrate element stiffness */
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		double scale = constK*(*Det++)*(*Weight++);
	
		/* strain displacement matrix */
/*		B(fShapes->CurrIP(), fB);
		fD.SetToScaled(scale, fCurrMaterial->c_ijkl());							
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
*/
	    fLHS = 0.0;
	}
}

/* calculate the internal force contribution ("-k*d") */
void FluidElementT::FormKd(double constK)
{
	/* integration parameters */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	
	int nsd = NumSD();
	dMatrixT grad;
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		/* set field gradient */

/*		B(fShapes->CurrIP(), fB);

		fB.MultTx(fCurrMaterial->s_ij(), fNEEvec);

		fRHS.AddScaled(-constK*(*Weight++)*(*Det++), fNEEvec);
*/
	    fRHS = 0.0;
	}	
}

/** describe the parameters needed by the interface */
void FluidElementT::DefineParameters(ParameterListT& list) const
{
  WriteCallLocation("DefineParameters"); //DEBUG
  /* inherited */
  ElementBaseT::DefineParameters(list);
}

/** information about subordinate parameter lists */
void FluidElementT::DefineSubs(SubListT& sub_list) const
{
  WriteCallLocation("DefineSubs"); //DEBUG
  /* inherited */
  ContinuumElementT::DefineSubs(sub_list);

  sub_list.AddSub("fluid_element_nodal_output", ParameterListT::ZeroOrOnce);
  sub_list.AddSub("fluid_element_element_output", ParameterListT::ZeroOrOnce);
  sub_list.AddSub("fluid_element_stab_param", ParameterListT::ZeroOrOnce);
  sub_list.AddSub("fluid_element_block", ParameterListT::OnePlus);
}

/** a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FluidElementT::NewSub(const StringT& name) const
{
  WriteCallLocation("NewSub"); //DEBUG
  if (name == "fluid_element_nodal_output")
  {
    ParameterContainerT* node_output = new ParameterContainerT(name);
    /* all false by default */
    for (int i = 0; i < NumNodalOutputCodes; i++)
    {
      ParameterT output(ParameterT::Integer, NodalOutputNames[i]);
      output.SetDefault(1);
      node_output->AddParameter(output, ParameterListT::ZeroOrOnce);
    }
    return node_output;
  }
  else if (name == "fluid_element_element_output")
  {
    ParameterContainerT* element_output = new ParameterContainerT(name);
    /* all false by default */
    for (int i = 0; i < NumElementOutputCodes; i++)
    {
      ParameterT output(ParameterT::Integer, ElementOutputNames[i]);
      output.SetDefault(1);
      element_output->AddParameter(output, ParameterListT::ZeroOrOnce);
    }
    return element_output;
  }
  else if (name == "fluid_element_stab_param")
  {
    ParameterContainerT* stab_param = new ParameterContainerT(name);
    /* all false by default */
    for (int i = 0; i < NumStabParamCodes; i++)
    {
      ParameterT output(ParameterT::Integer, StabParamNames[i]);
      output.SetDefault(1);
      stab_param->AddParameter(output, ParameterListT::ZeroOrOnce);
    }
    return stab_param;
  }  
	else if (name == "fluid_element_block")
	{
		ParameterContainerT* block = new ParameterContainerT(name);
		
		/* list of element block ID's (defined by ElementBaseT) */
		block->AddSub("block_ID_list", ParameterListT::Once);
	
		/* choice of materials lists (inline) */
		block->AddSub("fluid_material", ParameterListT::Once);
	
		/* set this as source of subs */
		block->SetSubSource(this);
		
		return block;
	}
  else /* inherited */
    return ContinuumElementT::NewSub(name);  
}

/** accept parameter list */
void FluidElementT::TakeParameterList(const ParameterListT& list)
{
  WriteCallLocation("TakeParameterList"); //DEBUG
  const char caller[] = "FluidElementT::TakeParameterList";

  /* inherited */
  ContinuumElementT::TakeParameterList(list);

	fPresIndex = NumDOF() - 1;
	/* allocate work space */
	fB.Dimension(dSymMatrixT::NumValues(NumSD()), NumSD()*NumElementNodes());
	fD.Dimension(dSymMatrixT::NumValues(NumSD()));

	/* allocate gradient lists */
	fVel_list.Dimension(NumIP());
	fPres_list.Dimension(NumIP());
	fGradVel_list.Dimension(NumIP());
	fGradPres_list.Dimension(NumIP());
	for (int i = 0; i < NumIP(); i++) {
		fGradVel_list[i].Dimension(NumSD());
		fGradPres_list[i].Dimension(NumSD());
		fVel_list[i].Dimension(NumSD());
	}

  /* nodal output codes */
  fNodalOutputCodes.Dimension(NumNodalOutputCodes);
  fNodalOutputCodes = IOBaseT::kAtNever;
  const ParameterListT* node_output = list.List("fluid_element_nodal_output");
  if (node_output)
  {
    /* set flags */
    for (int i = 0; i < NumNodalOutputCodes; i++)
    {
      /* look for entry */
      const ParameterT* nodal_value = node_output->Parameter(NodalOutputNames[i]);
      if (nodal_value)
      {
        int do_write = *nodal_value;        
        if (do_write == 0)
        {
          fNodalOutputCodes[i] = IOBaseT::kAtInc;
          WriteCallLocation("Picked: fNodalOutputCodes: coordinates"); //DEBUG
        }
        else if (do_write == 1)
        {
          fNodalOutputCodes[i] = IOBaseT::kAtInc;
          WriteCallLocation("Picked: fNodalOutputCodes: velocities"); //DEBUG
        }
        else if (do_write == 2)
        {
          /* check order of the time integrator */
          if (fIntegrator->Order() < 1)
            ExceptionT::GeneralFail(caller, "expecting time integrator order of 1 or higher");          
          fNodalOutputCodes[i] = IOBaseT::kAtInc;
          WriteCallLocation("Picked: fNodalOutputCodes: accelerations"); //DEBUG
        }
        else if (do_write == 3)
        {
          fNodalOutputCodes[i] = IOBaseT::kAtInc;
          WriteCallLocation("Picked: fNodalOutputCodes: pressures"); //DEBUG
        }
      }
    }
  }
  
  /* element output codes */
  fElementOutputCodes.Dimension(NumElementOutputCodes);
  fElementOutputCodes = IOBaseT::kAtNever;
  const ParameterListT* element_output = list.List("fluid_element_element_output");
  if (element_output)
  {
    /* set flags */
    for (int i = 0; i < NumElementOutputCodes; i++)
    {
      /* look for entry */
      const ParameterT* element_value = element_output->Parameter(ElementOutputNames[i]);
      if (element_value)
      {
        int do_write = *element_value;
        if (do_write == 0)
        {
          fElementOutputCodes[i] = IOBaseT::kAtInc;
          WriteCallLocation("Picked: fElementOutputCodes: NONE"); //DEBUG
        }
      }
    }
  }

  /* stabilization parameter codes */
  const ParameterListT* stab_param = list.List("fluid_element_stab_param");
  if (stab_param)
  {
    /* set flags */
    for (int i = 0; i < NumStabParamCodes; i++)
    {
      /* look for entry */
      const ParameterT* stab_param_value = stab_param->Parameter(StabParamNames[i]);
      if (stab_param_value)
      {
        int param_value = *stab_param_value;
        if (param_value == 0)
        {
          /* implement \tau_m = \tau_c = \tau_PSPG = \tau_SUPG \tau
           *
           * \tau = [ (2/(del t))^2 + (2/h * ||v||) + (4\mu / h^2) ]^(-1/2)
           * ** check reference
           */
           WriteCallLocation("Picked: \tau_m = \tau_c"); //DEBUG
        }
      }
    }
  }  
    
}

/* extract the list of material parameters */
void FluidElementT::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{
	const char caller[] = "FluidElementT::CollectMaterialInfo";
	
	/* initialize */
	mat_params.Clear();

	/* set materials list name */
	mat_params.SetName("fluid_material");
	
	/* collected material parameters */
	int num_blocks = all_params.NumLists("fluid_element_block");
	for (int i = 0; i < num_blocks; i++) {

		/* block information */	
		const ParameterListT& block = all_params.GetList("fluid_element_block", i);
		
		/* collect material parameters */
		const ParameterListT& mat_list = block.GetList(mat_params.Name());
		const ArrayT<ParameterListT>& mat = mat_list.Lists();
		mat_params.AddList(mat[0]);
	}
}

/* construct output labels array */
void FluidElementT::SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
  iArrayT& counts) const
{
  WriteCallLocation("SetNodalOutputCodes"); //DEBUG
  /* initialize */
  counts.Dimension(flags.Length());
  counts = 0;

  /* set output flags */
  if (flags[iNodalCrd] == mode) counts[iNodalCrd] = NumSD();
  if (flags[iNodalVel] == mode) counts[iNodalVel] = NumSD();
  if (flags[iNodalAcc] == mode) counts[iNodalAcc] = NumSD();
  if (flags[iNodalPrs] == mode) counts[iNodalPrs] = FluidElementT::kPressureNDOF;
  if (flags[iNodalPrs] == mode) counts[iNodalPrs] = 1;
}

void FluidElementT::SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
  iArrayT& counts) const
{
  WriteCallLocation("SetElementOutputCodes"); //DEBUG
  /* initialize */
  counts.Dimension(flags.Length());
  counts = 0;

  /* set output flags */
  if (fElementOutputCodes[iNONE] == mode) counts[iNONE] = 0;
}

void FluidElementT::GenerateOutputLabels(const iArrayT& n_codes,
    ArrayT<StringT>& n_labels, const iArrayT& e_codes, ArrayT<StringT>& e_labels) const
{
  WriteCallLocation("GenerateOutputLabels"); //DEBUG
  const char caller[] = "FluidElementT::GenerateOutputLabels";

  /* allocate node labels */
  n_labels.Dimension(n_codes.Sum());

  int count = 0;
  if (n_codes[iNodalCrd])
  {
    const char* xlabels[] = {"x1", "x2", "x3"};
    for (int i = 0; i < NumSD(); i++)
      n_labels[count++] = xlabels[i];
  }

  if (n_codes[iNodalVel])
  {
    const char* xlabels[] = {"v1", "v2", "v3"};
    for (int i = 0; i < NumSD(); i++)
      n_labels[count++] = xlabels[i];
  }

  if (n_codes[iNodalAcc])
  {
    const char* xlabels[] = {"a1", "a2", "a3"};
    for (int i = 0; i < NumSD(); i++)
      n_labels[count++] = xlabels[i];
  }

  if (n_codes[iNodalPrs])
  {
    const char* xlabels[] = {"p"};
    for (int i = 0; i < FluidElementT::kPressureNDOF; i++)
      n_labels[count++] = xlabels[i];
  }

  if (e_codes.Sum() != 0)
    ExceptionT::GeneralFail("FluidElementT::GenerateOutputLabels",
      "not expecting any element output codes");  
}

/***********************************************************************
 * Private
 ***********************************************************************/


/** FOR DEBUGGING PURPOSES ONLY */
void FluidElementT::WriteCallLocation( char* loc ) const {
cout << "Inside of FluidElementT::" << loc << endl;
}