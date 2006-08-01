/* $Header: /home/regueiro/tahoe_cloudforge_repo_snapshots/development/src/elements/fluid_element/FluidElementT.cpp,v 1.16 2006-08-01 23:34:26 a-kopacz Exp $ */
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

const int FluidElementT::NumStabParamCodes = 2;
static const char* StabParamNames[] = {
  "tau_m_is_tau_c", /* T.E. Tezduyar, Stabilized Finite Element Formulations for Incompressible Flow Computations, Adv. Appl. Mech. 28(1991) 1-44. */
  "NONE"}; /* Galerkin formulation; ie \tau_m and \tau_c = 0

/* parameters */
const int FluidElementT::kPressureNDOF = 1;

/* constructor */
FluidElementT::FluidElementT(const ElementSupportT& support):
  ContinuumElementT(support),
  /* dofs: vels and press */
  fLocDisp(LocalArrayT::kDisp),
  fLocLastDisp(LocalArrayT::kLastDisp),
  fLocVel(LocalArrayT::kVel),

  /* velocity */
  fLocCurVel(LocalArrayT::kUnspecified),
  fLocOldVel(LocalArrayT::kUnspecified),
  /* accelaration */
  fLocCurAcc(LocalArrayT::kUnspecified),
  /* pressure */
  fLocCurPrs(LocalArrayT::kUnspecified),
  fFluidMatSupport(NULL)
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
      flags[iNodalCrd] = NumSD();   WriteCallLocation("iNodalCrd");
    break;
    case iNodalVel:
      flags[iNodalVel] = NumSD();   WriteCallLocation("iNodalVel");
    break;
    case iNodalAcc:
      flags[iNodalAcc] = NumSD();   WriteCallLocation("iNodalAcc");
    break;
    case iNodalPrs:
      flags[iNodalPrs] = FluidElementT::kPressureNDOF;   WriteCallLocation("iNodalPrs");
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

/** initialize local arrays */
void FluidElementT::SetLocalArrays(void)
{
  WriteCallLocation("SetLocalArrays"); //DEBUG

  /* inherited */
  ContinuumElementT::SetLocalArrays();

  /* allocate */
  int nen = NumElementNodes();
  int ndof  = NumDOF();

  fLocDisp.Dimension(nen, ndof);
  fLocLastDisp.Dimension(nen, ndof);
  fLocVel.Dimension(nen, ndof);

  /* set source */
  Field().RegisterLocal(fLocDisp);
  Field().RegisterLocal(fLocLastDisp);

  if (fIntegrator->Order() > 0)
    Field().RegisterLocal(fLocVel);

  /* map */
  const double* p_cur_vel = fLocDisp(0);
  fLocCurVel.Alias(nen, ndof - 1, p_cur_vel);

  const double* p_cur_pres = fLocDisp(ndof-1);
  fLocCurPrs.Alias(nen, 1, p_cur_pres);

  const double* p_old_vel = fLocLastDisp(0);
  fLocOldVel.Alias(nen, ndof - 1, p_old_vel);

  const double* p_cur_acc = fLocVel(0);
  fLocCurAcc.Alias(nen, ndof - 1, p_cur_acc);
}

/** allocate and initialize shape function objects */
void FluidElementT::SetShape(void)
{
  WriteCallLocation("SetShape"); //DEBUG
  fShapes = new ShapeFunctionT(GeometryCode(), NumIP(), fLocInitCoords);
  if (!fShapes ) throw ExceptionT::kOutOfMemory;
  fShapes->Initialize();
}

/** construct a new material support and return a pointer */
MaterialSupportT* FluidElementT::NewMaterialSupport(MaterialSupportT* p) const
{
  WriteCallLocation("NewMaterialSupport"); //DEBUG
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
  WriteCallLocation("NewMaterialList"); //DEBUG
	/* no match */
	if (name != "fluid_material")
		return NULL;

	if (size > 0)
	{
    /* material support */
		if (!fFluidMatSupport)
    {
			fFluidMatSupport = TB_DYNAMIC_CAST(FluidMatSupportT*, NewMaterialSupport());
			if (!fFluidMatSupport)
        ExceptionT::GeneralFail("FluidElementT::NewMaterialList");
		}
    /* allocate */
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

#pragma unused(e_values)
  if (e_out > 0)
    ExceptionT::GeneralFail("FluidElementT::ComputeOutput", "element output not supported");

  /* dimensions */
  int nen = NumElementNodes();

  /* reset averaging workspace */
  ElementSupport().ResetAverage(n_out);

  /* work arrays */
  dArray2DT nodal_space(nen, n_out);
  dArray2DT nodal_all(nen, n_out);
  dArray2DT coords, vel, acc, prs;

  /* set shallow copies */
  double* pall = nodal_space.Pointer();
  coords.Set(nen, n_codes[iNodalCrd], pall);
  pall += coords.Length();
  vel.Set(nen, n_codes[iNodalVel], pall);
  pall += vel.Length();
  acc.Set(nen, n_codes[iNodalAcc], pall);
  pall += acc.Length();
  prs.Set(nen, n_codes[iNodalPrs], pall);
  pall += prs.Length();

  Top();
  while (NextElement())
  {
    /* initialize */
    nodal_space = 0.0;
    /* global shape function values */
    SetGlobalShape();
    SetLocalU(fLocInitCoords);

    /* coordinates and displacements all at once */
    if (n_codes[iNodalCrd])  fLocInitCoords.ReturnTranspose(coords);
    if (n_codes[iNodalVel])  fLocCurVel.ReturnTranspose(vel);
    if (n_codes[iNodalAcc])  fLocCurAcc.ReturnTranspose(acc);
    if (n_codes[iNodalPrs])  fLocCurPrs.ReturnTranspose(prs);

    /* copy in the cols (in sequence of output) */
    int colcount = 0;
    nodal_all.BlockColumnCopyAt(coords, colcount); colcount += coords.MinorDim();
    nodal_all.BlockColumnCopyAt(vel  , colcount); colcount += vel.MinorDim();
    nodal_all.BlockColumnCopyAt(acc, colcount); colcount += acc.MinorDim();
    nodal_all.BlockColumnCopyAt(prs, colcount);

    /* accumulate - extrapolation done from ip's to corners => X nodes */
    ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_all);
  }

  /* get nodally averaged values */
  ElementSupport().OutputUsedAverage(n_values);
}

/* current element operations */
bool FluidElementT::NextElement(void)
{
  WriteCallLocation("NextElement"); //DEBUG
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

/* form shape functions and derivatives */
void FluidElementT::SetGlobalShape(void)
{
  WriteCallLocation("SetGlobalShape"); //DEBUG
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
    dArrayT& ovel = fOldVel_list[i];
		dMatrixT& gradV = fGradVel_list[i];
		dMatrixT gradP;
		gradP.Set(1,NumSD(), fGradPres_list[i].Pointer());
    dArrayT pres;
    double* p = fPres_list.Pointer();
    pres.Set(1,p+i);
		fShapes->GradU(fLocCurVel, gradV, i);
		fShapes->GradU(fLocCurPrs, gradP, i);
		fShapes->InterpolateU(fLocCurPrs, pres, i);
		fShapes->InterpolateU(fLocCurVel, vel, i);
		fShapes->InterpolateU(fLocOldVel, ovel, i);
	}
}

/* set the \e B matrix at the specified integration point */
void FluidElementT::B(int ip, dMatrixT& B_matrix) const
{
    const dArray2DT& DNa = fShapes->Derivatives_U(ip);
    int nnd = DNa.MinorDim();
    double* pB = B_matrix.Pointer();
    /* 2D */
    if (DNa.MajorDim() == 2)
    {
        const double* pNax = DNa(0);
        const double* pNay = DNa(1);
        for (int i = 0; i < nnd; i++)
        {
            *pB++ = *pNax++;
            *pB++ = *pNay++;
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
            *pB++ = *pNax++;
            *pB++ = *pNay++;
            *pB++ = *pNaz++;
        }
    }
}

void FluidElementT::Set_B(const dArray2DT& DNa, dMatrixT& B) const
{
  WriteCallLocation("Set_B"); //DEBUG
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

/* set initial velocities */
void FluidElementT::InitialCondition(void)
{
  WriteCallLocation("InitialCondition"); //DEBUG
	/* inherited */
	ContinuumElementT::InitialCondition();
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
      tau_m = 0.0;
      tau_c = 0.0;

			/* global shape function values */
			SetGlobalShape();

      if(StabParamNames[fStabParam]=="tau_m_is_tau_c")
      {
        double viscosity = fCurrMaterial->Shear_Modulus();
        double h_nsum=0.0;
        double h=0.0;
        double OldVelMag=0.0;


        /* degrees of freedom */
			  int ndof = NumDOF();
        /* dimensions */
        int  nsd = NumSD();
			  int  nen = NumElementNodes();
			  int  nun = fLocDisp.NumberOfNodes();

			  const double* Det    = fShapes->IPDets();
			  const double* Weight = fShapes->IPWeights();
        
        /* loop over integration points */
 				fShapes->TopIP();
				while (fShapes->NextIP())
				{
          const double* Na          = fShapes->IPShapeU();                   /* [nun] */
	        const dArray2DT& GradNa   = fShapes->Derivatives_U();              /* [nsd x nun] */
          const dArrayT& OldVel     = fOldVel_list[fShapes->CurrIP()];       /* [nsd] */

          /* integration factor */
					double temp1 = (*Weight++)*(*Det++);

          /* calculate stabilization terms */
          for (int lnd = 0; lnd < nun; lnd++)
          {
            if ( nsd == 2 )
            {
              OldVelMag += temp1*Na[lnd]*sqrt( pow(OldVel[0],2)+pow(OldVel[1],2) );
              h_nsum += temp1*sqrt( pow(OldVel[0]*GradNa(0,lnd)+OldVel[0]*GradNa(1,lnd)+OldVel[0]*GradNa(2,lnd),2)
                                   +pow(OldVel[1]*GradNa(0,lnd)+OldVel[1]*GradNa(1,lnd)+OldVel[1]*GradNa(2,lnd),2) );
            }
            else /* 3D */
            {
              OldVelMag += temp1*Na[lnd]*sqrt( pow(OldVel[0],2)+pow(OldVel[1],2)+pow(OldVel[2],2) );
              h_nsum += temp1*sqrt( pow(OldVel[0]*GradNa(0,lnd)+OldVel[0]*GradNa(1,lnd)+OldVel[0]*GradNa(2,lnd),2)
                                   +pow(OldVel[1]*GradNa(0,lnd)+OldVel[1]*GradNa(1,lnd)+OldVel[1]*GradNa(2,lnd),2)
                                   +pow(OldVel[2]*GradNa(0,lnd)+OldVel[2]*GradNa(1,lnd)+OldVel[2]*GradNa(2,lnd),2) );                   
            }
          }
        }
   
        /* check for zero values */
        if ( h_nsum == 0.0 ) h_nsum = 1e-12;
        if ( OldVelMag == 0.0 ) OldVelMag = 1e-12;

        /* T.E. Tezduyar, Y.Osawa / Comput. Methods Appl. Mech. Engrg. 190 (2000) 411-430 {eq. 54} */
        h=2*OldVelMag*(1/h_nsum);
        /* T.E. Tezduyar, Y.Osawa / Comput. Methods Appl. Mech. Engrg. 190 (2000) 411-430 {eq. 58} */
        tau_m = 1/sqrt( (2*by_dt)*(2*by_dt) + (2*OldVelMag/h)*(2*OldVelMag/h) + (4*viscosity/(h*h))*(4*viscosity/(h*h)) );
        /* T.E. Tezduyar, Y.Osawa / Comput. Methods Appl. Mech. Engrg. 190 (2000) 411-430 {eq. 59} */
        tau_c=tau_m;
      }

			if (formKd)
      {
        SetLocalU(fLocDisp);
        FormKd(-constKd);
      }

			if (formBody) AddBodyForce(fLocVel);

			/* add internal contribution */
			double density = fCurrMaterial->Density();
			FormMa(kConsistentMass, -constCv*density, axisymmetric, &fLocVel, NULL, NULL);

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
  WriteCallLocation("FormMa"); //DEBUG
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
      /* degrees of freedom */
			int ndof = NumDOF();
      /* dimensions */
      int  nsd = NumSD();
			int  nen = NumElementNodes();
			int  nun = nodal_values->NumberOfNodes();

			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();

			if (axisymmetric)
			{
        ExceptionT::BadInputValue("FluidElementT::FormMa", "axisymmetric not applicable");
			}
			else /* not axisymmetric */
			{
				fShapes->TopIP();
				while (fShapes->NextIP())
				{
					/* interpolate nodal values to ip */
					if (nodal_values) fShapes->InterpolateU(*nodal_values, fDOFvec);

					/* ip sources */
					if (ip_values)  fDOFvec -= (*ip_values)(fShapes->CurrIP());

					/* accumulate in element force vector */
					double*	pfRHS             = fRHS.Pointer();
          const double* Na          = fShapes->IPShapeU();                   /* [nun] */
	        const dArray2DT& GradNa   = fShapes->Derivatives_U();              /* [nsd x nun] */
          const dArrayT& OldVel     = fOldVel_list[fShapes->CurrIP()];       /* [nsd] */
          double temp0;

					/* integration factor */
					double temp1 = constM*(*Weight++)*(*Det++);
					if (ip_weight) temp1 *= *ip_weight++;
    
					for (int lnd = 0; lnd < nun; lnd++)
					{
            /* temp0 = v_{k,old}*N_{A,k} */
            temp0 = 0.0;
            if ( nsd == 2 )
              temp0 += OldVel[0]*GradNa(0,lnd)+OldVel[1]*GradNa(1,lnd);
            else /* 3D */
              temp0 += OldVel[0]*GradNa(0,lnd)+OldVel[1]*GradNa(1,lnd)+OldVel[2]*GradNa(2,lnd);

            /* temp3 = N_{A,i}*\dot{v_{i}} */
            double temp3 = 0.0;
            double* pacc = fDOFvec.Pointer();
            for (int dof = 0; dof < nsd; dof++)
            {
              /* term : N_{A}*\rho*\dot{v_{i}} */
              *pfRHS += temp1*(*Na)*(*pacc);

              /* term : \tau^{m}*v_{k,old}*N_{A,k}*\rho*\dot{v_{i}} */
              *pfRHS += temp1*tau_m*temp0*(*pacc);
              *pfRHS++;
                
              /* temp3 = N_{A,i}*\dot{v_{i}} */
              temp3 += GradNa(dof,lnd)*(*pacc);
              *pacc++;  
            }
            /* term : \tau^{c}*N_{A,i}*\rho*\dot{v_{i}} */
            *pfRHS += temp1*tau_c*temp3;
            *pfRHS++;
            *pacc++;           
            *Na++;
					}
				}
			}
			break;
		}
		default:
			ExceptionT::BadInputValue("FluidElementT::FormMass", "unknown mass matrix code");
	}
}

/* calculate the internal force contribution ("-k*d") */
void FluidElementT::FormKd(double constK)
{
  WriteCallLocation("FormKd"); //DEBUG

  /* degrees of freedom */
  int ndof = NumDOF();
  /* dimensions */
  int  nsd = NumSD();
  int  nen = NumElementNodes();
  int  nun = fLocDisp.NumberOfNodes();
      
	/* integration parameters */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
  double Density = fCurrMaterial->Density();

	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
    double*	pfRHS             = fRHS.Pointer();							/* [nun] */
    const double* Na          = fShapes->IPShapeU();                   /* [nun] */
	  const dArray2DT& GradNa   = fShapes->Derivatives_U();              /* [nsd x nun] */

    const dArrayT& OldVel     = fOldVel_list[fShapes->CurrIP()];       /* [nsd] */

    const dArrayT& Vel        = fVel_list[fShapes->CurrIP()];          /* [nsd] */
    const dMatrixT& GradVel   = fGradVel_list[fShapes->CurrIP()];      /* [nsd x nsd] */

    const double& Pres        = fPres_list[fShapes->CurrIP()];         /* [1] */
    const dArrayT& GradPres   = fGradPres_list[fShapes->CurrIP()];     /* [nsd] */
    double temp0, temp5;
    dArrayT temp7(nsd), temp4(nsd);
    const dSymMatrixT& s_ij = fCurrMaterial->s_ij();					/*[nsd x nsd] or [numstress]
    
	  /* integration factor */
	  double temp1 = constK*(*Weight++)*(*Det++);

    /* temp4[i] = v_{j}*v{i,j} */
    temp4 = 0.0;
    /* temp5 = v_{j,j} */
    temp5 = 0.0;
    if ( nsd ==2 )
    {
      temp4[0] = OldVel[0]*GradVel(0,0)+OldVel[1]*GradVel(0,1);
      temp4[1] = OldVel[0]*GradVel(1,0)+OldVel[1]*GradVel(1,1);
      temp5    = GradVel(0,0)+GradVel(1,1);
    }
    else /* 3D */
    {
      temp4[0] = OldVel[0]*GradVel(0,0)+OldVel[1]*GradVel(0,1)+OldVel[2]*GradVel(0,2);
      temp4[1] = OldVel[0]*GradVel(1,0)+OldVel[1]*GradVel(1,1)+OldVel[2]*GradVel(1,2);
      temp4[2] = OldVel[0]*GradVel(2,0)+OldVel[1]*GradVel(2,1)+OldVel[2]*GradVel(2,2);
      temp5    = GradVel(0,0)+GradVel(1,1)+GradVel(2,2);
    }

    for (int lnd = 0; lnd < nen; lnd++)
    {
      /* temp0 = v_{k,old}*N_{A,k} */
      temp0 = 0.0;
      /* temp7[i] = \sigma_{ij}*N_{A,j}*/
      temp7 = 0.0;
      if ( nsd == 2 )
      {
        temp0 += OldVel[0]*GradNa(0,lnd)+OldVel[1]*GradNa(1,lnd);
        temp7[0] += s_ij[0]*GradNa(0,lnd) + s_ij[2]*GradNa(1,lnd);
        temp7[1] += s_ij[1]*GradNa(1,lnd) + s_ij[2]*GradNa(0,lnd);
      }
      else /* 3D */
      {
        temp0 += OldVel[0]*GradNa(0,lnd)+OldVel[1]*GradNa(1,lnd)+OldVel[2]*GradNa(2,lnd);
        temp7[0] += s_ij[0]*GradNa(0,lnd) + s_ij[4]*GradNa(2,lnd)+ s_ij[5]*GradNa(1,lnd);
        temp7[1] += s_ij[1]*GradNa(1,lnd) + s_ij[3]*GradNa(2,lnd)+ s_ij[5]*GradNa(0,lnd);
        temp7[2] += s_ij[2]*GradNa(2,lnd) + s_ij[3]*GradNa(1,lnd)+ s_ij[4]*GradNa(0,lnd);
      }

      /* temp3 = N_{A,i}*[ v_{j}*v_{i,j} ] */
      double temp3 = 0.0;
      /* temp6 = N_{A,i}*p_{,i} */
      double temp6 = 0.0;
      for (int dof = 0; dof < nsd; dof++)
      {
        /* term : N_{A}*\rho*v_{j}*v_{i,j} */
        *pfRHS += temp1*(*Na)*Density*temp4[dof];
          
        /* term : \tau^{m}*v_{k,old}*N_{A,k}*\rho*v_{j}*v_{i,j} */
        *pfRHS += temp1*tau_m*temp0*Density*temp4[dof];

        /* term : N_{A,j}*\sigma_{ij} */
        *pfRHS += temp1*temp7[dof];
          
        /* term : \tau^{m}*v_{k,old}*N_{A,k}*p_(,i) */
        *pfRHS += temp1*tau_m*temp0*GradPres[dof];
          
        /* term : \tau^{c}*N_{A,i}*v_{j,j} */
        *pfRHS += temp1*tau_c*GradNa(dof,lnd)*temp5;
          
        *pfRHS++;

        /* temp3 = N_{A,i}*[ v_{j}*v_{i,j} ] */
        temp3 += GradNa(dof,lnd)*temp4[dof];

        /* temp6 = N_{A,i}*p_{,i} */
        temp6 += GradNa(dof,lnd)*GradPres[dof];
      }
      /* term : \tau^{c}*N_{A,i}*\rho*v_{j}*v_{i,j} */
      *pfRHS += temp1*tau_c*Density*temp3;

      /* term : \tau^{c}*N_{A,i}*p_{,i} */
      *pfRHS += temp1*tau_c*temp6;
      
      /* term : N_{A}*v_{j,j} */
      *pfRHS += temp1*(*Na)*temp5;

      *pfRHS++;
      *Na++;
    }
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
	double dt = ElementSupport().TimeStep();
	double by_dt = (fabs(dt) > kSmall) ? 1.0/dt: 0.0; /* for dt -> 0 */

	/* loop over elements */
	bool axisymmetric = Axisymmetric();
	Top();
	while (NextElement())
	{
		/* initialize */
		fLHS = 0.0;
    tau_m = 0.0;
    tau_c = 0.0;

    /* global shape function values */
		SetGlobalShape();

    if(StabParamNames[fStabParam]=="tau_m_is_tau_c")
    {
      double viscosity = fCurrMaterial->Shear_Modulus();
      double h_nsum=0.0;
      double h=0.0;
      double OldVelMag=0.0;

      /* degrees of freedom */
		  int ndof = NumDOF();
      /* dimensions */
      int  nsd = NumSD();
		  int  nen = NumElementNodes();
		  int  nun = fLocDisp.NumberOfNodes();

		  const double* Det    = fShapes->IPDets();
		  const double* Weight = fShapes->IPWeights();

      /* loop over integration points */
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
        const double* Na          = fShapes->IPShapeU();                   /* [nun] */
        const dArray2DT& GradNa   = fShapes->Derivatives_U();              /* [nsd x nun] */
        const dArrayT& OldVel     = fOldVel_list[fShapes->CurrIP()];       /* [nsd] */

        /* integration factor */
				double temp1 = (*Weight++)*(*Det++);

        /* calculate stabilization terms */
        for (int lnd = 0; lnd < nun; lnd++)
        {
          if ( nsd == 2 )
          {
            OldVelMag += temp1*Na[lnd]*sqrt( pow(OldVel[0],2)+pow(OldVel[1],2) );
            h_nsum += temp1*sqrt( pow(OldVel[0]*GradNa(0,lnd)+OldVel[0]*GradNa(1,lnd)+OldVel[0]*GradNa(2,lnd),2)
                                 +pow(OldVel[1]*GradNa(0,lnd)+OldVel[1]*GradNa(1,lnd)+OldVel[1]*GradNa(2,lnd),2) );
          }
          else /* 3D */
          {
            OldVelMag += temp1*Na[lnd]*sqrt( pow(OldVel[0],2)+pow(OldVel[1],2)+pow(OldVel[2],2) );
            h_nsum += temp1*sqrt( pow(OldVel[0]*GradNa(0,lnd)+OldVel[0]*GradNa(1,lnd)+OldVel[0]*GradNa(2,lnd),2)
                                 +pow(OldVel[1]*GradNa(0,lnd)+OldVel[1]*GradNa(1,lnd)+OldVel[1]*GradNa(2,lnd),2)
                                 +pow(OldVel[2]*GradNa(0,lnd)+OldVel[2]*GradNa(1,lnd)+OldVel[2]*GradNa(2,lnd),2) );
          }
        }
      }

      /* check for zero values */
      if ( h_nsum == 0.0 ) h_nsum = 1e-12;
      if ( OldVelMag == 0.0 ) OldVelMag = 1e-12;

      /* T.E. Tezduyar, Y.Osawa / Comput. Methods Appl. Mech. Engrg. 190 (2000) 411-430 {eq. 54} */
      h=2*OldVelMag*(1/h_nsum);
      /* T.E. Tezduyar, Y.Osawa / Comput. Methods Appl. Mech. Engrg. 190 (2000) 411-430 {eq. 58} */
      tau_m = 1/sqrt( (2*by_dt)*(2*by_dt) + (2*OldVelMag/h)*(2*OldVelMag/h) + (4*viscosity/(h*h))*(4*viscosity/(h*h)) );
      /* T.E. Tezduyar, Y.Osawa / Comput. Methods Appl. Mech. Engrg. 190 (2000) 411-430 {eq. 59} */
      tau_c=tau_m;
    }

		/* element mass */
		double density =  fCurrMaterial->Density();
		if (formC) FormMass(kConsistentMass, constC*density, axisymmetric, NULL);

		/* element stiffness */
		if (formK) FormStiffness(constK);

		/* add to global equations */
		AssembleLHS();
	}
}

/* form the element mass matrix */
void FluidElementT::FormMass(MassTypeT mass_type, double constM, bool axisymmetric, const double* ip_weight)
{
  WriteCallLocation("FormMass"); //DEBUG
	const char caller[] = "FluidElementT::FormMass";
#if __option(extended_errorcheck)
	if (fLocDisp.Length() != fLHS.Rows()) ExceptionT::SizeMismatch(caller);
#endif

	switch (mass_type)
	{
		case kNoMass:			/* no mass matrix */
		break;

		/*Only consider consistent mass for now.  Do we need a case for lumped mass? */
		case kConsistentMass:	/* consistent mass	*/
		{ 
      /* degrees of freedom */
			int ndof = NumDOF();
      /* dimensions */
		  int  nsd = NumSD();
			int  nun = fLocDisp.NumberOfNodes();

			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();

			if (axisymmetric)
			{
        ExceptionT::BadInputValue("FluidElementT::FormMass", "axisymmetric not applicable");
			}
			else /* not axisymmetric */
			{
        const LocalArrayT& coords = fShapes->Coordinates();
				fShapes->TopIP();
				while ( fShapes->NextIP() )
				{
					/* accumulate in element force vector */
					const double* Na          = fShapes->IPShapeU();                   /* [nun] */
					const dArray2DT& GradNa   = fShapes->Derivatives_U();              /* [nsd x nun] */
          const dArrayT& OldVel     = fOldVel_list[fShapes->CurrIP()];       /* [nsd] */
          double temp0;
          dArrayT temp3(nsd);

					/* integration factor */
					double temp1 = constM*(*Weight++)*(*Det++);
					if (ip_weight) temp1 *= *ip_weight++;

          int a,i,b,j,p,q;
					for (a = 0; a < nun; a++)
					{
						/* temp0 = v_{k,old}*N_{A,k} */
						temp0 = 0.0;
						if ( nsd == 2 )
							temp0 += OldVel[0]*GradNa(0,a)+OldVel[1]*GradNa(1,a);
						else /* 3D */
							temp0 += OldVel[0]*GradNa(0,a)+OldVel[1]*GradNa(1,a)+OldVel[2]*GradNa(2,a);

            for (i = 0; i < nsd; i++)
            {
              p = a*ndof + i;
              for (b = 0; b < nun; b++)
                for (j = 0; j < nsd; j++)
                {
                  q = b*ndof + j;
                  if(i == j)
                  {
                    /* term : \rho*N_{A}*N_{B} \delta_{ij} */
                    fLHS(p,q) += temp1*Na[a]*Na[b];

                    /* term : \tau^{m}*v_{k,old}*N_{A,k}*\rho*N_{B} \delta_{ij} */
                    fLHS(p,q) += temp1*tau_m*temp0*Na[b];
                  }
                }
            }
            /* i4th term */
            p = a*ndof + 3;
            for (b = 0; b < nun; b++)
              for (j = 0; j < nsd; j++)
              {
                q = b*ndof + j;

                /* term : \tau^{c}*N_{A,j}*\rho*N_{B} */
                fLHS(p,q) += temp1*tau_c*Na[b]*GradNa[j,a];
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
  WriteCallLocation("FormStiffness"); //DEBUG
	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

  /* degrees of freedom */
  int ndof = NumDOF();
  /* dimensions */
  int  nsd = NumSD();
  int  nen = NumElementNodes();
  int  nun = fLocDisp.NumberOfNodes();

	/* integration parameters */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
  double Density = fCurrMaterial->Density();

	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
    double*	pfRHS             = fRHS.Pointer();							/* [nun] */
    const double* Na          = fShapes->IPShapeU();                   /* [nun] */
	  const dArray2DT& GradNa   = fShapes->Derivatives_U();              /* [nsd x nun] */

    const dArrayT& OldVel     = fOldVel_list[fShapes->CurrIP()];       /* [nsd] */

    const dArrayT& Vel        = fVel_list[fShapes->CurrIP()];          /* [nsd] */
    const dMatrixT& GradVel   = fGradVel_list[fShapes->CurrIP()];      /* [nsd x nsd] */

    const double& Pres        = fPres_list[fShapes->CurrIP()];         /* [1] */
    const dArrayT& GradPres   = fGradPres_list[fShapes->CurrIP()];     /* [nsd] */
    double temp0, temp2, temp3, temp4;        
    dMatrixT temp5(nsd);
    const dMatrixT& c = fCurrMaterial->c_ijkl();
    
    /* integration factor */
	  double temp1 = constK*(*Weight++)*(*Det++);
    
    int a,i,b,j,p,q;
		for (a = 0; a < nun; a++)
		{
  		/* temp0 = v_{k,old}*N_{A,k} */
			temp0 = 0.0;
      /* temp2 = v_{l}*N_{B,l} */
			temp2 = 0.0;
      /* temp4 = N_{A,i}*N_{B,i} */
      temp4 = 0.0;
      /* N_{A,k}*c_{ikjl}*N_{B,l} == > B^{T}*D*B */
      temp5 = 0.0;
			if ( nsd == 2 )
      {
					temp0 += OldVel[0]*GradNa(0,a)+OldVel[1]*GradNa(1,a);
					temp2 += Vel[0]*GradNa(0,a)+Vel[1]*GradNa(1,a);
					temp4 += GradNa(0,a)*GradNa(0,a)+GradNa(1,a)*GradNa(1,a);

          /* note: GradNa == GradNb */
          temp5(0,0) += ( c(0,0)*GradNa(0,a)+c(2,0)*GradNa(1,a) )*GradNa(0,a)+
                        ( c(0,2)*GradNa(0,a)+c(2,2)*GradNa(1,a) )*GradNa(1,a);
                        
          temp5(0,1) += ( c(0,1)*GradNa(0,a)+c(2,1)*GradNa(1,a) )*GradNa(1,a)+
                        ( c(0,2)*GradNa(0,a)+c(2,2)*GradNa(1,a) )*GradNa(0,a);
          /**/              
          temp5(1,0) += ( c(1,0)*GradNa(1,a)+c(2,0)*GradNa(0,a) )*GradNa(0,a)+
                        ( c(1,2)*GradNa(1,a)+c(2,2)*GradNa(0,a) )*GradNa(1,a);
                        
          temp5(1,1) += ( c(1,1)*GradNa(1,a)+c(2,1)*GradNa(0,a) )*GradNa(1,a)+
                        ( c(1,2)*GradNa(1,a)+c(2,2)*GradNa(0,a) )*GradNa(0,a);           
      }
			else /* 3D */
      {
					temp0 += OldVel[0]*GradNa(0,a)+OldVel[1]*GradNa(1,a)+OldVel[2]*GradNa(2,a);
					temp2 += Vel[0]*GradNa(0,a)+Vel[1]*GradNa(1,a)+Vel[2]*GradNa(2,a);
					temp4 += GradNa(0,a)*GradNa(0,a)+GradNa(1,a)*GradNa(1,a)+GradNa(2,a)*GradNa(2,a);

          /* note: GradNa == GradNb */
          temp5(0,0) += ( c(0,0)*GradNa(0,a)+c(4,0)*GradNa(2,a)+c(5,0)*GradNa(1,a) )*GradNa(0,a)+
                        ( c(0,4)*GradNa(0,a)+c(4,4)*GradNa(2,a)+c(5,4)*GradNa(1,a) )*GradNa(2,a)+
                        ( c(0,5)*GradNa(0,a)+c(4,5)*GradNa(2,a)+c(5,5)*GradNa(1,a) )*GradNa(1,a);
                        
          temp5(0,1) += ( c(0,1)*GradNa(0,a)+c(4,1)*GradNa(2,a)+c(5,1)*GradNa(1,a) )*GradNa(1,a)+
                        ( c(0,3)*GradNa(0,a)+c(4,3)*GradNa(2,a)+c(5,3)*GradNa(1,a) )*GradNa(2,a)+
                        ( c(0,5)*GradNa(0,a)+c(4,5)*GradNa(2,a)+c(5,5)*GradNa(1,a) )*GradNa(0,a);
                        
          temp5(0,2) += ( c(0,2)*GradNa(0,a)+c(4,2)*GradNa(2,a)+c(5,2)*GradNa(1,a) )*GradNa(2,a)+
                        ( c(0,3)*GradNa(0,a)+c(4,3)*GradNa(2,a)+c(5,3)*GradNa(1,a) )*GradNa(1,a)+
                        ( c(0,4)*GradNa(0,a)+c(4,4)*GradNa(2,a)+c(5,4)*GradNa(1,a) )*GradNa(0,a);
          /**/
          temp5(1,0) += ( c(1,0)*GradNa(1,a)+c(3,0)*GradNa(2,a)+c(5,0)*GradNa(0,a) )*GradNa(0,a)+
                        ( c(1,4)*GradNa(1,a)+c(3,4)*GradNa(2,a)+c(5,4)*GradNa(0,a) )*GradNa(2,a)+
                        ( c(1,5)*GradNa(1,a)+c(3,5)*GradNa(2,a)+c(5,5)*GradNa(0,a) )*GradNa(1,a);
                        
          temp5(1,1) += ( c(1,1)*GradNa(1,a)+c(3,1)*GradNa(2,a)+c(5,1)*GradNa(0,a) )*GradNa(1,a)+
                        ( c(1,3)*GradNa(1,a)+c(3,3)*GradNa(2,a)+c(5,3)*GradNa(0,a) )*GradNa(2,a)+
                        ( c(1,5)*GradNa(1,a)+c(3,5)*GradNa(2,a)+c(5,5)*GradNa(0,a) )*GradNa(0,a);
                        
          temp5(1,2) += ( c(1,2)*GradNa(1,a)+c(3,2)*GradNa(2,a)+c(5,2)*GradNa(0,a) )*GradNa(2,a)+
                        ( c(1,3)*GradNa(1,a)+c(3,3)*GradNa(2,a)+c(5,3)*GradNa(0,a) )*GradNa(1,a)+
                        ( c(1,4)*GradNa(1,a)+c(3,4)*GradNa(2,a)+c(5,4)*GradNa(0,a) )*GradNa(0,a);
          /**/
          temp5(2,0) += ( c(2,0)*GradNa(2,a)+c(3,0)*GradNa(1,a)+c(4,0)*GradNa(0,a) )*GradNa(0,a)+
                        ( c(2,4)*GradNa(2,a)+c(3,4)*GradNa(1,a)+c(4,4)*GradNa(0,a) )*GradNa(2,a)+
                        ( c(2,5)*GradNa(2,a)+c(3,5)*GradNa(1,a)+c(4,5)*GradNa(0,a) )*GradNa(1,a);
                        
          temp5(2,1) += ( c(2,1)*GradNa(2,a)+c(3,1)*GradNa(1,a)+c(4,1)*GradNa(0,a) )*GradNa(1,a)+
                        ( c(2,3)*GradNa(2,a)+c(3,3)*GradNa(1,a)+c(4,3)*GradNa(0,a) )*GradNa(2,a)+
                        ( c(2,5)*GradNa(2,a)+c(3,5)*GradNa(1,a)+c(4,5)*GradNa(0,a) )*GradNa(0,a);
                        
          temp5(2,2) += ( c(2,2)*GradNa(2,a)+c(3,2)*GradNa(1,a)+c(4,2)*GradNa(0,a) )*GradNa(2,a)+
                        ( c(2,3)*GradNa(2,a)+c(3,3)*GradNa(1,a)+c(4,3)*GradNa(0,a) )*GradNa(1,a)+
                        ( c(2,4)*GradNa(2,a)+c(3,4)*GradNa(1,a)+c(4,4)*GradNa(0,a) )*GradNa(0,a);
      }
      
      /* temp3 = N_{a,i}*v_{i,j} */
      temp3 = 0.0;
      for (i = 0; i < nsd; i++)
      {
        p = a*ndof + i;
        for (b = 0; b < nun; b++)
        {
          for (j = 0; j < nsd; j++)
          {
            temp3 += GradNa(i,a)*GradVel(i,j);

            q = b*ndof + j;

            /* term : N_{A}*\rho*N_{B}*v_{i,j} */
            fLHS(p,q) += temp1*Na[a]*Density*Na[b]*GradVel(i,j);

            /* term : \tau^{m}*v_{k,old}*N_{A,k}*\rho*N_{B}*v_{i,j} */
            fLHS(p,q) += temp1*tau_m*temp0*Density*Na[b]*GradVel(i,j);

            if(i == j)
            {
              /* term : N_{A}*\rho*v_{l}*N_{B,l} \delta_{ij} */
              fLHS(p,q) += temp1*Na[a]*Density*temp2;

              /* term : \tau^{m}*v_{k,old}*N_{A,k}*\rho*v_{l}*N_{B,l} \delta{ij} */
              fLHS(p,q) += temp1*tau_m*temp0*Density*temp2;
            }

            /* term : N_{A,k}*c_{ikjl}*N_{B,l} */
            fLHS(p,q) += temp5(i,j);
            
            /* term : \tau^{c}*N_{A,i}*N_{B,j} */
            fLHS(p,q) += temp1*tau_c*GradNa(i,a)*GradNa(j,b);
          }
          /* j4th term */
          q = b*ndof + 3;

          /* term : \tau^{m}*v_{k,old}*N_{A,k}*N_{B,i} */
          fLHS(p,q) += temp1*tau_m*temp0*GradNa(i,b);
        }
      }
      /* i4th term */
      p = a*ndof + 3;
      for (b = 0; b < nun; b++)
      {
        for (j = 0; j < nsd; j++)
        {
          q = b*ndof + j;

          /* term : \tau^{c}*N_{A,i}*\rho*N_{B}*v_{i,j} */
          fLHS(p,q) += temp1*tau_c*temp3*Density*Na[b];

          /* term : \tau^{c}*N_{A,i}*\rho*v_{l}*N_{B,l} \delta_{ij} */
          fLHS(p,q) += temp1*tau_c*GradNa(j,a)*Density*temp2;
          
          /* term : N_{A}*N_{B,j} */
          fLHS(p,q) += temp1*Na[a]*GradNa(j,b);
        }
        /* j4th term */
        q = b*ndof + 3;

        /* term : \tau^{c}*N_{A,i}*N_{B,i} */
        fLHS(p,q) += temp1*tau_c*temp4;
      }
    }
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
  sub_list.AddSub("fluid_element_stab_param");
  sub_list.AddSub("fluid_element_block", ParameterListT::OnePlus);
}

/** a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FluidElementT::NewSub(const StringT& name) const
{
  WriteCallLocation("NewSub"); //DEBUG

  /* try construct stabilization paramaters */
  if (name == "fluid_element_stab_param")
  {
    ParameterContainerT* stab_param = new ParameterContainerT(name);
    stab_param->SetListOrder(ParameterListT::Choice);

    /* choice of parameters */
    for (int i = 0; i < NumStabParamCodes; i++)
      stab_param->AddSub(ParameterContainerT(StabParamNames[i]));

    return stab_param;
  }
  else if (name == "fluid_element_nodal_output")
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
	fOldVel_list.Dimension(NumIP());
	fPres_list.Dimension(NumIP());
	fGradVel_list.Dimension(NumIP());
	fGradPres_list.Dimension(NumIP());
	for (int i = 0; i < NumIP(); i++)
  {
    fVel_list[i].Dimension(NumSD());
    fOldVel_list[i].Dimension(NumSD());
		fGradVel_list[i].Dimension(NumSD());
		fGradPres_list[i].Dimension(NumSD());
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
        fNodalOutputCodes[i] = IOBaseT::kAtInc;
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
        fElementOutputCodes[i] = IOBaseT::kAtInc;
    }
  }

  /* stabilization parameter codes */
  const ParameterListT& stab_param = list.GetListChoice(*this, "fluid_element_stab_param");
  for (int i = 0; i < NumStabParamCodes; i++)
    if (stab_param.Name() == StabParamNames[i])
      fStabParam = StabParamCodeT(i);
}

/* extract the list of material parameters */
void FluidElementT::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{
  WriteCallLocation("CollectMaterialInfo"); //DEBUG
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
  if (counts.Sum() == 0)
    return;

  /* initialize */
  counts.Dimension(flags.Length());
  counts = 0;

  /* set output flags */
  if (flags[iNodalCrd] == mode) counts[iNodalCrd] = NumSD();
  if (flags[iNodalVel] == mode) counts[iNodalVel] = NumSD();
  if (flags[iNodalAcc] == mode) counts[iNodalAcc] = NumSD();
  if (flags[iNodalPrs] == mode) counts[iNodalPrs] = FluidElementT::kPressureNDOF;
}

void FluidElementT::SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
  iArrayT& counts) const
{
  WriteCallLocation("SetElementOutputCodes"); //DEBUG
#pragma unused(mode)
#pragma unused(flags)
  if (counts.Sum() != 0)
    ExceptionT::BadInputValue("FluidElementT::SetElementOutputCodes", "not implemented");
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
cout << "\n Inside of FluidElementT::" << loc << endl;
}
