/* $Id: PMLT.cpp,v 1.13 2004-07-15 08:27:59 paklein Exp $ */
#include "PMLT.h"

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include "toolboxConstants.h"

#include "fstreamT.h"
#include "ElementCardT.h"
#include "ShapeFunctionT.h"
#include "eIntegratorT.h"
#include "SolidMaterialT.h"
#include "MaterialListT.h"
#include "iAutoArrayT.h"
#include "SSSolidMatT.h"

using namespace Tahoe;

/* constructor */
PMLT::PMLT(const ElementSupportT& support, const FieldT& field):
	SolidElementT(support),
	fNeedsOffset(-1),
	fGradU(NumSD()),
	fLHSa(fLHS.Format()),
	fLHSb(fLHS.Format()),
	fDa(dSymMatrixT::NumValues(NumSD())),
	fDb(dSymMatrixT::NumValues(NumSD())),
	fdt(ElementSupport().TimeStep())

{
ExceptionT::GeneralFail("PMLT::PMLT", "out of date");
#if 0	
#if __option(extended_errorcheck)
	if (NumSD() > 2)
	{
		cout << "\n PML is implemented only for 2D geometries";
		throw ExceptionT::kBadInputValue;
	} 
#endif
	if (fStrainDispOpt != kStandardB ) throw ExceptionT::kBadInputValue;

//NumDOF() set correctly from the nodes
//	NumDOF() = NumSD()*NumSD();             
//	fNumElemEqnos = NumElementNodes()*NumDOF();
	fNEESub = NumSD()*NumElementNodes();
#endif
}

/* called immediately after constructor */
void PMLT::Initialize(void)
{
ExceptionT::GeneralFail("PMLT::Initialize", "out of date");
#if 0	
	/* inherited */
	SolidElementT::Initialize();
	
//	fDOFvec.Dimension(NumDOF());
		
// Q:	fLocDisp.Dimension(NumElementNodes(), NumDOF()); Is this needed here?
	fTotDisp.Dimension(NumElementNodes(),NumSD());
	fTotLastDisp.Dimension(NumElementNodes(),NumSD());
	fTotVel.Dimension(NumElementNodes(),NumSD());
	fTotAcc.Dimension(NumElementNodes(), NumSD());

#if 0
	Field().RegisterLocal(fTotDisp);
	Field().RegisterLocal(fTotLastDisp);
	Field().RegisterLocal(fTotVel);
	Field().RegisterLocal(fTotAcc);
#endif
//only need to register arrays that need to collect values for the
//element from a global array

	/* allocates decomposition of strain-displacement matrix */
	fBa.Dimension(dSymMatrixT::NumValues(NumSD()), NumSD()*NumElementNodes());
	fBb.Dimension(dSymMatrixT::NumValues(NumSD()), NumSD()*NumElementNodes());

	/* resizes LHS and RHS components*/
	fLHSa.Dimension(fNEESub);
	fLHSb.Dimension(fNEESub);
	
	fRHSa.Set(fNEESub,fRHS.Pointer());
	fRHSb.Set(fNEESub,fRHS.Pointer()+fNEESub);
		
	fBody.Dimension(NumSD());

		/* what's needed */
	bool need_strain = false;
	bool need_strain_last = false;
	for (int i = 0; i < fMaterialNeeds.Length(); i++)
	{
		const ArrayT<bool>& needs = fMaterialNeeds[i];
		need_strain = need_strain || needs[fNeedsOffset + kstrain];
		need_strain_last = need_strain_last || needs[fNeedsOffset + kstrain_last];
	}

	/* allocate deformation gradient list */
	if (need_strain)
	{
		fGradU_List.Dimension(NumIP());
		for (int i = 0; i < NumIP(); i++)
			fGradU_List[i].Dimension(NumSD());
	}
	
	/* allocate "last" deformation gradient list */
	if (need_strain_last)
	{
		fGradU_last_List.Dimension(NumIP());
		for (int i = 0; i < NumIP(); i++)
			fGradU_last_List[i].Dimension(NumSD());
	}
#endif
}

/* construct the field */
void PMLT::NodalDOFs(const iArrayT& nodes, dArray2DT& DOFs) const
{
#pragma unused(nodes)
#pragma unused(DOFs)

	cout << "\n PMLT::NodalDOFs: not implemented" << endl;
	throw ExceptionT::kGeneralFail;

#if 0
	if (nodes.Length() == fNumElementNodes)
	{
		fLocDisp.SetLocal(nodes);
	
		//construct total displacement in DOFs
	}
	else
	{
	
	
	}


#endif	
	
}

/* construct the effective mass matrix */
void PMLT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
	/* inherited */
	SolidElementT::LHSDriver(sys_type);

	/* element contribution */
	ElementLHSDriver();
}
void PMLT::ElementLHSDriver(void)
{
	/* set components and weights */
	double constM = 0.0;
	double constC = 0.0;
	double constK = 0.0;
	
	int formM = fIntegrator->FormM(constM);
	int formC = fIntegrator->FormC(constC);
	int formK = fIntegrator->FormK(constK);

	/* override algorithm */
	if (fMassType == kNoMass) formM = 0;

	/* quick exit */
	if ((formM == 0 && formC == 0 && formK == 0) ||
	    (fabs(constM) < kSmall &&
	     fabs(constC) < kSmall &&
	     fabs(constK) < kSmall)) return;

	/* loop over elements */
	Top();
	while (NextElement())
	{
		double constKe = constK;
		double constMe = constM;
	
		/* initialize */
		fLHS = 0.0;
		
		/* set shape function derivatives */
		SetGlobalShape();

		if (fabs(constMe) > kSmall)
		{
			FormMass(fMassType, constMe*(fCurrMaterial->Density()));
//Q: Should modify to check if damping factor is non zero before forming Cv
			FormDamping(fMassType, constMe*(fCurrMaterial->Density()));
		}
		/* element stiffness */
		if (fabs(constKe) > kSmall)
			FormStiffness(constKe);
	
		/* add to global equations */
		AssembleLHS();		
	}
}

void PMLT::RHSDriver(void)
{
	/* inherited */
	SolidElementT::RHSDriver();

	/* element contribution */
	ElementRHSDriver();
}

/* form the residual force vector */
void PMLT::ElementRHSDriver(void)
{
	/* set components and weights */
	double constMa = 0.0;
	double constCv = 0.0;
	double constKd = 0.0;
	
	/* components dicated by the algorithm */
	int formMa = fIntegrator->FormMa(constMa);
	int formCv = fIntegrator->FormCv(constCv);
	int formKd = fIntegrator->FormKd(constKd);

	/* body forces */
	int formBody = 0;
	if (fMassType != kNoMass &&
	   (fBodySchedule && fBody.Magnitude() > kSmall))
	{	
		cout << "\nPML does not yet incorporate body forces";
		if (!formMa) constMa = 1.0; /* override */
	}

	/* override controller */
	if (fMassType == kNoMass) formMa = 0;

	Top();
	while (NextElement())
	{
		fRHS = 0.0;
		
		/* global shape function values */
		SetGlobalShape();
			
		/* internal force contribution */	
		if (formKd) FormKd(-1.0);
				
		/* damping */
		if (formCv)
		{
			SetLocalU(fLocVel);
			FormCv_PML(fMassType, -(fCurrMaterial->Density()), fLocAcc);			  		
		}

		/* inertia forces */
		if (formMa)
		{
			SetLocalU(fLocAcc);
			FormMa(fMassType, -(fCurrMaterial->Density()), fLocAcc);			  		
		}
								
		/* assemble */
		AssembleRHS();
	}
}

void PMLT::FormMass(int mass_type, double constM)
{
#if __option(extended_errorcheck)
	if (fLocDisp.Length() != fLHS.Rows()) throw ExceptionT::kSizeMismatch;
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
			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();
			
			int nen = NumElementNodes();
			
			/* matrix form */
			int a = 0, zero = 0;
			int& b_start = (fLHS.Format() == ElementMatrixT::kSymmetricUpper) ? a : zero;
			
			fShapes->TopIP();	
			while ( fShapes->NextIP() )
			{
				double temp = constM*(*Weight++)*(*Det++);
				const double* Na = fShapes->IPShapeU();
								
				for (a = 0; a < nen; a++)
					for (int i = 0; i < NumSD(); i++)
					{
						int p = a*NumSD() + i;
						
						/* upper triangle only */
						for (int b = b_start; b < nen; b++)
							for (int j = 0; j < NumSD(); j++)
								if(i == j)
								{									
									int q = b*NumSD() + j;
									fLHSa(p,q) += temp*Na[a]*Na[b];
									fLHSb(p,q) += temp*Na[a]*Na[b];
								}
					}
			}
	
			double* p1 = fLHS.Pointer();
			double* p2 = p1+((fLHS.Rows()+1)*fNEESub);
			double* pa = fLHSa.Pointer();
			double* pb = fLHSb.Pointer();
	
			for (int cols = 0; cols < fNEESub; cols++)
			{
				for (int rows = 0; rows< fNEESub; rows++)
				{
					*p1++ = *pa++;
					*p2++ = *pb++;
				}
				p1 += fNEESub;
				p2 += fNEESub;
			}
			break;
		}
			
		case kLumpedMass:	/* lumped mass */
		{
			int nen = fLocDisp.NumberOfNodes();

		    double dsum   = 0.0;
		    double totmas = 0.0;
		    fNEEvec = 0.0;

			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();

			/* total mass and diagonal sum */
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				double temp1     = constM*(*Weight++)*(*Det++);
				const double* Na = fShapes->IPShapeU();

				totmas += temp1;
				for (int lnd = 0; lnd < nen; lnd++)
				{
					double temp2 = temp1*Na[lnd]*Na[lnd];
					dsum += temp2;
					fNEEvec[lnd] += temp2;
				}
			}	
				
			/* scale diagonal to conserve total mass */
			double diagmass = totmas/dsum;
			
			/* lump mass onto diagonal */
			int inc = fLHS.Rows() + 1;

			double* p1 = fLHS.Pointer();
			double* p2 = p1+((fLHS.Rows()+1)*fNEESub);

			for (int lnd = 0; lnd < nen; lnd++)
			{
				double temp = diagmass*fNEEvec[lnd];
				for (int ed = 0; ed < NumSD(); ed++)
				{
					*p1 += temp;
					p1 += inc;	
					*p2 += temp;
					p2 += inc;
				}
			}
			break;
		}			
		default:
		
			cout << "\n Elastic::FormMass: unknown mass matrix code\n" << endl;
			throw ExceptionT::kBadInputValue;
	}
}

void PMLT::FormDamping(int mass_type, double constC)
{
#if __option(extended_errorcheck)
	if (fLocDisp.Length() != fLHS.Rows()) throw ExceptionT::kSizeMismatch;
#endif
	//To be replaced later by C1Functions
	double damp_a = 0.0;
	double damp_b = 0.0;
	
	switch (mass_type)
	{
		case kNoMass:			/* no mass matrix */
		
			break;
		
		case kConsistentMass:	/* consistent mass	*/
		{
			// integration of the element mass is done
			// in the reference configuration since density
			// is mass/(undeformed volume)
			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();
			
			int nen = NumElementNodes();
			
			/* matrix form */
			int a = 0, zero = 0;
			int& b_start = (fLHS.Format() == ElementMatrixT::kSymmetricUpper) ? a : zero;
			
			fShapes->TopIP();	
			while ( fShapes->NextIP() )
			{
				double temp = constC*(*Weight++)*(*Det++);
				const double* Na = fShapes->IPShapeU();
								
				for (a = 0; a < nen; a++)
					for (int i = 0; i < NumSD(); i++)
					{
						int p = a*NumSD() + i;
						
						/* upper triangle only */
						for (int b = b_start; b < nen; b++)
							for (int j = 0; j < NumSD(); j++)
								if(i == j)
								{									
									int q = b*NumSD() + j;
									fLHSa(p,q) += damp_a*temp*Na[a]*Na[b];
									fLHSb(p,q) += damp_b*temp*Na[a]*Na[b];
								}
					}
			}
	
			double* p1 = fLHS.Pointer();
			double* p2 = p1+((fLHS.Rows()+1)*fNEESub);
			double* pa = fLHSa.Pointer();
			double* pb = fLHSb.Pointer();
	
			for (int cols = 0; cols < fNEESub; cols++)
			{
				for (int rows = 0; rows< fNEESub; rows++)
				{
					*p1++ = *pa++;
					*p2++ = *pb++;
				}
				p1 += fNEESub;
				p2 += fNEESub;
			}
			break;
		}
			
		case kLumpedMass:	/* lumped mass */
		{
			int nen = fLocDisp.NumberOfNodes();

		    double dsum   = 0.0;
		    double totmas = 0.0;
		    fNEEvec = 0.0;

			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();

			/* total mass and diagonal sum */
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				double temp1     = constC*(*Weight++)*(*Det++);
				const double* Na = fShapes->IPShapeU();

				totmas += temp1;
				for (int lnd = 0; lnd < nen; lnd++)
				{
					double temp2 = temp1*Na[lnd]*Na[lnd];
					dsum += temp2;
					fNEEvec[lnd] += temp2;
				}
			}	
				
			/* scale diagonal to conserve total mass */
			double diagmass = totmas/dsum;
			
			/* lump mass onto diagonal */
			int inc = fLHS.Rows() + 1;

			double* p1 = fLHS.Pointer();
			double* p2 = p1+((fLHS.Rows()+1)*fNEESub);

			for (int lnd = 0; lnd < nen; lnd++)
			{
				double temp = diagmass*fNEEvec[lnd];
				for (int ed = 0; ed < NumSD(); ed++)
				{
					*p1 += temp * damp_a;
					p1 += inc;	
					*p2 += temp * damp_b;
					p2 += inc;
				}
			}
			break;
		}			
		default:
		
			cout << "\n Elastic::FormMass: unknown mass matrix code\n" << endl;
			throw ExceptionT::kBadInputValue;
	}
}
/* form the element stiffness matrix */
void PMLT::FormStiffness(double constK)
{
	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;
		
   //Q: How is symmetry accounted for?

	/* integrate element stiffness */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/* clear */
	fLHSa = 0.0;
	fLHSb = 0.0;
	
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		double scale = constK*(*Det++)*(*Weight++);
	
		/* strain displacement matrix */
		Set_B(fShapes->Derivatives_U() , fB);
		Ba(fBa,fB);
		Bb(fBb,fB);

		double dampa = 0.0;
		double dampb = 0.0;
		
		double scale_a = scale/(1+dampa*fdt*0.5);
		double scale_b = scale/(1+dampb*fdt*0.5);
		/* get D matrix */
		fDa.SetToScaled(scale_a, fCurrMaterial->c_ijkl());
		fDb.SetToScaled(scale_b, fCurrMaterial->c_ijkl());
							
		/* multiply b_alpha(transpose) * d*b */		
		fLHSa.MultATBC(fBa, fDa, fBa, dMatrixT::kWhole, dMatrixT::kAccumulate);	
		fLHSa.MultATBC(fBa, fDb, fBb, dMatrixT::kWhole, dMatrixT::kAccumulate);
	
		fLHSb.MultATBC(fBb, fDa, fBa, dMatrixT::kWhole, dMatrixT::kAccumulate);	
		fLHSb.MultATBC(fBb, fDb, fBb, dMatrixT::kWhole, dMatrixT::kAccumulate);
	/*Map into fLHS*/
	}
	int next = NumElementNodes()*NumDOF()*fNEESub;
	double* p1 = fLHS.Pointer();
	double* p2 = p1+fNEESub;
	double* pa = fLHSa.Pointer();
	double* pb = fLHSb.Pointer();
	
	for (int cols = 0; cols < fNEESub; cols++)
	{
		for (int rows = 0; rows< fNEESub; rows++)
		{
			*(p1+next)=*pa;
			*p1++ = *pa++;
			*(p2+next)=*pa;
			*p2++ = *pb++;
		}
		p1 += fNEESub;
		p2 += fNEESub;

	}
}

/* solution calls */
void PMLT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
	/* not my field */
	if (&field != &(Field())) return;

	/* quick exit */
	bool hasnode = false;
	for (int i=0; i < fBlockData.Length() && !hasnode; i++)
		if (fConnectivities[i]->HasValue(node)) hasnode = true;
	if (!hasnode) return;

	/* set components and weights */
	double constMa = 0.0;
	double constCv = 0.0;
	double constKd = 0.0;
	
	/* components dicated by the algorithm */
	int formMa = fIntegrator->FormMa(constMa);
	int formCv = fIntegrator->FormCv(constCv);
	int formKd = fIntegrator->FormKd(constKd);

	/* body forces */
	int formBody = 0;
	if (fMassType != kNoMass && 
	   (fBodySchedule && fBody.Magnitude() > kSmall))
	{	
		cout << "\nWarning: Body forces not yet implemented in PML";
		if (!formMa) constMa = 1.0; /* override */
	}

	/* override controller */
	if (fMassType == kNoMass) formMa = 0;
	
	/* temp for nodal force */
	dArrayT nodalforce;
	
	Top();
	while (NextElement())
	{
		int nodeposition;
		if (CurrentElement().NodesU().HasValue(node, nodeposition))
		{
			/* initialize */
			fRHS = 0.0;
	
			/* effective accelerations and displacements */
			//ComputeEffectiveDVA(formBody, formMa, constMa, formCv, constCv, formKd, constKd);
			//DEV - Rayleigh damping is poorly formulated

			/* global shape function values */
			SetGlobalShape();
	
			/* internal force contribution */	
			if (formKd) FormKd(1.0);
				
			if (formCv)
			{
				SetLocalU(fLocVel); 
				FormCv_PML(fMassType, fCurrMaterial->Density(), fLocVel);
			}
	
			/* inertia forces */
			if (formMa)
			{
				SetLocalU(fLocAcc);
				FormMa(fMassType, fCurrMaterial->Density(), fLocAcc);
			}

			/* components for node */
			nodalforce.Set(NumDOF(), &fRHS[NumDOF()*nodeposition]);
	
			/* accumulate */
			force += nodalforce;
		}
	}
}
			
/* calculate the body force contribution */
void PMLT::FormMa(int mass_type, double constM, const LocalArrayT& body_force)
{
	int nen = NumElementNodes();
	int next = nen*NumSD();
	LocalArrayT body_forcea;
	LocalArrayT body_forceb;
	body_forcea.Alias(nen,NumSD(),body_force.Pointer());
	body_forceb.Alias(nen,NumSD(),body_force.Pointer()+next);
	
	switch (mass_type)
	{
		case kConsistentMass:	
		{
#if __option(extended_errorcheck)
			if (fRHS.Length() != body_force.Length()) throw ExceptionT::kSizeMismatch;
#endif
						
			dArrayT fNSDveca(NumSD());
			dArrayT fNSDvecb(NumSD());
			
			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();

			fShapes->TopIP();
			while ( fShapes->NextIP() )
			{					
				/* integration point accelerations */
				fShapes->InterpolateU(body_forcea, fNSDveca);
				fShapes->InterpolateU(body_forceb, fNSDvecb);

				/* accumulate in element residual force vector */				
				double*	resa      = fRHSa.Pointer();
				double* resb      = fRHSb.Pointer();
				const double* Na = fShapes->IPShapeU();

				double temp = constM*(*Weight++)*(*Det++);				
				for (int lnd = 0; lnd < nen; lnd++)
				{
					double  temp2 = temp*(*Na++);
					double* pacca = fNSDveca.Pointer();
					double* paccb = fNSDvecb.Pointer();
					
					for (int nsd = 0; nsd < NumSD(); nsd++)			
						*resa++ += temp2*(*pacca++);
						*resb++ += temp2*(*paccb++);
				}
			}
			break;
		}
		case kLumpedMass:
		{
			//cout << "\n ContinuumElementT::FormMa: inertial forces with lumped mass not supported";
			//cout << endl;
			//throw ExceptionT::kGeneralFail;
			
			//for now, no inertial force for lumped mass
			//but should probably generalize the FormMass and
			//FormStiffness routines by passing in a target object
			//in which to place the data

#if __option(extended_errorcheck)
			if (fLHS.Rows() != body_force.Length()) throw ExceptionT::kSizeMismatch;
#endif
			fLHS = 0.0; //hope there's nothing in there!
			FormMass(kLumpedMass, constM);

			dArrayT fNEEveca(fNEESub);
			dArrayT fNEEvecb(fNEESub);
			
			body_forcea.ReturnTranspose(fNEEveca);
			body_forceb.ReturnTranspose(fNEEvecb);

			double* pAcca = fNEEveca.Pointer();
			double* pAccb = fNEEvecb.Pointer();
			double* pResa = fRHSa.Pointer();
			double* pResb = fRHSb.Pointer();
			int     massdexa = 0;
			for (int i = 0; i < fNEESub; i++)
			{
				int massdexb = massdexa + fNEESub;
				*pResa++ += (*pAcca++)*fLHS(massdexa,massdexa);
				*pResb++ += (*pAccb++)*fLHS(massdexb,massdexb);
				massdexa++;
			}
			
			break;
		}
	}
}

/* calculate the body force contribution */
void PMLT::FormCv_PML(int mass_type, double constC, const LocalArrayT& body_force)
{
	int nen = NumElementNodes();
	int next = nen*NumSD();
	LocalArrayT body_forcea;
	LocalArrayT body_forceb;
	body_forcea.Alias(nen,NumSD(),body_force.Pointer());
	body_forceb.Alias(nen,NumSD(),body_force.Pointer()+next);

	switch (mass_type)
	{
		case kConsistentMass:	
		{
#if __option(extended_errorcheck)
			if (fRHS.Length() != body_force.Length()) throw ExceptionT::kSizeMismatch;
#endif
						
			dArrayT fNSDveca(NumSD());
			dArrayT fNSDvecb(NumSD());
			
		//To be replaced later by C1Functions
			double damp_a = 0.0;
			double damp_b = 0.0;
		
			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();

			fShapes->TopIP();
			while ( fShapes->NextIP() )
			{					
				/* integration point velocities */
				fShapes->InterpolateU(body_forcea, fNSDveca);
				fShapes->InterpolateU(body_forceb, fNSDvecb);

				/* accumulate in element residual force vector */				
				double*	resa      = fRHSa.Pointer();
				double* resb      = fRHSb.Pointer();
				const double* Na = fShapes->IPShapeU();

				double temp = constC*(*Weight++)*(*Det++);				
				for (int lnd = 0; lnd < nen; lnd++)
				{
					double  temp2 = temp*(*Na++);
					double* pvela = fNSDveca.Pointer();
					double* pvelb = fNSDvecb.Pointer();
					
					for (int nsd = 0; nsd < NumSD(); nsd++)			
						*resa++ += damp_a*temp2*(*pvela++);
						*resb++ += damp_b*temp2*(*pvelb++);
				}
			}
			break;
		}	
		case kLumpedMass:
		{
			//cout << "\n ContinuumElementT::FormMa: inertial forces with lumped mass not supported";
			//cout << endl;
			//throw ExceptionT::kGeneralFail;
			
			//for now, no inertial force for lumped mass
			//but should probably generalize the FormMass and
			//FormStiffness routines by passing in a target object
			//in which to place the data

#if __option(extended_errorcheck)
			if (fLHS.Rows() != body_force.Length()) throw ExceptionT::kSizeMismatch;
#endif
			fLHS = 0.0; //hope there's nothing in there!
			FormDamping(kLumpedMass, constC);
			dArrayT fNEEveca(fNEESub);
			dArrayT fNEEvecb(fNEESub);
			
			body_forcea.ReturnTranspose(fNEEveca);
			body_forceb.ReturnTranspose(fNEEvecb);

			double* pVela = fNEEveca.Pointer();
			double* pVelb = fNEEvecb.Pointer();
			double* pResa = fRHSa.Pointer();
			double* pResb = fRHSb.Pointer();
			int     dampdexa = 0;
			for (int i = 0; i < fNEESub; i++)
			{
				int dampdexb = dampdexa + fNEESub;
				*pResa++ += (*pVela++)*fLHS(dampdexa,dampdexa);
				*pResb++ += (*pVelb++)*fLHS(dampdexb,dampdexb);
				dampdexa++;
			}
			
			break;
		}
	}
}
void PMLT::FormKd(double constK)
{
	dArrayT fNEEveca(fNEESub);
	dArrayT fNEEvecb(fNEESub);

	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
		
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		/* get strain-displacement matrix */
		Set_B(fShapes->Derivatives_U(), fB);
		Ba(fBa,fB);
		Bb(fBb,fB);
		/* B^T * Cauchy stress */
		fBa.MultTx(fCurrMaterial->s_ij(), fNEEveca);
		fBb.MultTx(fCurrMaterial->s_ij(), fNEEvecb);
		
		/* accumulate */
		fRHSa.AddScaled(constK*(*Weight++)*(*Det++), fNEEveca);
		fRHSb.AddScaled(constK*(*Weight++)*(*Det++), fNEEvecb);
	}	
}

void PMLT::AddLinearMomentum(dArrayT& momentum)
{
	/* check */
	if (momentum.Length() != NumSD()) throw ExceptionT::kSizeMismatch;
		
	/* loop over elements */
	dArrayT vec_ndof(NumDOF());
	Top();
	while (NextElement())
	{
		/* global shape function derivatives, jacobians, local coords */
		SetGlobalShape();
		
		/* get velocities */
		SetLocalU(fLocVel);
		LocalArrayT& LocVel = TotalVel(fLocVel);

		/* material density */
		double density = fCurrMaterial->Density();

		/* integration */
		const double* Det    = fShapes->IPDets();
		const double* Weight = fShapes->IPWeights();
	
		fShapes->TopIP();
		while ( fShapes->NextIP() )
		{					
			double temp  = density*(*Det++)*(*Weight++);

			/* integration point velocities */
			fShapes->InterpolateU(LocVel, vec_ndof);

			double* p    = momentum.Pointer();
			double* pvel = vec_ndof.Pointer();					
			for (int nsd = 0; nsd < NumSD(); nsd++)			
				*p++ += temp*(*pvel++);
		}
	}
}

void PMLT::SendOutput(int kincode)
{
	/* output flags */
	iArrayT flags(fNodalOutputCodes.Length());

	/* set flags to get desired output */
	flags = IOBaseT::kAtNever;
	switch (kincode)
	{
		case iNodalDisp:
		    flags[iNodalDisp] = NumSD();
			break;
		case iNodalStress:
		    flags[iNodalStress] = 1;
			break;
		case iEnergyDensity:
		    flags[iEnergyDensity] = 1;
			break;
		case iPrincipal:
			flags[iPrincipal] = NumSD();
			break;
		default:
			cout << "\n SolidElementT::SendKinematic: invalid output code: ";
			cout << kincode << endl;
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
/* construct output labels array */
void PMLT::SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
	/* initialize */
	counts.Dimension(flags.Length());
	counts = 0;

	/* set output flags */
	if (flags[iNodalCoord] == mode)
		counts[iNodalCoord] = NumSD();
	if (flags[iNodalDisp] == mode)
		counts[iNodalDisp] = NumSD();
	if (flags[iNodalStress] == mode)
		counts[iNodalStress] = dSymMatrixT::NumValues(NumSD());
	if (flags[iPrincipal] == mode)
		counts[iPrincipal] = NumSD();
	if (flags[iEnergyDensity] == mode)
		counts[iEnergyDensity] = 1;
	if (flags[iWaveSpeeds] == mode)
		counts[iWaveSpeeds] = NumSD();
	if (flags[iMaterialData] == mode)
		counts[iMaterialData] = (*fMaterialList)[0]->NumOutputVariables();
}
void PMLT::SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
	/* initialize */
	counts.Dimension(flags.Length());
	counts = 0;

	/* set output flags */
	if (fElementOutputCodes[iCentroid] == mode) counts[iCentroid] = NumSD();
	if (fElementOutputCodes[iMass] == mode) counts[iMass] = 1;
	if (fElementOutputCodes[iStrainEnergy] == mode) counts[iStrainEnergy] = 1;
	if (fElementOutputCodes[iKineticEnergy] == mode) counts[iKineticEnergy] = 1;
	if (fElementOutputCodes[iLinearMomentum] == mode) counts[iLinearMomentum] = NumSD();
	if (fElementOutputCodes[iIPStress] == mode) counts[iIPStress] = dSymMatrixT::NumValues(NumSD())*NumIP();
	if (fElementOutputCodes[iIPMaterialData] == mode) 
		counts[iIPMaterialData] = (*fMaterialList)[0]->NumOutputVariables()*NumIP();
}
void PMLT::GenerateOutputLabels(const iArrayT& n_codes, ArrayT<StringT>& n_labels, 
	const iArrayT& e_codes, ArrayT<StringT>& e_labels) const
{
	/* allocate */
	n_labels.Dimension(n_codes.Sum());

	int count = 0;
	if (n_codes[iNodalDisp])
	{
		/* labels from the field */
		const ArrayT<StringT>& labels = Field().Labels();
		for (int i = 0; i < labels.Length(); i++)
			n_labels[count++] = labels[i];
	}

	if (n_codes[iNodalCoord])
	{
		const char* xlabels[] = {"x1", "x2", "x3"};
		for (int i = 0; i < NumSD(); i++)
			n_labels[count++] = xlabels[i];
	}

	if (n_codes[iNodalStress])
	{
		const char* slabels2D[] = {"s11", "s22", "s12"};
		const char* slabels3D[] = {"s11", "s22", "s33", "s23", "s13", "s12"};
		const char**    slabels = (NumSD() == 2) ? slabels2D : slabels3D;
		for (int i = 0; i < dSymMatrixT::NumValues(NumSD()); i++)
			n_labels[count++] = slabels[i];
	}
		
	if (n_codes[iPrincipal])
	{
		const char* plabels[] = {"s1", "s2", "s3"};
		for (int i = 0; i < NumSD(); i++)
			n_labels[count++] = plabels[i];
	}
		
	if (n_codes[iEnergyDensity]) n_labels[count++] = "phi";
	if (n_codes[iWaveSpeeds])
	{
		const char* clabels2D[] = {"cd", "cs"};
		const char* clabels3D[] = {"cd", "cs_min", "cs_max"};
		const char**    clabels = (NumSD() == 2) ? clabels2D : clabels3D;
		for (int i = 0; i < NumSD(); i++)
			n_labels[count++] = clabels[i];		
	}

	/* material output labels */
	if (n_codes[iMaterialData])
	{
		ArrayT<StringT> matlabels;
		(*fMaterialList)[0]->OutputLabels(matlabels);	
		
		for (int i = 0; i < matlabels.Length(); i++)
			n_labels[count++] = matlabels[i];
	}

	/* allocate */
	e_labels.Dimension(e_codes.Sum());
	count = 0;
	if (e_codes[iCentroid])
	{
		const char* xlabels[] = {"xc_1", "xc_2", "xc_3"};
		for (int i = 0; i < NumSD(); i++)
			e_labels[count++] = xlabels[i];
	}
	if (e_codes[iMass]) e_labels[count++] = "mass";
	if (e_codes[iStrainEnergy]) e_labels[count++] = "U";
	if (e_codes[iKineticEnergy]) e_labels[count++] = "T";
	if (e_codes[iLinearMomentum])
	{
		const char* plabels[] = {"L_X", "L_Y", "L_Z"};
		for (int i = 0; i < NumSD(); i++)
			e_labels[count++] = plabels[i];
	}
	if (e_codes[iIPStress])
	{
		const char* slabels2D[] = {"s11", "s22", "s12"};
		const char* slabels3D[] = {"s11", "s22", "s33", "s23", "s13", "s12"};
		const char**    slabels = (NumSD() == 2) ? slabels2D : slabels3D;

		/* over integration points */
		for (int j = 0; j < NumIP(); j++)
		{
			StringT ip_label;
			ip_label.Append("ip", j+1);
			
			/* over stress components */
			for (int i = 0; i < dSymMatrixT::NumValues(NumSD()); i++)
			{
				e_labels[count].Clear();
				e_labels[count].Append(ip_label, ".", slabels[i]);
				count++;
			}
		}		
	}

	/* material output labels */
	if (e_codes[iIPMaterialData])
	{
		ArrayT<StringT> matlabels;
		(*fMaterialList)[0]->OutputLabels(matlabels);	

		/* over integration points */
		for (int j = 0; j < NumIP(); j++)
		{
			StringT ip_label;
			ip_label.Append("ip", j+1);
			
			/* over stress components */
			for (int i = 0; i < matlabels.Length(); i++)
			{
				e_labels[count].Clear();
				e_labels[count].Append(ip_label, ".", matlabels[i]);
				count++;
			}
		}		
	}
}
/* extrapolate the integration point stresses and strains and extrapolate */
void PMLT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	const iArrayT& e_codes, dArray2DT& e_values)
{
	/* number of output values */
	int n_out = n_codes.Sum();
	int e_out = e_codes.Sum();

	/* nothing to output */
	if (n_out == 0 && e_out == 0) return;

	/* reset averaging workspace */
	ElementSupport().ResetAverage(n_out);
	
	/* allocate element results space */
	e_values.Dimension(NumElements(), e_out);

	/* nodal work arrays */
	dArray2DT nodal_space(NumElementNodes(), n_out);
	dArray2DT nodal_all(NumElementNodes(), n_out);
	dArray2DT coords, disp;
	dArray2DT nodalstress, princstress, matdat;
	dArray2DT energy, speed;

	/* ip values */
	dSymMatrixT cauchy(NumSD());
	dArrayT ipmat(n_codes[iMaterialData]), ipenergy(1);
	dArrayT ipspeed(NumSD()), ipprincipal(NumSD());

	/* set shallow copies */
	double* pall = nodal_space.Pointer();
	coords.Set(NumElementNodes(), n_codes[iNodalCoord], pall)      ; pall += coords.Length();
	disp.Set(NumElementNodes(), n_codes[iNodalDisp], pall)         ; pall += disp.Length();
	nodalstress.Set(NumElementNodes(), n_codes[iNodalStress], pall); pall += nodalstress.Length();
	princstress.Set(NumElementNodes(), n_codes[iPrincipal], pall)  ; pall += princstress.Length();
	energy.Set(NumElementNodes(), n_codes[iEnergyDensity], pall)   ; pall += energy.Length();
	speed.Set(NumElementNodes(), n_codes[iWaveSpeeds], pall)       ; pall += speed.Length();
	matdat.Set(NumElementNodes(), n_codes[iMaterialData], pall);

	/* element work arrays */
	dArrayT element_values(e_values.MinorDim());
	pall = element_values.Pointer();
	dArrayT centroid, ip_centroid;
	if (e_codes[iCentroid])
	{
		centroid.Set(NumSD(), pall); pall += NumSD();
		ip_centroid.Dimension(NumSD());
	}
	double m_tmp, w_tmp, ke_tmp;
	double& mass = (e_codes[iMass]) ? *pall++ : m_tmp;
	double& strain_energy = (e_codes[iStrainEnergy]) ? *pall++ : w_tmp;
	double& kinetic_energy = (e_codes[iKineticEnergy]) ? *pall++ : ke_tmp;
	dArrayT linear_momentum, ip_velocity;
	if (e_codes[iLinearMomentum])
	{
		linear_momentum.Set(NumSD(), pall); pall += NumSD();
		ip_velocity.Dimension(NumSD());
	}
	dArray2DT ip_stress;
	if (e_codes[iIPStress])
	{
		ip_stress.Set(NumIP(), e_codes[iIPStress]/NumIP(), pall);
		pall += ip_stress.Length();
	}
	dArray2DT ip_material_data;
	if (e_codes[iIPMaterialData])
	{
		ip_material_data.Set(NumIP(), e_codes[iIPMaterialData]/NumIP(), pall);
		pall += ip_material_data.Length();
		ipmat.Dimension(ip_material_data.MinorDim());
	}

	/* check that degrees are displacements */
//Q: What are interpolant_DOFs
	int interpolant_DOF = InterpolantDOFs();

	Top();
	while (NextElement())
	{
		/* initialize */
	    nodal_space = 0.0;

		/* global shape function values */
		SetGlobalShape();
		
		/* collect nodal values */
		if (e_codes[iKineticEnergy] || e_codes[iLinearMomentum])
			SetLocalU(fLocVel);
			
		/* coordinates and displacements all at once */
		if (n_codes[iNodalCoord]) fLocInitCoords.ReturnTranspose(coords);
		if (n_codes[ iNodalDisp])
		{
			if (interpolant_DOF)
			{
				LocalArrayT& TotDisp = TotalDisp(fLocDisp);
				TotDisp.ReturnTranspose(disp);
			}
			else
			{
//Q:  Help
				NodalDOFs(CurrentElement().NodesX(), disp);
			}
		}
		
		/* initialize element values */
		mass = strain_energy = kinetic_energy = 0;
		if (e_codes[iCentroid]) centroid = 0.0;
		if (e_codes[iLinearMomentum]) linear_momentum = 0.0;
		const double* j = fShapes->IPDets();
		const double* w = fShapes->IPWeights();
		double density = fCurrMaterial->Density();

		/* integrate */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/* element integration weight */
			double ip_w = (*j++)*(*w++);
		
			/* get Cauchy stress */
			cauchy = fCurrMaterial->s_ij();

			/* stresses */
			if (n_codes[iNodalStress]) fShapes->Extrapolate(cauchy, nodalstress);
			if (e_codes[iIPStress]) ip_stress.SetRow(fShapes->CurrIP(), cauchy);

			/* wave speeds */
			if (n_codes[iWaveSpeeds])
			{
				/* acoustic wave speeds */
				fCurrMaterial->WaveSpeeds(fNormal, ipspeed);
				fShapes->Extrapolate(ipspeed, speed);
			}

			/* principal values - compute principal before smoothing */
			if (n_codes[iPrincipal])
			{
				/* compute eigenvalues */
				cauchy.PrincipalValues(ipprincipal);
				fShapes->Extrapolate(ipprincipal, princstress);	
			}

			/* strain energy density */
			if (n_codes[iEnergyDensity] || e_codes[iStrainEnergy])
			{
				double ip_strain_energy = fCurrMaterial->StrainEnergyDensity();
			
				/* nodal average */
				if (n_codes[iEnergyDensity])
				{
					ipenergy[0] = ip_strain_energy;
					fShapes->Extrapolate(ipenergy,energy);
				}
				
				/* integrate over element */
				if (e_codes[iStrainEnergy])
					strain_energy += ip_w*ip_strain_energy;
			}

			/* material stuff */
			if (n_codes[iMaterialData] || e_codes[iIPMaterialData])
			{
				/* compute material output */
				fCurrMaterial->ComputeOutput(ipmat);
				
				/* store nodal data */
				if (n_codes[iMaterialData]) fShapes->Extrapolate(ipmat, matdat);
				
				/* store element data */
				if (e_codes[iIPMaterialData]) ip_material_data.SetRow(fShapes->CurrIP(), ipmat);
			}
			
			/* mass averaged centroid */
			if (e_codes[iCentroid] || e_codes[iMass])
			{
				/* mass */
				mass += ip_w*density;
			
				/* moment */
				if (e_codes[iCentroid])
				{
					fShapes->IPCoords(ip_centroid);
					centroid.AddScaled(ip_w*density, ip_centroid);
				}
			}
			
			/* kinetic energy/linear momentum */
			if (e_codes[iKineticEnergy] || e_codes[iLinearMomentum])
			{
				/* velocity at integration point */
				fShapes->InterpolateU(fLocVel, ip_velocity);
				
				/* kinetic energy */
				if (e_codes[iKineticEnergy])
					kinetic_energy += 0.5*ip_w*density*dArrayT::Dot(ip_velocity, ip_velocity);
					
				/* linear momentum */
				if (e_codes[iLinearMomentum])
					linear_momentum.AddScaled(ip_w*density, ip_velocity);
			}
		}

		/* copy in the cols */
		int colcount = 0;
		nodal_all.BlockColumnCopyAt(disp       , colcount); colcount += disp.MinorDim();
		nodal_all.BlockColumnCopyAt(coords     , colcount); colcount += coords.MinorDim();
		nodal_all.BlockColumnCopyAt(nodalstress, colcount); colcount += nodalstress.MinorDim();
		nodal_all.BlockColumnCopyAt(princstress, colcount); colcount += princstress.MinorDim();
		nodal_all.BlockColumnCopyAt(energy     , colcount); colcount += energy.MinorDim();
		nodal_all.BlockColumnCopyAt(speed      , colcount); colcount += speed.MinorDim();
		nodal_all.BlockColumnCopyAt(matdat     , colcount);

		/* accumulate - extrapolation done from ip's to corners => X nodes */
		ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_all);
		
		/* element values */
		if (e_codes[iCentroid]) centroid /= mass;
		
		/* store results */
		e_values.SetRow(CurrElementNumber(), element_values);
	}

	/* get nodally averaged values */
	ElementSupport().OutputUsedAverage(n_values);
}

/*********************************************************************
*Private
**********************************************************************/
/* form of tangent matrix */
GlobalT::SystemTypeT PMLT::TangentType(void) const
{
	return GlobalT::kNonSymmetric;
}

LocalArrayT& PMLT::TotalDisp(LocalArrayT& disp)
{
	double* ptot = fTotDisp.Pointer();
	double* pdispa = disp.Pointer();
	double* pdispb = pdispa+fNEESub;

	for (int i = 0; i<fNEESub; i++)
		*ptot++ = (*pdispa++) + (*pdispb++);
			
	return (fTotDisp);
}

LocalArrayT& PMLT::TotalLastDisp(LocalArrayT& disp)
{
	double* ptot = fTotLastDisp.Pointer();
	double* pdispa = disp.Pointer();
	double* pdispb = pdispa+fNEESub;

	for (int i = 0; i<fNEESub; i++)
		*ptot++ = (*pdispa++) + (*pdispb++);
	
	return (fTotLastDisp);
}
LocalArrayT& PMLT::TotalAcc(LocalArrayT& acc)
{
	double* ptot = fTotAcc.Pointer();
	double* pacca = acc.Pointer();
	double* paccb = pacca+fNEESub;

	for (int i = 0; i<fNEESub; i++)
		*ptot++ = (*pacca++) + (*paccb++);
	
	return (fTotAcc);
}
LocalArrayT& PMLT::TotalVel(LocalArrayT& vel)
{
	double* ptot = fTotVel.Pointer();
	double* pvela = vel.Pointer();
	double* pvelb = pvela+fNEESub;

	for (int i = 0; i<fNEESub; i++)
		*ptot++ = (*pvela++) + (*pvelb++);
	
	return (fTotVel);
}
	
void PMLT::Ba(dMatrixT& Ba_matrix, dMatrixT& B_matrix)
{
#if __option(extended_errorcheck)
	if (B_matrix.Rows() != Ba_matrix.Rows() ||
	    B_matrix.Cols() != Ba_matrix.Rows())
	    throw ExceptionT::kSizeMismatch;
	if (NumSD() != 2)
	{
		cout << "\n PML has been implemented only for 2D geometry";
		throw ExceptionT::kBadInputValue;
	}
#endif

	double* pB = B_matrix.Pointer();
	double* pBa = Ba_matrix.Pointer();
	
	for (int i = 0; i< NumElementNodes(); i++)
	{	
		*pBa++ = *pB++;
		*pBa++ = *pB++;
		*pBa++ = 0.0; pB++; 
		
		*pBa++ = *pB++;
		*pBa++ = 0.0; pB++;
		*pBa++ = *pB++;
	}
}

void PMLT::Bb(dMatrixT& Bb_matrix, dMatrixT& B_matrix)
{
#if __option(extended_errorcheck)
	if (B_matrix.Rows() != Bb_matrix.Rows() ||
	    B_matrix.Cols() != Bb_matrix.Rows())
	    throw ExceptionT::kSizeMismatch;
	if (NumSD() != 2)
	{
		cout << "\n PML has been implemented only for 2D geometry";
		throw ExceptionT::kBadInputValue;
	}
#endif

	double* pB = B_matrix.Pointer();
	double* pBb = Bb_matrix.Pointer();
	
	for (int i = 0; i< NumElementNodes(); i++)
	{	
		*pBb++ = 0.0; pB++;
		*pBb++ = *pB++;
		*pBb++ = *pB++; 
		
		*pBb++ = *pB++;
		*pBb++ = *pB++;
		*pBb++ = 0.0; pB++;
	}
}

	
/***********************************************************************
* Protected
***********************************************************************/

/* construct list of materials from the input stream */
void PMLT::ReadMaterialData(ifstreamT& in)
{
ExceptionT::GeneralFail("PMLT::ReadMaterialData", "out of date");
#if 0	
	/* inherited */
	SolidElementT::ReadMaterialData(in);

	/* offset to class needs flags */
	fNeedsOffset = fMaterialNeeds[0].Length();
	
	/* set material needs */
	for (int i = 0; i < fMaterialNeeds.Length(); i++)
	{
		/* needs array */
		ArrayT<bool>& needs = fMaterialNeeds[i];

		/* resize array */
		needs.Resize(needs.Length() + 2, true);

		/* casts are safe since class contructs materials list */
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[i];
		SSSolidMatT* mat = (SSSolidMatT*) pcont_mat;

		/* collect needs */
		needs[fNeedsOffset + kstrain     ] = mat->Need_Strain();
		needs[fNeedsOffset + kstrain_last] = mat->Need_Strain_last();
		
		/* consistency */
		needs[kNeedDisp] = needs[kNeedDisp] || needs[fNeedsOffset + kstrain];
		needs[KNeedLastDisp] = needs[KNeedLastDisp] || needs[fNeedsOffset + kstrain_last];
	}
#endif
}

/* increment current element */
void PMLT::SetGlobalShape(void)
{
	/* inherited */
	SolidElementT::SetGlobalShape();
	/* material information */
	int material_number = CurrentElement().MaterialNumber();
	/* material needs */
	const ArrayT<bool>& needs = fMaterialNeeds[material_number];
	
	/* loop over integration points */
	for (int i = 0; i < NumIP(); i++)
	{
		/* deformation gradient */
		if (needs[fNeedsOffset + kstrain])
		{
			/* displacement gradient */
			fShapes->GradU(TotalDisp(fLocDisp), fGradU, i);

			/* symmetric part */
			 fGradU_List[i]=fGradU;
		}

		/* "last" deformation gradient */
		if (needs[fNeedsOffset + kstrain_last])
		{
			/* displacement gradient */
			fShapes->GradU(TotalLastDisp(fLocLastDisp), fGradU, i);

			/* symmetric part */
			 fGradU_last_List=fGradU;
		}
	}

}

