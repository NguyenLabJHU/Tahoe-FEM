/* $Id: C1_LineT.cpp,v 1.1 2004-09-02 18:25:08 rdorgan Exp $ */ 
#include "C1_LineT.h"
#include "ifstreamT.h"
#include "StringT.h"
#include "ExceptionT.h"
#include "toolboxConstants.h"

#include "LocalArrayT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"

using namespace Tahoe;

C1_LineT::C1_LineT(GeometryT::CodeT geometry_code, int numip, int numnd_u,  int numdof_u, const LocalArrayT& coords):
	ShapeTools     (),

	fGeometryCode  (geometry_code),

	fCoords        (coords),
	fNumSD         (GeometryT::GeometryToNumSD(fGeometryCode)),   // only 1D: 1
	fNumIP         (numip),       // only 1,2,3,4

	fNumNodes_X    (fCoords.NumberOfNodes()),
	fNumNodes_U    (numnd_u),     // only vertex nodes: 2 (Line)
	fNumDOF_X      (fNumSD),
	fNumDOF_U      (numdof_u),

	fLNaU          (fNumIP, fNumNodes_U*fNumDOF_U),
	fLDNaU         (fNumIP),
	fLDDNaU        (fNumIP),

	fGNaU          (fNumIP, fNumNodes_U*fNumDOF_U),
	fGDNaU         (fNumIP),
	fGDDNaU        (fNumIP),

	fLNaX          (fNumIP, fNumNodes_X),
	fLDNaX         (fNumIP),

	fWeights       (fNumIP),
	fArray1        (fNumSD),
	fArray2        (fNumSD),
	fJacobian      (fNumSD),
	fMatx1         (fNumSD)
{
	const char caller[] = "C1_LineT::Constructor";

	/* check spatial dimensions */
	if (fNumSD != 1)
		ExceptionT::BadInputValue(caller, "nsd!=1");

	/* check element type */
	if (fGeometryCode != GeometryT::kLine)
		ExceptionT::BadInputValue(caller, "GeometryCode!=1");

	/* additional dimension checks */
	if (fNumNodes_U != 2)
		ExceptionT::BadInputValue(caller, "numnode_u != 2 (nsd=1)");
	if (fNumIP != 1 &&
		fNumIP != 2 && 
		fNumIP != 3 && 
		fNumIP != 4)
		ExceptionT::BadInputValue(caller, "numint != 1, 2, 3 or 4 (nsd=1)");
	if (fNumDOF_U != 2)
		ExceptionT::BadInputValue(caller, "numdof_u != 2 (nsd=1)");

	/* addional allocation for shape function derivatives */
	for (int i = 0; i < fNumIP; i++)
	{
		fLDNaU[i].Dimension(fNumSD, fNumNodes_U*fNumDOF_U);
		fGDNaU[i].Dimension(fNumSD, fNumNodes_U*fNumDOF_U);

		fLDDNaU[i].Dimension(fNumSD, fNumNodes_U*fNumDOF_U);
		fGDDNaU[i].Dimension(fNumSD, fNumNodes_U*fNumDOF_U);

		fLDNaX[i].Dimension(fNumSD, fNumNodes_X);
	}
}

/* set all local parameters */
void C1_LineT::Initialize(void)
{
	/* initialize arrays */
	fLNaU = 0.0;
	fLNaX = 0.0;

	for (int i = 0; i < fNumIP; i++) fLDNaU[i] = 0.0;
	for (int i = 0; i < fNumIP; i++) fLDNaX[i] = 0.0;
	for (int i = 0; i < fNumIP; i++) fLDDNaU[i] = 0.0;

	/* set local shape functions and their derivatives */
	SetLocalShape();
}


/* chain rule jacobian of shapefunctions wrt coordinates that
 * are passed in, for all integration points at once
 *
 * fCoords (#nodes x #sd), fLDNaX[#ip x (#sd x #nodes)]
 * fLNaU[#ip x #nodes#ndof] -> fGNaU[#ip x #nodes#ndof]
 * fLDNaU[#ip x (#sd x #nodes#ndof)] -> fGDNaU[#ip x (#sd x #nodes#ndof)]
 * fLDDNaU[#ip x (#sd x #nodes#ndof)] -> fGDDNaU[#ip x (#sd x #nodes#ndof)] */
void C1_LineT::SetDerivatives()
{
	const char caller[] = "C1_LineT::ComputeGDNaU";

	/* loop over all integration points */
	for (int i = 0; i < fNumIP; i++)
    {
		/* calculate the Jacobian matrix */
		Jacobian(fCoords, fLDNaX[i], fJacobian); 
		double det = fJacobian.Det();
     
		/* element check */
		if (det <= 0.0) ExceptionT::GeneralFail(caller, "det(Jac) < ksmall");

		dMatrixT jac = fJacobian;
		dMatrixT& jac_inv = fJacobian.Inverse();
		
		/* calculate the global shape function derivatives */
		double* pLNa = fLNaU(i);
		double* pGNa = fGNaU(i);
			
		double* pLNax = fLDNaU[i](0);
		double* pGNax = fGDNaU[i](0);
			
		double* pLNaxx = fLDDNaU[i](0);
		double* pGNaxx = fGDDNaU[i](0);
			
		double* pjac = jac.Pointer();
		double* pjac_inv = jac_inv.Pointer();

		for (int j = 0; j < fNumNodes_U; j++)
		{
			*pGNa++   = *pLNa++;
			*pGNa++   = *pLNa++   * pjac[0];

			*pGNax++  = *pLNax++  * pjac_inv[0];
			*pGNax++  = *pLNax++;

			*pGNaxx++ = *pLNaxx++ * pjac_inv[0] * pjac_inv[0];
			*pGNaxx++ = *pLNaxx++ * pjac_inv[0];
		}
    }
}

/* PRIVATE MEMBER FUNCTIONS */

/* compute local shape functions and derivatives */
/* fLNaU[#IP x #nodes#dof), fLDNaU[#IP x (#sd x #nodes#dof)], fLDDNaU[#IP x (#sd x #nodes#dof) */
/* fLNaX[#IP x #nodes], fLDNaX[#IP x (#sd x #nodes)] */
void C1_LineT::SetLocalShape()
{
	const char caller[] = "C1_LineT::SetLocalShape";

	/* 1 point */
	double r1[1] = {0.0};
        
	/* 2 point */
	double x2    = sqrt(1.0/3.0);
	double r2[2] = {-x2, x2};

	/* 3 point */
	double    x3 = sqrt(3.0/5.0);
	double r3[3] = {-x3, 0.0, x3};

	/* 4 point */
	double x41 = sqrt((3.0 - 2.0*sqrt(6.0/5.0))/7.0);
	double x42 = sqrt((3.0 + 2.0*sqrt(6.0/5.0))/7.0);
	double r4[4] = {-x42,-x41, x41, x42};
        
	/* integration coordinates and weights */
	double *xa;
	switch (fNumIP)
		{
		case 1: {
			xa = r1;                        
			fWeights[0] = 2.0;
			break;
		}

		case 2: {
			xa = r2;                                                   
			fWeights[0] = 1.0;
			fWeights[1] = 1.0;
			break;
		}

		case 3: {
			xa = r3;
			double a = 5.0/9.0;
			double b = 8.0/9.0;                          
			fWeights[0] = a;
			fWeights[1] = b;
			fWeights[2] = a;
			break;
		}

		case 4: {
			xa = r4;                                                    
			
			double w1 = (18.0 + sqrt(30.0))/36.0;
			double w2 = (18.0 - sqrt(30.0))/36.0;
			fWeights[0] = w2;
			fWeights[1] = w1;
			fWeights[2] = w1;
			fWeights[3] = w2;
			break;
		}

		default:
			ExceptionT::GeneralFail(caller, "Bad numint (nsd=1)");
		}
  
	/* set field shape functions and derivatives at IPs */
	switch (fNumNodes_U)
		{
		case 2: {         
			for (int i = 0; i < fNumIP; i++)	
			{	
				double* na_field   = fLNaU(i);
				double* nax_field  = fLDNaU[i](0);
				double* naxx_field = fLDDNaU[i](0);

				/* Na */
				na_field[0] =  0.25*(xa[i] - 1)*(xa[i] - 1)*(xa[i] + 2);
				na_field[2] = -0.25*(xa[i] + 1)*(xa[i] + 1)*(xa[i] - 2);

				na_field[1] =  0.25*(xa[i] - 1)*(xa[i] - 1)*(xa[i] + 1);
				na_field[3] =  0.25*(xa[i] + 1)*(xa[i] + 1)*(xa[i] - 1);
                
				/* Na,x */
				nax_field[0] =  0.75*(xa[i] - 1)*(xa[i] + 1);
				nax_field[2] = -0.75*(xa[i] + 1)*(xa[i] - 1);

				nax_field[1] =  0.25*(xa[i] - 1)*(3*xa[i] + 1);
				nax_field[3] =  0.25*(xa[i] + 1)*(3*xa[i] - 1);

				/* Na,xx */
				naxx_field[0] =  1.5*xa[i];
				naxx_field[2] = -1.5*xa[i];

				naxx_field[1] =  0.5*(3*xa[i] - 1);
				naxx_field[3] =  0.5*(3*xa[i] + 1);
			}
			break;
		}
		
		default:
			ExceptionT::GeneralFail(caller, "Bad fNumNodes_U");
		}

	/* set geometry shape functions and derivatives at IPs */
	switch (fNumNodes_X)
		{
		case 2: {
			for (int i = 0; i < fNumIP; i++)	
			{	
				double* na_geom   = fLNaX(i);
				double* nax_geom  = fLDNaX[i](0);

				/* Na */
				na_geom[0] = 0.5*(1.0 - xa[i]);
				na_geom[1] = 0.5*(1.0 + xa[i]);
		
				/* Na,x */
				nax_geom[0] =-0.5;
				nax_geom[1] = 0.5;				
			}
			break;
		}
		case 3: {
			for (int i = 0; i < fNumIP; i++)	
			{
				double* na_geom  = fLNaX(i);
				double* nax_geom = fLDNaX[i](0);

		    	/* Na */
		    	na_geom[0] =-xa[i]*0.5*(1.0 - xa[i]);
		    	na_geom[1] = xa[i]*0.5*(1.0 + xa[i]);
		    	na_geom[2] = (1.0 - xa[i])*(1.0 + xa[i]);
		
		        /* Na,x */
		    	nax_geom[0] =-0.5 + xa[i];
		    	nax_geom[1] = 0.5 + xa[i];
		    	nax_geom[2] =-2.0*xa[i];
		    }
			break;
		}
		
		default:
			ExceptionT::GeneralFail(caller, "Bad fNumNodes_X");
		}
}

/* nodal (#nodes x #sd), LDNaX (#sd x #nodes), fJac (#sd x #sd) */
/* note: #sd=LDNaX.MajorDim(), #nodes=LDNaX.MinorDim() */
void C1_LineT::Jacobian(const LocalArrayT& nodal, const dArray2DT& LDNaX,
							   dMatrixT& jac)
{
	// jacobian
	jac = 0.;
	for (int i = 0; i < fNumNodes_X; i++)
	{
		LDNaX.ColumnCopy(i, fArray2);
		for (int j = 0; j < fNumSD; j++) fArray1[j] = nodal(i,j);
		fMatx1.Outer(fArray1, fArray2);
		jac += fMatx1;
	}
}
