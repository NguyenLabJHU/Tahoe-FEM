/* $Id: C0_QuadT.cpp,v 1.1 2004-09-02 18:25:08 rdorgan Exp $ */ 
#include "C0_QuadT.h"
#include "ifstreamT.h"
#include "StringT.h"
#include "ExceptionT.h"
#include "toolboxConstants.h"

#include "LocalArrayT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"

using namespace Tahoe;

/* parameters */
const double sqrt3 = sqrt(3.0);

C0_QuadT::C0_QuadT(GeometryT::CodeT geometry_code, int numip, int numnd_u,  int numdof_u, const LocalArrayT& coords):
	ShapeTools     (),

	fGeometryCode  (geometry_code),

	fCoords        (coords),
	fNumSD         (GeometryT::GeometryToNumSD(fGeometryCode)),   // only 1D: 1
	fNumIP         (numip), // only 1,4,9

	fNumNodes_X    (fCoords.NumberOfNodes()), // only 4
	fNumNodes_U    (numnd_u),                 // only 4
	fNumDOF_X      (fNumSD),                  // only 2
	fNumDOF_U      (numdof_u),                // only 1

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
	const char caller[] = "C0_QuadT::Constructor";

	/* check spatial dimensions */
	if (fNumSD != 2)
		ExceptionT::BadInputValue(caller, "nsd!=2");

	/* check element type */
	if (fGeometryCode != GeometryT::kQuadrilateral)
		ExceptionT::BadInputValue(caller, "GeometryCode!=kQuadrilateral");

	/* additional dimension checks */
	if (fNumNodes_U != 4)
		ExceptionT::BadInputValue(caller, "numnode_u != 4 (nsd=2)", fNumNodes_U);
	if (fNumIP != 1 &&
		fNumIP != 4 && 
		fNumIP != 9)
		ExceptionT::BadInputValue(caller, "numint != 1, 4, or 9 (nsd=2)");
	if (fNumDOF_U != 1)
		ExceptionT::BadInputValue(caller, "numdof_u != 1 (nsd=2)");

	/* addional allocation for shape function derivatives */
	for (int i = 0; i < fNumIP; i++)
	{
		fLDNaU[i].Dimension(fNumSD, fNumNodes_U*fNumDOF_U);
		fGDNaU[i].Dimension(fNumSD, fNumNodes_U*fNumDOF_U);

		fLDDNaU[i].Dimension(dSymMatrixT::NumValues(fNumSD), fNumNodes_U*fNumDOF_U);
		fGDDNaU[i].Dimension(fNumSD, fNumNodes_U*fNumDOF_U);

		fLDNaX[i].Dimension(fNumSD, fNumNodes_X);
	}
}

/* set all local parameters */
void C0_QuadT::Initialize(void)
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
 * fLDDNaU[#ip x (#dSymMatrixT::NumValues(fNumSD) x #nodes#ndof)] -> fGDDNaU[#ip x (#sd x #nodes#ndof)] */
void C0_QuadT::SetDerivatives()
{
	const char caller[] = "C0_QuadT::ComputeGDNaU";

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
		double* pLNay = fLDNaU[i](1);

		double* pGNax = fGDNaU[i](0);
		double* pGNay = fGDNaU[i](1);
			
		double* pLNaxx = fLDDNaU[i](0);
		double* pLNayy = fLDDNaU[i](1);
		double* pLNaxy = fLDDNaU[i](2);

		double* pGNaxx = fGDDNaU[i](0);
		double* pGNayy = fGDDNaU[i](1);
			
		double* pjac = jac.Pointer();
		double* pjac_inv = jac_inv.Pointer();

		for (int j = 0; j < fNumNodes_U; j++)
		{
			*pGNa++   = (*pLNa++);

			*pGNax++  = (*pLNax  )*pjac_inv[0]+(*pLNay  )* pjac_inv[1];
			*pGNay++  = (*pLNax++)*pjac_inv[2]+(*pLNay++)* pjac_inv[3];

			*pGNaxx    = (*pLNaxx)*pjac_inv[0]*pjac_inv[0]+(*pLNaxy)*pjac_inv[0]*pjac_inv[1];
			*pGNaxx++ += (*pLNaxy)*pjac_inv[1]*pjac_inv[0]+(*pLNayy)*pjac_inv[1]*pjac_inv[1];

			*pGNayy    = (*pLNaxx++)*pjac_inv[2]*pjac_inv[2]+(*pLNaxy  )*pjac_inv[2]*pjac_inv[3];
			*pGNayy++ += (*pLNaxy++)*pjac_inv[3]*pjac_inv[2]+(*pLNayy++)*pjac_inv[3]*pjac_inv[3];
		}
    }
}

/* PRIVATE MEMBER FUNCTIONS */

/* compute local shape functions and derivatives */
/* fLNaU[#IP x #nodes#dof), fLDNaU[#IP x (#sd x #nodes#dof)], fLDDNaU[#IP x (#dSymMatrixT::NumValues(fNumSD) x #nodes#dof) */
/* fLNaX[#IP x #nodes], fLDNaX[#IP x (#sd x #nodes)] */
void C0_QuadT::SetLocalShape()
{
	const char caller[] = "C0_QuadT::SetLocalShape";

	/* integration point coordinates */
    double  ra[9] = {-1.0, 1.0, 1.0,-1.0, 0.0, 1.0, 0.0,-1.0, 0.0};
    double  sa[9] = {-1.0,-1.0, 1.0, 1.0,-1.0, 0.0, 1.0, 0.0, 0.0};
    double *xa, *ya;
    double g;
  
	/* integration weights */
    switch (fNumIP)
		{
		case 1:	
			xa = ra;
			ya = sa;
			g = 0.0;
			fWeights[0] = 4.0;
			break;
	
		case 4:
			xa = ra;			  
			ya = sa;			  
			g = 1.0/sqrt3;
			fWeights[0] = 1.0;
			fWeights[1] = 1.0;
			fWeights[2] = 1.0;
			fWeights[3] = 1.0;
			break;
	
		case 9:
        {
			xa = ra;			  
			ya = sa;			  
			g = sqrt(3.0/5.0);
			double a = 25.0/81.0;
			double b = 40.0/81.0;
			double c = 64.0/81.0;
			fWeights[0] = a;
			fWeights[1] = a;
			fWeights[2] = a;
			fWeights[3] = a;
			fWeights[4] = b;
			fWeights[5] = b;
			fWeights[6] = b;
			fWeights[7] = b;
			fWeights[8] = c;
			break;
        }
	
		default:
			ExceptionT::GeneralFail(caller, "Bad numint (nsd=2)");
		}
    
	/* set field shape functions and derivatives at IPs */
	/* vertex nodes */
    for (int in = 0; in < fNumIP; in++)	
	{
		double* na_field  = fLNaU(in);
		
		double* nax_field = fLDNaU[in](0);
		double* nay_field = fLDNaU[in](1);

		double* naxx_field = fLDDNaU[in](0);
		double* nayy_field = fLDDNaU[in](1);
		double* naxy_field = fLDDNaU[in](2);

		double r = g*xa[in];
		double s = g*ya[in];
	
		for (int lnd = 0; lnd < fNumNodes_U; lnd++)
		{
			double tempr1 = 1.0 + ra[lnd]*r;
			double temps1 = 1.0 + sa[lnd]*s;
	    
			*na_field++  = 0.25*tempr1*temps1;

			*nax_field++ = 0.25*ra[lnd]*temps1;
			*nay_field++ = 0.25*tempr1*sa[lnd];

			*naxx_field++ = 0.0;
			*nayy_field++ = 0.0;

			*naxy_field++ = 0.25*ra[lnd]*sa[lnd];
		}
	}

	/* set geometry shape functions and derivatives at IPs */
	/* vertex nodes */
    for (int in = 0; in < fNumIP; in++)	
	{
		double* na_geom   = fLNaX(in);
		
		double* nax_geom = fLDNaX[in](0);
		double* nay_geom = fLDNaX[in](1);

		double r = g*xa[in];
		double s = g*ya[in];
	
		for (int lnd = 0; lnd < fNumNodes_X; lnd++)
		{
			double tempr1 = 1.0 + ra[lnd]*r;
			double temps1 = 1.0 + sa[lnd]*s;
	    
			*na_geom++  = 0.25*tempr1*temps1;

			*nax_geom++ = 0.25*ra[lnd]*temps1;
			*nay_geom++ = 0.25*tempr1*sa[lnd];
		}
	}
}

/* nodal (#nodes x #sd), LDNaX (#sd x #nodes), fJac (#sd x #sd) */
/* note: #sd=LDNaX.MajorDim(), #nodes=LDNaX.MinorDim() */
void C0_QuadT::Jacobian(const LocalArrayT& nodal, const dArray2DT& LDNaX,
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
