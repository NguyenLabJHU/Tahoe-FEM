/* $Id: C0_HexahedronT.cpp,v 1.1 2004-09-02 18:25:08 rdorgan Exp $ */ 
#include "C0_HexahedronT.h"
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

C0_HexahedronT::C0_HexahedronT(GeometryT::CodeT geometry_code, int numip, int numnd_u,  int numdof_u, const LocalArrayT& coords):
	ShapeTools     (),

	fGeometryCode  (geometry_code),

	fCoords        (coords),
	fNumSD         (GeometryT::GeometryToNumSD(fGeometryCode)),   // only 3D: 3
	fNumIP         (numip), // only 1,8,27

	fNumNodes_X    (fCoords.NumberOfNodes()), // only 8
	fNumNodes_U    (numnd_u),                 // only 8
	fNumDOF_X      (fNumSD),                  // only 3
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
	const char caller[] = "C0_HexahedronT::Constructor";

	/* check spatial dimensions */
	if (fNumSD != 3)
		ExceptionT::BadInputValue(caller, "nsd!=3");

	/* check element type */
	if (fGeometryCode != GeometryT::kHexahedron)
		ExceptionT::BadInputValue(caller, "GeometryCode!=kHexahedron");

	/* additional dimension checks */
	if (fNumNodes_U != 8)
		ExceptionT::BadInputValue(caller, "numnode_u != 8 (nsd=3)", fNumNodes_U);
	if (fNumIP != 1 &&
		fNumIP != 8 && 
		fNumIP != 27)
		ExceptionT::BadInputValue(caller, "numint != 1, 8, or 27 (nsd=3)");
	if (fNumDOF_U != 1)
		ExceptionT::BadInputValue(caller, "numdof_u != 1 (nsd=3)");

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
void C0_HexahedronT::Initialize(void)
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
 * fLDDNaU[#ip x (#sd x #nodes#ndof x 6)] -> fGDDNaU[#ip x (#sd x #nodes#ndof)] */
void C0_HexahedronT::SetDerivatives()
{
	const char caller[] = "C0_HexahedronT::ComputeGDNaU";

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
		double* pLNaz = fLDNaU[i](2);

		double* pGNax = fGDNaU[i](0);
		double* pGNay = fGDNaU[i](1);
		double* pGNaz = fGDNaU[i](2);
			
		double* pLNaxx = fLDDNaU[i](0);
		double* pLNayy = fLDDNaU[i](1);
		double* pLNazz = fLDDNaU[i](2);

		double* pLNayz = fLDDNaU[i](3);
		double* pLNaxz = fLDDNaU[i](4);
		double* pLNaxy = fLDDNaU[i](5);

		double* pGNaxx = fGDDNaU[i](0);
		double* pGNayy = fGDDNaU[i](1);
		double* pGNazz = fGDDNaU[i](2);
			
		double* pjac = jac.Pointer();
		double* pjac_inv = jac_inv.Pointer();

		for (int j = 0; j < fNumNodes_U; j++)
		{
			*pGNa++   = (*pLNa++);
			
			*pGNax++  = (*pLNax  )*pjac_inv[0]+(*pLNay   )* pjac_inv[1]+(*pLNaz  )*pjac_inv[2];
			*pGNay++  = (*pLNax  )*pjac_inv[3]+(*pLNay   )* pjac_inv[4]+(*pLNaz  )*pjac_inv[5];
			*pGNaz++  = (*pLNax++)*pjac_inv[6]+(*pLNay++ )* pjac_inv[7]+(*pLNaz++)*pjac_inv[8];

			*pGNaxx    = (*pLNaxx)*pjac_inv[0]*pjac_inv[0]+(*pLNaxy)*pjac_inv[0]*pjac_inv[1]+(*pLNaxz)*pjac_inv[0]*pjac_inv[2];
			*pGNaxx   += (*pLNaxy)*pjac_inv[1]*pjac_inv[0]+(*pLNayy)*pjac_inv[1]*pjac_inv[1]+(*pLNayz)*pjac_inv[1]*pjac_inv[2];
			*pGNaxx++ += (*pLNaxz)*pjac_inv[2]*pjac_inv[0]+(*pLNayz)*pjac_inv[2]*pjac_inv[1]+(*pLNazz)*pjac_inv[2]*pjac_inv[2];

			*pGNayy    = (*pLNaxx)*pjac_inv[3]*pjac_inv[3]+(*pLNaxy)*pjac_inv[3]*pjac_inv[4]+(*pLNaxz)*pjac_inv[3]*pjac_inv[5];
			*pGNayy   += (*pLNaxy)*pjac_inv[4]*pjac_inv[3]+(*pLNayy)*pjac_inv[4]*pjac_inv[4]+(*pLNayz)*pjac_inv[4]*pjac_inv[5];
			*pGNayy++ += (*pLNaxz)*pjac_inv[5]*pjac_inv[3]+(*pLNayz)*pjac_inv[5]*pjac_inv[4]+(*pLNazz)*pjac_inv[5]*pjac_inv[5];

			*pGNazz    = (*pLNaxx++)*pjac_inv[6]*pjac_inv[6]+(*pLNaxy  )*pjac_inv[6]*pjac_inv[7]+(*pLNaxz  )*pjac_inv[6]*pjac_inv[8];
			*pGNazz   += (*pLNaxy++)*pjac_inv[7]*pjac_inv[6]+(*pLNayy++)*pjac_inv[7]*pjac_inv[7]+(*pLNayz  )*pjac_inv[7]*pjac_inv[8];
			*pGNazz++ += (*pLNaxz++)*pjac_inv[8]*pjac_inv[6]+(*pLNayz++)*pjac_inv[8]*pjac_inv[7]+(*pLNazz++)*pjac_inv[8]*pjac_inv[8];
		}
    }
}

/* PRIVATE MEMBER FUNCTIONS */

/* compute local shape functions and derivatives */
/* fLNaU[#IP x #nodes#dof), fLDNaU[#IP x (#sd x #nodes#dof)], fLDDNaU[#IP x (#sd x #nodes#dof x 6) */
/* fLNaX[#IP x #nodes], fLDNaX[#IP x (#sd x #nodes)] */
void C0_HexahedronT::SetLocalShape()
{
	const char caller[] = "C0_HexahedronT::SetLocalShape3D";


	/* integration point coordinates */
    double  ra[] = {-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0};
    double  sa[] = {-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0};
    double  ta[] = {-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0};
    double  xa_27[27], ya_27[27], za_27[27];
    double  *xa, *ya, *za;
    double  g;
    
	/* integration weights */
    switch (fNumIP)
		{
		case 1:	
			g = 0.0;
			xa = ra;
			ya = sa;
			za = ta;
			fWeights[0] = 8.0;
			break;
	
		case 8:
			g = 1.0/sqrt3;
			xa = ra;
			ya = sa;
			za = ta;
			fWeights = 1.0;
			break;
	
		case 27:
		{
			/* coordinates */
			double b1 = sqrt(3.0/5.0);
			double b_1D[3] = {-b1, 0.0, b1}; 
	
			/* weights */
			double w1 = 5.0/9.0;
			double w2 = 8.0/9.0;
			double w_1D[3] = {w1, w2, w1};
			int x_i = 0;
			int y_i = 0;
			int z_i = 0;
			for (int i = 0; i < 27; i++)
			{
				xa_27[i]   = b_1D[x_i];
				ya_27[i]   = b_1D[y_i];
				za_27[i]   = b_1D[z_i];
				fWeights[i] = w_1D[x_i]*w_1D[y_i]*w_1D[z_i];

				if (++x_i == 3)
				{
					x_i = 0;
					if (++y_i == 3)
					{
						y_i = 0;
						z_i++;
					}
				}
			}						

			xa = xa_27;
			ya = ya_27;
			za = za_27;
			g  = 1.0;		
			break;
		}
	
		default:
			ExceptionT::GeneralFail(caller, "Bad numint (nsd=3)");
		}
    
	/* set field shape functions and derivatives at IPs */
	/* vertex nodes */
    for (int in = 0; in < fNumIP; in++)	
	{
		double* na_field  = fLNaU(in);
		
		double* nax_field = fLDNaU[in](0);
		double* nay_field = fLDNaU[in](1);
		double* naz_field = fLDNaU[in](2);

		double* naxx_field = fLDDNaU[in](0);
		double* nayy_field = fLDDNaU[in](1);
		double* nazz_field = fLDDNaU[in](2);

		double* nayz_field = fLDDNaU[in](3);
		double* naxz_field = fLDDNaU[in](4);
		double* naxy_field = fLDDNaU[in](5);

		double r = g*xa[in];
		double s = g*ya[in];
		double t = g*za[in];
	
		for (int lnd = 0; lnd < fNumNodes_U; lnd++)
		{
			double tempr1 = 1.0 + ra[lnd]*r;
			double temps1 = 1.0 + sa[lnd]*s;
			double tempt1 = 1.0 + ta[lnd]*t;
	    
			*na_field++  = 0.125*tempr1*temps1*tempt1;

			*nax_field++ = 0.125*ra[lnd]*temps1*tempt1;
			*nay_field++ = 0.125*tempr1*sa[lnd]*tempt1;
			*naz_field++ = 0.125*tempr1*temps1*ta[lnd];

			*naxx_field++ = 0.0;
			*nayy_field++ = 0.0;
			*nazz_field++ = 0.0;

			*naxy_field++ = 0.125*ra[lnd]*sa[lnd]*tempt1;
			*nayz_field++ = 0.125*tempr1*sa[lnd]*ta[lnd];
			*naxz_field++ = 0.125*ra[lnd]*temps1*ta[lnd];
		}
	}

	/* set geometry shape functions and derivatives at IPs */
	/* vertex nodes */
    for (int in = 0; in < fNumIP; in++)	
	{
		double* na_geom   = fLNaX(in);
		
		double* nax_geom = fLDNaX[in](0);
		double* nay_geom = fLDNaX[in](1);
		double* naz_geom = fLDNaX[in](2);

		double r = g*xa[in];
		double s = g*ya[in];
		double t = g*za[in];
	
		for (int lnd = 0; lnd < fNumNodes_X; lnd++)
		{
			double tempr1 = 1.0 + ra[lnd]*r;
			double temps1 = 1.0 + sa[lnd]*s;
			double tempt1 = 1.0 + ta[lnd]*t;
	    
			*na_geom++  = 0.125*tempr1*temps1*tempt1;

			*nax_geom++ = 0.125*ra[lnd]*temps1*tempt1;
			*nay_geom++ = 0.125*tempr1*sa[lnd]*tempt1;
			*naz_geom++ = 0.125*tempr1*temps1*ta[lnd];
		}
	}
}

/* nodal (#nodes x #sd), LDNaX (#sd x #nodes), fJac (#sd x #sd) */
/* note: #sd=LDNaX.MajorDim(), #nodes=LDNaX.MinorDim() */
void C0_HexahedronT::Jacobian(const LocalArrayT& nodal, const dArray2DT& LDNaX,
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
