/* $Id: ShapeFunctionT.cpp,v 1.9.2.1 2002-09-21 08:59:07 paklein Exp $ */
/* created: paklein (06/26/1996) */

#include "ShapeFunctionT.h"
#include "ParentDomainT.h"
#include "LocalArrayT.h"
#include "dSymMatrixT.h"

using namespace Tahoe;

/* constructor */
ShapeFunctionT::ShapeFunctionT(GeometryT::CodeT geometry_code, int numIP,
	const LocalArrayT& coords):
	DomainIntegrationT(geometry_code, numIP, coords.NumberOfNodes()),
	fCoords(coords),
	fGrad_x_temp(NULL)
{
	/* consistency */
	if (GeometryT::GeometryToNumSD(geometry_code) != fCoords.MinorDim())
	{
		cout << "\n ShapeFunctionT::ShapeFunctionT: geometry code "
		     << geometry_code << " does not match\n"
		     <<   "     the number of spatial dimensions of the coordinates "
		     << fCoords.MinorDim() << endl;
		throw eGeneralFail;
	}

	/* configure workspace */
	Construct();
}

ShapeFunctionT::ShapeFunctionT(const ShapeFunctionT& link, const LocalArrayT& coords):
	DomainIntegrationT(link),
	fCoords(coords),
	fGrad_x_temp(NULL)
{
	/* configure workspace */
	Construct();
}	

/* compute local shape functions and derivatives */ 	
void ShapeFunctionT::SetDerivatives(void)
{
	fDomain->ComputeDNa(fCoords, fDNaX, fDet);
}

/* field gradients at specific parent domain coordinates. */
void ShapeFunctionT::GradU(const LocalArrayT& nodal, dMatrixT& grad_U, 
	const dArrayT& coord) const
{
#pragma unused(nodal)
#pragma unused(grad_U)
#pragma unused(coord)
	cout << "\n ShapeFunctionT::GradU: not implemented yet" << endl;
	throw eGeneralFail;
}

/************************ for the current integration point *********************/
void ShapeFunctionT::InterpolateU(const LocalArrayT& nodal,
	dArrayT& u) const
{
#if __option(extended_errorcheck)
	if (nodal.MinorDim() != u.Length() ||
	    nodal.NumberOfNodes() != pNaU->MinorDim()) throw eSizeMismatch;
#endif

	int num_u = nodal.MinorDim();
	for (int i = 0; i < num_u; i++)
		u[i] = pNaU->DotRow(fCurrIP, nodal(i));
}

void ShapeFunctionT::InterpolateU(const LocalArrayT& nodal,
	dArrayT& u, int ip) const
{
#if __option(extended_errorcheck)
	if (nodal.MinorDim() != u.Length() ||
	    nodal.NumberOfNodes() != pNaU->MinorDim()) throw eSizeMismatch;
#endif

	int num_u = nodal.MinorDim();
	for (int i = 0; i < num_u; i++)
		u[i] = pNaU->DotRow(ip, nodal(i));
}

#if 0
/* strain displacement matrix B */
void ShapeFunctionT::B(const dArray2DT& DNa, dMatrixT& B_matrix) const
{
#if __option(extended_errorcheck)
	if (B_matrix.Rows() != dSymMatrixT::NumValues(DNa.MajorDim()) ||
	    B_matrix.Cols() != DNa.Length())
	    throw eSizeMismatch;
#endif

	int numnodes = DNa.MinorDim();
	double*   pB = B_matrix.Pointer();

	/* standard strain-displacement operator */
	if (fB_option == kStandardB)
	{
	        /* 1D */
                if (DNa.MajorDim() == 1)
		{
			double* pNax = DNa(0);

			for (int i = 0; i < numnodes; i++)
			{
			        /* currently assuming that DNa gets 1D shape functions
                                   correctly from LineT somehow...*/
				*pB++ = *pNax++;
			}
		}
		/* 2D */
		else if (DNa.MajorDim() == 2)
		{
			double* pNax = DNa(0);
			double* pNay = DNa(1);

			for (int i = 0; i < numnodes; i++)
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
			double* pNax = DNa(0);
			double* pNay = DNa(1);
			double* pNaz = DNa(2);
			
			for (int i = 0; i < numnodes; i++)
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
	/* B-bar: mean dilatation - using values calculated during the last
	 * call to SetMeanDilatation */
	else
	{
	        /* 1D */
                if (DNa.MajorDim() == 1)
		{
		  	cout << "\n ShapeFunctionT::B: not implemented yet for 1D B-bar" << endl;
	                throw eGeneralFail;
		}
	        /* 2D */
		else if (DNa.MajorDim() == 2)
		{
			double* pNax = DNa(0);
			double* pNay = DNa(1);
			
			double* pBmx = fB_workspace(0);
			double* pBmy = fB_workspace(1);
			
			for (int i = 0; i < numnodes; i++)
			{
				double factx = ((*pBmx++) - (*pNax))/3.0;
				double facty = ((*pBmy++) - (*pNay))/3.0;
			
				/* Hughes (4.5.11-16) */
				*pB++ = *pNax + factx;
				*pB++ = factx;
				*pB++ = *pNay;
	
				*pB++ = facty;
				*pB++ = *pNay + facty;
				*pB++ = *pNax;
				
				pNax++; pNay++;
			}
		}
		/* 3D */
		else		
		{
			double* pNax = DNa(0);
			double* pNay = DNa(1);
			double* pNaz = DNa(2);

			double* pBmx = fB_workspace(0);
			double* pBmy = fB_workspace(1);
			double* pBmz = fB_workspace(2);
			
			for (int i = 0; i < numnodes; i++)
			{
				double factx = ((*pBmx++) - (*pNax))/3.0;
				double facty = ((*pBmy++) - (*pNay))/3.0;
				double factz = ((*pBmz++) - (*pNaz))/3.0;

				/* Hughes (4.5.11-16) */
				*pB++ = *pNax + factx;
				*pB++ = factx;
				*pB++ = factx;
				*pB++ = 0.0;
				*pB++ = *pNaz;
				*pB++ = *pNay;
	
				*pB++ = facty;
				*pB++ = *pNay + facty;
				*pB++ = facty;
				*pB++ = *pNaz;
				*pB++ = 0.0;
				*pB++ = *pNax;
	
				*pB++ = factz;
				*pB++ = factz;
				*pB++ = *pNaz + factz;
				*pB++ = *pNay;
				*pB++ = *pNax;
				*pB++ = 0.0;
				
				pNax++; pNay++; pNaz++;
			}
		}
	}
}
#endif

/* shape function gradients matrix (Hughes,4.90) */
void ShapeFunctionT::GradNa(const dArray2DT& DNa, dMatrixT& grad_Na) const
{
#if __option(extended_errorcheck)
	if (DNa.MajorDim() != grad_Na.Rows() ||
	    DNa.MinorDim() != grad_Na.Cols())
	    throw eSizeMismatch;
#endif

	int numsd    = DNa.MajorDim();
	int numnodes = DNa.MinorDim();
	double* p    = grad_Na.Pointer();

	if (numsd == 2)
	{
		double* pNax = DNa(0);
		double* pNay = DNa(1);
		
		for (int i = 0; i < numnodes; i++)
		{
			*p++ = *pNax++;
			*p++ = *pNay++;
		}
	}
	else if (numsd == 3)
	{
		double* pNax = DNa(0);
		double* pNay = DNa(1);
		double* pNaz = DNa(2);
		
		for (int i = 0; i < numnodes; i++)
		{
			*p++ = *pNax++;
			*p++ = *pNay++;
			*p++ = *pNaz++;
		}
	}
	else /* 1D case defaults to this */
	{
		for (int i = 0; i < numsd; i++)	
			for (int a = 0; a < numnodes; a++)	
				grad_Na(i,a) = DNa(i,a);
	}
}

#if 0
/* strain displacement matrix B */
void ShapeFunctionT::B_q(const dArray2DT& DNa, dMatrixT& B_matrix) const
{
#if __option(extended_errorcheck)
	if (DNa.MajorDim() != (*pDNaU)[fCurrIP].MajorDim() ||
	    DNa.MinorDim() != (*pDNaU)[fCurrIP].MinorDim())
	    throw eGeneralFail;
#endif

	int numnodes = DNa.MinorDim();
double*   pB = B_matrix.Pointer();

	/* 2D */
	if (DNa.MajorDim() == 2)
	{
		double* pNax = DNa(0);
		double* pNay = DNa(1);

		for (int i = 0; i < numnodes; i++)
		{
			*pB++ = *pNax++;
			*pB++ = *pNay++;
		}
	}
	/* 3D */
	else		
	{
		double* pNax = DNa(0);
		double* pNay = DNa(1);
		double* pNaz = DNa(2);
		
		for (int i = 0; i < numnodes; i++)
		{
			*pB++ = *pNax++;
			*pB++ = *pNay++;
			*pB++ = *pNaz++;
		}
	}
}
#endif

/********************************************************************************/

/* print the shape function values to the output stream */
void ShapeFunctionT::Print(ostream& out) const
{
/* inherited */
	DomainIntegrationT::Print(out);

	out << "\n Domain shape function derivatives:\n";
	for (int i = 0; i < fDNaX.Length(); i++)
		fDNaX[i].WriteNumbered(out);

	out << "\n Field shape function derivatives:\n";
	if (fDNaX.Pointer() == pDNaU->Pointer())
	    out << " isoparametric \n";
	else	
	    for (int i = 0; i < pDNaU->Length(); i++)
			(*pDNaU)[i].WriteNumbered(out);
}

/***********************************************************************
* Protected
***********************************************************************/

void ShapeFunctionT::DoTransformDerivatives(const dMatrixT& changeofvar, 
	const dArray2DT& original, dArray2DT& transformed)
{
	int  numnodes = original.MinorDim();

	/* allocate memory */
	transformed.Allocate(original.MajorDim(),numnodes);

	/* apply chain rule derivative */
	for (int i = 0; i < numnodes; i++)
	{
		/* fetch values */
		original.ColumnCopy(i,fv1);

		/* transform */
		changeofvar.MultTx(fv1,fv2);
		
		/* write back */	
		transformed.SetColumn(i,fv2);
	}
}

/***********************************************************************
* Private
***********************************************************************/

/* configure work space arrays - initializes shape function to be
* isoparametric */
void ShapeFunctionT::Construct(void)
{
	/* check local array type (ambiguous for linear geometry) */
	if (fCoords.Type() != LocalArrayT::kInitCoords &&
	    fCoords.Type() != LocalArrayT::kCurrCoords) throw eGeneralFail;

	/* dimensions */
	int numXnodes = fCoords.NumberOfNodes();
	int numUnodes = numXnodes; // assume isoparametric
	int numsd     = fCoords.MinorDim();

	/* parent domain jacobian */
	fDet.Allocate(fNumIP),

	/* memory for the derivatives */
	fDNaX.Allocate(fNumIP);
	for (int i = 0; i < fNumIP; i++)
		fDNaX[i].Allocate(numsd, numXnodes);		

	/* initialize to isoparametric */
	pNaU  = &(fDomain->Na());
	pDNaU = &fDNaX;
	
	/* work space */
	fv1.Allocate(numsd);
	fv2.Allocate(numsd);
}
