/* $Id: ShapeFunctionT.cpp,v 1.15 2004-08-10 23:13:08 raregue Exp $ */
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
		throw ExceptionT::kGeneralFail;
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
void ShapeFunctionT::GradU(const LocalArrayT& nodal, dMatrixT& grad_U, const dArrayT& coord, 
	dArrayT& Na, dArray2DT& DNa) const
{
	/* dimensions */
	int nnd = fCoords.NumberOfNodes();
	int nsd = fCoords.MinorDim();
	Na.Dimension(nnd);
	DNa.Dimension(nsd, nnd);

	/* compute the shape function derivatives */
	EvaluateShapeFunctions(coord, Na, DNa);

	/* compute the gradient */
	fDomain->Jacobian(nodal, DNa, grad_U);
}

/************************ for the current integration point *********************/
void ShapeFunctionT::InterpolateU(const LocalArrayT& nodal,
	dArrayT& u) const
{
#if __option(extended_errorcheck)
	if (nodal.MinorDim() != u.Length() ||
	    nodal.NumberOfNodes() != pNaU->MinorDim()) throw ExceptionT::kSizeMismatch;
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
	    nodal.NumberOfNodes() != pNaU->MinorDim()) throw ExceptionT::kSizeMismatch;
#endif

	int num_u = nodal.MinorDim();
	for (int i = 0; i < num_u; i++)
		u[i] = pNaU->DotRow(ip, nodal(i));
}

/* shape function gradients matrix (Hughes,4.90) */
void ShapeFunctionT::GradNa(const dArray2DT& DNa, dMatrixT& grad_Na) const
{
#if __option(extended_errorcheck)
	if (DNa.MajorDim() != grad_Na.Rows() ||
	    DNa.MinorDim() != grad_Na.Cols())
	    throw ExceptionT::kSizeMismatch;
#endif

	int numsd    = DNa.MajorDim();
	int numnodes = DNa.MinorDim();
	double* p    = grad_Na.Pointer();

	if (numsd == 2)
	{
		const double* pNax = DNa(0);
		const double* pNay = DNa(1);
		
		for (int i = 0; i < numnodes; i++)
		{
			*p++ = *pNax++;
			*p++ = *pNay++;
		}
	}
	else if (numsd == 3)
	{
		const double* pNax = DNa(0);
		const double* pNay = DNa(1);
		const double* pNaz = DNa(2);
		
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


/* laplacian of field at current ip */
/*
void ShapeFunctionT::LaplaceU(const LocalArrayT& field, dArrayT& laplacian) const
{
	
}
*/

/* laplacian of strain at current ip */
/*
void ShapeFunctionT::LaplaceStrain(const dSymMatrixT& strain, dSymMatrixT& laplacian) const
{
	
}
*/
	

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
	transformed.Dimension(original.MajorDim(),numnodes);

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
	    fCoords.Type() != LocalArrayT::kCurrCoords) throw ExceptionT::kGeneralFail;

	/* dimensions */
	int numXnodes = fCoords.NumberOfNodes();
	int numUnodes = numXnodes; // assume isoparametric
	int numsd     = fCoords.MinorDim();

	/* parent domain jacobian */
	fDet.Dimension(fNumIP),

	/* memory for the derivatives */
	fDNaX.Dimension(fNumIP);
	for (int i = 0; i < fNumIP; i++)
		fDNaX[i].Dimension(numsd, numXnodes);		

	/* initialize to isoparametric */
	pNaU  = &(fDomain->Na());
	pDNaU = &fDNaX;
	
	/* work space */
	fv1.Dimension(numsd);
	fv2.Dimension(numsd);
}
