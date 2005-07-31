/* $Id: ShapeFunctionT.cpp,v 1.1.1.1 2001-01-29 08:20:31 paklein Exp $ */
/* created: paklein (06/26/1996)                                          */

#include "ShapeFunctionT.h"
#include "ParentDomainT.h"
#include "LocalArrayT.h"

/* constructor */
ShapeFunctionT::ShapeFunctionT(GeometryT::CodeT geometry_code, int numIP,
	const LocalArrayT& coords, int B_option):
DomainIntegrationT(geometry_code, numIP, coords.NumberOfNodes()),
	fCoords(coords),
	fB_option(B_option),
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

	/* check option */
	if (fB_option != kStandardB &&
	    fB_option != kMeanDilBbar) throw eBadInputValue;

	/* configure workspace */
	Construct();
}

ShapeFunctionT::ShapeFunctionT(const ShapeFunctionT& link, const LocalArrayT& coords):
	DomainIntegrationT(link),
	fCoords(coords),
	fB_option(link.fB_option),
	fGrad_x_temp(NULL)
{
	/* check */
	if (fB_option != kStandardB)
	{
		cout << "\n ShapeFunctionT::ShapeFunctionT: WARNING: non-standard B-option not\n";
		cout << "     expected for linked shape functions" << endl;
	}

	/* configure workspace */
	Construct();
}	

/* compute local shape functions and derivatives */ 	
void ShapeFunctionT::SetDerivatives(void)
{
	fDomain->ComputeDNa(fCoords, fDNaX, fDet);
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

/* strain displacement matrix B */
void ShapeFunctionT::B(const dArray2DT& DNa, dMatrixT& B_matrix) const
{
#if __option(extended_errorcheck)
	if (DNa.MajorDim() != (*pDNaU)[fCurrIP].MajorDim() ||
	    DNa.MinorDim() != (*pDNaU)[fCurrIP].MinorDim())
	    throw eGeneralFail;
#endif

	int numnodes = DNa.MinorDim();
double*   pB = B_matrix.Pointer();

	/* standard strain-displacement operator */
	if (fB_option == kStandardB)
	{
		/* 2D */
		if (DNa.MajorDim() == 2)
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
	/* B-bar: mean dilatation */
	else
	{
		/* non-const this */
		ShapeFunctionT* non_const_this = (ShapeFunctionT*) this;
	
		/* compute mean dilatation, Hughes (4.5.23) */
		non_const_this->SetMeanDilatation();
	
		/* 2D */
		if (DNa.MajorDim() == 2)
		{

			//fB_workspace.Allocate(numsd, numUnodes);

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

/* shape function gradients matrix (Hughes,4.90) */
void ShapeFunctionT::GradNa(const dArray2DT& DNa, dMatrixT& grad_Na) const
{
#if __option(extended_errorcheck)
	if (DNa.MajorDim() != (*pDNaU)[fCurrIP].MajorDim() ||
	    DNa.MinorDim() != (*pDNaU)[fCurrIP].MinorDim())
	    throw eGeneralFail;
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
	else
	{
		for (int i = 0; i < numsd; i++)	
			for (int a = 0; a < numnodes; a++)	
				grad_Na(i,a) = DNa(i,a);
	}
}

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

/* extrapolate integration point values to the nodes
*    IPvalues[numvals] : values from a single integration point
*    nodalvalues[fNumNodes x numvals] : extrapolated values */
void ShapeFunctionT::Extrapolate(const dArrayT& IPvalues,
	dArray2DT& nodalvalues) const
{
	fDomain->NodalValues(IPvalues, nodalvalues, CurrIP());
}	

/* convert shape function derivatives by applying a chain rule
* transformation:
*
*      d Na / d x_i = (d Na / d X_J) (d X_J/d x_i)
*/
void ShapeFunctionT::SetChainRule(const dMatrixT& changeofvar,
	dArray2DT& derivatives)
{
	dArray2DT& DNa = (*pDNaU)[fCurrIP];
	int   numnodes = DNa.MinorDim();

	/* allocate memory */
	derivatives.Allocate(DNa.MajorDim(),numnodes);

	/* apply chain rule derivative */
	for (int i = 0; i < numnodes; i++)
	{
		/* fetch values */
		DNa.ColumnCopy(i,fv1);

		/* transform */
		changeofvar.MultTx(fv1,fv2);
		
		/* write back */	
		derivatives.SetColumn(i,fv2);
	}
}

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

	/* B-option workspace */
	if (fB_option == kMeanDilBbar) fB_workspace.Allocate(numsd, numUnodes);

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

/* compute mean dilatation, Hughes (4.5.23) */
void ShapeFunctionT::SetMeanDilatation(void)
{
	/* volume */
	const double* w = IPWeights();
	double* det = fDet.Pointer();
	double  vol = 0.0;
	for (int i = 0; i < fNumIP; i++)
		vol += (*w++)*(*det++);

	/* initialize */
	fB_workspace = 0.0;			

	/* integrate */
	w   = IPWeights();
	det = fDet.Pointer();
	for (int l = 0; l < fNumIP; l++)
		fB_workspace.AddScaled((*w++)*(*det++)/vol, (*pDNaU)[l]);
}
