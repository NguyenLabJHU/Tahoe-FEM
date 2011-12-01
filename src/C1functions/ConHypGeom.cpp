/* $Id: ConHypGeom.cpp,v 1.7 2011-12-01 20:37:57 beichuan Exp $ */
// CONHYPGEOM is the confluent hypergeometric function of the 1st kind
// (a.k.a. Kummer function) Because of the rapid growth/decay of individual
// components, for sufficiently large domain values, an asymptotic approximation
// is utilized (see Abramawitz/Stegun, page 511.)
#include "ConHypGeom.h"
#include <cmath>
#include <iostream>
#include "ExceptionT.h"
#include "dArrayT.h"
#include "Gamma.h"

using namespace Tahoe;

/* constants */
const double EPS = 1.0e-12;
const double MaxAn = 10;
const double PI = 3.141592653589793; 

/*
* Constructor
*/
ConHypGeom::ConHypGeom(double A, double B): fA(A), fB(B) { }

/*
* Destructor
*/
ConHypGeom::~ConHypGeom()
{
}

/*
 * Methods
 */
double ConHypGeom::PochhammerRat(double numer, double denom, int index) const
// Modified Pochhammer function - determine ratio of Pochhammer functions
// specified by the index, n. Input [a,b,n], returns P[a,n]/(P[b,n] n!)
// where P=Pochhammer function.
{
	double value = 1.0;
	
	if (index>0)
	{
		for (int i=1; i<=index; i++)
			value *= ((numer+index-i)/(double(i)*(denom+index-i)));
	}
	else
	{
		cout << "\n Negative index sent to ConHypGeom::PochhammerRat.\n";
		throw ExceptionT::kBadInputValue;
	}
	
	return value;
}

double ConHypGeom::PochhammerProd(double lhs, double rhs, int index) const
// Modified Pochhammer function - determine product of Pochhammer functions
// specified by the index, n. Inout [a,b,n], returns P[a,n]*P[b,n]/(n!)
// where P=Pochhammer function.
{
	double value = 1.0;
	
	if (index>0)
	{
		for (int i=1; i<=index; i++)
			value *= (lhs+index-i)*(rhs+index-i)/double(i);
	}
	else
	{
		cout << "\n Negative index sent to ConHypGeom::PochhammerProd.\n";
		throw ExceptionT::kBadInputValue;
	}
	
	return value;
}

/*
* I/O
*/
void ConHypGeom::Print(ostream& out) const
{
	/* parameters */
	out << " ConHypGeom parameter a. . . . . . . . . . . . . . . . = " << fA << '\n';
	out << " ConHypGeom parameter b. . . . . . . . . . . . . . . . = " << fB << '\n';
}

void ConHypGeom::PrintName(ostream& out) const
{
	out << "ConHypGeom Function of 1st Kind (Confluent Hypergeometric Function)\n";
}

/*
* Returning values
*/
double ConHypGeom::Function(double x) const
{
	double value = 1.0;
	
	if ((fabs(x) < MaxAn) && (fabs(x) > 0.0))
	{
		double sum=1.0, next_term=0.0;
		int k=1;   
	// need better method for determining accuracy - run the risk of infinite loop!
		do {
			next_term = PochhammerRat(fA,fB,k)*pow(x,k);
			sum += next_term;
			k++;
		} while (fabs(next_term)>EPS);
		
		value = sum;
	}
	else if (fabs(x) >= MaxAn)
	// Real portion of the asymptotic approximation of CHG function
	// using 13.5.1 of Abramawitz/Stegun, pg. 508. 
	{
		double term[2];
		int i, maxterm = 8;	// number of terms in expansion
		Gamma g;
		
		term[0] = term[1] = 1.0;	// initialize entries with index=0 value
		for (i=1; i<maxterm; i++)
		{
			term[0] += (PochhammerProd(fA,1+fA-fB,i)*pow(-x,-i));
			term[1] += (PochhammerProd(fB-fA,1-fA,i)*pow(x,-i));
		}
		
		if (x<0)
		{
			term[0] *= (pow(fabs(x),-fA)/g.Function(fB-fA));
			term[1] *= (exp(x)*cos(PI*(fA-fB))*pow(fabs(x),fA-fB)/g.Function(fA));
		}
		else if (x>0)
		{
			term[0] *= (cos(PI*fA)/(pow(x,fA)*g.Function(fB-fA)));
			term[1] *= (exp(x)*pow(x,fA-fB)/g.Function(fA));
		}
		
		value = g.Function(fB)*(term[0]+term[1]);
	}
		
	return value;
}

double ConHypGeom::DFunction(double x) const
{
	cout << "\n***ERROR! DFunction of ConHypGeom incomplete!\n";
	return 0.0*x;
}

double ConHypGeom::DDFunction(double x) const
{
	cout << "\n***ERROR! DDFunction of ConHypGeom incomplete!\n";
	return 0.0*x;
}

/*
* Returning values in groups - derived classes should define
* their own non-virtual function called within this functon
* which maps in to out w/o requiring a virtual function call
* everytime. Default behavior is just to map the virtual functions
* above.
*/
dArrayT& ConHypGeom::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl = in.Pointer();
	double* pU = out.Pointer();
	double sum=0.0;
	Gamma g;
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = *pl++;
		double value = 1.0;
	
		if ((fabs(r) < MaxAn) && (fabs(r) > 0.0))
		{
			double sum=1.0, next_term=0.0;
			int k=1;   
		// need better method for determining accuracy - run the risk of infinite loop!
			do {
				next_term = PochhammerRat(fA,fB,k)*pow(r,k);
				sum += next_term;
				k++;
			} while (fabs(next_term)>EPS);
		
			value = sum;
		}
		else if (fabs(r) >= MaxAn)
		// Real portion of the asymptotic approximation of CHG function
		// using 13.5.1 of Abramawitz/Stegun, pg. 508. 
		{
			double term[2];
			int i, maxterm = 8;	// number of terms in expansion
			Gamma g;
		
			term[0] = term[1] = 1.0;	// initialize entries with index=0 value
			for (i=1; i<maxterm; i++)
			{
				term[0] += (PochhammerProd(fA,1+fA-fB,i)*pow(-r,-i));
				term[1] += (PochhammerProd(fB-fA,1-fA,i)*pow(r,-i));
			}
		
			if (r<0)
			{
				term[0] *= (pow(fabs(r),-fA)/g.Function(fB-fA));
				term[1] *= (exp(r)*cos(PI*(fA-fB))*pow(fabs(r),fA-fB)/g.Function(fA));
			}
			else if (r>0)
			{
				term[0] *= (cos(PI*fA)/(pow(r,fA)*g.Function(fB-fA)));
				term[1] *= (exp(r)*pow(r,fA-fB)/g.Function(fA));
			}	
			value = g.Function(fB)*(term[0]+term[1]);
		}				
		*pU++ = sum;
	}
	return out;
}

dArrayT& ConHypGeom::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl  = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = *pl++;					
		*pdU++ = 0.0;
	}
	return out;
}

dArrayT& ConHypGeom::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = *pl++;				
		*pddU++ = 0.0;
	}
	return out;
}

