/* $Id: LocalizeT.cpp,v 1.7 2004-07-15 08:28:31 paklein Exp $ */
/* created: paklein (09/11/1997) */

#include "LocalizeT.h"
#include "ElementSupportT.h"
//#include "CommunicatorT.h"
//#include "CommManagerT.h"
#include "LocalArrayT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "ifstreamT.h"

using namespace Tahoe;

/* constructor */
LocalizeT::LocalizeT(const ElementSupportT& support):
  fTime(support.Time()),
  fPi(acos(-1.0)){
ExceptionT::GeneralFail("LocalizeT::LocalizeT", "out of date");
#if 0	  
  /*Dimension member data*/
  int nsd = support.NumSD();
  fStress.Dimension(nsd);
  fModulus.Dimension(dSymMatrixT::NumValues(nsd));
  fNormal.Dimension(nsd);
  fNormal = 0.0;
  fElemCenter.Dimension(nsd);
  fElemCenter = 0.0;

  ifstreamT& in = support.Input();
  ostream& out = support.Output();

  in >> fCheck;
  if (fCheck == 1){
    out << "\nCheck Localization"<<endl;

    int numblocks;
    in >> numblocks;
    fBlockList.Dimension(numblocks);
    for (int i = 0; i < numblocks; i++)
    {
      int index;
      in >> index;
      fBlockList[i] = index-1;
    }
    out << "\nIn Blocks "<<fBlockList<<endl;

    const StringT& input_file = in.filename();
    flocalize_file.Root(input_file);
    flocalize_file.Append(".loc");
 
    fout.open(flocalize_file);
    fout << setw(8) << "Time"
	 << setw(16) << "Coord_X"
	 << setw(16) << "Coord_Y"
	 << setw(16) << "Normal_X"
	 << setw(16) << "Normal_Y"
	   << endl;
    fout.close();
  }
  else if (fCheck == 0) {
    out <<"\nNo Localization Check"<<endl;
  }
  else {
    cout << "\nLocalizeT::LocalizeT invalid input for check localization flag.";
    throw ExceptionT::kBadInputValue; 
  }
#endif  
}


void LocalizeT::WriteLocalize(const iArrayT& flags, const dArray2DT& elem_centers, const dArray2DT& normals)
{
  int nel = flags.Length();
  
  fout.open_append(flocalize_file);
  for (int i = 0; i < nel; i++)
  {
    if (flags[i] == 1)
      fout << setw(8) << fTime
	   << setw(16) << elem_centers(i,0)
	   << setw(16) << elem_centers(i,1)
	   << setw(16) << normals(i,0)
	   << setw(16) << normals(i,1)    
	   << endl;
  }
  fout.close();
}

void LocalizeT::WriteLocalize(void)
{
  fout.open_append(flocalize_file);
  fout << setw(8) <<fTime
       << setw(16) << fElemCenter[0]
       << setw(16) << fElemCenter[1]
       << setw(16) << fNormal[0]
       << setw(16) << fNormal[1]      
       << endl;
  fout.close();
}

int LocalizeT::CheckLocalizeFS(const dSymMatrixT& stress, const dMatrixT& modulus, const LocalArrayT& initcoords)
{
  fStress = stress; 
  fModulus = modulus;

  initcoords.Average(fElemCenter);
  if (fStress.Rows() == 2)
    return DetCheck2D(fNormal);
  else {
    cout << "LocalizeT::CheckLocalizationFS: finite strain localization check not implemented in 3D \n";
    throw ExceptionT::kGeneralFail;
  }
}

int LocalizeT::CheckLocalizeSS(const dSymMatrixT& stress, const dMatrixT& modulus, const LocalArrayT& initcoords)
{
  fStress = stress; 
  fModulus = modulus;

  initcoords.Average(fElemCenter);
  if (fStress.Rows() == 2){
    double theta = 0.0;
    int check = 0;
    SPINLOC_localize(fModulus.Pointer(), &theta, &check);
    if (check == 1){
      fNormal[0] = cos(theta);
      fNormal[1] = sin(theta);
    }
    return check;
  }
  else {
    cout << "DetCheck2D::IsLocalized: Small Strain localization check not implemented in 3D \n";
    throw ExceptionT::kGeneralFail;
  }
  
}

/**********************************************************************
* Private
**********************************************************************/

/* 2D determinant check function */
int LocalizeT::DetCheck2D(dArrayT& normal)
{
	/* set det function */
	ComputeCoefficients();

	/* quick exit */
	if ( fA0 > 0 && fA0 > (fabs(fA2) + fabs(fA4)) ) return 0;

	/* local minima */
	double min2  = ((fA2 > 0.0) ? 3.0*fPi/4.0 : fPi/4.0) - fphi2;

	double dat[4] = {-fPi/8.0, 3.0*fPi/8.0, 7.0*fPi/8.0, 11.0*fPi/8.0};

	dArrayT mins4(4,dat);
	if (fA4 > 0) mins4 -= fphi4;
	else mins4 -= (fphi4 + fPi/4.0);

	/* closest local minima in sin2t and sin4t*/
	int    dex;
	double min = fPi;
	for (int i = 0; i < 4; i++)
	{	
		double diff = fabs(mins4[i] - min2);

		if (diff < min)
		{	
			dex = i;
			min = diff;
		}
	}
	
	/* search starting value */
	double angle = 0.5*(min2 + mins4[dex]);
	double maxstep = fPi/18.0; //10 degrees max per step
	
	/* Newton search */
	double res  = ddet(angle);
	double res0 = res;
	int   count = 0;
	while ( fabs(res/res0) > 1.0e-10 && ++count < 15)
	{
		double dangle = res/dddet(angle);
		dangle = (fabs(dangle) > maxstep) ? maxstep*((dangle < 0) ? -1:1) : dangle;
	
		angle -= dangle;
		res  = ddet(angle);
	}

	/* check minima */
	if (dddet(angle) < 0.0)
	{
		cout << "\n LocalizeT::IsLocalized: ERROR:\n";
		cout << "f =  " << fA0 << " + " << fA2 << "*sin(2.0*(";
		cout << fphi2 << "+t)) + " << fA4 << "*sin(4.0*(";
		cout << fphi4 << "+t))" << '\n';
		cout << "root found at " << angle << " rad" << endl;
		cout << "starting at   " << 0.5*(min2 + mins4[dex]) << " rad" << endl;


		/* search starting value */
		double angle = 0.5*(min2 + mins4[dex]);
		
		/* Newton search */
		double res  = ddet(angle);
		double res0 = res;
		int   count = 0;
		while( fabs(res/res0) > 1.0e-10 && ++count < 10)
		{
			double dangle = res/dddet(angle);

			cout << count << '\t' << res << '\t' << dangle;
		
			angle -= dangle;
			res  = ddet(angle);
		}
			throw ExceptionT::kGeneralFail;
	}
	
	if (det(angle) > 0.0)
		return 0;
	else
	{
		normal.Dimension(2);

		/* angle in [0, fPi] */
		if (angle > fPi) angle -= fPi;
		
		/* compute normal */
		normal[0] = cos(angle);
		normal[1] = sin(angle);

		return 1;
	}
}

/* compute coefficients of det(theta) function */
void LocalizeT::ComputeCoefficients(void)
{
	/* moduli components */

	double c11 = fModulus(0,0);
	double c22 = fModulus(1,1);
	double c33 = fModulus(2,2);
	double c23 = fModulus(1,2);
	double c13 = fModulus(0,2);
	double c12 = fModulus(0,1);
	
	/* stress components */
	double s11 = fStress[0];
	double s22 = fStress[1];
	double s12 = fStress[2];

	/* intermediate values */
	double s2t = (-(c12*c13) + c13*c22 + c11*c23 - c12*c23 + c13*s11 +
			c23*s11 + c11*s12 + c22*s12 + 2*c33*s12 + 2*s11*s12 + c13*s22 +
			c23*s22 + 2*s12*s22)/2;
	double s4t = (-(c12*c13) - c13*c22 + c11*c23 + c12*c23 + c13*s11 +
			c23*s11 + c11*s12 - c22*s12 + 2*s11*s12 - c13*s22 - c23*s22 -
			2*s12*s22)/4;
			
	double c2t = (-c13*c13 + c23*c23 + c11*c33 - c22*c33 + c11*s11 + c33*s11
			+ s11*s11 - c22*s22 - c33*s22 - s22*s22)/2;
      	double c4t = (c12*c12 - c13*c13 - c11*c22 - 2*c13*c23 - c23*c23 +
			c11*c33 + 2*c12*c33 + c22*c33 + c11*s11 - c22*s11 + s11*s11 -
			4*c13*s12 - 4*c23*s12 - 4*s12*s12 - c11*s22 + c22*s22 - 2*s11*s22
			+ s22*s22)/8;

	/* phase shifts */
       	if (fabs(c2t) < 1.0e-12 && fabs(s2t) < 1.0e-12)
	  fphi2 = fPi/4.0;
	else fphi2 = atan2(c2t,s2t)/2.0;
	
	if (fabs(c4t) < 1.0e-12 && fabs(s4t) < 1.0e-12)
	  fphi4 = fPi/8.0;
	  else fphi4 = atan2(c4t,s4t)/4.0;

	if (fabs(fphi4) == 0.0)
	  cout << "\nDivide by zero: c4t = "<<c4t<<" s4t = "<<s4t
	       << "\nfabs c4t = "<<fabs(c4t)<<" fabs s4t = " <<fabs(s4t);
	if (fabs(fphi2) == 0.0)
	  cout << "\nDivide by zero: c2t = "<<c2t<<" s2t = "<<s2t 
	       << "\nfabs c2t = "<<fabs(c2t)<<" fabs s2t = " <<fabs(s2t);

	/* amplitudes */
	fA0 = (-c12*c12 - 3*c13*c13 + c11*c22 + 2*c13*c23 - 3*c23*c23 +
			3*c11*c33 - 2*c12*c33 + 3*c22*c33 + 3*c11*s11 + c22*s11 +
			4*c33*s11 + 3*s11*s11 + 4*c13*s12 + 4*c23*s12 + 4*s12*s12 +
			c11*s22 + 3*c22*s22 + 4*c33*s22 + 2*s11*s22 + 3*s22*s22)/8;

	fA2 = c2t/sin(2.0*fphi2);
	fA4 = c4t/sin(4.0*fphi4);

	/*	cout << "\nstress: "<<fStress;
	cout << "\nmodulus: " <<fModulus;

	cout << "\ns2t: "<<s2t
	     << "  s4t: "<<s4t
	     << "  c2t: "<<c2t
	     << "  c4t: "<<c4t<<endl;

	cout << "\nA0: "<<fA0;
	cout << "\nA2: "<<fA2;
	cout << "\nA4: "<<fA4;
	cout << "\nphi2: "<<fphi2;
	cout << "\nphi4: "<<fphi4;*/
}


int LocalizeT::SPINLOC_localize(double *c__, double *thetan, int *loccheck)
{
    /* Initialized data */
    double zero = 0.;
    double one = 1.;
    double two = 2.;
    double three = 3.;
    double four = 4.;
    double tol = .01;

    /* System generated locals */
    int i__1;
    double d__1, d__2, d__3, d__4;

    /* Local variables */
    double capa, capb, half, fmin, temp, xmin, temp2, temp3, a, b, 
	    f;
    int i__, n;
    double p, q, r__, x[3], theta, third, a0, a1, a2, a3, a4, qq, rad;

    /* Parameter adjustments */
    c__ -= 4;

    /* Function Body */
    half = one / two;
    third = one / three;
    rad = four * atan(one) / 180.;


    //  cout << "c__=\n";
    //  cout << c__[4] << ' ' <<  c__[7] << ' ' << c__[10] << '\n';
    //  cout << c__[5] << ' ' <<  c__[8] << ' ' << c__[11] << '\n';
    //  cout << c__[6] << ' ' <<  c__[9] << ' ' << c__[12] << '\n';


    a0 = c__[4] * c__[12] - c__[10] * c__[6];
    a1 = c__[4] * (c__[9] + c__[11]) - c__[10] * c__[5] - c__[6] * c__[7];
    a2 = c__[4] * c__[8] + c__[10] * c__[9] + c__[6] * c__[11] - c__[7] * (
	    c__[12] + c__[5]) - c__[12] * c__[5];
    a3 = c__[8] * (c__[10] + c__[6]) - c__[11] * c__[7] - c__[9] * c__[5];
    a4 = c__[12] * c__[8] - c__[9] * c__[11];

    p = three / four * (a3 / a4);
    q = a2 / a4 * (one / two);
    r__ = a1 / a4 * (one / four);

	/* Computing 2nd power */
    d__1 = p;
    a = (three * q - d__1 * d__1) * third;
/* Computing 3rd power */
    d__1 = p, d__2 = d__1;
    b = (two * (d__2 * (d__1 * d__1)) - p * 9. * q + r__ * 27.) / 27.;

/* Computing 2nd power */
    d__1 = b;
/* Computing 3rd power */
    d__2 = a, d__3 = d__2;
    qq = d__1 * d__1 / four + d__3 * (d__2 * d__2) / 27.;
    if (fabs(qq) < 1e-8) {
	qq = zero;
    }

    temp = p * third;

    if (qq > zero || qq == 0.f) {
	temp2 = one;
	temp3 = -half * b + sqrt(qq);
	if (temp3 < zero) {
	    temp2 = -one;
	}
	d__1 = fabs(temp3);
	capa = temp2 * pow(d__1, third);
	temp2 = one;
	temp3 = -half * b - sqrt(qq);
	if (temp3 < zero) {
	    temp2 = -one;
	}
	d__1 = fabs(temp3);
	capb = temp2 * pow(d__1, third);
	x[0] = capa + capb - temp;
	x[1] = -(capa + capb) * half - temp;
	x[2] = x[1];
    } else {
	if (a < zero) {
/* Computing 3rd power */
	    d__2 = a * third, d__3 = d__2;
	    theta = acos(-half * b / sqrt((d__1 = -(d__3 * (d__2 * d__2)), 
		    fabs(d__1))));
	    temp2 = two * sqrt((d__1 = -a * third, fabs(d__1)));
	    x[0] = temp2 * cos(theta * third) - temp;
	    x[1] = -temp2 * cos(theta * third + rad * 60.) - temp;
	    x[2] = -temp2 * cos(theta * third - rad * 60.) - temp;
	} else {
		cout << "\n LocalizeT::SPINLOC_localize: a is positive when it should be negative" << endl;
	}
    }

    fmin = 1e50;
    n = 3;
    if (fabs(qq) < 1e-8) {
	n = 2;
    }
    if (qq > zero) {
	n = 1;
    }
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 4th power */
	d__1 = x[i__ - 1], d__1 *= d__1;
/* Computing 3rd power */
	d__2 = x[i__ - 1], d__3 = d__2;
/* Computing 2nd power */
	d__4 = x[i__ - 1];
	f = a4 * (d__1 * d__1) + a3 * (d__3 * (d__2 * d__2)) + a2 * (d__4 * 
		d__4) + a1 * x[i__ - 1] + a0;
	if (f <= fmin) {
	    fmin = f;
	    xmin = x[i__ - 1];
	}
/* L5: */
    }

/* .. output */

    *thetan = atan(xmin);

/* 	if(fmin.lt.tol) then */
    if (fmin / c__[4] < tol) {
/* localized */
	*loccheck = 1;
    } else {
/* not localized */
	*loccheck = 0;
    }
    return 0;
}

