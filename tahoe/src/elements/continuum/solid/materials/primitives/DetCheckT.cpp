/* $Id: DetCheckT.cpp,v 1.7 2001-08-17 23:32:55 cfoster Exp $ */
/* created: paklein (09/11/1997) */

#include "DetCheckT.h"
#include <math.h>
#include "ExceptionCodes.h"
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dMatrixEXT.h"
#include "dArrayT.h"

/* constants */
const double Pi = acos(-1.0);

/* constructor */
DetCheckT::DetCheckT(const dSymMatrixT& s_jl, const dMatrixT& c_ijkl):
	fs_jl(s_jl),
	fc_ijkl(c_ijkl)
{

}

/* inline functions */

/* determinant function and derivatives */
inline double DetCheckT::det(double t) const
{
	return A0 + A2*sin(2.0*(phi2+t)) + A4*sin(4.0*(phi4+t));
}

inline double DetCheckT::ddet(double t) const
{
	return 2.0*A2*cos(2.0*(phi2+t)) + 4.0*A4*cos(4.0*(phi4+t));
}

inline double DetCheckT::dddet(double t) const
{
	return -4.0*A2*sin(2.0*(phi2+t)) - 16.0*A4*sin(4.0*(phi4+t));
}

/*
* Returns 1 if acoustic tensor isn't positive definite,
* and returns the normal to the surface of localization.
* Returns 0, otherwise.
*/
int DetCheckT::IsLocalized(dArrayT& normal)
{
	if (fs_jl.Rows() == 2)
		return DetCheck2D(normal);
	else
	  {
	//TEMP - not implemented
	  cout << "Finite Strain Localization Check not implemented in 3D \n";
	  return 5;
	  }
}


/* check ellipticity of tangent modulus for small strain formulation 
* 2D version is closed form solution for plane strain taken from 
* R.A.Regueiro's SPINLOC.
* 3d is a numerical search algorithm after Ortiz, et. al. (1987) */

int DetCheckT::IsLocalized_SS(dArrayT& normal)
{
	if (fs_jl.Rows() == 2)
	{
		/* call SPINLOC routine */
		double theta = 0.0;
		int check = 0; 
		SPINLOC_localize(fc_ijkl.Pointer(), &theta, &check);
		if (check == 0)
		{
			normal[0] = cos(theta);
			normal[1] = sin(theta);		
		}
		return check;
	}
	else
          // not ready yet. remove to run problem to debug
	  //return DetCheck3D_SS(normal);
	 return 7;    
}


/**********************************************************************
* Private
**********************************************************************/

/* 2D determinant check function */
int DetCheckT::DetCheck2D(dArrayT& normal)
{
	/* set det function */
	ComputeCoefficients();

	/* quick exit */
	if ( A0 > 0 && A0 > (fabs(A2) + fabs(A4)) ) return 0;

	/* local minima */
	double min2  = ((A2 > 0.0) ? 3.0*Pi/4.0 : Pi/4.0) - phi2;

	double dat[4] = {-Pi/8.0, 3.0*Pi/8.0, 7.0*Pi/8.0, 11.0*Pi/8.0};

	dArrayT mins4(4,dat);
	if (A4 > 0) mins4 -= phi4;
	else mins4 -= (phi4 + Pi/4.0);

	/* closest local minima in sin2t and sin4t*/
	int    dex;
	double min = Pi;
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
	double maxstep = Pi/18.0; //10 degrees max per step
	
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
		cout << "\n DetCheckT::IsLocalized: ERROR:\n";
		cout << "f =  " << A0 << " + " << A2 << "*sin(2.0*(";
		cout << phi2 << "+t)) + " << A4 << "*sin(4.0*(";
		cout << phi4 << "+t))" << '\n';
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
			throw eGeneralFail;
	}
	
	if (det(angle) > 0.0)
		return 0;
	else
	{
		normal.Allocate(2);

		/* angle in [0, Pi] */
		if (angle > Pi) angle -= Pi;
		
		/* compute normal */
		normal[0] = cos(angle);
		normal[1] = sin(angle);
		
		return 1;
	}
}

/* 3D determinant check function */
/* assumes small strain formulation */
int DetCheckT::DetCheck3D_SS(dArrayT& normal)
{

  cout << "In DetCheck3D_SS \n";

  double theta, phi; //horizontal plane angle and polar angle for normal
  int i,j,k,l,m,n, control;        // counters and control variable
  double tol=.0000000001; 
  double leastmin=tol; 
  double guess; 
  int newtoncounter; //makes sure Newton iteration doesn't take too long
  double resid;

  double detA [72] [19]; //determinant of acoustic tensor at each increment
  int localmin [72] [19]; // 1 for local minimum, 0 if not
  dMatrixT A(3,3), Ainverse(3,3); //acoustic tensor
  dMatrixEXT J(3); // det(A)*C*Ainverse
  dArrayT prevnormal(3), finalnormal (3), eigs(3);


  double C [3] [3] [3] [3]; // rank 4 tangent modulus

 // convert fc_ijkl from rank 2 to rank 4
           for (k=0;k<3;k++){
	     for (l=0;l<3;l++){ 
            C [k] [k] [l] [l] = fc_ijkl (k,l);
	      }
	   }
	    
           for (k=0;k<3;k++){
             C [k] [k] [1] [2] = fc_ijkl (k,3);
             C [k] [k] [2] [1] = fc_ijkl (k,3);
             C [k] [k] [0] [2] = fc_ijkl (k,4);
             C [k] [k] [2] [0] = fc_ijkl (k,4);
             C [k] [k] [0] [1] = fc_ijkl (k,5);
             C [k] [k] [1] [0] = fc_ijkl (k,5);
           }

           for (k=0;k<3;k++){
             C [1] [2] [k] [k] = fc_ijkl (3,k);
             C [2] [1] [k] [k] = fc_ijkl (3,k);
             C [0] [2] [k] [k] = fc_ijkl (4,k);
             C [2] [0] [k] [k] = fc_ijkl (4,k);
             C [0] [1] [k] [k] = fc_ijkl (5,k);
             C [1] [0] [k] [k] = fc_ijkl (5,k);
           }

          for (k=0;k<3;k++){
             C [1] [2] [1] [2] = fc_ijkl (3,3);
             C [1] [2] [2] [1] = fc_ijkl (3,3);
             C [2] [1] [1] [2] = fc_ijkl (3,3);
             C [2] [1] [2] [1] = fc_ijkl (3,3);
             
             C [0] [1] [0] [1] = fc_ijkl (5,5);
             C [0] [1] [1] [0] = fc_ijkl (5,5);
             C [1] [0] [0] [1] = fc_ijkl (5,5);
             C [1] [0] [1] [0] = fc_ijkl (5,5);
             
             C [0] [2] [0] [2] = fc_ijkl (4,4);
             C [0] [2] [2] [0] = fc_ijkl (4,4);
             C [2] [0] [0] [2] = fc_ijkl (4,4);
             C [2] [0] [2] [0] = fc_ijkl (4,4);
             
             C [1] [2] [0] [2] = fc_ijkl (3,4);
             C [1] [2] [2] [0] = fc_ijkl (3,4);
             C [2] [1] [0] [2] = fc_ijkl (3,4);
             C [2] [1] [2] [0] = fc_ijkl (3,4);
           
             C [0] [2] [1] [2] = fc_ijkl (4,3);
             C [0] [2] [2] [1] = fc_ijkl (4,3);
             C [2] [0] [1] [2] = fc_ijkl (4,3);
             C [2] [0] [2] [1] = fc_ijkl (4,3);

             C [0] [2] [0] [1] = fc_ijkl (4,5);
             C [0] [2] [1] [0] = fc_ijkl (4,5);
             C [2] [0] [0] [1] = fc_ijkl (4,5);
             C [2] [0] [1] [0] = fc_ijkl (4,5);
             
             C [0] [2] [0] [1] = fc_ijkl (5,4);
             C [0] [2] [1] [0] = fc_ijkl (5,4);
             C [2] [0] [0] [1] = fc_ijkl (5,4);
             C [2] [0] [1] [0] = fc_ijkl (5,4);

             C [1] [2] [0] [1] = fc_ijkl (3,5);
             C [1] [2] [1] [0] = fc_ijkl (3,5);
             C [2] [1] [0] [1] = fc_ijkl (3,5);
             C [2] [1] [1] [0] = fc_ijkl (3,5);
             
             C [0] [2] [2] [1] = fc_ijkl (5,3);
             C [0] [2] [1] [2] = fc_ijkl (5,3);
             C [2] [0] [2] [1] = fc_ijkl (5,3);
             C [2] [0] [1] [2] = fc_ijkl (5,3);

 
           }


 /* initial sweep to determine approximate local minima */

  for (i=0;i<72;i++){
    for (j=0;j<19; j++){

      // cout << "i= " << i << "; j= " << j << '\n';

      theta=Pi/36*i;
      phi=Pi/36*j;

      normal[0]=cos(theta)*cos(phi);
      normal[1]=sin(theta)*cos(phi);
      normal[2]=sin(phi);

	     // initialize acoustic tensor A
   
             for (m=0;m<3;m++){
	       for (n=0;n<3;n++){
		 A (m,n) =0;
	       }
	     }
	    


	     // A=normal*C*normal
              
             for (m=0;m<3;m++){
	       for (n=0;n<3;n++){
		 for (k=0;k<3;k++){
		   for (l=0;l<3;l++){
                     A (m,n)+=normal[k]*C [k] [m] [n] [l]*normal[l]; 
		   }
		 }
	       }
	     }

             detA [i] [j]= A.Det();
    
	     /* if and case statement find approximations to local minima *
              * by checking if it's less than neighbors on 4 sides. *
              * Wraps around in theta and phi */
if ((i==0) && (j==0))
 control =1;
else 
  if ((i==0) && (j==18))
    control = 3;
  else
     if (i==0)
       control = 2;   
     else 
        if ((i==71) && (j==0))
           control = 7;
        else
          if ((i==71) && (j==18))
           control = 9;
          else
             if (i==71)
                control = 8;
             else
               if (j==0)
                 control = 4;
               else
                  if (j==17)
                     control = 6;
                  else
                     control = 5;

 switch(control){

    case 1:
          localmin [i] [j] = 1;
          break;

    case 2:
          if (detA [i] [j] < detA [i] [j-1])
            {
                localmin [i] [j] = 1;
                localmin [i] [j-1] =0;
            }
          else 
                localmin [i] [j] = 0;
          break;

    case 3:
	  if ((detA [i] [j] < detA [i] [j]) && (detA [i] [j] < detA [i] [0]))
	    {
                localmin [i] [j] = 1;
                localmin [i] [j-1] = 0;
                localmin [i] [0] = 0;
	    }
          else 
                localmin [i] [j] =0;
          break;

     case 4:
          if (detA [i] [j] < detA [i-1] [j])
            {
                localmin [i] [j] = 1;
                localmin [i-1] [j] = 0;
            }
          else 
                localmin [i] [j] = 0;
          break;

	
    case 5:
          if ((detA [i] [j] < detA [i-1] [j]) && (detA [i] [j] < detA [i] [j-1]))
            {
                localmin [i] [j] = 1;
                localmin [i-1] [j] =0;
                localmin [i] [j-1] =0;
            }
          else 
                localmin [i] [j] = 0;
          break;
  

    case 6:
        	  if ((detA [i] [j] < detA [i-1] [j]) && (detA [i] [j] < detA [i] [0]) && (detA [i] [j] < detA [i] [j-1]))
	    {
                localmin [i] [j] = 1;
                localmin [i-1] [j] = 0;
                localmin [i] [0] = 0;
                localmin [i] [j-1] =0;
	    }
          else
                localmin [i] [j] = 0;
          break;

        

    case 7:
          if ((detA [i] [j] < detA [i-1] [j]) && (detA [i] [j] < detA [0] [j]))
            {
                localmin [i] [j] = 1;
                localmin [i-1] [j] = 0;
                localmin [0] [j] = 0;
            }
          else 
                localmin [i] [j] = 0;
          break;


     case 8:
          if ((detA [i] [j] < detA [i-1] [j]) && (detA [i] [j] < detA [i] [j-1]) && (detA [i] [j] < detA [0] [j]))
            {
                localmin [i] [j] = 1;
                localmin [i-1] [j] =0;
                localmin [i] [j-1] =0;
                localmin [0] [j] =0;
            }
          else 
                localmin [i] [j] = 0;
          break;

     case 9: 
	  if ((detA [i] [j] < detA [i-1] [j]) && (detA [i] [j] < detA [i] [0]) && (detA [i] [j] < detA [i] [j-1]) && (detA [i] [j] < detA [0] [j]))
	    {
                localmin [i] [j] = 1;
                localmin [i-1] [j] = 0;
                localmin [i] [0] = 0;
                localmin [i] [j-1] =0;
                localmin [0] [j] =0;
	    }
          else
                localmin [i] [j] = 0;
          break; 
 }
    }   
  }
     /* Newton iteration to refine minima*/

 for (i=0;i<72;i++){
    for (j=0;j<19; j++){

       if (localmin [i] [j] ==1)
	 {
            theta=Pi/36*i;
            phi=Pi/36*j;

             normal[0]=cos(theta)*cos(phi);
             normal[1]=sin(theta)*cos(phi);
             normal[2]=sin(phi);

             prevnormal[0]=1;
             prevnormal[1]=1;
             prevnormal[2]=1;


             newtoncounter=0;
	    
	     resid=1.0;

	     // while (sqrt((prevnormal[0]-normal[0])*(prevnormal[0]-normal[0])+(prevnormal[1]-normal[1])*(prevnormal[1]-normal[1])+(prevnormal[2]-normal[2])*(prevnormal[2]-normal[2]))>tol){  
	     while (resid > tol){


	       cout << " normal = \n";
               cout << normal << '\n';

	       cout << " prevnormal = \n";
               cout << prevnormal << '\n';

	       cout << " norm = \n";
               cout << resid << '\n';
	       // cout << sqrt((prevnormal[0]-normal[0])*(prevnormal[0]-normal[0])+(prevnormal[1]-normal[1])*(prevnormal[1]-normal[1])+(prevnormal[2]-normal[2])*(prevnormal[2]-normal[2])) << ' \n';
             


               newtoncounter++;
               if (newtoncounter > 20)
		 {
		 cout << "Newton refinement did not converge after 20 iterations-Localization check failed \n"; 
                 return 7;
                 }


	       for (k=0;k<3;k++){
		 prevnormal[k]=normal[k];
	       }
 // initialize acoustic tensor A
              
             for (m=0;m<3;m++){
	       for (n=0;n<3;n++){
		 A (m,n) =0;
	       }
	     }
	    


	     // A=normal*C*normal
              
             for (m=0;m<3;m++){
	       for (n=0;n<3;n++){
		 for (k=0;k<3;k++){
		   for (l=0;l<3;l++){
                     A (m,n)+=normal[k]*C [k] [m] [n] [l]*normal[l]; 
		   }
		 }
	       }
	     }

             detA [i] [j]= A.Det();  
 	     Ainverse.Inverse(A);

   // initialize J
   
             for (m=0;m<3;m++){
	       for (n=0;n<3;n++){
		 J (m,n) =0;
	       }
	     }

	     //form Jmn=det(A)*Cmkjn*(A^-1)jk

 for (m=0;m<3;m++){
	       for (n=0;n<3;n++){
		 for (k=0;k<3;k++){
		   for (l=0;l<3;l++){
                     J (m,n)+=detA [i] [j]*C [m] [k] [j] [n]*Ainverse (j,k); 
		   }
		 }
	       }
	     }
 // find least eigenvector of J


 cout << "check1 \n";

eigs[0]=1;
eigs[1]=0;
eigs[2]=0;


//J.Diagonalize(eigs);

//if (J(0,0)<J(1,1))
//  guess = J(0,0);
//else
//  guess = J(1,1);

//if (J(2,2)<guess)
//  guess = J(2,2);

guess=0.0;

J.Eigenvector(guess, normal);

 cout << "check2 \n";

resid= (prevnormal[0]-normal[0])*(prevnormal[0]-normal[0])+(prevnormal[1]-normal[1])*(prevnormal[1]-normal[1])+(prevnormal[2]-normal[2])*(prevnormal[2]-normal[2]);

 cout << "resid= \n";
 cout << resid << '\n';

	     }

if (leastmin > detA [i] [j])
  {
   leastmin = detA [i] [j];
   finalnormal[0]=normal[0];
   finalnormal[1]=normal[1];
   finalnormal[2]=normal[2];

  }


	 }
    
    }
 }


  /* output of function */

  if (leastmin > 0)
      return 0;
  else
    {

             normal[0]=finalnormal[0];
             normal[1]=finalnormal[1];
             normal[2]=finalnormal[2];

             return 1;
    }

}

/* compute coefficients of det(theta) function */
void DetCheckT::ComputeCoefficients(void)
{
	/* moduli components */
	double c11 = fc_ijkl(0,0);
	double c22 = fc_ijkl(1,1);
	double c33 = fc_ijkl(2,2);
	double c23 = fc_ijkl(1,2);
	double c13 = fc_ijkl(0,2);
	double c12 = fc_ijkl(0,1);
	
	/* stress components */
	double s11 = fs_jl[0];
	double s22 = fs_jl[1];
	double s12 = fs_jl[2];

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
	phi2 = atan2(c2t,s2t)/2.0;
	phi4 = atan2(c4t,s4t)/4.0;
	
	/* amplitudes */
	A0 = (-c12*c12 - 3*c13*c13 + c11*c22 + 2*c13*c23 - 3*c23*c23 +
			3*c11*c33 - 2*c12*c33 + 3*c22*c33 + 3*c11*s11 + c22*s11 +
			4*c33*s11 + 3*s11*s11 + 4*c13*s12 + 4*c23*s12 + 4*s12*s12 +
			c11*s22 + 3*c22*s22 + 4*c33*s22 + 2*s11*s22 + 3*s22*s22)/8;

	A2 = c2t/sin(2.0*phi2);
	A4 = c4t/sin(4.0*phi4);
}

/* ***************************************************************** */
/* closed-form check for localization, assuming plane strain condition */
/* 1 is 11 */
/* 2 is 22 */
/* 3 is 12 */
/* angle theta subtends from the x1 axis to the band normal */
int DetCheckT::SPINLOC_localize(double *c__, double *thetan, int *loccheck)
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
    double p, q, r__, x[3], theta, third, a0, a1, a2, a3, a4, pi, 
	    qq, rad;

    /* Parameter adjustments */
    c__ -= 4;

    /* Function Body */
    half = one / two;
    third = one / three;
    rad = four * atan(one) / 180.;
    pi = four * atan(one);



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
		cout << "\n DetCheckT::SPINLOC_localize: a is positive when it should be negative" << endl;
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
