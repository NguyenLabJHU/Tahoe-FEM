/* $Id: DetCheckT.cpp,v 1.39 2005-03-04 04:06:07 cfoster Exp $ */
/* created: paklein (09/11/1997) */
#include "DetCheckT.h"
#include <math.h>
#include "ofstreamT.h"
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dMatrixEXT.h"
#include "dArrayT.h"
#include "dTensor4DT.h"

/* needed to access element information */
#include "SolidMatSupportT.h"

using namespace Tahoe;

/* initialize static variables */
bool DetCheckT::fFirstPass = true;

/* constants */
const double Pi = acos(-1.0);

/* constructor */
DetCheckT::DetCheckT(const dSymMatrixT& s_jl, const dMatrixT& c_ijkl, const dMatrixT& ce_ijkl):
	fs_jl(s_jl),
	fc_ijkl(c_ijkl),
	fce_ijkl(ce_ijkl),
	fStructuralMatSupport(NULL)
//	,fElement(NULL)
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


/* Tests whether two normals are equivalent, checking to see whether they
 * are the same within a given tolerance tol, or opposites within that same
 * tolerance (opposite normals are considered equilvalent for this purpose).
 * Returns true if the normals are distinct, false if not. 
 * Both vectors should be normalized before being passed into this function. 
 */
bool NormalCompare(const dArrayT& normal1, const dArrayT& normal2, double tol);
bool NormalCompare(const dArrayT& normal1, const dArrayT& normal2, double tol)
{
	dArrayT resid(3), negResid(3);
	
	resid.DiffOf(normal1, normal2);
	negResid.SumOf(normal1, normal2);
	
	double residMag = resid.Magnitude();
	double negResidMag = negResid.Magnitude();

	return ( residMag*residMag < tol || negResidMag*negResidMag < tol);
};

/* As above, but fits form for comparator function in AutoArrayT. 
 * 10e-10 seems to be a good tolerance for this comparison
 */
bool NormalCompare(const dArrayT& normal1, const dArrayT& normal2);
bool NormalCompare(const dArrayT& normal1, const dArrayT& normal2)
{
	return NormalCompare(normal1, normal2, 10e-10);
};

/*
* Returns 1 if acoustic tensor isn't positive definite,
* and returns the normal to the surface of localization.
* Returns 0, otherwise.
*/
int DetCheckT::IsLocalized(dArrayT& normal)
{
	if (fs_jl.Rows() == 2) return DetCheck2D(normal);
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

/*
bool DetCheckT::IsLocalized_SS(AutoArrayT <dArrayT> &normals,
							AutoArrayT <dArrayT> &slipdirs)
{
  double dummy =1.0;
  return IsLocalized_SS( normals, slipdirs, dummy);
}
*/

bool DetCheckT::IsLocalized_SS(AutoArrayT <dArrayT> &normals,
							AutoArrayT <dArrayT> &slipdirs, double detA)
{
  cout << "in DetCheckT::IsLocalized_SS\n"; 

  int nsd =fs_jl.Rows();
	dArrayT normal(nsd), slipdir(nsd);
	dTensor4DT C(nsd, nsd, nsd, nsd);
	dMatrixEXT A(nsd); //acoustic tensor 
	/* for eigen analysis */
	dArrayT realev(3), imev(3), altnormal_i(3), altnormal_ii(3);
	//double detA;
	
	if (fs_jl.Rows() == 2)
	{
		//this doesn't work for 2D, comment out slip dir calc for now
		//C.ConvertTangentFrom2DTo4D(C, fc_ijkl);
		/* call SPINLOC routine */
		double eigVal;
		AutoArrayT <double> theta;
		int numev = 0; 
		bool check = false; 
		/* clear normals and slipdirs */
		normals.Free();
		slipdirs.Free();
		theta.Free();
		SPINLOC_localize(fc_ijkl.Pointer(), theta, &check);
		if (check)
		{
		  theta.Top();
		  while (theta.Next())
		    {
		      normal[0] = cos(theta.Current());
		      normal[1] = sin(theta.Current());
		      normals.Append(normal);
		      
			  //form acoustic tensor
			  A = 0.0;
			  A(0,0) = normal[0] * fc_ijkl(0,0) * normal[0] + normal[0] * fc_ijkl(0,2) * normal[1] 
				+ normal[1] * fc_ijkl(2,0) * normal[0] + normal[1] * fc_ijkl(2,2) * normal[1];
			  A(0,1) = normal[0] * fc_ijkl(0,2) * normal[0] + normal[0] * fc_ijkl(0,1) * normal[1] 
				+ normal[1] * fc_ijkl(2,2) * normal[0] + normal[1] * fc_ijkl(2,1) * normal[1];
			   A(1,0) = normal[0] * fc_ijkl(2,0) * normal[0] + normal[0] * fc_ijkl(2,2) * normal[1] 
				+ normal[1] * fc_ijkl(1,0) * normal[0] + normal[1] * fc_ijkl(1,2) * normal[1];
			  A(1,1) = normal[0] * fc_ijkl(2,2) * normal[0] + normal[0] * fc_ijkl(2,1) * normal[1] 
				+ normal[1] * fc_ijkl(1,2) * normal[0] + normal[1] * fc_ijkl(1,1) * normal[1];
			  
			  //find mininum eigenvalue by quadratic formula
			  double minus_b = A(0,0) + A(1,1);
			  detA = A(0,0) * A(1,1) - A(0,1)*A(1,0);
			  
			  double discriminant = minus_b*minus_b - 4 * detA;
			  if (discriminant < 0)
				{
					cout << "Complex eigenvalue detected in DetCheckT::IsLocalized_SS" << flush;
					throw ExceptionT::kGeneralFail;
				}
			  
			  double lambda_min =  .5*(minus_b - sqrt(discriminant));
			  
			  if (lambda_min > 0)
				{
					cout << "Warning: Minimun eigenvalue is positive even though localization was detected. DetCheckT::IsLocalized_SS" << flush;
				}
				 
			  //find slip direction
			  A.Eigenvector(lambda_min, slipdir);
			  slipdirs.Append(slipdir);
		    }

		} // if(check)
		return check;
	}   //if(fs_jl.Row()s==2)
	else
		//return DetCheck3D_SS(normals,slipdirs);
		return DetCheck3D_SS(normals,slipdirs, detA);
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
		throw ExceptionT::kGeneralFail;
	}
	
	if (det(angle) > 0.0)
		return 0;
	else
	{
		normal.Dimension(2);

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
bool DetCheckT::DetCheck3D_SS(AutoArrayT <dArrayT> &normals,
//							AutoArrayT <dArrayT> &slipdirs)
							AutoArrayT <dArrayT> &slipdirs, double detAmin)
{
	int i,j,k,l,m,n; // counters 
	
	//double detAmin;
	
	/* calculated normal at particular angle increment */
	dArrayT normal(3);

	/* principal tensors under analysis */
	dMatrixEXT A(3), Ae(3), Ainverse(3); //acoustic tensor and inverse
	dTensor4DT C(3,3,3,3), Ce(3,3,3,3); // Rank 4 Tensor version of Tangent Modulus
	dMatrixEXT J(3); // det(A)*C*Ainverse

	/* for eigen analysis */
	dArrayT prevnormal(3), realev(3), imev(3), altnormal_i(3), altnormal_ii(3);

	/* for initial sweep */
	double theta, phi; //horizontal plane angle and polar angle for normal
	double detA [numThetaChecks] [numPhiChecks]; //determinant of acoustic tensor at each increment
	double detAe [numThetaChecks] [numPhiChecks]; //determinant of elastic acoustic tensor at each increment
	int localmin [numThetaChecks] [numPhiChecks]; //1 for local minimum, 0 if not

	/* for refinement of normals */
	double normalTol = 10e-20; // tolerance for convergence of normals, based on square of norm of difference
	int newtoncounter=0; //makes sure Newton iteration doesn't take too long
  
	/* for choosing normals w/ least determinant */
	double setTol=1.0e-5; //setTol=1.0e-7 // tolerance for if normals should be in normal set 
	double leastmin=2.0*setTol; 
	double leastdetAe;

	/* determination of slip direction */
	dArrayT slipdir(3); //normal in direction of slip plane
	double eigVal; //final eigenvalue for determining slip direction
	int numev = 0; //number of eigenvectors for given eigenvalue

	/* variables for normalSet output */
	const int outputPrecision = 10;
	const int outputFileWidth = outputPrecision + 8;
	AutoArrayT <dArrayT> normalSet, slipdirSet;

	//const ElementCardT* element = (fElement) ? &(fElement->CurrentElement()) : NULL;
	//const ElementSupportT* support = (fElement) ? &(fElement->ElementSupport()) : NULL;
	//??SetfStructuralMatSupport(support);
  
	/* Set up output file */
	normalSet.Free();
	ofstreamT normal_out;
	if (fFirstPass) 
	{
		normal_out.open("normal.info");
		fFirstPass = false;
	}
	else normal_out.open_append("normal.info");
  
	normal_out << "\n\n time step    element #    ip# \n"
			<< ((fStructuralMatSupport) ? fStructuralMatSupport->StepNumber() : 0) 
			<< setw(outputFileWidth) 
			<< ((fStructuralMatSupport) ? fStructuralMatSupport->StepNumber() : 0) 
			<< setw(outputFileWidth) 
			<< ((fStructuralMatSupport) ? fStructuralMatSupport->StepNumber() : 0) 
			<< "\n\n"
			<< setw(outputFileWidth) << "approx normal0" << setw(outputFileWidth) << "approx normal1" 
			<< setw(outputFileWidth) << "approx normal2" <<  setw(outputFileWidth) << "normal0" 
			<< setw(outputFileWidth) << "normal1" <<  setw(outputFileWidth) << "normal2" 
			<< setw(outputFileWidth) << "detA" <<  setw(outputFileWidth) << "in normalSet?" << endl;
  
	// initialize variables

	normal=0.0;
	A = 0.0;
	Ae = 0.0;
	leastdetAe = 1.0;
	C = 0.0;
	Ce = 0.0;
	Ainverse = 0.0;
	J = 0.0;
	prevnormal = 0.0;
	realev = 0.0;
	imev = 0.0;
	altnormal_i = 0.0;
	altnormal_ii = 0.0;
	slipdir=0.0;	
	for (m=0; m<numThetaChecks; m++)
		for (n=0; n<numPhiChecks; n++)
		{
			detA [m] [n] = 0.0;
			detAe [m] [n] = 1.0;
			localmin [m] [n] = 0;
		}
  
  
	/* Get C from fc_ijkl */
	C.ConvertTangentFrom2DTo4D(C, fc_ijkl);
	Ce.ConvertTangentFrom2DTo4D(Ce, fce_ijkl);
  
	/* Sweep through angles and find approximate local minima */
	FindApproxLocalMins(detA, localmin, C);
  
	/* sweep angle increments and use Newton iteration to refine minima */
	int maxcount = 100;
	for (i=0; i<numThetaChecks; i++)
	{
		for (j=0 ;j<numPhiChecks; j++)
		{
			if (localmin [i] [j] ==1)
			{
				/* reform starting approximate normal */
				theta=Pi/180.0*sweepIncr*i;
				phi=Pi/180.0*sweepIncr*j;

				normal[0]=cos(theta)*cos(phi);
				normal[1]=sin(theta)*cos(phi);
				normal[2]=sin(phi);
	     
				/* output to normal.info */
				normal_out << setprecision(outputPrecision) << endl << setw(outputFileWidth) 
						<< normal[0] << setw(outputFileWidth) << normal[1] << setw(outputFileWidth) << normal[2];

				/* Iteration to refine normal */
				newtoncounter=0;	    
				prevnormal = 0.0;

				while ( !NormalCompare(normal, prevnormal, normalTol) && newtoncounter <= maxcount )
				{
					newtoncounter++;

					/* if too many iterations */
					if (newtoncounter > maxcount)
					{
						cout << "Warning: Bifurcation check failed. Newton refinement did not converge after 100 iterations. \n"; 
						normal_out << setw(2*outputFileWidth) << "Did not converge";
						//return 8;
						//normal=0.0;
						//return 0;
					}

					prevnormal=normal;
		  
					// initialize acoustic tensor A
					A = 0.0;
					A.formacoustictensor(A, C, normal);
					detA [i] [j] = A.Det();
					Ainverse.Inverse(A);
		  
					// form detAe for relative tolerance check
					Ae = 0.0;
					Ae.formacoustictensor(Ae, Ce, normal);
					detAe [i] [j] = Ae.Det();
		  
					// initialize J
					J = 0.0;

					//form Jmn=det(A)*Cmkjn*(A^-1)jk
					for (m=0;m<3;m++)
						for (n=0;n<3;n++)
							for (k=0;k<3;k++)
								for (l=0;l<3;l++)
									J (m,n) += detA [i] [j]*C(m,k,l,n)*Ainverse(l,k); 

					// find least eigenvector of J, the next approx normal
					normal = ChooseNewNormal(prevnormal, J);		  
				} //end while statement

				/* output info to normal.info */
				normal_out << setw(outputFileWidth) <<  normal[0] << setw(outputFileWidth) << normal [1] 
						<< setw(outputFileWidth) << normal[2] << setw(outputFileWidth) << detA [i] [j];

				/* Determine which normals have least value of DetA and record
				 * them. Typically, there are two distinct normals which produce
				 * the same minimum value. Choose between these later*/
				//if (detA [i] [j] - leastmin < -setTol && (detA [i] [j] - leastmin)/fabs(leastmin) < -setTol )
				if ( fabs(leastmin) < 3.0*setTol )
				{
					if ( detA [i] [j] < setTol || fabs( (detA [i] [j])/(detAe [i] [j]) ) < setTol )
					{
						//clear auto array
						normalSet.Free();

						// add normal to auto array and reset leastmin
						normalSet.Append(normal);
						leastmin = detA [i] [j];
						leastdetAe = detAe [i] [j];

						/* output to normal.info */
						normal_out << setw(outputFileWidth) << "Yes - 1st";	
					}
					else
					{
						/* output to normal.info */
						normal_out << setw(outputFileWidth) << "No";
					}
				}
				else if ( fabs(leastmin) > 3.0*setTol )
				{
					//else if (fabs(detA [i] [j] - leastmin) < setTol || fabs((detA [i] [j] - leastmin)/leastmin) < setTol )
					if ( detA [i] [j] - leastmin < -setTol )
					{
						// add normal to auto array and reset leastmin
						if (normalSet.AppendUnique(normal, NormalCompare))
						{
							leastmin = detA [i] [j];
							leastdetAe = detAe [i] [j];
							normal_out << setw(outputFileWidth) << " Yes - not 1st - but is now least min detA";							
						}
						else
							normal_out << setw(outputFileWidth) << "Already Exists";
					}
					else if ( detA [i] [j] < setTol || fabs( (detA [i] [j])/(detAe [i] [j]) ) < setTol )
					{
						// add normal to auto array and output to normal.info
						if (normalSet.AppendUnique(normal, NormalCompare))
							normal_out << setw(outputFileWidth) << "Yes - Added";
						else
							normal_out << setw(outputFileWidth) << "Already Exists";
					}
				}

			} //end if localmin  
		} //end j		  
	} //end i
	
	slipdirSet.Free();
	normals.Free();
	slipdirs.Free();
	
	if (leastmin/leastdetAe > setTol)  //no bifurcation occured
	{
		detAmin = 1.0;
		return false;
	}
	else //bifurcation occured
	{	
		detAmin = leastmin;
		//choose normal from set of normals producing least detA
		//normal = ChooseNormalFromNormalSet(normalSet, C);
		
		normalSet.Top();
		int num_normals = normalSet.Length();
		int count_normals = 0;
		dArrayT normal_curr(3);
		
		/* transfer normals from normalSet to normals (more than 3 unique?, should not be */
		while (normalSet.Next())
		{
			normal_curr = normalSet.Current();
			normals.Append(normal_curr);
			normal_out << endl << setw(outputFileWidth) <<  normal_curr[0] << setw(outputFileWidth) << normal_curr[1] 
					<< setw(outputFileWidth) << normal_curr[2];
			
			A = 0.0;
			A.formacoustictensor(A, C, normal_curr);
			A.eigvalfinder(A, realev, imev);
			eigVal = realev[0];
			if (realev[1] < eigVal) eigVal = realev[1];
			if (realev[2] < eigVal) eigVal = realev[2];
			A.eigenvector3x3(A, eigVal, numev, slipdir, altnormal_i, altnormal_ii);
			slipdirs.Append(slipdir);
			normal_out << setw(outputFileWidth) <<  slipdir[0] << setw(outputFileWidth) << slipdir[1] 
					<< setw(outputFileWidth) << slipdir[2];

		}	
		return true;
	}
	
} // end DetCheckT::DetCheck3D_SS



/* initial sweep in sweepIncr-degree increments to determine approximate 
 *local minima of determine of acoustic tensor A as fn of normal. Sweeps
 * only over half sphere since A(n) = A(-n)*/
void DetCheckT::FindApproxLocalMins(double detA [numThetaChecks] [numPhiChecks],
		int localmin [numThetaChecks] [numPhiChecks], dTensor4DT& C)
{
	int i,j,k;
	double theta, phi; //horizontal plane angle and polar angle for normal, resp
	dMatrixEXT A(3), Ainverse(3); //acoustic tensor and its inverse
	dArrayT normal (3);

	/* Find determinant of A as function of angle */
	//case where j = numPhiChecks - 1, i.e. pole
	normal = 0.0;
	normal [2] = 1.0;
  
	A = 0.0;
	A.formacoustictensor(A, C, normal);

	detA [0] [numPhiChecks - 1] = A.Det();

	for (i = 1; i < numThetaChecks; i++)
		detA [i] [numPhiChecks - 1] = detA [0] [numPhiChecks - 1];
  
	// general case
	for (i=0; i<numThetaChecks; i++)
		for (j=0; j<numPhiChecks - 1; j++)
		{
			theta=Pi/180.0*sweepIncr*i;
			phi=Pi/180.0*sweepIncr*j;

			normal[0]=cos(theta)*cos(phi);
			normal[1]=sin(theta)*cos(phi);
			normal[2]=sin(phi);

			A.formacoustictensor(A, C, normal);	

			detA [i] [j]= A.Det();
		}
	
	/* Check detA against 4 values around it to see if it is local min */
	for (i=0; i<numThetaChecks; i++)
		for (j=0; j<numPhiChecks; j++)
		{
			if (j == numPhiChecks - 1) //pole, checks only vs. points around it
			{
				if (i == 0)
				{
					localmin [i] [j] = 1;
					for (k=0; k<numThetaChecks; k++)
						if (detA [i] [j] > detA [k] [j - 1])
						{
							localmin [i] [j] = 0;
							break;
						}
				}
			}
			else
			{
				if (detA [i] [j] < detA [(i-1)%numThetaChecks] [j]) 
					if (detA [i] [j] <= detA [(i+1)%numThetaChecks] [j])
						if (detA [i] [j] < detA [i] [j+1])
							if ( j == 0)
							{
								// check against normal across sphere
								/* only need to check to 180 degrees, rest are opposite
								 * of already tested, speeds up plane strain problems */
								if ( i < numThetaChecks/2 ) 
									if (detA [i] [j] <= detA [i + numThetaChecks/2] [j+1])
										localmin [i] [j] = 1;
							}
							else
							{
								//standard case
								if (detA [i] [j] <= detA [i] [j-1])
									localmin [i] [j] = 1;
							}
			}
			
		} //end j, end i  
	 
} // end FindApproxLocalMins


/* Finds next iteration on a normal by finding Eigenvalues of Matrix
 * J and choosing the best one, i.e. closest in norm to previous vector */
dArrayT DetCheckT::ChooseNewNormal(dArrayT& prevnormal, dMatrixEXT& J)
{
	double tol = 10e-10;
	dArrayT normal(3), trialNormal(3);
	dArrayT altnormal_i(3), altnormal_ii(3);
	dArrayT realev(3), imev(3);
	int numev = 0, i;
	dMatrixEXT Atrial(3); //trial acoustic tensor
	double inprod, maxInprod;
	Atrial = 0.0;
	trialNormal = 0.0;

	// chooses eigvector by closest approx to previous normal
	J.eigvalfinder(J, realev, imev);

	for (i=0; i<3; i++)
	{  
		if ( fabs(imev[i])<tol || fabs(0.001*imev[i]/realev[i])<tol )
		{
			J.eigenvector3x3(J, realev[i], numev, trialNormal, altnormal_i, altnormal_ii);
			//J.Eigenvector(realev[i], trialNormal);
			inprod = fabs(trialNormal.Dot(trialNormal, prevnormal));
			if ( i==0 )
			{
				maxInprod = inprod;
				normal = trialNormal;
			}
			else if (inprod > maxInprod)
			{
				maxInprod = inprod;
				normal = trialNormal;
			}	    
		}
	}
	
	return normal;
	
}


/* Chooses normal from a set that have essentially same detA that is least
 * of all the normals. Typically there are two
 */
dArrayT DetCheckT::ChooseNormalFromNormalSet(AutoArrayT <dArrayT> &normalSet, dTensor4DT &C)
{
	dMatrixEXT A(3);
	dArrayT trialNormal(3), bestNormal(3);
	double leastmin, detA;

	bestNormal = 0.0;

	normalSet.Top();

	while (normalSet.Next())
	{
		/* Currently chooses normal by least minimum value of the determinant
		* of the acoustic tensor A. Experience suggests that two unique 
		* normals are created and that this choice is based simply on 
		* numerical error. Better to choose another criterion. Wells and Sluys
		* (2001) suggest maximum plastic dissipation, which may be best
		* handled in the function calling IsLocalized_SS, after it is called.  
		*/

		trialNormal = normalSet.Current();
		A.formacoustictensor(A, C, trialNormal);
		detA = A.Det();
		if (detA < leastmin || normalSet.Position() == 0)
		{
			leastmin = detA;
			bestNormal = trialNormal;
		}
	}      
	
	return bestNormal;

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
bool DetCheckT::SPINLOC_localize(const double *c__, AutoArrayT <double>
&thetan, bool *loccheck)
{
	/* Initialized data */
	double zero = 0.;
	double one = 1.;
	double two = 2.;
	double three = 3.;
	double four = 4.;
	double tol = 1.0e-3;
	double detTol = 1.0e-1; 
	//want very loose tolernace for detTol, algorithm is not exact and can
	// easily be sorted out at element level

	/* System generated locals */
	int i__1;
	double d__1, d__2, d__3, d__4;

	/* Local variables */
	double capa, capb, half, fmin, temp, temp2, temp3, a, b, f;
	AutoArrayT <double> xmin;
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
	if (fabs(qq) < 1e-8) qq = zero;

	temp = p * third;

	if (qq > zero || qq == 0.f) 
	{
		temp2 = one;
		temp3 = -half * b + sqrt(qq);
		if (temp3 < zero) temp2 = -one;
		d__1 = fabs(temp3);
		capa = temp2 * pow(d__1, third);
		temp2 = one;
		temp3 = -half * b - sqrt(qq);
		if (temp3 < zero) temp2 = -one;
		d__1 = fabs(temp3);
		capb = temp2 * pow(d__1, third);
		x[0] = capa + capb - temp;
		x[1] = -(capa + capb) * half - temp;
		x[2] = x[1];
	} 
	else 
	{
		if (a < zero) 
		{
			/* Computing 3rd power */
			d__2 = a * third, d__3 = d__2;
			theta = acos(-half * b / sqrt((d__1 = -(d__3 * (d__2 * d__2)), fabs(d__1))));
			temp2 = two * sqrt((d__1 = -a * third, fabs(d__1)));
			x[0] = temp2 * cos(theta * third) - temp;
			x[1] = -temp2 * cos(theta * third + rad * 60.) - temp;
			x[2] = -temp2 * cos(theta * third - rad * 60.) - temp;
		} 
		else 
		{
			cout << "\n DetCheckT::SPINLOC_localize: a is positive when it should be negative" << endl;
		}
	}
	
	fmin = 1e50;
	n = 3;
	if (fabs(qq) < 1e-8) n = 2;
	if (qq > zero) n = 1;
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) 
	{
		/* Computing 4th power */
		d__1 = x[i__ - 1], d__1 *= d__1;
		/* Computing 3rd power */
		d__2 = x[i__ - 1], d__3 = d__2;
		/* Computing 2nd power */
		d__4 = x[i__ - 1];
		f = a4 * (d__1 * d__1) + a3 * (d__3 * (d__2 * d__2)) + a2 * (d__4 * d__4) + a1 * x[i__ - 1] + a0;

		if ( fabs(f - fmin) < detTol * fabs(f) || fabs(f - fmin) < detTol * fabs(fmin) )
		  {
		    xmin.Append(x[i__ - 1]);
		  }
		else if (f <= fmin) 
		  {
		    xmin.Free();
		    fmin = f;
		    xmin.Append(x[i__ - 1]);
		}
		/* L5: */
	}

	/* .. output */

	xmin.Top();
	while(xmin.Next())
	  thetan.Append( atan(xmin.Current()));

	/*if(fmin.lt.tol) then */
	if (fmin / c__[4] < tol) 
	{
		/* localized */
		*loccheck = true;
	} 
	else 
	{
		/* not localized */
		*loccheck = false;
	}
	
	return false;
	
}

