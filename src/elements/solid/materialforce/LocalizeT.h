/* $Id: LocalizeT.h,v 1.3 2003-11-19 06:09:46 thao Exp $ */
/* created: paklein (09/11/1997) */

#ifndef _LOCALIZET_H_
#define _LOCALIZET_H_

#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "ofstreamT.h"

namespace Tahoe {

class ElementSupportT;
class LocalArrayT;
class iArrayT;
class dArray2DT;

class LocalizeT
{
public:

       /**constructor*/
	LocalizeT(const ElementSupportT& support);

	/*print localization results*/
	void WriteLocalize(const iArrayT& locflags, const dArray2DT& elem_center, const dArray2DT& normal);
	void WriteLocalize (void);

	/*check localization*/
	int CheckLocalizeFS(const dSymMatrixT& stress, const dMatrixT& modulus,
			    const LocalArrayT& initcoords);
	int CheckLocalizeSS(const dSymMatrixT& stress, const dMatrixT& modulus,
			    const LocalArrayT& initcoords);

	/*accessors*/
	const dArrayT& LocalizedNormal(void) const;
	const dArrayT& LocalizedElemCenter(void) const;

 private:
	/*Coppied straight from DetCheckT*/

	/* 2D determinant check function */
	int DetCheck2D(dArrayT& normal);

	/*compute coefficients of det(theta) function */
	void ComputeCoefficients(void);

	/*closed-form check for localization, assuming plane strain conditions.
	 * Taken from R.A.Regueiro's SPINLOC.
	 * 1 is 11
	 * 2 is 22
	 * 3 is 12
	 * angle theta subtends from the x1 axis to the band normal */
	int SPINLOC_localize(double *c__, double *thetan, int *loccheck);

	/* determinant function and derivatives */
	double det(double t) const;
	double ddet(double t) const;
	double dddet(double t) const;
	
 protected:
	/*check localization flag*/
	int fCheck;	
	iArrayT fBlockList;
	iArrayT fip_loc;
 private:
	
	/*member data*/
	dSymMatrixT fStress;    /*Cauchy stress tensor*/
	dMatrixT fModulus;      /*Spatial tangent modulus*/
	dArrayT fNormal;        /*Normal to localization plane*/
	dArrayT fElemCenter;    /*Center coordinates of localized element*/
	const double& fTime;    /*Calculation time*/

	double fphi2, fphi4;	/* phase shifts */
	double fA0, fA2, fA4;	/* amplitudes   */
	const double fPi;       /* pi*/
       
	/*output*/ 
	ofstreamT fout;         /*output stream*/
	StringT flocalize_file; /*output file name*/


};

 inline const dArrayT& LocalizeT::LocalizedNormal(void) const {return fNormal;}
 inline const dArrayT& LocalizeT::LocalizedElemCenter(void)const
 {
   return fElemCenter;
 }

 /********private***************/
/* inline functions */

 /* determinant function and derivatives */
 inline double LocalizeT::det(double t) const
 {
     return fA0 + fA2*sin(2.0*(fphi2+t)) + fA4*sin(4.0*(fphi4+t));
 }
 
 inline double LocalizeT::ddet(double t) const
 {
   return 2.0*fA2*cos(2.0*(fphi2+t)) + 4.0*fA4*cos(4.0*(fphi4+t));
 }
 
 inline double LocalizeT::dddet(double t) const
 {
   return -4.0*fA2*sin(2.0*(fphi2+t)) - 16.0*fA4*sin(4.0*(fphi4+t));
 }

} // namespace Tahoe 
#endif /* _LocalLizeT_H_ */
