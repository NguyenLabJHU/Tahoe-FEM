/* $Id: RandomNumberT.h,v 1.4 2003-04-14 17:26:21 cjkimme Exp $ */
#ifndef _RANDOM_NUMBER_T_H_
#define _RANDOM_NUMBER_T_H_

/* Environmental */
//#include "Environment.h"

//#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class dArrayT;
class ifstreamT;

class RandomNumberT
{
public:
	enum DistributionT {
		kUniform = 0,
		kGaussian = 1,
		kParadynUniform = 2,
		kParadynGaussian = 3
	};

	/** \name constructors */
	/*@{*/
	RandomNumberT(DistributionT type = kUniform);
	/*@}*/

	/*@{*/
	RandomNumberT(ifstreamT& in);
	/*@}*/
	
	void sRand(long seed, long a = 16807, long rm = 2147483647);

	double Rand(void);

	
	dArrayT& RandomArray(dArrayT& fillArray);
	
	long RandSeed(void);

    // Need some way to write seed to restart. 

 private:

	double UniformRandom(void);

	double GaussianRandom(void);

	double ParadynUniformRandom(void);

	long fseed, fa, frm; 
	double dseed, da, drm;
	DistributionT randomType, uniformType;

	double (RandomNumberT::*randFunc)(void);
	double (RandomNumberT::*uniformFunc)(void);

};

inline double RandomNumberT::Rand(void)
{
	return (this->*randFunc)();
}


}//namespace Tahoe
#endif /* _RANDOM_NUMBER_T_H_ */
