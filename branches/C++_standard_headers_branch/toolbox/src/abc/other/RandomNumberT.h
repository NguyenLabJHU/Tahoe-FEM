/* $Id: RandomNumberT.h,v 1.6 2003-11-21 22:41:39 paklein Exp $ */
#ifndef _RANDOM_NUMBER_T_H_
#define _RANDOM_NUMBER_T_H_

namespace Tahoe {

/* forward declarations */
class dArrayT;
class dArray2DT;
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
	
	/** set the seed */
	void sRand(long seed, long a = 16807, long rm = 2147483647);

	/** return a random number of the appropriate type */
	double Rand(void);

	/* Fill arrays with random numbers */
	dArrayT& RandomArray(dArrayT& fillArray);
	dArray2DT& RandomArray(dArray2DT& fillArray);
	
	/** Accessor for the seed */
	long RandSeed(void);

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
