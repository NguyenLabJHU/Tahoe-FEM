/* $Id: MaterialSupportT.h,v 1.2.8.1 2002-10-26 16:24:21 paklein Exp $ */
#ifndef _MATERIAL_SUPPORT_T_H_
#define _MATERIAL_SUPPORT_T_H_

namespace Tahoe {

/* forward declarations */
class ContinuumElementT;
class ElementCardT;

/** support for the Tahoe materials classes. */
class MaterialSupportT
{
public:

	/** constructor */
	MaterialSupportT(int nsd, int ndof, int nip);

	/** destructor */
	~MaterialSupportT(void);

	/** \name dimensions */
	/*@{*/
	/** number of spatial dimensions */
	int NumSD(void) const { return fNumSD; };
	
	/** number of degrees of freedom (per node) */
	int NumDOF(void) const { return fNumDOF; };

	/** stress evaluation points per element */
	int NumIP(void) const { return fNumIP; };
	/*@}*/

	/** \name run time status */
	/*@{*/
	/** current stress evaluation point within the element. If
	 * no source for the current point is set using 
	 * MaterialSupportT::SetCurrIP, will return -1. */
	int CurrIP(void) const;

	/** the iteration number for the current time increment. If
	 * no source for the iteration number is set using 
	 * MaterialSupportT::SetIterationNumber, will return -1. */
	int IterationNumber(void) const;

	/** set the source for the current stress evaluation point */
	void SetCurrIP(const int& curr_ip);

	/** set the source for the iteration number */
	void SetIterationNumber(const int& iter);
	/*@}*/
	
	/** \name host code information */
	/*@{*/
	/** return a pointer to the host element. Returns NULL if no
	 * no element information in available. The ContinuumElementT
	 * pointer is set using MaterialSupportT::SetContinuumElement. */
	const ContinuumElementT* ContinuumElement(void) const;

	int NumElements(void) const;
	int CurrElementNumber(void) const;
	ElementCardT& ElementCard(int card) const;
	ElementCardT& CurrentElement(void) const;
	/*@}*/
	
	/** set the element group pointer */
	void SetContinuumElement(ContinuumElementT* p);

  private:
  
  	/** \name dimensions */
  	/*@{*/
	/** number of degrees of spatial dimensions */
	int fNumSD;

	/** number of degrees of freedom */
	int fNumDOF;
	
	/** number of integration points */
	int fNumIP;
  	/*@}*/
  	
  	/** pointer to run time information */
  	/*@{*/
	const int* fCurrIP;
	const int* fIterationNumber;
  	/*@}*/
  
  	/** pointer to the continuum element */
  	ContinuumElementT* fContinuumElement;

};

/* inlines functions */
inline const ContinuumElementT* MaterialSupportT::ContinuumElement(void) const
{
	return fContinuumElement;
}

inline void MaterialSupportT::SetContinuumElement(ContinuumElementT* p)
{
	fContinuumElement = p;
}

/* run time status */
inline int MaterialSupportT::CurrIP(void) const
{
	if (fCurrIP) return *fCurrIP;
	else return -1;
}

inline int MaterialSupportT::IterationNumber(void) const
{
	if (fIterationNumber) return *fIterationNumber;
	else return -1;
}

inline void MaterialSupportT::SetCurrIP(const int& curr_ip)
{
	fCurrIP = &curr_ip;
}

inline void MaterialSupportT::SetIterationNumber(const int& iter)
{
	fIterationNumber = &iter;
}

} /* namespace Tahoe */
#endif /* _SS_HOOKEAN_MAT_H_ */
