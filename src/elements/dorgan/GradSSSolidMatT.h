/* $Id: GradSSSolidMatT.h,v 1.11 2004-07-15 08:28:12 paklein Exp $ */
#ifndef _GRAD_SS_SOLID_MAT_T_H_
#define _GRAD_SS_SOLID_MAT_T_H_

/* base class */
#include "SSSolidMatT.h"

/* direct members */
#include "dMatrixT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class GradSSMatSupportT;

/** defines the interface for gradient dependent small strain continuum materials */
class GradSSSolidMatT: public SSSolidMatT
{
public:

	/** constructor */
	GradSSSolidMatT(ifstreamT& in, const GradSSMatSupportT& support);
	GradSSSolidMatT(void);

	/** destructor */
	~GradSSSolidMatT(void);

	/* apply pre-conditions at the current time step */
	virtual void InitStep(void);
	
	/** \name field */
	/*@{*/
	const double& Field(void) const;
	const double& Field(int ip) const;
	/*@}*/
	
	/** \name field from the end of the previous time step */
	/*@{*/
	const double& Field_last(void) const;
	const double& Field_last(int ip) const;
	/*@}*/
	
	/** \name gradient field */
	/*@{*/
	const double& GradField(void) const;
	const double& GradField(int ip) const;
	/*@}*/
	
	/** \name gradient field from the end of the previous time step */
	/*@{*/
	const double& GradField_last(void) const;
	const double& GradField_last(int ip) const;
	/*@}*/
	
	/** \name Laplacian field */
	/*@{*/
	const double& LapField(void) const;
	const double& LapField(int ip) const;
	/*@}*/
	
	/** \name Laplacian field from the end of the previous time step */
	/*@{*/
	const double& LapField_last(void) const;
	const double& LapField_last(int ip) const;
	/*@}*/
	
	/** number of degrees of freedom for field field (per node) in the host
	 * element group. */
	int NumDOF_Field(void) const;
	
	/** number of total degrees of freedom (per node) in the host
	 * element group. */
	int NumDOF_Total(void) const;
	
	/** the total number of integration points per element in the
	 * host element group for field field. */
	int NumIP_Field(void) const;
	
	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus for Kaa_bb */
	virtual const dMatrixT& dm_bb_ijkl(void) = 0;
	
	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void) = 0;
	
	/** off diagonal moduli for Kar_bh */
	virtual const dMatrixT& om_bh_ij(void) = 0;
	
	/** off diagonal moduli for Kra_hb */
	virtual const dMatrixT& om_hb_ij(void) = 0;
	
	/** moduli for local term in Krr_hh */
	virtual const dMatrixT& gm_hh(void) = 0;

	/** moduli for gradient term in Krr_hp */
	virtual const dMatrixT& gm_hp(void) = 0;

	/** moduli for gradient term in Krr_hp */
	virtual const dMatrixT& gm_hq(void) = 0;

	/** yield criteria moduli */
	virtual double yc(void) = 0;
	
	/** incremental change in Field_bar */
	virtual double del_Field(void) = 0;
	/*@}*/
	
	/** incremental change in GradField */
	virtual double del_GradField(void) = 0;
	
	/** incremental change in LapField */
	virtual double del_LapField(void) = 0;
	/*@}*/
	
protected:

	/** number of degrees of freedom for field */
	int fNumDOF_Field;
	
	/** total number of degrees of freedom */
	int fNumDOF_Total;
	
	/** number of integration points for field */
	int fNumIP_Field;
	
	/** small strain material support */
	const GradSSMatSupportT* fGradSSMatSupport;
};

/* inlines */
inline int GradSSSolidMatT::NumDOF_Field(void) const { return fNumDOF_Field; }
inline int GradSSSolidMatT::NumDOF_Total(void) const { return fNumDOF_Total; }
inline int GradSSSolidMatT::NumIP_Field(void) const { return fNumIP_Field; }

} // namespace Tahoe 
#endif /* _GRAD_SS_SOLID_MAT_T_H_ */
