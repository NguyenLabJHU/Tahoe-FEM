/* $Id: SimoQ1P0.h,v 1.3 2002-10-05 20:10:45 paklein Exp $ */
#ifndef _SIMO_Q1_P0_H_
#define _SIMO_Q1_P0_H_

/* base classes */
#include "UpdatedLagrangianT.h"

namespace Tahoe {

/** finite strain, mixed element formulation.
 * Formulation due to Simo, Taylor, and Pister, CMAME \b 51, 
 * 177-208, 1985. The element is formulated in three
 * dimensions; however, both hexahedral and quadrilateral
 * element geometries are allowed. In two dimensions, plane
 * strain is assumed, i.e., the out-of-plane stretch is assumed
 * to be unity. \note current implementation does not
 * implement the full consistent tangent. */
class SimoQ1P0: public UpdatedLagrangianT
{
public:

	/** constructor */
	SimoQ1P0(const ElementSupportT& support, const FieldT& field);

	/** destructor */
	~SimoQ1P0(void);

	/** data initialization */
	virtual void Initialize(void);

	/** finalize current step - step is solved */
	virtual void CloseStep(void);
	
	/** restore last converged state */
	virtual void ResetStep(void);

	/** read restart information from stream */
	virtual void ReadRestart(istream& in);

	/** write restart information from stream */
	virtual void WriteRestart(ostream& out) const;

protected:

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/** form the element stiffness matrix. \note this is not the
	 * consistent linearization of SimoQ1P0::FormKd. Only B is
	 * replaced with B-bar. */
	virtual void FormStiffness(double constK);

	/** calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);

	/** read materials data. For two dimensional problems, checks
	 * that all materials are plane strain. */
	virtual void ReadMaterialData(ifstreamT& in);	
	
private:

	/** compute mean shape function gradient, H (reference volume), and
	 * current element volume, equation (2.20) */
	void SetMeanGradient(dArray2DT& mean_gradient, double& H, double& v) const;
	
	/** */
	void bSp_bRq_to_KSqRp(const dMatrixT& b, dMatrixT& K) const;
	
protected:

	/** \name element volume */
	/*@{*/
	/** deformed element volume */
	dArrayT fElementVolume;

	/** deformed element volume from the last time step */
	dArrayT fElementVolume_last;
	
	/** flag to indicate SimoQ1P0::fElementVolume_last has been initialized */
	bool fLastVolumeInit;
	/*@}*/
	
	/** element pressure. Calculated during SimoQ1P0::FormKd. */
	dArrayT fPressure;

	/** \name work space */
	/*@{*/
	dArray2DT fMeanGradient; /**< mean gradient over element */
	dMatrixT fF_tmp; /**< F workspace */
	dMatrixT fNEEmat; /**< dimension of stiffness matrix */
	dMatrixT fdiff_b;
	dMatrixT fb_bar;
	dMatrixT fb_sig;
	dMatrixT fStressStiff2;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _SIMO_Q1_P0_H_ */
