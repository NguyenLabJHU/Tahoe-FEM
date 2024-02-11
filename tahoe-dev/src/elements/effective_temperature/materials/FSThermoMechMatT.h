
#if !defined(_FSThermoMechMatT_)
#define _FSThermoMechMatT_

#include <cassert>

#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "FSSolidMatT.h"
#include "FSThermoMechSupportT.h"
#include "RGSplitT2.h"

namespace Tahoe {

  class FSThermoMechMatT: public RGSplitT2
  {

  public:

    /* constructors */
    FSThermoMechMatT(void);

    /* information about subordinate parameter lists*/
    virtual void DefineSubs(SubListT& sub_list) const;

    /*dS/dT,return fCoupledModulus */
    virtual const dSymMatrixT& d_ij(void);

    /* spatial dphi_visc/dC, where phi =Ce Dev(Sneq):Lv is the viscous dissipation */
    virtual const dSymMatrixT& h_ij(void);

    /*spatial heat flux*/
    virtual const dArrayT& q_i(void);

    /*b1 is the d(themrmal body force)/dT as (d(rho*deltac*dotTf-Wp)/dT) */
    virtual  double b1(void);
      virtual double dtime(void);

    /*this is the residual for the thermal part represented as:rho*deltac*dotTf-Wp */
    virtual   double heatres(void);
 	
    const double Density(void) const;
    const double SpecificHeat(void) const;
    virtual  const double Capacity(void) const;
    virtual const dMatrixT& Conductivity(void);

    void SetFSThermoMechSupport(const FSThermoMechSupportT* support);
    
      const dArrayT TemperatureGradient(void);
      
      /*temperature gradient at given integration point*/
      const dArrayT TemperatureGradient(int ip);
      
      /* Number of stress and structural relaxation process*/
      const int NumStressRelaxation(void) const;
      const int NumStructuralRelaxation(void) const;
      
      /*retrieve effective temeperatures at all ip*/
      virtual const dArrayT& Effective_Temperature(void) = 0;
      virtual const dArrayT& Effective_Temperature_last(void) = 0;

  protected:
      virtual void Initialize(void);

    const FSThermoMechSupportT* fFSThermoMechSupport;
	/** \name parameters */
	/*@{*/
	double   frho;
      double fcg;
    double   fCapacity;
	/*@}*/
	  
	/** \name return values */
	/*@{*/
	/** heat flux return value */
	dArrayT fHeatFlux;
	  
	/** conductivity */
      double fk;
      dMatrixT fkij;
      dSymMatrixT fhij;

	  /** thermomechanical coupling*/
	  dSymMatrixT fCoupledModulus;
      dArrayT fTemperatureGradient;
	/*@}*/

      /*number of structural relaxation process*/
      /*number of stress relaxation process is given by fNumProcess defined in RGVsicoelasticity*/
      int fNumR;

      /*for spectral decomposition of spatial tensors into eigevnvalues and eigventors*/
      SpectralDecompT fSpectralDecompSpat;
      

  };

	/* returns the density */
	inline const double FSThermoMechMatT::Density(void) const { return frho; }

	/* returns the glassy specific heat capacity */
	inline const double FSThermoMechMatT::SpecificHeat(void) const { return fcg; }

	/* returns glassy heat capacity */
	inline const double FSThermoMechMatT::Capacity(void) const { return frho*fcg; }

	/*returns number of stress relaxation processes*/
	inline const int FSThermoMechMatT::NumStressRelaxation(void) const {return fNumProcess;}

	/*returns number of structural relaxation processes*/
	inline const int FSThermoMechMatT::NumStructuralRelaxation(void) const {return fNumR;}

	
	/* conductivity */
    //Rui:should be the same as the conducutivity
//	inline const dMatrixT& FSThermoMechMatT::Conductivity(void) const {return fkij; }
} // namespace Tahoe

#endif // _FSThermoMechMatT_
