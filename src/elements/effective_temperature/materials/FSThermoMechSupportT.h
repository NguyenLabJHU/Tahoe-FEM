#if !defined(_FSThermoMechSupportT_)
#define _FSThermoMechSupportT_

#include <cassert>

#include "FSMatSupportT.h"

namespace Tahoe {
    
    /* forward declarations*/
    class FSThermoMechT;
    
    /* support for effective temperature class. interface to pass strains and temperature gradients from element class to material class*/
    class FSThermoMechSupportT: public FSMatSupportT
    {
        
    public:
        /* constructor*/
        FSThermoMechSupportT(int ndof, int nip);
        /* \name host code information
           @{
           return a pointer to the host element. Returns 0 if no
            element information in available.
         */
        const FSThermoMechT* FSThermoMech() const;

        /* set the element group pointer.  this allows material to access element data*/
        virtual void SetContinuumElement(const ContinuumElementT* p);

        /* temperature gradient at current integration point. accesor for material class*/
        const dArrayT& TemperatureGradient() const;
        
        /* temperature gradient at given integration point*/
        const dArray2DT& TeGradient_last(int ip) const;
   
        /* temperature gradient at current integration point. accesor for material class*/
        const dArray2DT& TeGradient_last() const;
        
        /* temperature gradient at given integration point*/
        const dArrayT& TemperatureGradient(int ip) const;
        /* Set source for temperature gradient*/
        void SetTemperatureGradient(const ArrayT<dArrayT>* TempGrad_List);
        
        /* Set source for effective temperature gradient*/
        void SetTeGradient_last(const ArrayT<dArray2DT>* GradTe_last_List);
        
    private:
        
        /* \name Sources for electric field and pointers to localarrays
         @{
        Temperature field
        */
        const ArrayT<dArrayT>* fTempGrad_List;
        const ArrayT<dArray2DT>* fGradTe_last_List;
        
        /* pointer to the host element*/
        const FSThermoMechT* fFSThermoMech;
        
    };
    
    
	inline const FSThermoMechT*	FSThermoMechSupportT::FSThermoMech() const
	{
		return fFSThermoMech;
	}
    
	inline const dArrayT& FSThermoMechSupportT::TemperatureGradient() const
	{
		if (!fTempGrad_List) throw ExceptionT::kGeneralFail;
		return (*fTempGrad_List)[CurrIP()];
	}
	
    /* temperature gradient at given integration point*/
    inline const dArrayT& FSThermoMechSupportT::TemperatureGradient(int ip) const
	{
		if (!fTempGrad_List) throw ExceptionT::kGeneralFail;
		return (*fTempGrad_List)[ip];
	}

    inline const dArray2DT& FSThermoMechSupportT::TeGradient_last() const
	{
		if (!fGradTe_last_List) throw ExceptionT::kGeneralFail;
		return (*fGradTe_last_List)[CurrIP()];
	}
	
    /* temperature gradient at given integration point*/
    inline const dArray2DT& FSThermoMechSupportT::TeGradient_last(int ip) const
	{
		if (!fGradTe_last_List) throw ExceptionT::kGeneralFail;
		return (*fGradTe_last_List)[ip];
	}

    
	
} // namespace Tahoe

#endif // _FSThermoMechSupportT_
