#include "FSThermoMechSupportT.h"
#include "FSThermoMechT.h"

namespace Tahoe{
        
    FSThermoMechSupportT::FSThermoMechSupportT(int ndof, int nip):
        FSMatSupportT(ndof, nip),
        fTempGrad_List(NULL),
        fGradTe_last_List(NULL),
        fFSThermoMech(NULL)
	{
		
	}

    void FSThermoMechSupportT::SetContinuumElement(const ContinuumElementT* p)
    {
    /* inherited*/
        FSMatSupportT::SetContinuumElement(p);
        fFSThermoMech = dynamic_cast<const FSThermoMechT*>(p);
    }
    
    /* set source for the deformation gradient */
    void FSThermoMechSupportT::SetTemperatureGradient(const ArrayT<dArrayT>* temperaturegradient)
    {
        fTempGrad_List = temperaturegradient;
    }
    
    void FSThermoMechSupportT::SetTeGradient_last(const ArrayT<dArray2DT>* GradTe_last_List)
    {
        /* keep pointer */
        fGradTe_last_List = GradTe_last_List;
    }


} //namespace Tahoe
