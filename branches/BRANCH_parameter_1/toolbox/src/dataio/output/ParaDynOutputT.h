#ifndef _PARADYNOUTPUT_T_H_
#define _PARADYNOUTPUT_T_H_

#include "OutputBaseT.h"
#include "AutoArrayT.h"
#include "StringT.h"
#include "ParaDynT.h"

#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"

namespace Tahoe {

class ParaDynOutputT : public OutputBaseT
  {
  public:
    ParaDynOutputT (ostream& out, 
		    const ArrayT<StringT>& out_strings);
    void WriteGeometry (void);
    
  private:
    
    void OpenGeometryFile (ParaDynT& par, ofstream& geo) const;
    StringT CreateFileName (const StringT& label) const;
    
    void WriteBounds (ostream& geo, const ParaDynT& par) const;

    void WriteCoordinates (ostream& geo, ParaDynT& par, 
			   const iArrayT& nodes_used) const;
    void WritePart (ostream& geo, ParaDynT& par, int index) const;

  };
 
} // namespace Tahoe 
#endif

