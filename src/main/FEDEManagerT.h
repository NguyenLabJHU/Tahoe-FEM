// This FEDEManagerT class is supposed to be inherited from FEManagerT, 
// and includes an object of DEManagerT class, i.e., the Discrete Element
// manager.
#ifndef _FE_DE_MANAGER_H_
#define _FE_DE_MANAGER_H_

/* base class  */
#include "FEManagerT.h"

namespace Tahoe {

class FEDEManagerT: public FEManagerT
{
public:

    /** constructor */
    FEDEManagerT(const StringT& input_file, ofstreamT& output, CommunicatorT& comm,
		 const ArrayT<StringT>& argv, TaskT task);
    
    /** destructor */
    virtual ~FEDEManagerT(void);
    
    /** solve all the time sequences */
    virtual void Solve(void);

    void ApplyDEForce(void);

    void ApplyNodeForce(int node_num, int dof, double force);

protected:
    int ijk;

};

} /* namespace Tahoe */

#endif /*_FE_DE_MANAGER_H_ */
