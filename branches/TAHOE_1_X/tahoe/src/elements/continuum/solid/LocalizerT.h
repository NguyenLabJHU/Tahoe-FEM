/* $Id: LocalizerT.h,v 1.7 2004-06-17 07:40:07 paklein Exp $ */
/* created: paklein (02/19/1998) */
#ifndef _LOCALIZER_T_H_
#define _LOCALIZER_T_H_

/* base class */
#include "UpdatedLagrangianT.h"

/* direct members */
#include "ofstreamT.h"
#include "MonitorT.h"
#include "AutoArrayT.h"
#include "dMatrixEXT.h"

namespace Tahoe {

class LocalizerT: public UpdatedLagrangianT
{
public:

	/* constructors */
	LocalizerT(const ElementSupportT& support, const FieldT& field);

	/* set work space */
	virtual void Initialize(void);

	/* finalize time increment */
	virtual void CloseStep(void);

	/* writing results */
	virtual void WriteOutput(void);

	/* returns true if the internal force has been changed since
	 * the last time step */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/* initial condition/restart functions
	 *
	 * Set to initial conditions.  The restart functions
	 * should read/write any data that overrides the default
	 * values */
	virtual void InitialCondition(void);
	virtual void ReadRestart(istream& in);
	virtual void WriteRestart(ostream& out) const;
	
	/* resets to the last converged solution */
	virtual GlobalT::RelaxCodeT ResetStep(void);
		
protected:

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;

	/* increment current element */
	virtual bool NextElement(void);

	/* element data */
	virtual void ReadMaterialData(ifstreamT& in);

	/* add localization checking during solution finding */
	virtual void FormStiffness(double constK);
	
private:

	/* element specific input */
	void EchoData(ifstreamT& in, ostream& out);

	/* element localization check driver - loop over all/list */
	void Localization(void);

	/* check localization in current element returns the number of
	 * localized ip's and writes their coordinates to out */
	int CheckLocalization(ostream& out);

	/* output list of localized element numbers */
	void PrintLocalized(ostream& out);

	/* check strain magnitude in current element */
	void CheckStrain(void);
		//this needs to be a material call.

protected:

	/* flags */
	int    fStrainCheckFlag;
	double fCriticalStretch;
	
	int fLocCheckInc;
	int fLocCheckCount;
	// NOTE: currently checking strain magnitude are localization
	//       conditions and not implemented in a  unified way. use
	//       one or the other.

	/* output stream for localization data */
	ofstreamT fLocOut;
	ofstreamT fLocTOC;

	/* Element status monitor */
	MonitorT fElementMonitor;
	// need better plan for storage size fNumElements for
	// variable size element groups
	
	/* list of elements checked for localization */
	AutoArrayT<int> fLocCheckList;
	// list of zero length implies ALL are checked
	
	/* element neighborlists */
	iArray2DT fNeighborList;
	
	/* work space */
	dArrayT     fNormal;
	dSymMatrixT fAvgStretch;

private: /* localization workspace */

	iArrayT  fnoRBeqs; // element geometry dependent
	dMatrixT fKloc;
	
	/* stiffness matrix eigenvalues */
	dMatrixEXT fEigSolver;
	dArrayT    fEigs;

//NOTE - want to send localized elements as separate "changing geometry"
//       group or just send flags as "element values"???
//	iArray2DT fOutputConn;
//	int fLocalizerConnects;
};

} // namespace Tahoe 
#endif /* _LOCALIZER_H_ */
