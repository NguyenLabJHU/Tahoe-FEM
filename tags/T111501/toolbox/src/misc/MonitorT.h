/* $Id: MonitorT.h,v 1.1.1.1 2001-01-25 20:56:25 paklein Exp $ */
/* created: paklein (09/30/1996)                                          */

#ifndef _MONITOR_T_H_
#define _MONITOR_T_H_

/* direct members */
#include "iArrayT.h"
#include "VariArrayT.h"

/* forward declarations */
#include "ios_fwd_decl.h"

class MonitorT
{
public:

	/* predefined status flags */
	enum StatusFlagT {kOFF = 0,
                       kON = 1,
                   kMarked = 2};

	/* constructor */
	MonitorT(void);
	MonitorT(int dim);

	/* resize with the given status */
	void Resize(int length, bool copy_in);
		
	/* sets status flag for index to Marked */
	void Mark(int index);
	
	/* reset functions */
	void MarkedToOFF(void);
	void MarkedToON(void);

	/* checking flags */
	int Count(int value) const;
	
	/* using other status flags */
	void SetFlag(int index, int to);
	void Reset(void);                  //resets changed flag
	void ResetFlags(int from, int to); //resets changed flag
	
	/* reset all status flags to ON */
	void AllToON(void);
	
	/* returns status of the given position */
	int Status(int index) const;
	const iArrayT& Status(void) const;

	/* returns 1 if any flags have changed since the last Reset */
	int Changed(void) const;

	/* print indexes (real numbering) flagged as kMonitorOFF
	 * to the output stream with a header */
	void PrintOFF(ostream& out) const;	
	
	/* print indexes (real numbering) with the prescribed flag value
	 * to the output stream */
	void PrintValued(ostream& out, int value, int wrapat) const;

	/* read/write all status flags */
	void ReadStatus(istream& in);
	void WriteStatus(ostream& out) const;

private:

	int		        fChanged;
	iArrayT	        fStatus;
	VariArrayT<int> fStatusWrapper; //dynamic resize wrapper for fStatus
	
};

/* inlines */

/* checking flags */
inline int MonitorT::Count(int value) const
{
	return(fStatus.Count(value));
}

/* returns 0 if Off else returns 1 */
inline int MonitorT::Status(int index) const { return fStatus[index]; }
inline const iArrayT& MonitorT::Status(void) const {  return fStatus; }

/* returns 1 if any flags have changed since the last Reset */
inline int MonitorT::Changed(void) const { return(fChanged); }

#endif /* _MONITOR_T_H_ */
