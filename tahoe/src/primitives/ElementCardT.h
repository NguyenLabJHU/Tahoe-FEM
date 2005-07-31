/* $Id: ElementCardT.h,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (05/24/1996)                                          */
/* Empty organizer class - needs manager class to control data.           */

#ifndef _ELEMENT_CARD_T_H_
#define _ELEMENT_CARD_T_H_

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class dMatrixT;
class ElementStorageT;

class ElementCardT
{
public:

	/* constructors */
	ElementCardT(void);
	ElementCardT(const ElementCardT& source);
	
	/* destructor */
	~ElementCardT(void);

	/* assignment operator */
	ElementCardT& operator=(const ElementCardT& rhs);
	
	/* set material number */
	void SetMaterialNumber(int matnum);

	/* setting/getting the activity flags */
	int& Flag(void);
						
	/* accessors */
	int MaterialNumber(void) const;
	iArrayT& NodesX(void);             // geometry nodes
	const iArrayT& NodesX(void) const; // geometry nodes
	iArrayT& NodesU(void);             // field nodes
	iArrayT& Equations(void);

	/* reset field nodes array pointer (non-isoparametric) */
	void SetNodesU(iArrayT& nodesU);
		
	/* restart operations */
	void ReadRestart(istream& in);
	void WriteRestart(ostream& out) const;

	/* element storage accessors/modifiers */
	int IsAllocated(void) const;
	void Allocate(int i_size, int d_size);
	iArrayT& IntegerData(void) const;
	dArrayT& DoubleData(void) const;
	
private:
	
	int fMatNum;
	int fFlag;

	/* geometry nodes */
	iArrayT fNodesX;

	/* field nodes */
	iArrayT* fNodesU; // &fNodesX by default
	iArrayT	 fEqnos;
	
	/* element storage */
	ElementStorageT* fData;

	/* junk - return values if not allocated */
	static iArrayT i_junk;
	static dArrayT d_junk;
};

class ElementStorageT
{
	friend class ElementCardT;

private:

	/* constructor */
	ElementStorageT(int i_size, int d_size);
	ElementStorageT(const ElementStorageT& source);
	
	/* I/O operators */
	friend istream& operator>>(istream& in, ElementStorageT& data);
	friend ostream& operator<<(ostream& out, const ElementStorageT& data);

	/* assignment operator */
	ElementStorageT& operator=(const ElementStorageT& rhs);
	
private:

	/* data */
	iArrayT fIntegerData;
	dArrayT fDoubleData;
};

/* in-lines */

/* setting/getting the activity flags */
inline int& ElementCardT::Flag(void) { return fFlag; }

/* accessors */
inline int ElementCardT::MaterialNumber(void) const { return fMatNum; }

inline iArrayT& ElementCardT::NodesX(void) { return fNodesX;  }
inline const iArrayT& ElementCardT::NodesX(void) const { return fNodesX;  }
inline iArrayT& ElementCardT::NodesU(void) { return *fNodesU; }
inline iArrayT& ElementCardT::Equations(void) { return fEqnos;   }

/* reset field nodes array pointer (non-isoparametric) */
inline void ElementCardT::SetNodesU(iArrayT& nodesU)
{
	fNodesU = &nodesU;
}

/* element storage accessors/modifiers */
inline int ElementCardT::IsAllocated(void) const { return (fData != NULL); }
inline iArrayT& ElementCardT::IntegerData(void) const
{
	return (!fData) ? i_junk : fData->fIntegerData;
}

inline dArrayT& ElementCardT::DoubleData(void) const
{
	return (!fData) ? d_junk : fData->fDoubleData;
}

/* constructor */
inline ElementStorageT::ElementStorageT(int i_size, int d_size):
	fIntegerData(i_size),
	fDoubleData(d_size) { }

/* (deep) copy constructor */
inline ElementStorageT::ElementStorageT(const ElementStorageT& source):
	fIntegerData(source.fIntegerData),
	fDoubleData(source.fDoubleData) { }

/* assignment operator */
inline ElementStorageT& ElementStorageT::operator=(const ElementStorageT& rhs)
{
	fIntegerData = rhs.fIntegerData;
	fDoubleData  = rhs.fDoubleData;

	return *this;
}

#endif /* _ELEMENT_CARD_T_H_ */
