/* $Id: MemoryGroupT.h,v 1.1.1.1 2001-01-25 20:56:22 paklein Exp $ */
/* created: paklein (04/17/1998)                                          */
/* Base class to handle memory (re-/de-) allocation for                   */
/* derived classes managing grouped arrays with memory                    */
/* divided into equally-sized blocks                                      */
/* NOTE: derived class need to define a function which                    */
/* (1) determines the new block size based on derived                     */
/* class parameters.                                                      */
/* (2) resets the block size with a call to SetBlockSize                  */
/* (3) resets the pointers and size parameters for all                    */
/* members of fArrays.                                                    */
/* Since the number of arguments for (1) and (3) can be just              */
/* about anything, no prototypes are provided in this base                */
/* class.                                                                 */

#ifndef _MEMORYGROUP_T_H_
#define _MEMORYGROUP_T_H_

/* direct members */
#include "ArrayT.h"
#include "AutoArrayT.h"

template <class TYPE>
class MemoryGroupT
{
public:

	/* constructor */
	MemoryGroupT(int headroom);

	/* destructor */
	~MemoryGroupT(void);

	/* over-allocation parameter */
	int HeadRoom(void) const;
	void SetHeadRoom(int headroom);

	/* add array to list of managed */
	void Register(ArrayT<TYPE>& array);
	bool IsRegistered(const ArrayT<TYPE>& array) const;

protected:

	/* current block size */
	int BlockSize(void) const;

	/* return a pointer to the specified block */
	TYPE* BlockPointer(int block) const;

	/* memory (re-) allocation and copy old data if specified */
	void SetBlockSize(int newblocksize, bool copy_in);
	
private:

	/* NOT DEFINED - no copy construction */
	MemoryGroupT(const MemoryGroupT& source);

	/* NOT DEFINED - no assignment operator */
	MemoryGroupT& operator=(MemoryGroupT&);

protected:

	/* list of managed */
	AutoArrayT<ArrayT<TYPE>*> fArrays;

private:

	/* oversize parameter */
	int fHeadRoom;
	
	/* grouped data */
	TYPE* fData;
	int   fBlockSize;
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructor */
template <class TYPE>
MemoryGroupT<TYPE>::MemoryGroupT(int headroom):
	fHeadRoom(headroom),
	fData(NULL),
	fBlockSize(0)
{
	/* error check */
	if (fHeadRoom < 0) throw eGeneralFail;
}

/* destructor */
template <class TYPE>
MemoryGroupT<TYPE>::~MemoryGroupT(void)
{
	delete[] fData;
	fData = NULL;
}

/* over-allocation parameter */
template <class TYPE>
inline int MemoryGroupT<TYPE>::HeadRoom(void) const
{
	return fHeadRoom;
}

template <class TYPE>
inline void MemoryGroupT<TYPE>::SetHeadRoom(int headroom)
{
	fHeadRoom = headroom;

	/* check */
	if (fHeadRoom < 0) throw eGeneralFail;
}


/* add array to list of managed */
template <class TYPE>
void MemoryGroupT<TYPE>::Register(ArrayT<TYPE>& array)
{
	/* only until memory is allocated */
	if (fBlockSize > 0)
	{
		cout << "\n MemoryGroupT<TYPE>::Register: all arrays must be registered\n";
		cout <<   "     before initial allocation\n" << endl;
		throw eGeneralFail;
	}

	/* add to list */
	fArrays.Append(&array);
}

template <class TYPE>
bool MemoryGroupT<TYPE>::IsRegistered(const ArrayT<TYPE>& array) const
{
	for (int i = 0; i < fArrays.Length(); i++)
		if (fArrays[i] == &array)
			return true;
		
	return false;
}

/**********************************************************************
*  Private
**********************************************************************/

/* current block size */
template <class TYPE>
inline int MemoryGroupT<TYPE>::BlockSize(void) const { return(fBlockSize); }

/* return a pointer to the specified block */
template <class TYPE>
TYPE* MemoryGroupT<TYPE>::BlockPointer(int block) const
{
	/* range check */
	if (block < 0 || block >= fArrays.Length()) throw(eOutOfRange);
	
	return(fData + fBlockSize*block);
}

/* memory (re-) allocation */
template <class TYPE>
void MemoryGroupT<TYPE>::SetBlockSize(int newblocksize, bool copy_in)
{
	/* take the smaller */
	int copysize = (fBlockSize < newblocksize) ? fBlockSize:newblocksize;

	/* new memory (with extra space) */
	newblocksize += newblocksize*fHeadRoom/100;
	TYPE* newdata;

#ifdef __NEW_THROWS__
	try { newdata = new TYPE[fArrays.Length()*newblocksize]; }
	catch (bad_alloc) { newdata = NULL; }
#else
	newdata = new TYPE[fArrays.Length()*newblocksize];
#endif

	if (!newdata) throw(eOutOfMemory);

	/* copy data in */
	if (copy_in)
	{
		TYPE* pold = fData;
		TYPE* pnew = newdata;
		for (int i = 0; i < fArrays.Length(); i++)
		{
			memcpy(pnew,pold,sizeof(TYPE)*copysize);
			
			pold += fBlockSize;
			pnew += newblocksize;
		}
	}

	/* reset grouped pointer */
	delete[] fData;
	fData      = newdata;
	fBlockSize = newblocksize;
}

#endif /* _MEMORYGROUP_T_H_ */
