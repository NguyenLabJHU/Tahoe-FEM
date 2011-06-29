/* $Id: MemoryGroupT.h,v 1.6 2005-06-05 06:21:07 paklein Exp $ */
/* created: paklein (04/17/1998) */
#ifndef _MEMORYGROUP_T_H_
#define _MEMORYGROUP_T_H_

/* direct members */
#include "ArrayT.h"
#include "AutoArrayT.h"

namespace Tahoe {

/** Base class to handle memory (re-/de-) allocation for
 * derived classes managing grouped arrays with memory
 * divided into equally-sized blocks
 * NOTE: derived class need to define a function which
 * (1) determines the new block size based on derived
 * class parameters.
 * (2) resets the block size with a call to SetBlockSize
 * (3) resets the pointers and size parameters for all
 * members of fArrays.
 * Since the number of arguments for (1) and (3) can be just
 * about anything, no prototypes are provided in this base
 * class. Memory for the arrays in the group may be either
 * pooled or separated. Pooled memory is preferable for a
 * large number of small arrays while separate memory is
 * preferable for a small number of large arrays.
 */
template <class TYPE>
class MemoryGroupT
{
public:

	/** constructor */
	MemoryGroupT(int headroom, bool pool_memory);

	/** destructor */
	~MemoryGroupT(void);

	/** \name over-allocation parameter */
	/*@{*/
	int HeadRoom(void) const;
	void SetHeadRoom(int headroom);
	/*@}*/

	/** \name add array to list of managed */
	/*@{*/
	/** register an array with the group and returns the index of the array */
	int Register(ArrayT<TYPE>& array);
	bool IsRegistered(const ArrayT<TYPE>& array) const;
	int NumRegistered(void) const { return fArrays.Length(); };
	/*@}*/

protected:

	/** current block size */
	int BlockSize(void) const;

	/** return a pointer to the specified block */
	TYPE* BlockPointer(int block) const;

	/** memory (re-) allocation and copy old data if specified */
	void SetBlockSize(int newblocksize, bool copy_in);
	
private:

	/** \name not allowed */
	/*@{*/
	/** copy construction */
	MemoryGroupT(const MemoryGroupT& source);

	/* assignment operator */
	MemoryGroupT& operator=(MemoryGroupT&);
	/*@}*/

protected:

	/** list of managed */
	AutoArrayT<ArrayT<TYPE>*> fArrays;

private:

	/** oversize parameter */
	int fHeadRoom;
	
	/** true if memory for the group is pooled, false if it is separate */
	bool fPoolMemory;
	
	/** \name memory  */
	/*@{*/
	/** memory for the registered arrays. List will be length 1 if 
	 * MemoryGroupT::fPoolMemory is true; otherwise, will be the same
	 * length as MemoryGroupT::fArrays */
	AutoArrayT<TYPE*> fData;

	/** current size of the memory per array */
	int fBlockSize;
	/*@}*/
};

/*************************************************************************
 * Implementation
 *************************************************************************/

/* constructor */
template <class TYPE>
MemoryGroupT<TYPE>::MemoryGroupT(int headroom, bool pool_memory):
	fHeadRoom(headroom),
	fPoolMemory(pool_memory),
	fBlockSize(0)
{
	/* error check */
	if (fHeadRoom < 0) ExceptionT::GeneralFail();

	/* pooled memory */
	if (fPoolMemory) fData.Append(NULL);
}

/* destructor */
template <class TYPE>
MemoryGroupT<TYPE>::~MemoryGroupT(void)
{
	for (int i = 0; i < fData.Length(); i++)
		delete[] fData[i];
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
	if (fHeadRoom < 0) ExceptionT::GeneralFail();
}

/* add array to list of managed */
template <class TYPE>
int MemoryGroupT<TYPE>::Register(ArrayT<TYPE>& array)
{
	/* only until memory is allocated */
	if (fBlockSize > 0 && fPoolMemory)
		ExceptionT::GeneralFail("MemoryGroupT<TYPE>::Register", 
			"no registration after dimensioning with pooled memory");

	/* add to list */
	fArrays.Append(&array);

	/* memory */
	if (!fPoolMemory)
	{
		if (fBlockSize > 0) /* allocate space for new array */ {
			TYPE* new_data = ArrayT<TYPE>::New(fBlockSize);
			fData.Append(new_data);
		}
		else
	 		fData.Append(NULL);
	}

	return fData.Length() - 1;
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
inline int MemoryGroupT<TYPE>::BlockSize(void) const { return fBlockSize; }

/* return a pointer to the specified block */
template <class TYPE>
TYPE* MemoryGroupT<TYPE>::BlockPointer(int block) const
{
	if (fPoolMemory) {	
		if (block < 0 || block >= fArrays.Length()) ExceptionT::OutOfRange();	
		return fData[0] + fBlockSize*block;
	}
	else
		return fData[block];
}

/* memory (re-) allocation */
template <class TYPE>
void MemoryGroupT<TYPE>::SetBlockSize(int newblocksize, bool copy_in)
{
	/* take the smaller */
	int copysize = (fBlockSize < newblocksize) ? fBlockSize : newblocksize;

	/* new memory (with extra space) */
	newblocksize += newblocksize*fHeadRoom/100;

	/* memory for all arrays pooled */
	if (fPoolMemory)
	{
		/* allocate new array */
		TYPE* newdata = ArrayT<TYPE>::New(fArrays.Length()*newblocksize);

		/* copy data in */
		if (copy_in) {
			TYPE* pold = fData[0];
			TYPE* pnew = newdata;
			for (int i = 0; i < fArrays.Length(); i++) {
				memcpy(pnew, pold, sizeof(TYPE)*copysize);
			
				pold += fBlockSize;
				pnew += newblocksize;
			}
		}
	
		/* reset grouped pointer */
		delete[] fData[0];
		fData[0] = newdata;
		fBlockSize = newblocksize;
	}
	else
	{
		fBlockSize = newblocksize;
		for (int i = 0; i < fArrays.Length(); i++) {
		
			/* allocate */
			TYPE* newdata = ArrayT<TYPE>::New(fBlockSize);
			
			/* copy data */
			if (copy_in) memcpy(newdata, fData[i], sizeof(TYPE)*copysize);
	
			/* reset grouped pointer */
			delete[] fData[i];
			fData[i] = newdata;
		}
	}
}

} // namespace Tahoe 
#endif /* _MEMORYGROUP_T_H_ */
