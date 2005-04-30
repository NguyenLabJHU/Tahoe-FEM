/* $Id: AutoFill2DT.h,v 1.16 2005-04-30 21:14:31 paklein Exp $ */
/* created: paklein (01/19/1999) */
#ifndef _AUTO_ARRAY2D_T_H_
#define _AUTO_ARRAY2D_T_H_

/* direct members */
#include "iArrayT.h"

namespace Tahoe {

/** class to facilitate dynamically building a two-dimensional
 * data structure where every row can have its own length.
 * Data is stored in chunks, which may be as small as a single
 * row, or as large as the whole array, i.e., a single chunk.
 * Values are added to the rows of the structure with the
 * methods AutoFill2DT::Append and AutoFill2DT::AppendUnique.
 * Currently, the structure just grows in memory.
 * \note one possible future development would be to have the
 *       array adaptively break into chunks based on the
 *       distribution of row lengths.
 */
template <class TYPE>
class AutoFill2DT
{
public:

	/** \name constructors */
	/*@{*/
	/** default constructor. Requires a subsequent call to AutoFill2DT::SetSize
	 * before any values can be added to the array. */
	AutoFill2DT(void);
	AutoFill2DT(int majordim, int numchunks, int headroom, int maxminordim);
//	AutoFill2DT(int majordim, int headroom); //add later
//	AutoFill2DT(const AutoFill2DT<TYPE>& source);
	/*@}*/

	/** destructor */
	~AutoFill2DT(void);

	/** \name dimensions */
	/*@{*/
	int MajorDim(void) const;
	int HeadRoom(void) const;
	int MinorDim(int row) const;
	int MaxMinorDim(void) const;
	int MinMinorDim(void) const;	
	int LogicalSize(void) const; // size of data (<= total memory allocation)
	int MinorDimCount(int dim) const; // returns number of rows with the specified size
	
	int NumChunks(void) const { return fChunks.Length(); };
	/*@}*/

	/** re-sizing */
	/*@{*/
	/** set the major dimension and number of chunks. Erases previous contents
	 * of the array. */
	void SetSize(int majordim, int numchunks);
	void SetMaxMinorDim(int maxminordim, bool make_headroom = false);
	void SetHeadRoom(int headroom);
	
	/** free memory */
	void Free(void);
	/*@}*/
	
	/** \name accessors
	 * Because the array is stored in chunks, each access requires an integer
	 * division. Also, one cannot assume that the values from one row continue
	 * to the values on the next row since the rows may be in separate chunks. */
	/*@{*/
	/** single element accessor */
	TYPE& operator()(int row, int col);

	/** const single element accessor */
	const TYPE& operator()(int row, int col) const;

	/** row accessor */
	TYPE* operator()(int row);

	/** row accessor */
	const TYPE* operator()(int row) const;
	/*@}*/

	/** \name (re-)set logical size to 0 */
	/*@{*/
	/** reset all rows */
	void Reset(void);
	
	/** reset just the specified row */
	void Reset(int row);
	/*@}*/

	/** copy of source - dimensions must match */
//	void Copy(const AutoFill2DT<TYPE>& source);

	/** returns the index of the first occurence of value in the given row or -1
	 * if not present - NOTE: requires "==" */
	int PositionInRow(int row, const TYPE& value) const;
	
	/** \name adding values to rows */
	/*@{*/
	void Append(int row, const TYPE& value);
	void Append(int row, const ArrayT<TYPE>& source);

	/** insert value into the given position within the row. The column
	 * must be within the length of the row or 1 beyond the end to append. */
	void Insert(int row, const TYPE& value, int col);
	/*@}*/

	/** \name add values only if not already present */
	/*@{*/
	/** append one value. Returns 1 if the value is added or 0 if it not */	
	int AppendUnique(int row, const TYPE& value);

	/** append an array of values. Returns the number of values that were appended */	
	int AppendUnique(int row, const ArrayT<TYPE>& source);
	/*@}*/	
	
private:

	/** dimension all chunks to the same minor dim */
	void Dimension(int minordim, int headroom);

	/** resize data and copy in previous values. Minor dimension of the chunk grows by 
	 * at least 1 with every call */
	void Dimension(int chunk, int maxminordim, int headroom);
	
	/** compute chunk dimensions */
	void ChunkDimensions(int major_dim, int num_chunks, int& chunk_major_dim, int& last_major_dim);

	/** map the global row number to a chunk and chunk row offset */
	void Chunk(int row, int& chunk, int& chunk_row) const;

	/** \name disallowed */
	/*@{*/
	AutoFill2DT(const AutoFill2DT<TYPE>& source);
	AutoFill2DT& operator=(const AutoFill2DT<TYPE>& rhs);
	/*@}*/

private:

	/** \name dimensions */
	/*@{*/
	int fMajorDim;
	int fHeadRoom; /**< percent over-allocation */
	int fMaxMinorDim;
	/*@}*/

	/** number of values per row */
	iArrayT fCounts;

	/** \name chunk dimensions */
	/*@{*/
	/** minor dimension of every chunk */
	ArrayT<int> fChunkMinorDim;

	/** number of rows in each chunk except the last. Needed because the number of rows
	 * is probably not a multiple of the number of chunks */
	int fChunkMajorDim;

	/** number of rows in last chunk. Needed because the number of rows
	 * is probably not a multiple of the number of chunks */
	int fLastChunkMajorDim;
	/*@}*/

	/** the data as array of pointers to chunks */
	ArrayT<TYPE*> fChunks;	
};

/*************************************************************************
 * Implementation
 *************************************************************************/

/* map the global row number to a chunk and chunk row offset */
template <class TYPE>
inline void AutoFill2DT<TYPE>::Chunk(int row, int& chunk, int& chunk_row) const
{
	chunk = row/fChunkMajorDim;
	chunk = (chunk == fChunkMinorDim.Length()) ? chunk-1 : chunk; /* last chunk catches overflow */
	chunk_row = row - chunk*fChunkMajorDim;
}

template <class TYPE>
AutoFill2DT<TYPE>::AutoFill2DT(void):
	fMajorDim(0),
	fHeadRoom(0),
	fMaxMinorDim(0),
	fChunkMajorDim(0),
	fLastChunkMajorDim(0)
{

}

template <class TYPE>
AutoFill2DT<TYPE>::AutoFill2DT(int majordim, int numchunks, int headroom, int maxminordim):
	fMajorDim(majordim),
	fHeadRoom(headroom),
	fMaxMinorDim(0),
	fCounts(fMajorDim),
	fChunkMinorDim(numchunks),
	fChunkMajorDim(0),
	fLastChunkMajorDim(0),
	fChunks(numchunks)
{
	/* check */
	if (fMajorDim < 0 ||
		fMaxMinorDim < 0 ||
	    fHeadRoom < 0) ExceptionT::GeneralFail("AutoFill2DT<TYPE>::AutoFill2DT");	

	/* set major dimensions of each chunk */
	ChunkDimensions(fMajorDim, numchunks, fChunkMajorDim, fLastChunkMajorDim);

	/* initialize arrays */
	fCounts = 0;
	fChunkMinorDim = 0;
	fChunks = NULL;

	/* initialize memory */
	Dimension(maxminordim, 0);
}

template <class TYPE>
AutoFill2DT<TYPE>::~AutoFill2DT(void)
{
	/* free memory */
	for (int i = 0; i < fChunks.Length(); i++)
		delete[] fChunks[i];
}

/* dimensions */
template <class TYPE>
inline int AutoFill2DT<TYPE>::MajorDim(void) const { return fMajorDim; }

template <class TYPE>
inline int AutoFill2DT<TYPE>::HeadRoom(void) const { return fHeadRoom; }

template <class TYPE>
inline int AutoFill2DT<TYPE>::MinorDim(int row) const
{
#if __option(extended_errorcheck)
	/* range check */
	if (row < 0 || row >= fMajorDim) ExceptionT::OutOfRange();
#endif
	return fCounts[row];
}

template <class TYPE>
int AutoFill2DT<TYPE>::MaxMinorDim(void) const { return fMaxMinorDim; }

template <class TYPE>
inline int AutoFill2DT<TYPE>::MinMinorDim(void) const { return fCounts.Min(); }

template <class TYPE>
int AutoFill2DT<TYPE>::LogicalSize(void) const { return fCounts.Sum(); }

template <class TYPE>
inline int AutoFill2DT<TYPE>::MinorDimCount(int dim) const { return fCounts.Count(dim); }

/* set the major dimension and number of chunks */
template <class TYPE>
void AutoFill2DT<TYPE>::SetSize(int majordim, int numchunks)
{
	/* free memory */
	for (int i = 0; i < fChunks.Length(); i++)
		delete[] fChunks[i];

	fMajorDim = majordim;
	fCounts.Dimension(fMajorDim);
	fCounts = 0;
	fChunks.Dimension(numchunks);
	fChunks = NULL;
	fChunkMinorDim.Dimension(numchunks);
	fChunkMinorDim = 0;

	/* set major dimensions of each chunk */
	fChunkMajorDim = 0;
	fLastChunkMajorDim = 0;
	ChunkDimensions(fMajorDim, numchunks, fChunkMajorDim, fLastChunkMajorDim);
}

/* re-sizing */
template <class TYPE>
inline void AutoFill2DT<TYPE>::SetMaxMinorDim(int maxminordim, bool make_headroom)
{
	//TEMP - size can only grow for now
	if (maxminordim < fMaxMinorDim) ExceptionT::GeneralFail();

	/* set memory */
	Dimension(maxminordim, (make_headroom == 1) ? fHeadRoom : 0);
}

template <class TYPE>
inline void AutoFill2DT<TYPE>::SetHeadRoom(int headroom)
{
	if (headroom < 0) ExceptionT::GeneralFail();
	fHeadRoom = headroom;
}

/* free memory */
template <class TYPE>
void AutoFill2DT<TYPE>::Free(void)
{
	fCounts = 0;
	fMaxMinorDim = 0;
	fChunkMinorDim = 0;
	for (int i = 0; i < fChunks.Length(); i++)
	{
		delete[] fChunks[i];
		fChunks[i] = NULL;
	}
}

/* accessors */
template <class TYPE>
inline const TYPE& AutoFill2DT<TYPE>::operator()(int row, int col) const
{
#if __option(extended_errorcheck)
	/* checks */
	if (row < 0 || row >= fMajorDim) ExceptionT::OutOfRange();
	if (col < 0 || col >= fCounts[row]) ExceptionT::OutOfRange();
#endif

	/* resolve chunk */
	int chunk = 0;
	int chunk_row = 0;
	Chunk(row, chunk, chunk_row);

	/* map to single element */
	return *(fChunks[chunk] + chunk_row*fChunkMinorDim[chunk] + col);
}

template <class TYPE>
inline TYPE& AutoFill2DT<TYPE>::operator()(int row, int col)
{
#if __option(extended_errorcheck)
	/* checks */
	if (row < 0 || row >= fMajorDim) ExceptionT::OutOfRange();
	if (col < 0 || col >= fCounts[row]) ExceptionT::OutOfRange();
#endif

	/* resolve chunk */
	int chunk = 0;
	int chunk_row = 0;
	Chunk(row, chunk, chunk_row);

	/* map to single element */
	return *(fChunks[chunk] + chunk_row*fChunkMinorDim[chunk] + col);
}

template <class TYPE>
inline const TYPE* AutoFill2DT<TYPE>::operator()(int row) const
{
#if __option(extended_errorcheck)
	/* checks */
	if (row < 0 || row >= fMajorDim) ExceptionT::OutOfRange();
#endif

	/* resolve chunk */
	int chunk = 0;
	int chunk_row = 0;
	Chunk(row, chunk, chunk_row);
	
	/* map to row pointer */
	return fChunks[chunk] + chunk_row*fChunkMinorDim[chunk];
}

template <class TYPE>
inline TYPE* AutoFill2DT<TYPE>::operator()(int row)
{
#if __option(extended_errorcheck)
	/* checks */
	if (row < 0 || row >= fMajorDim) ExceptionT::OutOfRange();
#endif

	/* resolve chunk */
	int chunk = 0;
	int chunk_row = 0;
	Chunk(row, chunk, chunk_row);
	
	/* map to row pointer */
	return fChunks[chunk] + chunk_row*fChunkMinorDim[chunk];
}

/* set logical size to 0 */
template <class TYPE>
inline void AutoFill2DT<TYPE>::Reset(void)
{
	memset(fCounts.Pointer(), 0, sizeof(int)*fMajorDim);
}

template <class TYPE>
inline void AutoFill2DT<TYPE>::Reset(int row)
{
#if __option(extended_errorcheck)
	/* checks */
	if (row < 0 || row >= fMajorDim) ExceptionT::OutOfRange();
#endif
	fCounts[row] = 0;
}

#if 0
/* assignment operator - dimensions must match exactly */
template <class TYPE>
void AutoFill2DT<TYPE>::Copy(const AutoFill2DT<TYPE>& source)
{
	/* no copies to self */
	if (this != &source)
	{
		/* dimensions must match */
		if (fMajorDim    != source.fMajorDim ||
		    fMaxMinorDim != source.fMaxMinorDim) throw ExceptionT::kSizeMismatch;
		
		/* copy counts */
		memcpy(fCounts, source.fCounts, sizeof(int)*fMajorDim);
		
		/* copy data (whole block) */
		memcpy(fArray, source.fArray, sizeof(TYPE)*fMajorDim*fMaxMinorDim);
	}
}
#endif

/* adding values to rows */
template <class TYPE>
inline 
void AutoFill2DT<TYPE>::Append(int row, const TYPE& value)
{
#if __option(extended_errorcheck)
	if (row < 0 || row >= fMajorDim) ExceptionT::OutOfRange();
#endif

	int& count = fCounts[row];
	
	/* resolve chunk */
	int chunk, chunk_row;
	Chunk(row, chunk, chunk_row);
	
	/* need more memory */
	int& minor_dim = fChunkMinorDim[chunk];
	if (count == minor_dim) Dimension(chunk, minor_dim + 1, fHeadRoom);

	*(fChunks[chunk] + chunk_row*minor_dim + count) = value;
	count++;	
}

/* returns the index of the first occurence of value in the given row or -1 
 * if not present - NOTE: requires "==" */
template <class TYPE>
int  AutoFill2DT<TYPE>::PositionInRow(int row, const TYPE& value) const
{
	const TYPE* prow = (*this)(row);
	int length = MinorDim(row);
	for (int i = 0; i < length; i++)
		if (*prow++ == value)
			return i; /* found */

	/* not found */
	return -1;
}

template <class TYPE>
void AutoFill2DT<TYPE>::Append(int row, const ArrayT<TYPE>& source)
{
#if __option(extended_errorcheck)
	if (row < 0 || row >= fMajorDim) ExceptionT::OutOfRange();
#endif

	int length = source.Length();
	if (length > 0)
	{
		int& count = fCounts[row];

		/* resolve chunk */
		int chunk, chunk_row;
		Chunk(row, chunk, chunk_row);
		
		/* need more memory */
		int& minor_dim = fChunkMinorDim[chunk];
		if (count + length >= minor_dim) Dimension(chunk, minor_dim + length, fHeadRoom);
		
		/* copy data */
		TYPE* start = fChunks[chunk] + chunk_row*minor_dim + count;
		memcpy(start, source.Pointer(), sizeof(TYPE)*length);
		count += length;
	}
}

/* sorted append---inserts element at ith position */
template <class TYPE>
void AutoFill2DT<TYPE>::Insert(int row, const TYPE& value, int col)
{
	/* check range */
	int last_row = MinorDim(row);
	if (col < 0 || col > last_row) ExceptionT::OutOfRange();

	/* expand the row */
	Append(row, value);
	
	/* shift elements */
	if (row != last_row) {
	
		/* move data */
		TYPE* prow = (*this)(row);
		TYPE* from = prow + col;
		TYPE*   to = from + 1;
		int length = last_row - (col + 1);
		ArrayT<TYPE>::MemMove(to, from, length);
	
		/* insert value at position */
		prow[col] = value;
	}	
}

/* add only if not already in the list - returns 1 if */
template <class TYPE>
int AutoFill2DT<TYPE>::AppendUnique(int row, const TYPE& value)
{
#if __option(extended_errorcheck)
	if (row < 0 || row >= fMajorDim) throw ExceptionT::kOutOfRange;
#endif

	/* chunk information */
	int chunk, chunk_row;
	Chunk(row, chunk, chunk_row);

	/* scan logical size for duplicates */
	TYPE* pthis = fChunks[chunk] + chunk_row*fChunkMinorDim[chunk];
	int count = fCounts[row];
	for (int i = 0; i < count; i++)
		if (*pthis++ == value)
			return 0;
			
	/* append value on fall through */
	Append(row, value);			
	return 1;
}

template <class TYPE>
inline int AutoFill2DT<TYPE>::AppendUnique(int row, const ArrayT<TYPE>& source)
{	
	const TYPE* psrc = source.Pointer();
	int length = source.Length();
	int count = 0;
	for (int i = 0; i < length; i++)
		count += AppendUnique(row, *psrc++);
	
	return count;
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* dimension all chunks */
template <class TYPE>
inline void AutoFill2DT<TYPE>::Dimension(int minordim, int headroom)
{
	for (int i = 0; i < fChunks.Length(); i++)
		Dimension(i, minordim, headroom);
}

template <class TYPE>
void AutoFill2DT<TYPE>::Dimension(int chunk, int maxminordim, int headroom)
{
	/* determine new memory size */
	int memsize;
	if (headroom > 0)
		memsize  = int(maxminordim*(1.0 + double(headroom)/100.0));
	else
		memsize = maxminordim;
	
	/* row in the chunk */
	int major_dim = (chunk == fChunks.Length() - 1) ? fLastChunkMajorDim : fChunkMajorDim;

	/* allocate new array */
	TYPE* newfArray = ArrayT<TYPE>::New(memsize*major_dim);
	
	/* copy current data */
	TYPE* pdata = fChunks[chunk];
	if (pdata)
	{
		int chunk_minor_dim = fChunkMinorDim[chunk];
		int* pcount = fCounts.Pointer(chunk*fChunkMajorDim);
		TYPE* pnew  = newfArray;
		for (int i = 0; i < major_dim; i++)
		{
			int dim = *pcount++;
			if (dim > 0) memcpy(pnew, pdata, sizeof(TYPE)*dim);
			
			pdata += chunk_minor_dim;
			pnew  += memsize;
		}
	}

	/* swap memory */
	delete[] fChunks[chunk];
	fChunks[chunk] = newfArray;

	/* save dimensions */
	fChunkMinorDim[chunk] = memsize;
}

/* compute chunk dimensions */
template <class TYPE>
void AutoFill2DT<TYPE>::ChunkDimensions(int major_dim, int num_chunks, int& chunk_major_dim, int& last_major_dim)
{
	/* checks */
	if (num_chunks < 0) 
		ExceptionT::SizeMismatch("AutoFill2DT<TYPE>::ChunkDimensions", "bad number of chunks %d", num_chunks);

	/* integer division truncates */
	chunk_major_dim = major_dim/num_chunks;

	/* collect remainder */
	last_major_dim = major_dim - (num_chunks - 1)*chunk_major_dim;
}

} // namespace Tahoe 
#endif /* _AUTO_ARRAY2D_T_H_ */
