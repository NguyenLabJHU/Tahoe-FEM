/* $Id: RaggedArray2DT.h,v 1.21 2004-05-06 18:47:37 cjkimme Exp $ */
/* created: paklein (09/10/1998) */
#ifndef _RAGGED_ARRAY_2D_T_H_
#define _RAGGED_ARRAY_2D_T_H_

/* language support */
#include <iostream.h>

/* direct members */
#include "AutoArrayT.h"
#include "nArray2DT.h"
#include "AutoFill2DT.h"
#include "RowAutoFill2DT.h"
#include "LinkedListT.h"

namespace Tahoe {

/** 2D array with arbitrary row lengths. operator()'s provided for 
 * data retrieval. Memory is managed using AutoArrayT's to allow
 * the array to be reconfigured without thrashing memory too much.
 */
template <class TYPE>
class RaggedArray2DT
{
public:

	/** \name constructors */
	/*@{*/
	/** construct empty array */
	RaggedArray2DT(int headroom = 0);
	
	/** construct array with the same dimensions for every row. */
	RaggedArray2DT(int majordim, int minordim, int headroom = 0, int blocksize = 1);

	/** copy king */
	RaggedArray2DT(const RaggedArray2DT& source);
	/*@}*/

	/** set the array size with fixed row dimensions. No change occurs if the array
	 * is already the specified size. The previous contents of the array is
	 * not preserved. */
	void Dimension(int majordim, int minordim);

	/** \deprecated replaced by RaggedArray2DT::Dimension on 02/13/2002 */
	void Allocate(int majordim, int minordim) { Dimension(majordim, minordim); };

	/** return the dimension of the data block */		
	int Length(void) const;
	
	/** return the number of rows */
	int MajorDim(void) const;
	
	/** return the length of the specified row */
	int MinorDim(int row) const;
	
	/** retrieve the dimensions of all the rows */
	void MinorDim(ArrayT<int>& minordim) const;

	/** return the length of the longest row */
	int MaxMinorDim(void) const;

	/** return the smallest minor dimension.
	 * \param floor values smaller than floor are ignored when determining the smallest
	 *        row length */
	int MinMinorDim(int floor) const;

	/** return the smallest minor dimension and its location.
	 * \param dex returns with the location of the longest row
	 * \param floor values smaller than floor are ignored when determining the smallest
	 *        row length */
	int MinMinorDim(int& dex, int floor) const; // smallest minor dimension ( > floor)

	/** \name configuration methods */
	/*@{*/
	/** configuration the array using a list of row dimensions
	 * \param rowcounts array containing the length of each row
	 * \param blocksize allocated space is rowcounts[i]*blocksize
	 *        for each row. */
	void Configure(const ArrayT<int>& rowcounts, int blocksize = 1);

	/** configuration to be the scaled shape of the source.
	 * Only the shape of the array is duplicated. No data is copied from the source array
	 * \param source array providing its shape to this one.
	 * \param blocksize allocated space is rowcounts[i]*blocksize
	 *        for each row. */
	void Configure(const RaggedArray2DT& source, int blocksize = 1);
	
	/** configuration to be the scaled shape of the source.
	 * Only the shape of the linked list is duplicated. No data is copied from the source array
	 * \param source linkedlist providing its shape to this one.
	 * \param blocksize allocated space is rowcounts[i]*blocksize
	 *        for each row. */
	void Configure(const ArrayT<LinkedListT<TYPE> >& source, int blocksize = 1);
	/*@}*/

	/** shallow copy */
	void Alias(const RaggedArray2DT& source);

	/** shallow copy from nArray2DT */
	void Alias(const nArray2DT<TYPE>& source);

	/** \name assignment operators */
	/*@{*/
	/** assigment operator from another RaggedArray2DT */
	RaggedArray2DT<TYPE>& operator=(const RaggedArray2DT& source);

	/** set entire array to the same value */
	RaggedArray2DT<TYPE>& operator=(const TYPE& value);
	/*@}*/

	/** configure and construct using a AutoFill2DT */
	void Copy(const AutoFill2DT<TYPE>& source);

	/** configure and construct using a RowAutoFill2DT */
	void Copy(const RowAutoFill2DT<TYPE>& source);

	/** configure and construct using raw data.
	 * \param rowcounts array of number of entries per row
	 * \param data array of pointers to the row data */
	void Copy(const ArrayT<int>& rowcounts, const ArrayT<TYPE*>& data);

	/** configure and construct using a AutoFill2DT. Empty
	 * rows are removed. */
	void CopyCompressed(const AutoFill2DT<TYPE>& source);

	/** generate adjacency offset vector
	 * \param offsets returns with the offset of the start of every
	 *        of every row from the base address of the data array. */
	void GenerateOffsetVector(ArrayT<int>& offsets) const;

	/** write data to a row of the array. Dimension of the source
	 * must match the allocated size of the specified row.
	 * \param row destination row
	 * \param array source of data */
	void SetRow(int row, const ArrayT<TYPE>& array);

	/** write data to a row of the array.
	 * \param row destination row
	 * \param array source of data. The length of this array is
	 *        assumed to be at least as long as the size of the
	 *        specified row. */
	void SetRow(int row, const TYPE* array);

	/** write value to a row of the array */
	void SetRow(int row, const TYPE& value);

	/* write rows with numbers to the output stream */
	void WriteNumbered(ostream& out) const;

	/** shallow copy of a row.
	 * \param row row number to alias
	 * \param array array to alias to the row data */
	void RowAlias(int row, ArrayT<TYPE>& rowdata) const;

	/** return a pointer to first element in the specified row */
	TYPE* operator()(int row);

	/** return a pointer to first element in the specified row */
	const TYPE* operator()(int row) const;

	/** return a pointer to first element in the specified row */
	TYPE& operator()(int row, int col);

	/** return a pointer to first element in the specified row */
	const TYPE& operator()(int row, int col) const;
	
	/** \name return a pointer to the first element in the data array */
	/*@{*/
	TYPE* Pointer(void);
	const TYPE* Pointer(void) const;
	const ArrayT<int>& Data(void) const { return fData; };
	/*@}*/

	/** free memory */
	void Free(void);

	/** configure from stream. Read the data for the array from the stream.
	 * This data includes information about the internal structure of the
	 * array. This is the complement to RaggedArray2DT::Write. */
	void Read(istream& in);

	/** write to stream. Write the data for the array to the stream.
	 * This data includes information about the internal structure of the
	 * array. This is the complement to RaggedArray2DT::Read. */
	void Write(ostream& out) const;
	
	/** binary input of the data array. Read only the data portion of
	 * the array from the input stream. The structure of the array
	 * must be set previously. This is the complement to 
	 * RaggedArray2DT::WriteDataBinary. */
	void ReadDataBinary(istream& in) { fData.ReadBinary(in); };

	/** binary output of the data array. Write only the data portion of
	 * the array from the input stream. Information about the structure 
	 * of the array is not included in the output. Hence, arrays that read
	 * this data must be configured previously. This is the complement 
	 * to RaggedArray2DT::ReadDataBinary. */
	void WriteDataBinary(ostream& out) const { fData.WriteBinary(out); };

	/** text input of the data array. Read only the data portion of
	 * the array from the input stream. The structure of the array
	 * must be set previously. This is the complement to 
	 * RaggedArray2DT::WriteData. */
	void ReadData(istream& in);

	/** text output of the data array. Write only the data portion of
	 * the array from the input stream. Information about the structure 
	 * of the array is not included in the output. Hence, arrays that read
	 * this data must be configured previously. This is the complement 
	 * to RaggedArray2DT::ReadDataBinary. */
	void WriteData(ostream& out) const;

private:

	/** setting pointers for equal row sizes */
	void SetEvenPointers(int minordim);
	  	
protected:

	/** number of rows */
	int fMajorDim;

	/** smallest row dimension */
	int fMinMinorDim;

	/** largest row dimension */
	int fMaxMinorDim;

	/** pointers to the data array. Length is majordim + 1 */
	AutoArrayT<TYPE*> fPtrs;

	/** data array */
	AutoArrayT<TYPE>  fData;
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructors */
template <class TYPE>
inline RaggedArray2DT<TYPE>::RaggedArray2DT(int headroom):
	fMajorDim(0),
	fMinMinorDim(0),
	fMaxMinorDim(0),
	fPtrs(headroom),
	fData(headroom)
{
	/* length is majordim + 1 */
	fPtrs.Append(NULL);
}

template <class TYPE>
RaggedArray2DT<TYPE>::RaggedArray2DT(int majordim, int minordim, int headroom, 
	int blocksize):
	fPtrs(headroom),
	fData(headroom)
{
	/* configure */
	Dimension(majordim, minordim*blocksize);
}

/** copy king */
template <class TYPE>
RaggedArray2DT<TYPE>::RaggedArray2DT(const RaggedArray2DT& source):
	fMajorDim(0),
	fMinMinorDim(0),
	fMaxMinorDim(0),
	fPtrs(0),
	fData(0)
{
	operator=(source);
}
	
/* allocate array with fixed row dimensions */
template <class TYPE>
void RaggedArray2DT<TYPE>::Dimension(int majordim, int minordim)
{
	/* set dimensions */
	fMajorDim = majordim;
	fMinMinorDim = fMaxMinorDim = minordim;

	/* allocate space */
	fPtrs.Dimension(fMajorDim + 1);
	fData.Dimension(fMajorDim*minordim);

	/* set pointers */
	SetEvenPointers(minordim);
}

/* dimensions */
template <class TYPE>
inline int RaggedArray2DT<TYPE>::Length(void) const
{
	return fData.Length();
}

template <class TYPE>
inline int RaggedArray2DT<TYPE>::MajorDim(void) const
{
	return fMajorDim;
}

template <class TYPE>
inline int RaggedArray2DT<TYPE>::MinorDim(int row) const
{
#if __option(extended_errorcheck)
	if (row < 0 || row >= fMajorDim) ExceptionT::OutOfRange();
#endif

	TYPE** p = (TYPE**) fPtrs.Pointer() + row;
	return *(p + 1) - *p;
}

template <class TYPE>
inline void RaggedArray2DT<TYPE>::MinorDim(ArrayT<int>& minordim) const
{
#if __option(extended_errorcheck)
	if (minordim.Length() != fMajorDim) ExceptionT::SizeMismatch();
#endif

	TYPE** p = (TYPE**) fPtrs.Pointer();
	int* pdim = minordim.Pointer();
	for (int i = 0; i < fMajorDim; i++)
	{
		*pdim++ = *(p + 1) - *p;
		p++;
	}
}

template <class TYPE>
inline int RaggedArray2DT<TYPE>::MaxMinorDim(void) const
{
	return fMaxMinorDim;
}

/* smallest minor dimension ( > lower_limit) */
template <class TYPE>
int RaggedArray2DT<TYPE>::MinMinorDim(int floor) const
{
	int min = fMaxMinorDim;
	int len = MajorDim();
	for (int i = 0; i < len; i++)
	{
		int dim = MinorDim(i);
		min = (dim > floor && dim < min) ? dim : min;
	}

	return min;
}

template <class TYPE>
int RaggedArray2DT<TYPE>::MinMinorDim(int& dex, int floor) const
{
	dex = 0;
	int min = fMaxMinorDim;
	int len = MajorDim();
	for (int i = 0; i < len; i++)
	{
		int dim = MinorDim(i);
		if (dim > floor && dim < min)
		{
			dex = i;
			min = dim;
		}
	}

	return min;
}

/* configuration functions */
template <class TYPE>
void RaggedArray2DT<TYPE>::Configure(const ArrayT<int>& rowcounts, int blocksize)
{
	/* number of rows */
	fMajorDim = rowcounts.Length();

	/* count total entries */
	int size = 0;
	const int* pcount = rowcounts.Pointer();
	for (int j = 0; j < fMajorDim; j++)
		size += *pcount++;
		
	/* allocate memory */
	fPtrs.Dimension(fMajorDim + 1);
	fData.Dimension(size*blocksize);		

	/* set pointers */
	pcount = rowcounts.Pointer();
	fMinMinorDim = (rowcounts.Length() > 0) ? rowcounts[0] : 0;
	fMaxMinorDim = 0;
	
	TYPE*  pdata  = fData.Pointer();
	TYPE** p      = fPtrs.Pointer();	
	for (int i = 0; i < fMajorDim; i++)
	{
		/* track dimensions */
		if (*pcount < fMinMinorDim)
			fMinMinorDim = *pcount;
		else if (*pcount > fMaxMinorDim)
			fMaxMinorDim = *pcount;
	
		*p++   = pdata;
		pdata += blocksize*(*pcount++);
	}
	fPtrs[fMajorDim] = pdata;

	/* adjust for block size */
	fMinMinorDim *= blocksize;
	fMaxMinorDim *= blocksize;
}

/* configuration to be the scaled shape of the source */
template <class TYPE>
void RaggedArray2DT<TYPE>::Configure(const RaggedArray2DT& source, int blocksize)
{
	/* (scaled) dimensions */
	fMajorDim    = source.fMajorDim;
	fMinMinorDim = source.fMinMinorDim*blocksize;
	fMaxMinorDim = source.fMaxMinorDim*blocksize;

	/* allocate memory */
	fPtrs.Dimension(fMajorDim + 1);
	fData.Dimension(source.fData.Length()*blocksize);		

	/* set pointers */
	TYPE*  pdata = fData.Pointer();
	TYPE** p     = fPtrs.Pointer();	
	for (int i = 0; i < fMajorDim; i++)
	{	
		*p++   = pdata;
		pdata += blocksize*source.MinorDim(i);
	}
	fPtrs[fMajorDim] = pdata;
}

/* configuration to be the scaled shape of the source */
template <class TYPE>
void RaggedArray2DT<TYPE>::Configure(const ArrayT<LinkedListT<TYPE> >& source, int blocksize)
{
	ArrayT<int> rowCounts(source.Length());
	
	for (int i = 0; i < source.Length(); i++)
		rowCounts[i] = source[i].Length();

	Configure(rowCounts, blocksize);
}

/* shallow copy/conversion */
template <class TYPE>
inline void RaggedArray2DT<TYPE>::Alias(const RaggedArray2DT& source)
{
	/* dimensions */
	fMajorDim = source.fMajorDim;
	fMinMinorDim = source.fMinMinorDim;
	fMaxMinorDim = source.fMaxMinorDim;
	
	/* data */
	fPtrs.Alias(source.fPtrs);
	fData.Alias(source.fData);
}

template <class TYPE>
void RaggedArray2DT<TYPE>::Alias(const nArray2DT<TYPE>& source)
{
	/* dimensions */
	fMajorDim    = source.MajorDim();
	fMinMinorDim = fMaxMinorDim = source.MinorDim;
	
	/* configure memory */
	fPtrs.Dimension(fMajorDim + 1);
	fData.Set(source.Length(), source.Pointer());
	
	/* set pointers */
	SetEvenPointers(source.MinorDim());
}

/* assigment operator */
template <class TYPE>
RaggedArray2DT<TYPE>& RaggedArray2DT<TYPE>::operator=(const RaggedArray2DT& source)
{
	/* copy parameters */
	fPtrs.SetHeadRoom(source.fPtrs.HeadRoom());
	fData.SetHeadRoom(source.fData.HeadRoom());

	/* quick copy for repeated, even-spaced copies */
	if (source.fMinMinorDim == source.fMaxMinorDim &&
	    fMinMinorDim == fMaxMinorDim &&
	    fMinMinorDim == source.fMinMinorDim &&
	    fMajorDim == source.fMajorDim)
	{
		/* copy data */
		fData = source.fData;
	}
	else
	{
		/* (re-)dimension */
		fMajorDim    = source.fMajorDim;
		fMaxMinorDim = source.fMaxMinorDim;
		fMinMinorDim = source.fMinMinorDim;
	
		/* data */
		fData = source.fData;
		
		/* allocate space for data pointers */
		fPtrs.Dimension(source.fPtrs.Length());
	
		/* quick set */
		if (fMinMinorDim == fMaxMinorDim)
			SetEvenPointers(fMinMinorDim);
		else
		{
			/* set pointers */
			TYPE* ptr = fData.Pointer();
			for (int i = 0; i < fMajorDim; i++)
			{
				fPtrs[i] = ptr;
				ptr += source.MinorDim(i);
			}
			fPtrs[fMajorDim] = ptr;
		}
	}
	return *this;
}

/* assigment operator */
template <class TYPE>
RaggedArray2DT<TYPE>& RaggedArray2DT<TYPE>::operator=(const TYPE& value)
{
	/* set array */
	fData = value;
	return *this;
}

template <class TYPE>
void RaggedArray2DT<TYPE>::Copy(const AutoFill2DT<TYPE>& source)
{
	/* total memory size */
	fMajorDim = source.MajorDim();		
	fPtrs.Dimension(fMajorDim + 1);
	fData.Dimension(source.LogicalSize());

	fMinMinorDim = (fMajorDim > 0) ? source.MinorDim(0) : 0;
	fMaxMinorDim = 0;
	TYPE** pptrs = fPtrs.Pointer();
	TYPE*  pdata = fData.Pointer();
	for (int i = 0; i < fMajorDim; i++)
	{
		/* copy data */
		int length = source.MinorDim(i);
		if (length > 0) memcpy(pdata, source(i), sizeof(TYPE)*length);

		/* track dimensions */
		if (length < fMinMinorDim)
			fMinMinorDim = length;
		else if (length > fMaxMinorDim)
			fMaxMinorDim = length;
		
		/* set pointer */
		*pptrs = pdata;			
		pdata += length;
		
		pptrs++;
	}
	
	/* set trailing pointer */
	*pptrs = pdata;
	
	/* memory check */
	if (pdata - fData.Pointer() != fData.Length())
		ExceptionT::GeneralFail("RaggedArray2DT<TYPE>::Copy", "memory partitioning error");
}

template <class TYPE>
void RaggedArray2DT<TYPE>::Copy(const RowAutoFill2DT<TYPE>& source)
{
	/* total memory size */
	fMajorDim = source.MajorDim();		
	fPtrs.Dimension(fMajorDim + 1);
	fData.Dimension(source.LogicalSize());

	fMinMinorDim = (fMajorDim > 0) ? source.MinorDim(0) : 0;
	fMaxMinorDim = 0;
	TYPE** pptrs = fPtrs.Pointer();
	TYPE*  pdata = fData.Pointer();
	for (int i = 0; i < fMajorDim; i++)
	{
		/* copy data */
		int length = source.MinorDim(i);
		if (length > 0) memcpy(pdata, source(i), sizeof(TYPE)*length);

		/* track dimensions */
		if (length < fMinMinorDim)
			fMinMinorDim = length;
		else if (length > fMaxMinorDim)
			fMaxMinorDim = length;
		
		/* set pointer */
		*pptrs = pdata;			
		pdata += length;
		
		pptrs++;
	}
	
	/* set trailing pointer */
	*pptrs = pdata;
	
	/* memory check */
	if (pdata - fData.Pointer() != fData.Length())
		ExceptionT::GeneralFail("RaggedArray2DT<TYPE>::Copy", "memory partitioning error");
}

template <class TYPE>
void RaggedArray2DT<TYPE>::Copy(const ArrayT<int>& rowcounts,
	const ArrayT<TYPE*>& data)
{
	/* allocate */
	Configure(rowcounts);

	/* copy data */
	for (int i = 0; i < rowcounts.Length(); i++)
	{
		int length = rowcounts[i];
		if (length > 0) memcpy(fPtrs[i], data[i], sizeof(TYPE)*length);
	}
}

template <class TYPE>
void RaggedArray2DT<TYPE>::CopyCompressed(const AutoFill2DT<TYPE>& source) // removes empty rows
{
	/* total memory size */
	fMajorDim = source.MajorDim() - source.MinorDimCount(0); // number of non-empty rows		
	fPtrs.Dimension(fMajorDim + 1); // have trailing pointer
	fData.Dimension(source.LogicalSize());

	fMaxMinorDim = 0;
	TYPE** pptrs = fPtrs.Pointer();
	TYPE*  pdata = fData.Pointer();
	for (int i = 0; i < source.MajorDim(); i++)
	{
		/* copy data */
		int length = source.MinorDim(i);
		if (length > 0)
		{
			memcpy(pdata, source(i), sizeof(TYPE)*length);
	
			/* set pointer */
			*pptrs = pdata;

			/* next */		
			pptrs++;
			pdata += length;
		}
		
		/* track dimensions */
		if (length < fMinMinorDim)
			fMinMinorDim = length;
		else if (length > fMaxMinorDim)
			fMaxMinorDim = length;
	}
	
	/* set trailing pointer */
	*pptrs = pdata;
	
	/* memory check */
	if (pdata - fData.Pointer() != fData.Length())
		ExceptionT::GeneralFail("RaggedArray2DT<TYPE>::CopyCompressed", "memory partitioning error");
}

/* generate adjacency offset vector */
template <class TYPE>
void RaggedArray2DT<TYPE>::GenerateOffsetVector(ArrayT<int>& offsets) const
{
	/* same length as pointers */
	offsets.Dimension(fPtrs.Length());

	/* pointer to base address */
	const TYPE* base = fData.Pointer();
	int length = offsets.Length();
	for (int i = 0; i < length; i++)
		offsets[i] = fPtrs[i] - base;
}

/* write data */
template <class TYPE>
inline void RaggedArray2DT<TYPE>::SetRow(int row, const TYPE* array)
{
	memcpy(fPtrs[row], array, sizeof(TYPE)*MinorDim(row));
}

template <class TYPE>
inline void RaggedArray2DT<TYPE>::SetRow(int row, const ArrayT<TYPE>& array)
{
#if __option(extended_errorcheck)
	if (array.Length() != MinorDim(row)) throw ExceptionT::kSizeMismatch;
#endif

	SetRow(row, array.Pointer());
}

template <class TYPE>
inline void RaggedArray2DT<TYPE>::SetRow(int row, const TYPE& value)
{
	TYPE* prow = fPtrs[row];
	int dim = MinorDim(row);
	for (int i = 0; i < dim; i++)
		*prow++ = value;
}			

/* write rows with numbers to the output stream */
template <class TYPE>
void RaggedArray2DT<TYPE>::WriteNumbered(ostream& out) const
{
	int width = OutputWidth(out, fPtrs[0]);
	
	for (int i = 0; i < fMajorDim; i++)
	{
		out << setw(kIntWidth) << i+1;
		
		int minordim = MinorDim(i);
		const TYPE* p = (*this)(i);
		for (int j = 0; j < minordim; j++)
			out << setw(width) << *p++;
		out << '\n';
	}
}

/* data retrieval */
template <class TYPE>
inline void RaggedArray2DT<TYPE>::RowAlias(int row,
	ArrayT<TYPE>& rowdata) const
{
	rowdata.Set(MinorDim(row),fPtrs[row]);
}

/* data retrieval */
template <class TYPE>
inline TYPE* RaggedArray2DT<TYPE>::operator()(int row)
{
	return fPtrs[row];
}

template <class TYPE>
inline const TYPE* RaggedArray2DT<TYPE>::operator()(int row) const
{
	return fPtrs[row];
}

/* return a pointer to first element in the specified row */
template <class TYPE>
inline TYPE& RaggedArray2DT<TYPE>::operator()(int row, int col)
{
	TYPE* prow = fPtrs[row];
	return prow[col];
}

template <class TYPE>
inline const TYPE& RaggedArray2DT<TYPE>::operator()(int row, int col) const
{
	const TYPE* prow = fPtrs[row];
	return prow[col];
}

template <class TYPE>
inline TYPE* RaggedArray2DT<TYPE>::Pointer(void)
{
	return fData.Pointer();
}

template <class TYPE>
const inline TYPE* RaggedArray2DT<TYPE>::Pointer(void) const
{
	return fData.Pointer();
}

/* free memory (if allocated) */
template <class TYPE>
inline void RaggedArray2DT<TYPE>::Free(void)
{
	fMajorDim = 0;
	fMinMinorDim = 0;
	fMaxMinorDim = 0;
	fPtrs.Free();
	fData.Free();
}

/* I/O */
template <class TYPE>
void RaggedArray2DT<TYPE>::Read(istream& in)
{
	in >> fMajorDim;

	ArrayT<int> counts(fMajorDim);
	for (int i = 0; i < fMajorDim; i++)
		in >> counts[i];
		
	Configure(counts);

	for (int j = 0; j < fData.Length(); j++)
		 in >> fData[j];
}

template <class TYPE>
void RaggedArray2DT<TYPE>::Write(ostream& out) const
{
	out << fMajorDim << '\n';
	
	for (int i = 0; i < fMajorDim; i++)
		out << MinorDim(i) << '\n';
		
	for (int j = 0; j < fData.Length(); j++)
		out << fData[j] << '\n';
}

template <class TYPE>
void RaggedArray2DT<TYPE>::ReadData(istream& in)
{ 
	/* read */
	for (int j = 0; j < fData.Length(); j++)
		 in >> fData[j];
}

template <class TYPE>
void RaggedArray2DT<TYPE>::WriteData(ostream& out) const
{
	const int wrap = 5;
	int count = 0;
	for (int i = 0; i < fData.Length(); i++)
	{
		/* wrap */
		if (count == wrap)
		{
			out << '\n';
			count = 0;
		}
		else if (count > 0)
			out << " ";
		count++;

		out << fData[i];		
	}
}

/*************************************************************************
* Private
*************************************************************************/

/* setting pointers for equal row sizes */
template <class TYPE>
void RaggedArray2DT<TYPE>::SetEvenPointers(int minordim)
{
	/* set pointers */
	TYPE*  pdata = fData.Pointer();
	TYPE** p     = fPtrs.Pointer();
	for (int i = 0; i <= fMajorDim; i++)
	{
		*p++   = pdata;
		pdata += minordim;
	}
}

} // namespace Tahoe 
#endif /* _RAGGED_ARRAY_2D_T_H_ */
