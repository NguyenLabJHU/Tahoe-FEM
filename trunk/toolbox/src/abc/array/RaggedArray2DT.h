/* $Id: RaggedArray2DT.h,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (09/10/1998)                                          */
/* 2D array with arbitrary "row" lengths. NO functions are provided       */
/* for data retrieval. derived classes should use RowAlias()              */
/* to retrieve data and return an appropriate array type.                 */

#ifndef _RAGGED_ARRAY_2D_T_H_
#define _RAGGED_ARRAY_2D_T_H_

/* language support */
#include <iostream.h>

/* direct members */
#include "ArrayT.h"
#include "nArray2DT.h"
#include "AutoFill2DT.h"

template <class TYPE>
class RaggedArray2DT
{
public:

	/* constructors */
	RaggedArray2DT(void);
	RaggedArray2DT(int majordim, int minordim, int blocksize = 1);
		// equal row sizes
		
	/* dimensions */
	int Length(void) const;      // length of the data block
	int MajorDim(void) const;
	int MinorDim(int row) const;
	void MinorDim(ArrayT<int>& minordim) const; // all minor dim's
	int MaxMinorDim(void) const;          // largest  minor dimension
	int MinMinorDim(int floor) const; // smallest minor dimension ( > floor)
	int MinMinorDim(int& dex, int floor) const; // smallest minor dimension ( > floor)

	/* configuration functions */
	void Configure(const ArrayT<int>& rowcounts, int blocksize = 1);

	/* shallow copy */
	void Alias(const RaggedArray2DT& source);
	void Alias(const nArray2DT<TYPE>& source); // equal row sizes

	/* assigment */
	RaggedArray2DT<TYPE>& operator=(const RaggedArray2DT& source);
	void Copy(const AutoFill2DT<TYPE>& source);
	void Copy(const ArrayT<int>& rowcounts, const ArrayT<TYPE*>& data);
	void CopyCompressed(const AutoFill2DT<TYPE>& source); // removes empty rows

	/* write data */
	void SetRow(int row, const ArrayT<TYPE>& array);
	void SetRow(int row, const TYPE* array);

	/* write rows with numbers to the output stream */
	void WriteNumbered(ostream& out) const;

	/* data retrieval */
	void RowAlias(int row, ArrayT<TYPE>& rowdata) const;
	TYPE* operator()(int row) const;
	TYPE* Pointer(void) const; // pointer to the first data element

	/* free memory (if allocated) */
	void Free(void);

	/* I/O */
	void Read(istream& in);
	void Write(ostream& out) const;

private:

	/* setting pointers for equal row sizes */
	void SetEvenPointers(int minordim);

	/* disallowed */
	RaggedArray2DT(const RaggedArray2DT& source);
	  	
protected:

	/* number of rows */
	int fMajorDim;
	int fMaxMinorDim;

	/* data pointers */
	ArrayT<TYPE*> fPtrs; // length is majordim + 1
	ArrayT<TYPE>  fData;
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructors */
template <class TYPE>
inline RaggedArray2DT<TYPE>::RaggedArray2DT(void):
	fMajorDim(0),
	fMaxMinorDim(0)
{

}

template <class TYPE>
RaggedArray2DT<TYPE>::RaggedArray2DT(int majordim, int minordim, int blocksize):
	fMajorDim(majordim),
	fMaxMinorDim(fMajorDim),
	fPtrs(fMajorDim + 1),
	fData(fMajorDim*minordim*blocksize)
{
	/* set pointers */
	SetEvenPointers(minordim*blocksize);
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
	if (row < 0 || row >= fMajorDim) throw eOutOfRange;
#endif

	TYPE** p = fPtrs.Pointer() + row;
	return *(p + 1) - *p;
}

template <class TYPE>
inline void RaggedArray2DT<TYPE>::MinorDim(ArrayT<int>& minordim) const
{
#if __option(extended_errorcheck)
	if (minordim.Length() != fMajorDim) throw eSizeMismatch;
#endif

	TYPE**  p = fPtrs.Pointer();
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
	int    size = 0;
	int* pcount = rowcounts.Pointer();
	for (int j = 0; j < fMajorDim; j++)
		size += *pcount++;
		
	/* allocate memory */
	fPtrs.Allocate(fMajorDim + 1);
	fData.Allocate(size*blocksize);		

	/* set pointers */
	pcount = rowcounts.Pointer();
	fMaxMinorDim = 0;
	
	TYPE*  pdata  = fData.Pointer();
	TYPE** p      = fPtrs.Pointer();	
	for (int i = 0; i < fMajorDim; i++)
	{
		/* keep max dimension */
		fMaxMinorDim = (*pcount > fMaxMinorDim) ? *pcount : fMaxMinorDim;
	
		*p++   = pdata;
		pdata += blocksize*(*pcount++);
	}
	fPtrs[fMajorDim] = pdata; // set outside to keep pcount in range

	/* adjust for block size */
	fMaxMinorDim *= blocksize;
}

/* shallow copy/conversion */
template <class TYPE>
inline void RaggedArray2DT<TYPE>::Alias(const RaggedArray2DT& source)
{
	/* dimensions */
	fMajorDim = source.fMajorDim;
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
	fMaxMinorDim = source.MinorDim;
	
	/* configure memory */
	fPtrs.Allocate(fMajorDim + 1);
	fData.Set(source.Length(), source.Pointer());
	
	/* set pointers */
	SetEvenPointers(source.MinorDim());
}

/* assigment operator */
template <class TYPE>
RaggedArray2DT<TYPE>& RaggedArray2DT<TYPE>::operator=(const RaggedArray2DT& source)
{
	/* dimension */
	fMajorDim    = source.fMajorDim;
	fMaxMinorDim = source.fMaxMinorDim;
	
	/* data */
	fData = source.fData;
	fPtrs = source.fPtrs;

	/* shift */
	TYPE* shift = fData.Pointer() - (source.fData).Pointer();
	TYPE** p    = fPtrs.Pointer();	
	for (int i = 0; i < fMajorDim; i++)
		*p++ += shift;
	
	return *this;
}


template <class TYPE>
void RaggedArray2DT<TYPE>::Copy(const AutoFill2DT<TYPE>& source)
{
	/* total memory size */
	fMajorDim = source.MajorDim();		
	fPtrs.Allocate(fMajorDim + 1);
	fData.Allocate(source.LogicalSize());

	fMaxMinorDim = 0;
	TYPE** pptrs = fPtrs.Pointer();
	TYPE*  pdata = fData.Pointer();
	for (int i = 0; i < fMajorDim; i++)
	{
		/* copy data */
		int length = source.MinorDim(i);
		if (length > 0)
		{
			memcpy(pdata, source(i), sizeof(TYPE)*length);
	
			/* track size */
			fMaxMinorDim = (length > fMaxMinorDim) ? length : fMaxMinorDim;
		}		
		
		/* set pointer */
		*pptrs = pdata;			
		pdata += length;
		
		pptrs++;
	}
	
	/* set trailing pointer */
	*pptrs = pdata;
	
	/* memory check */
	if (pdata - fData.Pointer() != fData.Length())
	{
		cout << "\n RaggedArray2DT<TYPE>::Copy: memory partitioning error" << endl;
		throw eGeneralFail;
	}
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
	fPtrs.Allocate(fMajorDim + 1); // have trailing pointer
	fData.Allocate(source.LogicalSize());

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

			/* track size */
			fMaxMinorDim = (length > fMaxMinorDim) ? length : fMaxMinorDim;
		
			pptrs++;
			pdata += length;
		}
	}
	
	/* set trailing pointer */
	*pptrs = pdata;
	
	/* memory check */
	if (pdata - fData.Pointer() != fData.Length())
	{
		cout << "\n RaggedArray2DT<TYPE>::CopyCompressed: memory partitioning error" << endl;
		throw eGeneralFail;
	}
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
	if (array.Length() != MinorDim(row)) throw eSizeMismatch;
#endif

	SetRow(row, array.Pointer());
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
		TYPE* p = (*this)(i);
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
inline TYPE* RaggedArray2DT<TYPE>::operator()(int row) const
{
	return fPtrs[row];
}

template <class TYPE>
inline TYPE* RaggedArray2DT<TYPE>::Pointer(void) const
{
	return fData.Pointer();
}

/* free memory (if allocated) */
template <class TYPE>
inline void RaggedArray2DT<TYPE>::Free(void)
{
	fMajorDim = 0;
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

#endif /* _RAGGED_ARRAY_2D_T_H_ */
