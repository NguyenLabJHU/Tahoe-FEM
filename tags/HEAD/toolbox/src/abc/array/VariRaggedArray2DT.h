/* $Id: VariRaggedArray2DT.h,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (02/17/2000)                                          */
/* VariRaggedArray2DT with dynamics resizing functions                    */

#ifndef _VARI_RAGGED_2D_T_H_
#define _VARI_RAGGED_2D_T_H_

/* base class */
#include "RaggedArray2DT.h"

/* direct members */
#include "VariArrayT.h"

template <class TYPE>
class VariRaggedArray2DT: public RaggedArray2DT<TYPE>
{
public:

	/* constructors */
	VariRaggedArray2DT(int headroom);
	VariRaggedArray2DT(int majordim, int minordim, int blocksize = 1);

	/* overidden base class functions */
	void Configure(const ArrayT<int>& rowcounts, int blocksize = 1);
	VariRaggedArray2DT<TYPE>& operator=(const VariRaggedArray2DT& source);
	void Copy(const AutoFill2DT<TYPE>& source);
	void Copy(const ArrayT<int>& rowcounts, const ArrayT<TYPE*>& data);
	void CopyCompressed(const AutoFill2DT<TYPE>& source);
	void Free(void);

	/* adding rows */
	void AddRow(int row, const ArrayT<TYPE>& row_data);
	//void AddRowsAt(int row, const ArrayT<TYPE>& row_data);

	/* deleting rows */
	//void DeleteRow(int row);
	//void DeleteRows(const ArrayT<TYPE>& rows);

	/* write data (resizes if needed) */
	void SetRow(int row, const ArrayT<TYPE>& array);
	void SetRow(int row, const TYPE* array);

private:

	/* resize wrappers */
	VariArrayT<TYPE*> fPtrs_wrap;
	VariArrayT<TYPE>  fData_wrap;
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructors */
template <class TYPE>
VariRaggedArray2DT<TYPE>::VariRaggedArray2DT(int headroom):
	fPtrs_wrap(headroom, fPtrs),
	fData_wrap(headroom, fData)
{

}

template <class TYPE>
VariRaggedArray2DT<TYPE>::VariRaggedArray2DT(int majordim, int minordim, int blocksize):
	RaggedArray2DT<TYPE>(majordim, minordim, blocksize),
	fPtrs_wrap(15, fPtrs),
	fData_wrap(15, fData)
{

}

/* configuration functions */
template <class TYPE>
inline void VariRaggedArray2DT<TYPE>::Configure(const ArrayT<int>& rowcounts, int blocksize)
{
	/* free all memory */
	fPtrs_wrap.Free();
	fData_wrap.Free();

	/* inherited */
	RaggedArray2DT<TYPE>::Configure(rowcounts, blocksize);

}

/* assigment operator */
template <class TYPE>
inline VariRaggedArray2DT<TYPE>& VariRaggedArray2DT<TYPE>::operator=(const VariRaggedArray2DT& source)
{
	/* free all memory */
	fPtrs_wrap.Free();
	fData_wrap.Free();

	/* inherited */
	RaggedArray2DT<TYPE>::operator=(source);
}

template <class TYPE>
inline void VariRaggedArray2DT<TYPE>::Copy(const AutoFill2DT<TYPE>& source)
{
	/* free all memory */
	fPtrs_wrap.Free();
	fData_wrap.Free();

	/* inherited */
	RaggedArray2DT<TYPE>::Copy(source);
}

template <class TYPE>
inline void VariRaggedArray2DT<TYPE>::Copy(const ArrayT<int>& rowcounts,
	const ArrayT<TYPE*>& data)
{
	/* free all memory */
	fPtrs_wrap.Free();
	fData_wrap.Free();

	/* inherited */
	RaggedArray2DT<TYPE>::Copy(rowcounts, data);
}

template <class TYPE>
void VariRaggedArray2DT<TYPE>::CopyCompressed(const AutoFill2DT<TYPE>& source) // removes empty rows
{
	/* free all memory */
	fPtrs_wrap.Free();
	fData_wrap.Free();

	/* inherited */
	RaggedArray2DT<TYPE>::CopyCompressed(rowcounts, data);
}

/* free memory (if allocated) */
template <class TYPE>
inline void VariRaggedArray2DT<TYPE>::Free(void)
{	
	/* free all memory */
	fPtrs_wrap.Free();
	fData_wrap.Free();

	/* inherited */
	RaggedArray2DT<TYPE>::Free();
}

/* adding rows */
template <class TYPE>
void VariRaggedArray2DT<TYPE>::AddRow(int row, const ArrayT<TYPE>& row_data)
{
	/* check */
	if (row < 0 || row > MajorDim()) throw eOutOfRange;

	/* dimensions */
	int size = row_data.Length();
	int TYPE_size  = sizeof(TYPE);

	if (MajorDim() == 0)
	{
		fData_wrap.SetLength(size, false);
		memcpy(fData.Pointer(), row_data.Pointer(), TYPE_size*size);
	
		fPtrs_wrap.SetLength(2, false);
		fPtrs[0] = fData.Pointer();
		fPtrs[1] = fPtrs[0] + size;
	}
	else
	{
		/* resize data array */
		int old_data_length = fData.Length();
		fData_wrap.SetLength(old_data_length + size, true); //fPtrs is now stale

		/* reset pointers */
		if (fData.Pointer() != fPtrs[0])
		{
			int last_row_length = MinorDim(0);
			fPtrs[0] = fData.Pointer();
			for (int i = 1; i < fMajorDim; i++)
			{
				int row_length = MinorDim(i);
				fPtrs[i] = fPtrs[i-1] + last_row_length;
				last_row_length = row_length;
			}
			fPtrs[fMajorDim] = fPtrs[fMajorDim - 1] + last_row_length;
		}

		/* resize pointers array */
		int old_ptrs_length = fPtrs.Length();
		fPtrs_wrap.SetLength(old_ptrs_length + 1, true);
		memmove(&fPtrs[row + 1], &fPtrs[row], sizeof(TYPE*)*(old_ptrs_length - row));

		/* shift pointers */
		TYPE** p = fPtrs.Pointer(row + 1);
		for (int j = row + 1; j < fPtrs.Length(); j++)
			*p++ += size;
	
		/* push up data and insert row */
		int offset = fPtrs[row] - fPtrs[0];
		int copy_size = old_data_length - offset;
		TYPE* data_src = fPtrs[row];
		memmove(data_src + size, data_src, TYPE_size*copy_size);
		memcpy(data_src, row_data.Pointer(), TYPE_size*size);
	}
	
	/* correct dimensions */
	fMajorDim++;
	fMaxMinorDim = (size > fMaxMinorDim) ? size : fMaxMinorDim;

#if __option(extended_errorcheck)
	/* check pointer bounds */
	if (fPtrs.First() != fData.Pointer() ||
	    fPtrs.Last() - fPtrs.First() != fData.Length())
	{
		cout << "\n VariRaggedArray2DT<TYPE>::AddRow: pointer reassignment error" << endl;
		throw eGeneralFail;
	}
#endif
}

//void AddRowsAt(int row, const ArrayT<TYPE>& row_data);

/* deleting rows */
//void DeleteRow(int row);
//void DeleteRows(const ArrayT<TYPE>& rows);

/* write data */
template <class TYPE>
void VariRaggedArray2DT<TYPE>::SetRow(int row, const ArrayT<TYPE>& array) // with resizing
{
	int old_minor_dim = MinorDim(row);
	if (array.Length() == old_minor_dim)
		/* inherited */
		RaggedArray2DT<TYPE>::SetRow(row, array.Pointer());
	else
	{
		/* dimensions */
		int shift = array.Length() - old_minor_dim;
		int TYPE_size  = sizeof(TYPE);

		/* resize memory and reset pointers */
		if (shift != 0)
		{
			/* resize data array */
			int old_data_length = fData.Length();
			fData_wrap.SetLength(old_data_length + shift, true); // fPtrs is now stale
	
			/* reset pointers */
			if (fData.Pointer() != fPtrs[0])
			{
				int last_row_length = MinorDim(0);
				fPtrs[0] = fData.Pointer();
				for (int i = 1; i < fMajorDim; i++)
				{
					int row_length = MinorDim(i);
					fPtrs[i] = fPtrs[i-1] + last_row_length;
					last_row_length = row_length;
				}
				fPtrs[fMajorDim] = fPtrs[fMajorDim - 1] + last_row_length;
			}

			/* push data up and insert row */
			if (row != fMajorDim - 1)
			{
				int offset = fPtrs[row + 1] - fPtrs[0];
				int copy_size = old_data_length - offset;
				TYPE* data_src = fPtrs[row + 1];
				memmove(data_src + shift, data_src, TYPE_size*copy_size);
			}
	
			/* shift pointers */
			TYPE** p = fPtrs.Pointer(row + 1);
			for (int j = row + 1; j < fPtrs.Length(); j++)
				*p++ += shift;
		}
			
		/* write in data */
		memcpy(fPtrs[row], array.Pointer(), TYPE_size*array.Length());
	
		/* correct dimensions */
		if (array.Length() > fMaxMinorDim)
			fMaxMinorDim = array.Length();
		else if (fMaxMinorDim == old_minor_dim)
		{
			int new_max_minor_dim = 0;
			for (int i = 0; i < fMajorDim && new_max_minor_dim != fMaxMinorDim; i++)
			{
				int minor_dim = MinorDim(i);
				if (minor_dim > new_max_minor_dim)
					new_max_minor_dim = minor_dim;
			}
			
			/* set */
			fMaxMinorDim = new_max_minor_dim;
		}
	}
#if __option(extended_errorcheck)
	/* check pointer bounds */
	if (fPtrs.First() != fData.Pointer() ||
	    fPtrs.Last() - fPtrs.First() != fData.Length())
	{
		cout << "\n VariRaggedArray2DT<TYPE>::SetRow: pointer reassignment error" << endl;
		throw eGeneralFail;
	}
#endif
}

template <class TYPE>
inline void VariRaggedArray2DT<TYPE>::SetRow(int row, const TYPE* array) // no resizing
{
	/* inherited */
	RaggedArray2DT<TYPE>::SetRow(row, array);
}

#endif /* _VARI_RAGGED_2D_T_H_ */
