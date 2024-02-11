/* ***********************************************
 *
 * Copyright (c) 2010 Boris Petrov

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * ***********************************************
 */


#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <vector>
#include "parameters.h"
#include "realtypes.h"


namespace memFluid{

class SparseMatrix
{
public:

	SparseMatrix(int rows, int columns, REAL* vals, int* rowInds, int* colPtr);
	~SparseMatrix();

	// sets the value at that row and column - this modifies the matrix and could be O(number_of_nonzeros) in the worst case!
	void set_cell_value(int row, int col, REAL value);

	// Sparse matrix transposition.
	SparseMatrix* transposed() const;
	// Sparse matrix inversion - the matrix should be symmetric, lower triangular, positive definite. Cholesky decomposition is used.
	SparseMatrix* symmetric_lower_inverse() const;

	// Multiplies the matrix with a column vector.
	REAL* operator *(const REAL* right) const;
	// Returns the element at this row and column.
	REAL operator ()(int row, int column) const;

	// Equality and inequality.
	bool operator ==(const SparseMatrix&) const;
	bool operator !=(const SparseMatrix&) const;

	// Multiplies all elements in the matrix by a number IN PLACE!
	void multiply_by_number(REAL);

	SparseMatrix* cholesky_decompose_lower_triangular_returning_lower_triangular() const;
	SparseMatrix* ldlt_decompose_lower_triangular_returning_lower_triangular(REAL*& D) const;
	SparseMatrix* invert_diagonal_matrix() const;
	SparseMatrix* invert_lower_triangular_cholesky_decomposed() const;
	SparseMatrix* invert_lower_triangular_cholesky_decomposed_returning_lower_triangular() const;
	SparseMatrix* invert_lower_triangular_ldlt_decomposed(const REAL* D) const;
	SparseMatrix* invert_lower_triangular_ldlt_decomposed_returning_lower_triangular(const REAL* D) const;
	REAL* solve_eqn(std::vector<REAL>::iterator, std::vector<REAL>::iterator) const;
	REAL* solve_ldlt(const REAL* D, const REAL* b) const;
	void solve_ldlt_in_place(const SparseMatrix* LT, const REAL* D, REAL* b, int up_to = 0) const;

	// Adds two sparse matrices.
	SparseMatrix* add(const SparseMatrix& second) const;

	SparseMatrix* multiply_F(const SparseMatrix& second) const;
	SparseMatrix* multiply_LM(const SparseMatrix& second) const;
	SparseMatrix* multiply_diagonal_dense(const REAL* second, int secondCols) const;
	SparseMatrix* multiply_diagonal_sparse(const REAL* second, int secondCols) const;
	SparseMatrix* multiply_returning_unordered_F(const SparseMatrix& second) const;
	REAL* multiply_returning_diagonal(const SparseMatrix& second) const;
	SparseMatrix* multiply_returning_lower_triangular_F(const SparseMatrix& second) const;
	SparseMatrix* multiply_returning_lower_triangular_LM(const SparseMatrix& second) const;

	SparseMatrix* multiply_three_F(const SparseMatrix& second, const SparseMatrix& third) const;
	SparseMatrix* multiply_three_LM(const SparseMatrix& second, const SparseMatrix& third) const;
	REAL* multiply_three_returning_diagonal(const SparseMatrix& second, const SparseMatrix& third) const;
	SparseMatrix* multiply_three_returning_lower_triangular_F(const SparseMatrix& second, const SparseMatrix& third) const;
	SparseMatrix* multiply_three_returning_lower_triangular_LM(const SparseMatrix& second, const SparseMatrix& third) const;

	inline int columns_count() const { return cols; }
	inline int rows_count() const { return rows; }

	inline REAL* values() const { return vals; }
	inline int* row_indices() const { return rowind; }
	inline int* column_pointers() const { return colptr; }

	void write_matrix_file(const char *) const;
	static SparseMatrix* read_matrix_file(const char *);

	static SparseMatrix* deep_copy(const SparseMatrix* matrix);

	static SparseMatrix* generate_random(int rows, int cols, REAL XMin, REAL XMax, double fillPercent);
	static SparseMatrix* generate_random_lower(int size, REAL XMin, REAL XMax, double fillPercent);

private:

	SparseMatrix(int rows, int columns, int nnz)
	{
		this->rows = rows;
		cols = columns;

		vals = new REAL[nnz];
		colptr = new int[cols + 1];
		rowind = new int[nnz];
	}

	SparseMatrix(int rows, int columns)
	{
		this->rows = rows;
		cols = columns;

		vals = NULL;
		rowind = NULL;
		colptr = new int[cols + 1];
	}

	inline void setNNZ(int nnz)
	{
		vals = new REAL[nnz];
		rowind = new int[nnz];
	}

	int cols;
	int rows;

	int* colptr;	// length cols + 1
	int* rowind;	// length colptr[cols]
	REAL* vals;		// length colptr[cols]

	template <typename T>
	class Vector
	{
	public:
		inline Vector(int initialSize)
		{
			allocated = initialSize;
			count = 0;
			values = new T[allocated];
		}

		inline ~Vector()
		{
			delete[] values;
		}

		void add(const T& elem)
		{
			if (count == allocated)
			{
				T *temp = new T[allocated];
				for (int i = 0; i < allocated; i++)
				{
					temp[i] = values[i];
				}
				delete[] values;
				values = new T[allocated * 2];
				for (int i = 0; i < allocated; i++)
				{
					values[i] = temp[i];
				}
				delete[] temp;
				allocated *= 2;
			}
			values[count++] = elem;
		}

		inline const T& operator[] (int index) const
		{
			return values[index];
		}

		inline T& operator[] (int index)
		{
			return values[index];
		}

		inline int size() const
		{
			return count;
		}

	private:
		T *values;
		int count;
		int allocated;
	};
};


} // end memFluid

#endif // SPARSE_MATRIX_H
