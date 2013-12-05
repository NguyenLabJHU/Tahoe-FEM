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


#include "sparse.h"
#include <iostream>

namespace memFluid{

SparseMatrix::SparseMatrix(int rows, int columns, REAL* vals, int* rowInds, int* colPtr)
{
	this->rows = rows;
	cols = columns;

	this->vals = vals;
	colptr = colPtr;
	rowind = rowInds;
}

SparseMatrix::~SparseMatrix()
{
	if (colptr != NULL)
	{
		if (colptr == rowind)
		{
			rowind = NULL;
		}
		delete[] colptr;
	}

	if (vals != NULL)
	{
		delete[] vals;
	}

	if (rowind != NULL)
	{
		delete[] rowind;
	}
}

void SparseMatrix::set_cell_value(int row, int col, REAL value)
{
	int i, j;
	for (i = colptr[col]; i < colptr[col + 1]; i++)
	{
		if (rowind[i] == row)
		{
			vals[i] = value;
			return ;
		}
	}

	int* newRowind = new int[colptr[cols] + 1];
	REAL* values = new REAL[colptr[cols] + 1];

	for (i = 0; i < col; i++)
	{
		for (j = colptr[i]; j < colptr[i + 1]; j++)
		{
			newRowind[j] = rowind[j];
			values[j] = vals[j];
		}
	}
	for (j = colptr[i]; j < colptr[i + 1] && rowind[j] < row; j++)
	{
		newRowind[j] = rowind[j];
		values[j] = vals[j];
	}
	newRowind[j] = row;
	values[j] = value;
	for ( ; j < colptr[cols]; j++)
	{
		newRowind[j + 1] = rowind[j];
		values[j + 1] = vals[j];
	}
	for (i++; i <= cols; i++)
	{
		colptr[i]++;
	}

	if (vals != NULL)
	{
		delete[] vals;
	}
	if (rowind != NULL)
	{
		delete[] rowind;
	}

	vals = values;
	rowind = newRowind;
}

SparseMatrix* SparseMatrix::transposed() const
{
	int* temp = new int[rows];
	int i, j, q;

	REAL* values = new REAL[colptr[cols]];
	int* newColptr = new int[rows + 1];
	int* newRowind = new int[colptr[cols]];
	for (i = 0; i < rows; i++)
	{
		temp[i] = 0;
	}

	for (i = 0; i < colptr[cols]; i++)
	{
		temp[rowind[i]]++;
	}

	int total = 0;
	for (i = 0; i < rows; i++)
	{
		newColptr[i] = total;
		total += temp[i];
		temp[i] = newColptr[i];
	}
	newColptr[rows] = total;

	for (i = 0; i < cols; i++)
	{
		for (j = colptr[i]; j < colptr[i + 1]; j++)
		{
			q = temp[rowind[j]]++;
			newRowind[q] = i;
			values[q] = vals[j];
		}
	}

	delete[] temp;

	return new SparseMatrix(cols, rows, values, newRowind, newColptr);
}

SparseMatrix* SparseMatrix::symmetric_lower_inverse() const
{
	SparseMatrix* decomposed = this->cholesky_decompose_lower_triangular_returning_lower_triangular();
	if (decomposed != NULL)
	{
		SparseMatrix* inverse = decomposed->invert_lower_triangular_cholesky_decomposed_returning_lower_triangular();
		delete decomposed;
		SparseMatrix* inverseTrans = inverse->transposed();
		SparseMatrix* result = inverseTrans->multiply_F(*inverse);
		delete inverse;
		delete inverseTrans;
		return result;
	}
	else
	{
		return NULL;
	}
}

SparseMatrix* SparseMatrix::multiply_F(const SparseMatrix& second) const
{
	REAL *resultColumn = new REAL[rows];
	Vector<REAL> results(rows * second.cols / 2);
	Vector<int> rowInds(rows * second.cols / 2);
	int i, k, l;
	SparseMatrix *result = new SparseMatrix(rows, second.cols);
	result->colptr[0] = 0;
	for (i = 0; i < rows; i++)
	{
		resultColumn[i] = 0;
	}
	for (i = 0; i < second.cols; i++)
	{
		for (k = second.colptr[i]; k < second.colptr[i + 1]; k++)
		{
			for (l = colptr[second.rowind[k]]; l < colptr[second.rowind[k] + 1]; l++)
			{
				resultColumn[rowind[l]] += second.vals[k] * vals[l];
			}
		}
		for (k = 0; k < rows; k++)
		{
			if (resultColumn[k] != 0)
			{
				results.add(resultColumn[k]);
				rowInds.add(k);
				resultColumn[k] = 0;
			}
		}
		result->colptr[i + 1] = results.size();
	}

	delete[] resultColumn;

	result->setNNZ(results.size());
	for (i = 0; i < results.size(); i++)
	{
		result->vals[i] = results[i];
		result->rowind[i] = rowInds[i];
	}

	return result;
}

SparseMatrix* SparseMatrix::multiply_LM(const SparseMatrix& second) const
{
	REAL *resultColumn = new REAL[rows];
	int i, k, l;
	SparseMatrix *result;
	int nnz = 0;
	for (i = 0; i < rows; i++)
	{
		resultColumn[i] = 0;
	}
	for (i = 0; i < second.cols; i++)
	{
		for (k = second.colptr[i]; k < second.colptr[i + 1]; k++)
		{
			for (l = colptr[second.rowind[k]]; l < colptr[second.rowind[k] + 1]; l++)
			{
				resultColumn[rowind[l]] += second.vals[k] * vals[l];
			}
		}
		for (k = 0; k < rows; k++)
		{
			if (resultColumn[k] != 0)
			{
				nnz++;
				resultColumn[k] = 0;
			}
		}
	}
	result = new SparseMatrix(rows, second.cols, nnz);
	int count = 0;
	result->colptr[0] = 0;
	for (i = 0; i < second.cols; i++)
	{
		for (k = second.colptr[i]; k < second.colptr[i + 1]; k++)
		{
			for (l = colptr[second.rowind[k]]; l < colptr[second.rowind[k] + 1]; l++)
			{
				resultColumn[rowind[l]] += second.vals[k] * vals[l];
			}
		}
		for (k = 0; k < rows; k++)
		{
			if (resultColumn[k] != 0)
			{
				result->vals[count] = resultColumn[k];
				result->rowind[count++] = k;
				resultColumn[k] = 0;
			}
		}
		result->colptr[i + 1] = count;
	}

	delete[] resultColumn;

	return result;
}

SparseMatrix* SparseMatrix::multiply_diagonal_dense(const REAL* second, int secondCols) const
{
	int i, l;
	SparseMatrix *result = new SparseMatrix(rows, secondCols, colptr[cols]);
	result->colptr[0] = 0;
	int count = 0;
	for (i = 0; i < secondCols; i++)
	{
		for (l = colptr[i]; l < colptr[i + 1]; l++)
		{
			result->vals[count] = second[i] * vals[l];
			result->rowind[count++] = rowind[l];
		}
		result->colptr[i + 1] = count;
	}

	return result;
}

SparseMatrix* SparseMatrix::multiply_diagonal_sparse(const REAL* second, int secondCols) const
{
	int i, l, nnz = 0;
	for (int i = 0; i < secondCols; i++)
	{
		if (second[i] != 0)
		{
			nnz += colptr[i + 1] - colptr[i];
		}
	}
	SparseMatrix *result = new SparseMatrix(rows, secondCols, nnz);
	result->colptr[0] = 0;
	int count = 0;
	for (i = 0; i < secondCols; i++)
	{
		if (second[i] != 0)
		{
			for (l = colptr[i]; l < colptr[i + 1]; l++)
			{
				result->vals[count] = second[i] * vals[l];
				result->rowind[count++] = rowind[l];
			}
		}
		result->colptr[i + 1] = count;
	}

	return result;
}

SparseMatrix* SparseMatrix::multiply_returning_unordered_F(const SparseMatrix& second) const
{
	int* cache = new int[rows];
	Vector<REAL> results(rows * second.cols / 2);
	Vector<int> rowInds(rows * second.cols / 2);
	int i, k, l;
	SparseMatrix *result = new SparseMatrix(rows, second.cols);
	result->colptr[0] = 0;
	for (i = 0; i < rows; i++)
	{
		cache[i] = 0;
	}
	for (i = 0; i < second.cols; i++)
	{
		for (k = second.colptr[i]; k < second.colptr[i + 1]; k++)
		{
			for (l = colptr[second.rowind[k]]; l < colptr[second.rowind[k] + 1]; l++)
			{
				if (cache[rowind[l]] == 0)
				{
					cache[rowind[l]] = results.size() + 1;
					results.add(second.vals[k] * vals[l]);
					rowInds.add(rowind[l]);
				}
				else
				{
					results[cache[rowind[l]] - 1] += second.vals[k] * vals[l];
				}
			}
		}
		result->colptr[i + 1] = results.size();
		for (k = result->colptr[i]; k < result->colptr[i + 1]; k++)
		{
			cache[rowInds[k]] = 0;
		}
	}

	delete[] cache;

	result->setNNZ(results.size());
	for (i = 0; i < results.size(); i++)
	{
		result->vals[i] = results[i];
		result->rowind[i] = rowInds[i];
	}

	return result;
}

REAL SparseMatrix::operator ()(int row, int column) const
{
	if (column <= cols / 2)
	{
		int i = colptr[column];
		while (i < colptr[column + 1] && rowind[i] < row)
		{
			i++;
		}
		if (i < colptr[column + 1] && rowind[i] == row)
		{
			return vals[i];
		}
	}
	else
	{
		int i = colptr[column + 1] - 1;
		while (i >= colptr[column] && rowind[i] > row)
		{
			i--;
		}
		if (i >= colptr[column] && rowind[i] == row)
		{
			return vals[i];
		}
	}
	return 0;
}

REAL* SparseMatrix::operator *(const REAL* right) const
{
	int i, j;

	REAL* result = new REAL[rows];

	for (i = 0; i < rows; i++)
	{
		result[i] = 0;
	}

	for (i = 0; i < cols; i++)
	{
		if (right[i] != 0)
		{
			for (j = colptr[i]; j < colptr[i + 1]; j++)
			{
				result[rowind[j]] += right[i] * vals[j];
			}
		}
	}

	return result;
}

bool SparseMatrix::operator ==(const SparseMatrix& second) const
{
	if (cols != second.cols || rows != second.rows ||
		colptr[cols] != second.colptr[second.cols])
	{
		return false;
	}
	for (int i = 0; i < cols; i++)
	{
		if (colptr[i] != second.colptr[i])
		{
			return false;
		}
	}

	for (int i = 0; i < colptr[cols]; i++)
	{
		REAL temp = vals[i] / second.vals[i];
		if (rowind[i] != second.rowind[i] ||
			temp < 1 - sparseTHRESHOLD || temp > 1 + sparseTHRESHOLD)
		{
			return false;
		}
	}

	return true;
}

bool SparseMatrix::operator !=(const SparseMatrix& second) const
{
	return !(*this == second);
}

SparseMatrix* SparseMatrix::deep_copy(const SparseMatrix* matrix)
{
	SparseMatrix* m = new SparseMatrix(matrix->rows, matrix->cols, matrix->colptr[matrix->cols]);

	for (int i = 0; i <= matrix->cols; i++)
	{
		m->colptr[i] = matrix->colptr[i];
	}
	for (int i = 0; i < matrix->colptr[matrix->cols]; i++)
	{
		m->vals[i] = matrix->vals[i];
		m->rowind[i] = matrix->rowind[i];
	}

	return m;
}

SparseMatrix* SparseMatrix::cholesky_decompose_lower_triangular_returning_lower_triangular() const
{
	REAL* resultColumn = new REAL[cols];
	Vector<REAL> results(cols * cols / 2);
	Vector<int> rowInds(cols * cols / 2);

	int* colWhereTo = new int[cols];

	SparseMatrix *result = new SparseMatrix(cols, cols);

	result->colptr[0] = 0;

	for (int i = 0; i < cols; i++)
	{
		resultColumn[i] = 0;
	}

	for (int i = 0; i < cols; i++)
	{
		colWhereTo[i] = results.size() + 1;

		int column = 0;
		int col = result->colptr[column + 1];

		int count;

		while (column < i) // going through all the computed columns
		{
			count = colWhereTo[column];

			// if this element in the column is not zero
			if (count != -1 && rowInds[count] == i)
			{
				REAL columnMainValue = results[count++];
				resultColumn[i] -= columnMainValue * columnMainValue;

				// we multiply each element of the rest of the column with the "main value" which is
				// on the i -th row and then subtract that from the current row's sum
				while (count < col)
				{
					resultColumn[rowInds[count]] -= results[count] * columnMainValue;
					count++;
				}

				if (++colWhereTo[column] >= result->colptr[column + 1])
				{
					colWhereTo[column] = -1;
				}
			}
			col = result->colptr[++column + 1];
		}

		col = colptr[i];

		resultColumn[i] += vals[col];
		if (resultColumn[i] <= 0)
		{
			delete[] resultColumn;
			delete[] colWhereTo;
			delete result;
			return NULL;
		}
		resultColumn[i] = sqrt(resultColumn[i]);
		rowInds.add(i);
		results.add(resultColumn[i]);
		col++;

		for (int j = i + 1; j < cols; j++)
		{
			REAL val = resultColumn[j];
			if (col < colptr[i + 1] && j == rowind[col])
			{
				val += vals[col++];
			}
			if (val < -sparseTHRESHOLD || val > sparseTHRESHOLD)
			{
				rowInds.add(j);
				results.add(val / resultColumn[i]);
			}
			resultColumn[j] = 0;
		}
		resultColumn[i] = 0;

		result->colptr[i + 1] = results.size();
	}

	delete[] colWhereTo;
	delete[] resultColumn;

	result->setNNZ(results.size());
	for (int i = 0; i < results.size(); i++)
	{
		result->vals[i] = results[i];
		result->rowind[i] = rowInds[i];
	}

	return result;
}

SparseMatrix* SparseMatrix::ldlt_decompose_lower_triangular_returning_lower_triangular(REAL*& D) const
{
	REAL* resultColumn = new REAL[cols];
	Vector<REAL> results(cols * cols / 4);
	Vector<int> rowInds(cols * cols / 4);

	int* colWhereTo = new int[cols];

	SparseMatrix *result = new SparseMatrix(cols, cols);
	result->colptr[0] = 0;
	D = new REAL[cols];

	for (int i = 0; i < cols; i++)
	{
		resultColumn[i] = 0;
	}

	for (int i = 0; i < cols; i++)
	{
		colWhereTo[i] = results.size() + 1;

		// computing the diagonal element
		REAL diag;
		if (colptr[i] != colptr[cols])
		{
			diag = vals[colptr[i]];
		}
		else
		{
			diag = 0;
		}
		for (int j = 0; j < i; j++)
		{
			if (colWhereTo[j] != -1 && rowInds[colWhereTo[j]] == i)
			{
				REAL multiplier = D[j] * results[colWhereTo[j]];
				diag -= multiplier * results[colWhereTo[j]];

				if (++colWhereTo[j] >= result->colptr[j + 1])
				{
					colWhereTo[j] = -1;
				}
				else
				{
					for (int k = colWhereTo[j]; k < result->colptr[j + 1]; k++)
					{
						resultColumn[rowInds[k]] += multiplier * results[k];
					}
				}
			}
		}
		if (diag == 0)
		{
			delete[] resultColumn;
			delete[] colWhereTo;
			delete[] D;
			delete result;
			return NULL;
		}
		rowInds.add(i);
		results.add(1);
		D[i] = diag;

		diag = 1 / diag;
		// computing the rest of the column...
		for (int j = i + 1, k = colptr[i] + 1; j < cols; j++)
		{
			REAL element = diag * ((k < colptr[i + 1] && rowind[k] == j ? vals[k++] : 0) - resultColumn[j]);
			if (element != 0)
			{
				rowInds.add(j);
				results.add(element);
			}
			resultColumn[j] = 0;
		}

		result->colptr[i + 1] = results.size();
	}

	delete[] colWhereTo;
	delete[] resultColumn;

	result->setNNZ(results.size());
	for (int i = 0; i < results.size(); i++)
	{
		result->vals[i] = results[i];
		result->rowind[i] = rowInds[i];
	}

	return result;
}

REAL* SparseMatrix::solve_eqn(std::vector<REAL>::iterator b_beg, 
			      std::vector<REAL>::iterator b_end) const
{

    	std::vector<REAL> b(b_beg, b_end);

	REAL* result = new REAL[rows];

	SparseMatrix* LT = this->transposed();

	int col;
	REAL sum;
	for (int i = 0; i < cols; i++)
	{
		sum = b[i];
		for (col = LT->colptr[i]; LT->rowind[col] < i; col++)
		{
			sum -= LT->vals[col] * result[LT->rowind[col]];
		}
		result[i] = sum / LT->vals[col];
	}

	delete LT;

	for (int i = cols - 1; i >= 0; i--)
	{
		sum = result[i];
		for (col = colptr[i + 1] - 1; rowind[col] > i; col--)
		{
			sum -= vals[col] * result[rowind[col]];
		}
		result[i] = sum / vals[col];
	}

	return result;
}

REAL* SparseMatrix::solve_ldlt(const REAL* D, const REAL* b) const
{
	REAL* result = new REAL[rows];

	SparseMatrix* LT = this->transposed();

	int col;
	REAL sum;
	for (int i = 0; i < cols; i++)
	{
		sum = b[i];
		for (col = LT->colptr[i]; LT->rowind[col] < i; col++)
		{
			sum -= LT->vals[col] * result[LT->rowind[col]];
		}
		result[i] = sum;
	}

	delete LT;

	for (int i = cols - 1; i >= 0; i--)
	{
		sum = result[i] / D[i];
		for (col = colptr[i + 1] - 1; rowind[col] > i; col--)
		{
			sum -= vals[col] * result[rowind[col]];
		}
		result[i] = sum;
	}

	return result;
}

void SparseMatrix::solve_ldlt_in_place(const SparseMatrix* LT, const REAL* D, REAL* b, int up_to) const
{
	int col;
	REAL sum;
	for (int i = 0; i < cols; i++)
	{
		sum = b[i];
		for (col = LT->colptr[i]; LT->rowind[col] < i; col++)
		{
			sum -= LT->vals[col] * b[LT->rowind[col]];
		}
		b[i] = sum;
	}
	for (int i = cols - 1; i >= up_to; i--)
	{
		sum = b[i] / D[i];
		for (col = colptr[i + 1] - 1; rowind[col] > i; col--)
		{
			sum -= vals[col] * b[rowind[col]];
		}
		b[i] = sum;
	}
}

SparseMatrix* SparseMatrix::invert_diagonal_matrix() const
{
	SparseMatrix* inverse = new SparseMatrix(cols, cols, cols);
	for (int i = 0; i < cols; i++)
	{
		inverse->vals[i] = 1 / vals[colptr[i]];
		inverse->rowind[i] = i;
	}
	for (int i = 0; i <= cols; i++)
	{
		inverse->colptr[i] = i;
	}
	return inverse;
}

SparseMatrix* SparseMatrix::invert_lower_triangular_cholesky_decomposed() const
{
	SparseMatrix* transposed = this->transposed();

	Vector<REAL> results(colptr[cols] * 4);
	Vector<int> rowInds(colptr[cols] * 4);

	int i, j, k, l;
	REAL sum = 0;

	SparseMatrix *result = new SparseMatrix(cols, cols);
	result->colptr[0] = 0;

	for (i = 0; i < cols; i++)
	{
		for (j = 0; j < i; j++)
		{
			k = result->colptr[j];
			while (k < result->colptr[j + 1] && rowInds[k] < j)
			{
				k++;
			}
			for (l = transposed->colptr[i];
				k < result->colptr[j + 1] &&
				l < transposed->colptr[i + 1] &&
				rowInds[k] < i &&
				transposed->rowind[l] < i; )
			{
				if (transposed->rowind[l] < rowInds[k])
				{
					l++;
				}
				else if (rowInds[k] == transposed->rowind[l])
				{
					sum -= results[k] * transposed->vals[l];
					k++;
					l++;
				}
				else
				{
					k++;
				}
			}

			if (sum != 0)
			{
				rowInds.add(j);
				results.add(sum / vals[colptr[i]]);
				sum = 0;
			}
		}

		results.add(1 / vals[colptr[i]]);
		rowInds.add(i);

		for (j = i + 1; j < cols; j++)
		{
			int size = results.size();
			k = result->colptr[i];
			while (k < size && rowInds[k] < i)
			{
				k++;
			}
			for (l = transposed->colptr[j];
				k < size &&
				l < transposed->colptr[j + 1] &&
				rowInds[k] < j &&
				transposed->rowind[l] < j; )
			{
				if (transposed->rowind[l] < rowInds[k])
				{
					l++;
				}
				else if (rowInds[k] == transposed->rowind[l])
				{
					sum -= results[k] * transposed->vals[l];
					k++;
					l++;
				}
				else
				{
					k++;
				}
			}

			if (sum != 0)
			{
				rowInds.add(j);
				results.add(sum / vals[colptr[j]]);
				sum = 0;
			}
		}
		result->colptr[i + 1] = results.size();
	}

	delete transposed;

	result->setNNZ(results.size());
	for (i = 0; i < results.size(); i++)
	{
		result->rowind[i] = rowInds[i];
		result->vals[i] = results[i];
	}

	return result;
}

SparseMatrix* SparseMatrix::invert_lower_triangular_cholesky_decomposed_returning_lower_triangular() const
{
	SparseMatrix* transposed = this->transposed();

	Vector<REAL> results(colptr[cols] * 4);
	Vector<int> rowInds(colptr[cols] * 4);

	int i, j, k, l;
	REAL sum = 0;

	SparseMatrix *result = new SparseMatrix(cols, cols);
	result->colptr[0] = 0;

	for (i = 0; i < cols; i++)
	{
		results.add(1 / vals[colptr[i]]);
		rowInds.add(i);

		for (j = i + 1; j < cols; j++)
		{
			int size = results.size();
			for (k = result->colptr[i], l = transposed->colptr[j];
				k < size &&
				l < transposed->colptr[j + 1] &&
				rowInds[k] < j &&
				transposed->rowind[l] < j; )
			{
				if (transposed->rowind[l] < rowInds[k])
				{
					l++;
				}
				else if (rowInds[k] == transposed->rowind[l])
				{
					sum -= results[k] * transposed->vals[l];
					k++;
					l++;
				}
				else
				{
					k++;
				}
			}

			if (sum != 0)
			{
				rowInds.add(j);
				results.add(sum / vals[colptr[j]]);
				sum = 0;
			}
		}
		result->colptr[i + 1] = results.size();
	}

	delete transposed;

	result->setNNZ(results.size());
	for (i = 0; i < results.size(); i++)
	{
		result->rowind[i] = rowInds[i];
		result->vals[i] = results[i];
	}

	return result;
}

SparseMatrix* SparseMatrix::invert_lower_triangular_ldlt_decomposed(const REAL* D) const
{
	SparseMatrix* Ainv = this->invert_lower_triangular_ldlt_decomposed_returning_lower_triangular(D);
	SparseMatrix* trans = Ainv->transposed();
	for (int i = 0; i < cols; i++)
	{
		Ainv->vals[Ainv->colptr[i]] = 0;
	}
	SparseMatrix* result = Ainv->add(*trans);
	delete trans;
	delete Ainv;
	return result;
}

SparseMatrix* SparseMatrix::invert_lower_triangular_ldlt_decomposed_returning_lower_triangular(const REAL* D) const
{
	SparseMatrix* LT = this->transposed();

	REAL* temp = new REAL[cols];

	Vector<REAL> results(colptr[cols]);
	Vector<int> rowInds(colptr[cols]);

	SparseMatrix* result = new SparseMatrix(cols, cols);
	result->colptr[0] = 0;

	for (int i = 0; i < cols; i++)
	{
		temp[i] = 0;
	}

	for (int i = 0; i < cols; i++)
	{
		temp[i] = 1;

		solve_ldlt_in_place(LT, D, temp, i);

		for (int j = i; j < cols; j++)
		{
			if (temp[j] != 0)
			{
				results.add(temp[j]);
				rowInds.add(j);
				temp[j] = 0;
			}
		}

		result->colptr[i + 1] = results.size();
	}

	delete LT;
	delete[] temp;

	result->setNNZ(results.size());
	for (int i = 0; i < results.size(); i++)
	{
		result->rowind[i] = rowInds[i];
		result->vals[i] = results[i];
	}

	return result;
}

REAL* SparseMatrix::multiply_returning_diagonal(const SparseMatrix& second) const
{
	REAL* result = new REAL[second.cols];
	int i, k, l, secondN = second.cols / 2;
	for (i = 0; i < secondN; i++)
	{
		result[i] = 0;
		for (k = second.colptr[i]; k < second.colptr[i + 1]; k++)
		{
			l = colptr[second.rowind[k]];
			while (rowind[l] < i && l < colptr[second.rowind[k] + 1])
			{
				l++;
			}
			if (rowind[l] == i && l < colptr[second.rowind[k] + 1])
			{
				result[i] += vals[l] * second.vals[k];
			}
		}
	}
	for ( ; i < second.cols; i++)
	{
		result[i] = 0;
		for (k = second.colptr[i]; k < second.colptr[i + 1]; k++)
		{
			l = colptr[second.rowind[k] + 1] - 1;
			while (rowind[l] > i && l >= colptr[second.rowind[k]])
			{
				l--;
			}
			if (rowind[l] == i && l >= colptr[second.rowind[k]])
			{
				result[i] += vals[l] * second.vals[k];
			}
		}
	}

	return result;
}

REAL* SparseMatrix::multiply_three_returning_diagonal(const SparseMatrix& second, const SparseMatrix& third) const
{
	REAL* result = new REAL[third.cols];
	int i, k, j, l, thirdN = third.cols / 2;
	for (i = 0; i < thirdN; i++)
	{
		result[i] = 0;
		for (k = third.colptr[i]; k < third.colptr[i + 1]; k++)
		{
			for (l = second.colptr[third.rowind[k]]; l < second.colptr[third.rowind[k] + 1]; l++)
			{
				j = colptr[second.rowind[l]];
				while (rowind[j] < i && j < colptr[second.rowind[l] + 1])
				{
					j++;
				}
				if (rowind[j] == i && j < colptr[second.rowind[l] + 1])
				{
					result[i] += third.vals[k] * second.vals[l] * vals[j];
				}
			}
		}
	}
	for ( ; i < third.cols; i++)
	{
		result[i] = 0;
		for (k = third.colptr[i]; k < third.colptr[i + 1]; k++)
		{
			for (l = second.colptr[third.rowind[k]]; l < second.colptr[third.rowind[k] + 1]; l++)
			{
				j = colptr[second.rowind[l] + 1] - 1;
				while (rowind[j] > i && j >= colptr[second.rowind[l]])
				{
					j--;
				}
				if (rowind[j] == i && j >= colptr[second.rowind[l]])
				{
					result[i] += third.vals[k] * second.vals[l] * vals[j];
				}
			}
		}
	}
	return result;
}

SparseMatrix* SparseMatrix::multiply_three_F(const SparseMatrix& second, const SparseMatrix& third) const
{
	REAL *resultColumn = new REAL[rows];
	Vector<REAL> results(rows * third.cols / 2);
	Vector<int> rowInds(rows * third.cols / 2);
	int i, k, l, j;
	SparseMatrix *result = new SparseMatrix(rows, third.cols);
	result->colptr[0] = 0;
	for (i = 0; i < rows; i++)
	{
		resultColumn[i] = 0;
	}
	for (i = 0; i < third.cols; i++)
	{
		for (k = third.colptr[i]; k < third.colptr[i + 1]; k++)
		{
			for (l = second.colptr[third.rowind[k]]; l < second.colptr[third.rowind[k] + 1]; l++)
			{
				for (j = colptr[second.rowind[l]]; j < colptr[second.rowind[l] + 1]; j++)
				{
					resultColumn[rowind[j]] += third.vals[k] * second.vals[l] * vals[j];
				}
			}
		}
		for (k = 0; k < rows; k++)
		{
			if (resultColumn[k] != 0)
			{
				results.add(resultColumn[k]);
				rowInds.add(k);
				resultColumn[k] = 0;
			}
		}
		result->colptr[i + 1] = results.size();
	}

	delete[] resultColumn;

	result->setNNZ(results.size());
	for (i = 0; i < results.size(); i++)
	{
		result->vals[i] = results[i];
		result->rowind[i] = rowInds[i];
	}

	return result;
}

SparseMatrix* SparseMatrix::multiply_three_LM(const SparseMatrix& second, const SparseMatrix& third) const
{
	REAL *resultColumn = new REAL[rows];
	int i, k, l, j;
	SparseMatrix *result;
	int nnz = 0;
	for (i = 0; i < rows; i++)
	{
		resultColumn[i] = 0;
	}
	for (i = 0; i < third.cols; i++)
	{
		for (k = third.colptr[i]; k < third.colptr[i + 1]; k++)
		{
			for (l = second.colptr[third.rowind[k]]; l < second.colptr[third.rowind[k] + 1]; l++)
			{
				for (j = colptr[second.rowind[l]]; j < colptr[second.rowind[l] + 1]; j++)
				{
					resultColumn[rowind[j]] += third.vals[k] * second.vals[l] * vals[j];
				}
			}
		}
		for (k = 0; k < rows; k++)
		{
			if (resultColumn[k] != 0)
			{
				nnz++;
				resultColumn[k] = 0;
			}
		}
	}
	result = new SparseMatrix(rows, third.cols, nnz);
	result->colptr[0] = 0;
	int count = 0;
	for (i = 0; i < third.cols; i++)
	{
		for (k = third.colptr[i]; k < third.colptr[i + 1]; k++)
		{
			for (l = second.colptr[third.rowind[k]]; l < second.colptr[third.rowind[k] + 1]; l++)
			{
				for (j = colptr[second.rowind[l]]; j < colptr[second.rowind[l] + 1]; j++)
				{
					resultColumn[rowind[j]] += third.vals[k] * second.vals[l] * vals[j];
				}
			}
		}
		for (k = 0; k < rows; k++)
		{
			if (resultColumn[k] != 0)
			{
				result->vals[count] = resultColumn[k];
				result->rowind[count++] = k;
				resultColumn[k] = 0;
			}
		}
		result->colptr[i + 1] = count;
	}

	delete[] resultColumn;

	return result;
}

SparseMatrix* SparseMatrix::multiply_three_returning_lower_triangular_F(const SparseMatrix& second, const SparseMatrix& third) const
{
	REAL *resultColumn = new REAL[rows];
	Vector<REAL> results((rows * third.cols) / 2 + rows);
	Vector<int> rowInds((rows * third.cols) / 2 + rows);
	int i, k, l, j;
	SparseMatrix *result = new SparseMatrix(rows, third.cols);
	result->colptr[0] = 0;
	for (i = 0; i < rows; i++)
	{
		resultColumn[i] = 0;
	}
	for (i = 0; i < third.cols; i++)
	{
		for (k = third.colptr[i]; k < third.colptr[i + 1]; k++)
		{
			for (l = second.colptr[third.rowind[k]]; l < second.colptr[third.rowind[k] + 1]; l++)
			{
				for (j = colptr[second.rowind[l] + 1] - 1; j >= colptr[second.rowind[l]] && rowind[j] >= i; j--)
				{
					resultColumn[rowind[j]] += third.vals[k] * second.vals[l] * vals[j];
				}
			}
		}
		for (k = 0; k < rows; k++)
		{
			if (resultColumn[k] != 0)
			{
				results.add(resultColumn[k]);
				rowInds.add(k);
				resultColumn[k] = 0;
			}
		}
		result->colptr[i + 1] = results.size();
	}

	delete[] resultColumn;

	result->setNNZ(results.size());
	for (i = 0; i < results.size(); i++)
	{
		result->vals[i] = results[i];
		result->rowind[i] = rowInds[i];
	}

	return result;
}

SparseMatrix* SparseMatrix::multiply_three_returning_lower_triangular_LM(const SparseMatrix& second, const SparseMatrix& third) const
{
	REAL *resultColumn = new REAL[rows];
	int i, k, l, j;
	SparseMatrix *result;
	int nnz = 0;
	for (i = 0; i < rows; i++)
	{
		resultColumn[i] = 0;
	}
	for (i = 0; i < third.cols; i++)
	{
		for (k = third.colptr[i]; k < third.colptr[i + 1]; k++)
		{
			for (l = second.colptr[third.rowind[k]]; l < second.colptr[third.rowind[k] + 1]; l++)
			{
				for (j = colptr[second.rowind[l] + 1] - 1; j >= colptr[second.rowind[l]] && rowind[j] >= i; j--)
				{
					resultColumn[rowind[j]] += third.vals[k] * second.vals[l] * vals[j];
				}
			}
		}
		for (k = 0; k < rows; k++)
		{
			if (resultColumn[k] != 0)
			{
				nnz++;
				resultColumn[k] = 0;
			}
		}
	}
	result = new SparseMatrix(rows, third.cols, nnz);
	result->colptr[0] = 0;
	int count = 0;
	for (i = 0; i < third.cols; i++)
	{
		for (k = third.colptr[i]; k < third.colptr[i + 1]; k++)
		{
			for (l = second.colptr[third.rowind[k]]; l < second.colptr[third.rowind[k] + 1]; l++)
			{
				for (j = colptr[second.rowind[l] + 1] - 1; j >= colptr[second.rowind[l]] && rowind[j] >= i; j--)
				{
					resultColumn[rowind[j]] += third.vals[k] * second.vals[l] * vals[j];
				}
			}
		}
		for (k = 0; k < rows; k++)
		{
			if (resultColumn[k] != 0)
			{
				result->vals[count] = resultColumn[k];
				result->rowind[count++] = k;
				resultColumn[k] = 0;
			}
		}
		result->colptr[i + 1] = count;
	}

	delete[] resultColumn;

	return result;
}

SparseMatrix* SparseMatrix::multiply_returning_lower_triangular_F(const SparseMatrix& second) const
{
	REAL *resultColumn = new REAL[rows];
	Vector<REAL> results((rows * second.cols) / 2 + rows);
	Vector<int> rowInds((rows * second.cols) / 2 + rows);
	int i, k, l;
	SparseMatrix *result = new SparseMatrix(rows, second.cols);
	result->colptr[0] = 0;
	for (i = 0; i < rows; i++)
	{
		resultColumn[i] = 0;
	}
	for (i = 0; i < second.cols; i++)
	{
		for (k = second.colptr[i]; k < second.colptr[i + 1]; k++)
		{
			for (l = colptr[second.rowind[k] + 1] - 1; l >= colptr[second.rowind[k]] && rowind[l] >= i; l--)
			{
				resultColumn[rowind[l]] += second.vals[k] * vals[l];
			}
		}
		for (k = 0; k < rows; k++)
		{
			if (resultColumn[k] != 0)
			{
				results.add(resultColumn[k]);
				rowInds.add(k);
				resultColumn[k] = 0;
			}
		}
		result->colptr[i + 1] = results.size();
	}

	delete[] resultColumn;

	result->setNNZ(results.size());
	for (i = 0; i < results.size(); i++)
	{
		result->vals[i] = results[i];
		result->rowind[i] = rowInds[i];
	}

	return result;
}

SparseMatrix* SparseMatrix::multiply_returning_lower_triangular_LM(const SparseMatrix& second) const
{
	REAL *resultColumn = new REAL[rows];
	int i, k, l;
	SparseMatrix *result;
	int nnz = 0;
	for (i = 0; i < rows; i++)
	{
		resultColumn[i] = 0;
	}
	for (i = 0; i < second.cols; i++)
	{
		for (k = second.colptr[i]; k < second.colptr[i + 1]; k++)
		{
			for (l = colptr[second.rowind[k] + 1] - 1; l >= colptr[second.rowind[k]] && rowind[l] >= i; l--)
			{
				resultColumn[rowind[l]] += second.vals[k] * vals[l];
			}
		}
		for (k = 0; k < rows; k++)
		{
			if (resultColumn[k] != 0)
			{
				nnz++;
				resultColumn[k] = 0;
			}
		}
	}
	result = new SparseMatrix(rows, second.cols, nnz);
	result->colptr[0] = 0;
	int count = 0;
	for (i = 0; i < second.cols; i++)
	{
		for (k = second.colptr[i]; k < second.colptr[i + 1]; k++)
		{
			for (l = colptr[second.rowind[k] + 1] - 1; l >= colptr[second.rowind[k]] && rowind[l] >= i; l--)
			{
				resultColumn[rowind[l]] += second.vals[k] * vals[l];
			}
		}
		for (k = 0; k < rows; k++)
		{
			if (resultColumn[k] != 0)
			{
				result->vals[count] = resultColumn[k];
				result->rowind[count++] = k;
				resultColumn[k] = 0;
			}
		}
		result->colptr[i + 1] = count;
	}

	delete[] resultColumn;

	return result;
}

SparseMatrix* SparseMatrix::add(const SparseMatrix& second) const
{
	int first = 0, col = 1, sec = 0;
	int total = 0;

	while (col <= cols)
	{
		while (first < colptr[col] && sec < second.colptr[col])
		{
			if (rowind[first] == second.rowind[sec] && vals[first++] + second.vals[sec++] != 0)
			{
				total++;
			}
			else if (rowind[first] < second.rowind[sec])
			{
				first++;
				total++;
			}
			else
			{
				sec++;
				total++;
			}
		}
		total += colptr[col] - first + second.colptr[col] - sec;
		first = colptr[col];
		sec = second.colptr[col++];
	}

	SparseMatrix* result = new SparseMatrix(rows, cols, total);
	first = 0;
	col = 1;
	sec = 0;
	total = 0;

	int colp = 1;

	result->colptr[0] = 0;
	while (col <= cols)
	{
		while (first < colptr[col] && sec < second.colptr[col])
		{
			REAL res;
			if (rowind[first] == second.rowind[sec])
			{
				res = vals[first] + second.vals[sec++];
				if (res != 0)
				{
					result->rowind[total] = rowind[first];
					result->vals[total++] = res;
				}
				first++;
			}
			else if (rowind[first] < second.rowind[sec])
			{
				result->rowind[total] = rowind[first];
				result->vals[total++] = vals[first++];
			}
			else
			{
				result->rowind[total] = second.rowind[sec];
				result->vals[total++] = second.vals[sec++];
			}
		}
		while (first < colptr[col])
		{
			result->rowind[total] = rowind[first];
			result->vals[total++] = vals[first++];
		}
		while (sec < second.colptr[col])
		{
			result->rowind[total] = second.rowind[sec];
			result->vals[total++] = second.vals[sec++];
		}
		col++;
		result->colptr[colp++] = total;
	}

	return result;
}

void SparseMatrix::multiply_by_number(REAL n)
{
	for (int i = 0; i < colptr[cols]; i++)
	{
		vals[i] *= n;
	}
}

void SparseMatrix::write_matrix_file(const char *filename) const
{
	FILE *f = fopen(filename, "wt");
	fprintf(f, "%d %d %d\n", rows, cols, colptr[cols]);
	for (int i = 0; i < cols; i++)
	{
		for (int j = colptr[i]; j < colptr[i + 1]; j++)
		{
			fprintf(f, "%d %d %.20e\n", rowind[j] + 1, i + 1, (double) vals[j]);
		}
	}
	fclose(f);
}

SparseMatrix* SparseMatrix::read_matrix_file(const char* filename)
{
	FILE *f = fopen(filename, "rt");

	int row = 0, column = 0, nnz = 0;
	fscanf(f, "%d %d %d\n", &row, &column, &nnz);
	SparseMatrix* matrix = new SparseMatrix(row, column, nnz);

	int count = 0;
	int colptr = 1;
	matrix->colptr[0] = 0;
	int oldcol = 1;
	while (!feof(f))
	{
		double temp;
		fscanf(f, "%d %d %lf\n", &row, &column, &temp);

		while (oldcol != column)
		{
			matrix->colptr[colptr++] = count;
			oldcol++;
		}

		matrix->rowind[count] = row - 1;
		matrix->vals[count++] = temp;
	}
	while (colptr <= matrix->cols)
	{
		matrix->colptr[colptr++] = count;
	}
	fclose(f);

	return matrix;
}

SparseMatrix* SparseMatrix::generate_random(int rows, int cols, REAL XMin, REAL XMax, double fillPercent)
{
	int i;

	Vector<int> rowind(int(rows * cols * (fillPercent + 0.05)));
	int* colptr = new int[cols + 1];

	int colp = 1;
	colptr[0] = 0;
	for (i = 0; i < cols; i++)
	{
		for (int j = 0; j < rows; j++)
		{
			double r = double (rand()) / (double (RAND_MAX) + 1.0);
			if (r >= 1.0 - fillPercent)
			{
				rowind.add(j);
			}
		}
		colptr[colp++] = rowind.size();
	}

	colp--;
	REAL* values = new REAL[colptr[colp]];
	int* REALRowind = new int[colptr[colp]];
	for (int i = 0; i < colptr[colp]; i++)
	{
		REALRowind[i] = rowind[i];
		REAL X;
		do
		{
			REAL r = REAL (rand()) / (REAL (RAND_MAX) + 1.0);
			X = XMin + r * (XMax - XMin);
		}
		while (X == 0);
		values[i] = X;
	}

	return new SparseMatrix(rows, cols, values, REALRowind, colptr);
}

SparseMatrix* SparseMatrix::generate_random_lower(int size, REAL XMin, REAL XMax, double fillPercent)
{
	int i;

	Vector<int> rowind(int(size * size * (fillPercent + 0.05)));
	int* colptr = new int[size + 1];

	int colp = 1;
	colptr[0] = 0;
	for (i = 0; i < size; i++)
	{
		rowind.add(i);
		for (int j = i + 1; j < size; j++)
		{
			double r = double (rand()) / (double (RAND_MAX) + 1.0);
			if (r >= 1.0 - fillPercent)
			{
				rowind.add(j);
			}
		}
		colptr[colp++] = rowind.size();
	}

	colp--;
	REAL* values = new REAL[colptr[colp]];
	int* REALRowind = new int[colptr[colp]];
	int whichCol = 0;
	for (int i = 0; i < colptr[colp]; i++)
	{
		REALRowind[i] = rowind[i];
		REAL X;
		do
		{
			REAL r = REAL (rand()) / (REAL (RAND_MAX) + 1.0);
			X = XMin + r * (XMax - XMin);
		}
		while (X == 0 || (colptr[whichCol] == i && X < 0));
		if (colptr[whichCol] == i)
		{
			whichCol++;
		}
		values[i] = X;
	}

	return new SparseMatrix(size, size, values, REALRowind, colptr);
}

} // end memFluid

