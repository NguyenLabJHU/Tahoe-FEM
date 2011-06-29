/*
 * File: dTensor4DT.h 
 *
 */

/*
 * created      : PAK (05/25/97)
 * last modified: PAK (11/11/97)
 */

#include "dTensor4DT.h"
#include "dMatrixT.h"

/*
 * Constructor
 */

using namespace Tahoe;

dTensor4DT::dTensor4DT(void) { }
dTensor4DT::dTensor4DT(int dim0, int dim1, int dim2, int dim3):
	Tensor4DT<double>(dim0,dim1,dim2,dim3) { }
dTensor4DT::dTensor4DT(const dTensor4DT& source):
	Tensor4DT<double>(source) { }

/* Converts 3x3x3x3 Tensor form of the Tangent Modulus to the 6x6 matrix
 * form of the Tangent modulus */
/* Has not been tested !!!*/
void dTensor4DT::ConvertTangentFrom4DTo2D(dTensor4DT& C, dMatrixT fc_ijkl)
{
int i,j;
int HughesIndexArray [6] [2];

if (fc_ijkl.Rows() != 6 || fc_ijkl.Cols() != 6)
    ExceptionT::SizeMismatch("dTensor4DT::ConvertTangentFrom4DTo2D", "cannot convert matrix to rank 4 tensor");

/*if (fDim0 != 3 || fDim1 != 3 || fDim2 != 3 || fDim3 !=3 )
 * {
 *   cout << "Cannot convert matrix fc_ijkl to rank 4 tensor - tensor wrong size";
 *   ExceptionT::SizeMismatch(caller);
 * }
 */

/* Set values for Hughes index array. In this version
 * HughesIndexArray [i] is an array of two numbers which return
 * the 2D matrix version of a vectorized (Voigt) symmetric tensor,
 * e.g. s_ij(k) = sigma( HughesIndexArray [k] [0], HughesIndexArray [k] [1])
 */
for (i = 0; i < 3; i++)
  for (j = 0; j < 2; j++)
    HughesIndexArray[i] [j] = i;
HughesIndexArray [3] [0] = 1;
HughesIndexArray [3] [1] = 2;
HughesIndexArray [4] [0] = 0;
HughesIndexArray [4] [1] = 2;
HughesIndexArray [5] [0] = 1;
HughesIndexArray [5] [1] = 2; 

// convert tangent modulus from rank 2 to rank 4 using dTensor4DT
 for (i=0; i<6; i++)
   for (j=0; j<6; j++) 
     fc_ijkl(i, j) = C ( HughesIndexArray [i] [0], HughesIndexArray [i] [1],
			 HughesIndexArray [j] [0], HughesIndexArray [j] [1]);

} // end ConvertTangentFrom4DTo2D

/* Converts 6x6 Matrix form of the Tangent Modulus to the 3x3x3x3
 * form of the Tangent modulus */
void dTensor4DT::ConvertTangentFrom2DTo4D(dTensor4DT& C, dMatrixT fc_ijkl)
{
	
int i,j,k,l;
int HughesIndexArray [3] [3];

if (fc_ijkl.Rows() != 6 || fc_ijkl.Cols() != 6)
    ExceptionT::SizeMismatch("dTensor4DT::ConvertTangentFrom2DTo4D", "cannot convert matrix to rank 4 tensor");

/*if (fDim0 != 3 || fDim1 != 3 || fDim2 != 3 || fDim3 !=3 )
 * {
 *   cout << "Cannot convert matrix fc_ijkl to rank 4 tensor - tensor wrong size";
 *   ExceptionT::SizeMismatch(caller);
 * }
 */

/* Set values for Hughes index array. HughesIndexArray(i,j) = k,
 * where k is the location in the vectorized (Voigt) version of
 * a symmetric 2D tensor, e.g. sigma(i,j) = s_ij(k)
 */
for (i = 0; i < 3; i++)
  HughesIndexArray[i] [i] = i;
HughesIndexArray [0] [1] = 5;
HughesIndexArray [1] [0] = 5;
HughesIndexArray [0] [2] = 4;
HughesIndexArray [2] [0] = 4;
HughesIndexArray [1] [2] = 3;
HughesIndexArray [2] [1] = 3; 

// convert tangent modulus from rank 2 to rank 4 using dTensor4DT
 for (i=0; i<3; i++)
   for (j=0; j<3; j++)
     for (k=0; k<3; k++)
       for (l=0; l<3; l++)
	 C (i,j,k,l) = 
	   fc_ijkl(HughesIndexArray[i] [j], HughesIndexArray[k] [l]);

} // end ConvertTangentFrom2DTo4D
