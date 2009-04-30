/*=========================================================================

   Ellip3d output file reader for VTK -- Definition

   Written by Tom Buzbee

=========================================================================*/
#ifndef __vtkEllip3dReader_h
#define __vtkEllip3dReader_h

#include "vtkDataReader.h"

class vtkPolyData;

class VTK_IO_EXPORT vtkEllip3dReader : public vtkDataReader
{
   public:
      static vtkEllip3dReader *New();
      vtkTypeRevisionMacro(vtkEllip3dReader, vtkDataReader);
      void PrintSelf(ostream& os, vtkIndent indent);

      vtkPolyData *GetOutput();
      vtkPolyData *GetOutput(int idx);
      void SetOutput(vtkPolyData *output);

      vtkSetClampMacro(Resolution, int, 3, 64);
      vtkGetMacro(Resolution, int);

   protected:
      vtkEllip3dReader();
      ~vtkEllip3dReader();

      int RequestData(vtkInformation *, vtkInformationVector **,
                      vtkInformationVector *);

      int RequestUpdateExtent(vtkInformation *, vtkInformationVector **,
                              vtkInformationVector *);

      int FillOutputPortInformation(int, vtkInformation*);

   private:
      // Creates an ellipsoid vtkPolyData object defined by position, radii and
      // orthonormal basis vectors
      vtkPolyData *Ellipsoid(double position[3], double radius[3],
                             double axleA[3], double axleB[3], double axleC[3]);

      // Helper method to calculate the signed angle in degrees between vec1 and 
      // vec2 with respect to orthogonal vector vec3
      double angle(double vec1[3], double vec2[3], double vec3[3]);

      int Resolution;
};

#endif
