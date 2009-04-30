#ifndef __vtkPolyEllipsoidSource_h
#define __vtkPolyEllipsoidSource_h

#include "vtkPolyDataAlgorithm.h"

class VTK_GRAPHICS_EXPORT vtkPolyEllipsoidSource : public vtkPolyDataAlgorithm 
{
public:
  vtkTypeRevisionMacro(vtkPolyEllipsoidSource,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkPolyEllipsoidSource *New();

//  vtkSetClampMacro(Resolution, int, 3, 64);
//  vtkGetMacro(Resolution, int);

  vtkSetMacro(RadiusXPlus, double);
  vtkGetMacro(RadiusXPlus, double);
  vtkSetMacro(RadiusXMinus, double);
  vtkGetMacro(RadiusXMinus, double);
  vtkSetMacro(RadiusYPlus, double);
  vtkGetMacro(RadiusYPlus, double);
  vtkSetMacro(RadiusYMinus, double);
  vtkGetMacro(RadiusYMinus, double);
  vtkSetMacro(RadiusZPlus, double);
  vtkGetMacro(RadiusZPlus, double);
  vtkSetMacro(RadiusZMinus, double);
  vtkGetMacro(RadiusZMinus, double);

protected:
  vtkPolyEllipsoidSource();
  ~vtkPolyEllipsoidSource() {}

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  double RadiusXPlus;
  double RadiusXMinus;
  double RadiusYPlus;
  double RadiusYMinus;
  double RadiusZPlus;
  double RadiusZMinus;
  double Center[3];
  int Resolution;
};

#endif
