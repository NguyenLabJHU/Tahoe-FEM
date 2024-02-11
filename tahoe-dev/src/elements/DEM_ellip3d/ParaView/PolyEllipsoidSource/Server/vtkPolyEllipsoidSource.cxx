#include "vtkPolyEllipsoidSource.h"

#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkAppendPolyData.h"
#include "vtkSphereSource.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"

#include <math.h>

vtkCxxRevisionMacro(vtkPolyEllipsoidSource, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkPolyEllipsoidSource);

//----------------------------------------------------------------------------
vtkPolyEllipsoidSource::vtkPolyEllipsoidSource()
{
  this->Center[0] = 0.0;
  this->Center[1] = 0.0;
  this->Center[2] = 0.0;
  this->RadiusXPlus = 1;
  this->RadiusXMinus = 1;
  this->RadiusYPlus = 1;
  this->RadiusYMinus = 1;
  this->RadiusZPlus = 1;
  this->RadiusZMinus = 1;
  this->Resolution = 16;

  this->SetNumberOfInputPorts(0);
}

//----------------------------------------------------------------------------
int vtkPolyEllipsoidSource::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the output
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Organize the parameters for each octant so we can easily loop through
  // The octants are numbered as follows:
  // Octant 0: +x +y +z
  // Octant 1: -x +y +z
  // Octant 2: -x -y +z
  // Octant 3: +x -y +z
  // Octant 4: +x +y -z
  // Octant 5: -x +y -z
  // Octant 6: -x -y -z
  // Octant 7: +x -y -z
  double StartPhi[]  = {   0,   0,   0,   0,  90,  90,  90,  90};
  double EndPhi[]    = {  90,  90,  90,  90, 180, 180, 180, 180};
  double StartTheta[]= {   0,  90, 180, 270,   0,  90, 180, 270};
  double EndTheta[]  = {  90, 180, 270,   0,  90, 180, 270,   0};
  double RadiusX[]   = {RadiusXPlus,  RadiusXMinus, RadiusXMinus, RadiusXPlus,
                        RadiusXPlus,  RadiusXMinus, RadiusXMinus, RadiusXPlus};
  double RadiusY[]   = {RadiusYPlus,  RadiusYPlus,  RadiusYMinus, RadiusYMinus,
                        RadiusYPlus,  RadiusYPlus,  RadiusYMinus, RadiusYMinus};
  double RadiusZ[]   = {RadiusZPlus,  RadiusZPlus,  RadiusZPlus,  RadiusZPlus,
                        RadiusZMinus, RadiusZMinus, RadiusZMinus, RadiusZMinus};

  // PolyData objects for each of the 8 octants
  vtkPolyData *octant[8];

  vtkAppendPolyData *polyEllip = vtkAppendPolyData::New();

  // A sphere source, a transform and a filter we'll reuse for each octant
  vtkSphereSource *sphere = vtkSphereSource::New();
  sphere->SetPhiResolution(Resolution);
  sphere->SetThetaResolution(Resolution);
  vtkTransform *transform = vtkTransform::New();
  vtkTransformPolyDataFilter *filter = vtkTransformPolyDataFilter::New();

  // Create each octant from a partial sphere
  for (int i = 0; i < 8; i++)
  {
    // Set the angles for the sphere
    sphere->SetStartPhi(StartPhi[i]);
    sphere->SetEndPhi(EndPhi[i]);
    sphere->SetStartTheta(StartTheta[i]);
    sphere->SetEndTheta(EndTheta[i]);

    // Create the transformation matrix
    transform->Identity();
    transform->PreMultiply();
    transform->Scale(2 * RadiusX[i], 2 * RadiusY[i], 2 * RadiusZ[i]);
    transform->PostMultiply();
    transform->Translate(Center[0], Center[1], Center[2]);

    // Apply the transformation
    filter->SetTransform(transform);
    filter->SetInput(sphere->GetOutput());
    filter->Update();

    // Copy the results to the PolyData object and add to the vtkAppendPolyData
    octant[i] = vtkPolyData::New();
    octant[i]->DeepCopy(filter->GetOutput());
    polyEllip->AddInput(octant[i]);
  }


  // Copy the results into the output object
  polyEllip->Update();
  output->DeepCopy(polyEllip->GetOutput());

  // Delete leftovers
  filter->Delete();
  transform->Delete();
  sphere->Delete();
  polyEllip->Delete();
  for (int i = 0; i < 8; i++)
  {
    octant[i]->Delete();
  }

  return 1;
}

//----------------------------------------------------------------------------
void vtkPolyEllipsoidSource::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkPolyEllipsoidSource::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),
               -1);

  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_BOUNDING_BOX(),
               this->Center[0] - this->RadiusXMinus,
               this->Center[0] + this->RadiusXPlus,
               this->Center[1] - this->RadiusYMinus,
               this->Center[1] + this->RadiusYPlus,
               this->Center[2] - this->RadiusZMinus,
               this->Center[2] + this->RadiusZPlus);

  return 1;
}
