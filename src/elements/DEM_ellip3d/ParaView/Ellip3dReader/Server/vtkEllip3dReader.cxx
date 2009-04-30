/*=========================================================================

   Ellip3d output file reader for VTK -- Implementation

   Written by Tom Buzbee

=========================================================================*/
#include "vtkEllip3dReader.h"

#include "vtkCellArray.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkAppendPolyData.h"
#include "vtkMatrix4x4.h"
#include "vtkSphereSource.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"

#include <fstream>

#define DEFAULT_RESOLUTION 16
#define MAX_LINE_LENGTH 700
#define DEGREES 180.0 / 3.14159265

#define tryRead(command, lineNo, emsg)   command;\
if (file.rdstate())\
{\
   vtkErrorMacro(<< "Line " << lineNo << ": " << emsg);\
   file.close();\
   return 1;\
}

vtkCxxRevisionMacro(vtkEllip3dReader, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkEllip3dReader);

//----------------------------------------------------------------------------
vtkEllip3dReader::vtkEllip3dReader()
{
   vtkPolyData *output = vtkPolyData::New();
   this->SetOutput(output);
   // Releasing data for pipeline parallism.
   // Filters will know it is empty. 
   output->ReleaseData();
   output->Delete();

   Resolution = DEFAULT_RESOLUTION;
}

//----------------------------------------------------------------------------
vtkEllip3dReader::~vtkEllip3dReader()
{
}

//----------------------------------------------------------------------------
vtkPolyData* vtkEllip3dReader::GetOutput()
{
   return this->GetOutput(0);
}

//----------------------------------------------------------------------------
vtkPolyData* vtkEllip3dReader::GetOutput(int idx)
{
   return vtkPolyData::SafeDownCast(this->GetOutputDataObject(idx));
}

//----------------------------------------------------------------------------
void vtkEllip3dReader::SetOutput(vtkPolyData *output)
{
   this->GetExecutive()->SetOutputData(0, output);
}

//----------------------------------------------------------------------------
int vtkEllip3dReader::RequestUpdateExtent(vtkInformation *,
                                          vtkInformationVector **,
                                          vtkInformationVector *outputVector)
{
   vtkInformation *outInfo = outputVector->GetInformationObject(0);

   int piece, numPieces, ghostLevel;

   piece = outInfo->Get(
         vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
   numPieces = outInfo->Get(
         vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
   ghostLevel = outInfo->Get(
         vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
  
   // make sure piece is valid
   if (piece < 0 || piece >= numPieces)
   {
      return 1;
   }
  
   if (ghostLevel < 0)
   {
      return 1;
   }
  
   return 1;
}

//----------------------------------------------------------------------------
int vtkEllip3dReader::RequestData(vtkInformation *,
                                  vtkInformationVector **,
                                  vtkInformationVector *outputVector)
{
   // Special vtkPolyData object we will return our results in
   vtkInformation *outInfo = outputVector->GetInformationObject(0);
   vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(
                                                vtkDataObject::DATA_OBJECT()));

   // Composite vtkPolyData object we will use to group the ellipsoids
   vtkAppendPolyData *elements = vtkAppendPolyData::New();

   // Variables used for parsing
   ifstream file(this->FileName);
   int numParticles = 0;
   double ignore;
   double position[3];                      // {x, y, z}
   double radius[3];                        // Radii
   double axleA[3], axleB[3], axleC[3];     // Orientation vectors

   // Make sure the file is open the file for reading
   if (!file)
   {
      vtkErrorMacro(<<"Could not read file");
      return 1;
   }

   // Extract the number of particles from the first line
   tryRead(file >> numParticles, 1, "Invalid number of particles");
   tryRead(file.ignore(MAX_LINE_LENGTH, '\n'), 1, "Input error");

   // We need to store the particles separately so we can delete them later
   vtkPolyData *particles[numParticles];

   // Discard the dimensions and the column headers
   tryRead(file.ignore(MAX_LINE_LENGTH, '\n'), 2, "Input error");
   tryRead(file.ignore(MAX_LINE_LENGTH, '\n'), 3, "Input error");

   // Read and create each individual particle
   for (int i = 0; i < numParticles; i++)
   {
      // Get the fields we're interested in
      tryRead(file >> ignore, i + 4, "Invalid ID");                // ID
      tryRead(file >> ignore, i + 4, "Invalid type");              // type
      tryRead(file >> radius[0], i + 4, "Invalid radius_a");       // radius_a
      tryRead(file >> radius[1], i + 4, "Invalid radius_b");       // radius_b
      tryRead(file >> radius[2], i + 4, "Invalid radius_c");       // radius_c
      tryRead(file >> position[0], i + 4, "Invalid position_x");   // position_x
      tryRead(file >> position[1], i + 4, "Invalid position_y");   // position_y
      tryRead(file >> position[2], i + 4, "Invalid position_z");   // position_z
      tryRead(file >> axleA[0], i + 4, "Invalid axle_a_x");        // axle_a_x
      tryRead(file >> axleA[1], i + 4, "Invalid axle_a_y");        // axle_a_y
      tryRead(file >> axleA[2], i + 4, "Invalid axle_a_z");        // axle_a_z
      tryRead(file >> axleB[0], i + 4, "Invalid axle_b_x");        // axle_b_x
      tryRead(file >> axleB[1], i + 4, "Invalid axle_b_y");        // axle_b_y
      tryRead(file >> axleB[2], i + 4, "Invalid axle_b_z");        // axle_b_z
      tryRead(file >> axleC[0], i + 4, "Invalid axle_c_x");        // axle_c_x
      tryRead(file >> axleC[1], i + 4, "Invalid axle_c_y");        // axle_c_y
      tryRead(file >> axleC[2], i + 4, "Invalid axle_c_z");        // axle_c_z

      // Ignore the rest of the line
      tryRead(file.ignore(MAX_LINE_LENGTH, '\n'), i + 4, "Input error");

      // Create and add the particle to the set
      particles[i] = Ellipsoid(position, radius, axleA, axleB, axleC);
      elements->AddInput(particles[i]);

      // Update the progress bar every ten particles
      if (i % 10 == 0)
      {
         this->UpdateProgress((i + 1) / numParticles);
      }
   }

   // Copy the results to the output object
   elements->Update();
   output->DeepCopy(elements->GetOutput());

   // Final bookkeeping
   file.close();
   elements->Delete();
   for (int i = 0; i < numParticles; i++)
   {
      particles[i]->Delete();
   }

   return 1;
}

//----------------------------------------------------------------------------
int vtkEllip3dReader::FillOutputPortInformation(int, vtkInformation* info)
{
   info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
   return 1;
}

//----------------------------------------------------------------------------
void vtkEllip3dReader::PrintSelf(ostream& os, vtkIndent indent)
{
   this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
vtkPolyData* vtkEllip3dReader::Ellipsoid(double position[3], double radius[3],
                              double axleA[3], double axleB[3], double axleC[3])
{
   // Create a sphere and set the resolution
   vtkSphereSource *sphere = vtkSphereSource::New();
   sphere->SetPhiResolution(Resolution);
   sphere->SetThetaResolution(Resolution);

   // Create a transformation matrix
   vtkTransform *transform = vtkTransform::New();

   // Set up the rotations first
   double vec1[3], vec2[3], vec3[3];
   transform->PostMultiply();

   // We need to take the cos() of each element in the axle vectors to get
   // proper Cartesian lengths
   axleA[0] = cos(axleA[0]); axleA[1] = cos(axleA[1]); axleA[2] = cos(axleA[2]);
   axleB[0] = cos(axleB[0]); axleB[1] = cos(axleB[1]); axleB[2] = cos(axleB[2]);
   axleC[0] = cos(axleC[0]); axleC[1] = cos(axleC[1]); axleC[2] = cos(axleC[2]);

   // Rotate particle's X axis into place around the Z axis
   vec1[0] = 1;         vec1[1] = 0;         vec1[2] = 0;
   vec2[0] = axleA[0];  vec2[1] = axleA[1];  vec2[2] = 0;
   vec3[0] = 0;         vec3[1] = 0;         vec3[2] = 1;
   transform->RotateWXYZ(angle(vec1, vec2, vec3), vec3);

   // Rotate particle's Z axis into place around its Y axis
   vec2[0] = 1;         vec2[1] = 0;         vec2[2] = 0;
   transform->TransformPoint(vec2, vec1);
   vec2[0] = 0;         vec2[1] = 1;         vec2[2] = 0;
   transform->TransformPoint(vec2, vec3);
   transform->RotateWXYZ(angle(vec1, axleA, vec3), vec3);

   // Rotate particle's Z axis into place
   vec2[0] = 0;         vec2[1] = 0;         vec2[2] = 1;
   transform->TransformPoint(vec2, vec1);
   vec2[0] = 1;         vec2[1] = 0;         vec2[2] = 0;
   transform->TransformPoint(vec2, vec3);
   transform->RotateWXYZ(angle(vec1, axleC, vec3), vec3);

   // Pre-multiply a matrix for scaling the sphere into an ellipsoid
   transform->PreMultiply();
   transform->Scale(2 * radius[0], 2 * radius[1], 2 * radius[2]);

   // Post-multiply a matrix to translate the ellipsoid into its proper position
   transform->PostMultiply();
   transform->Translate(position[0], position[1], position[2]);

   // Encapsulate the transformation in a filter
   vtkTransformPolyDataFilter *filter = vtkTransformPolyDataFilter::New();
   filter->SetTransform(transform);
   filter->SetInput(sphere->GetOutput());
   filter->Update();

   // Copy the results into the output object
   vtkPolyData *output = vtkPolyData::New();
   output->DeepCopy(filter->GetOutput());

   // Delete the objects we no longer need
   transform->Delete();
   sphere->Delete();
   filter->Delete();

   return output;
}

//----------------------------------------------------------------------------
double vtkEllip3dReader::angle(double vec1[3], double vec2[3], double vec3[3])
{
   double angle, norm1, norm2, norm3;
   vtkMatrix4x4 *matrix = vtkMatrix4x4::New();

   // Normalize the vectors in case they aren't already
   norm1 = sqrt(vec1[0] * vec1[0] + vec1[1] * vec1[1] + vec1[2] * vec1[2]);
   norm2 = sqrt(vec2[0] * vec2[0] + vec2[1] * vec2[1] + vec2[2] * vec2[2]);
   norm3 = sqrt(vec3[0] * vec3[0] + vec3[1] * vec3[1] + vec3[2] * vec3[2]);
   vec1[0] = vec1[0] / norm1;
   vec1[1] = vec1[1] / norm1;
   vec1[2] = vec1[2] / norm1;
   vec2[0] = vec2[0] / norm2;
   vec2[1] = vec2[1] / norm2;
   vec2[2] = vec2[2] / norm2;
   vec3[0] = vec3[0] / norm3;
   vec3[1] = vec3[1] / norm3;
   vec3[2] = vec3[2] / norm3;

   // Calculate the angle
   angle = acos(vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2]);

   // Set angle to 0 if acos returned NaN
   if (isnan(angle))
      angle = 0;

   // Set up the matrix to determine the sign
   matrix->SetElement(0, 0, vec1[0]);
   matrix->SetElement(1, 0, vec1[1]);
   matrix->SetElement(2, 0, vec1[2]);
   matrix->SetElement(3, 0, 1);
   matrix->SetElement(0, 1, 0);
   matrix->SetElement(1, 1, 0);
   matrix->SetElement(2, 1, 0);
   matrix->SetElement(3, 1, 1);
   matrix->SetElement(0, 2, vec2[0]);
   matrix->SetElement(1, 2, vec2[1]);
   matrix->SetElement(2, 2, vec2[2]);
   matrix->SetElement(3, 2, 1);
   matrix->SetElement(0, 3, vec3[0]);
   matrix->SetElement(1, 3, vec3[1]);
   matrix->SetElement(2, 3, vec3[2]);
   matrix->SetElement(3, 3, 1);

   // Adjust the sign
   angle = matrix->Determinant() < 0 ? -angle : angle;
   matrix->Delete();

   // Return the angle in degrees
   return angle * DEGREES;
}
