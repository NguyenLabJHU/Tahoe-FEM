/* $Id: VTKBodyT.cpp,v 1.1 2001-10-24 18:20:36 paklein Exp $ */

#include "VTKBodyT.h"

#include "vtkPoints.h"

/* constructor */
VTKBodyT::VTKBodyT(const StringT& file_name): 
  inFile(file_name)
{
  /* allocate points */
  points = vtkPoints::New();

  //read model data
  //set up grid
  //set up mapper
  //set up actor
}

/* destructor */
VTKBodyT::~VTKBodyT(void)
{
  /* clean up */
  points->Delete();
}
