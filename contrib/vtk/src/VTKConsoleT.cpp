/* $Id: VTKConsoleT.cpp,v 1.2 2001-09-20 23:11:42 recampb Exp $ */

#include "VTKConsoleT.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkScalarBarActor.h"
#include "vtkCubeAxesActor2D.h"
#include "vtkRendererSource.h"
#include "vtkTIFFWriter.h"
#include "vtkWindowToImageFilter.h"
#include "vtkLookupTable.h"


VTKConsoleT::VTKConsoleT(void)
{
  /* set console name */
  iSetName("vtk");

  /* add variables to the console */
  iAddVariable("an_integer", fInteger);
  iAddVariable("a_double", fDouble);  
  iAddVariable("a_string", fString);
  iAddVariable("source_file", source_file);
  iAddVariable("output_file", output_file);
  iAddVariable("ugr", ugr);
  iAddVariable("renderer", renderer);
  iAddVariable("renWin", renWin);
  iAddVariable("iren", iren);
  iAddVariable("lut", lut);
  iAddVariable("ugridMapper", ugridMapper);
  iAddVariable("ugridActor", ugridActor);
  iAddVariable("wireActor", wireActor);
  iAddVariable("scalarBar", scalarBar);
  iAddVariable("renSrc", renSrc);

  ugr = vtkUnstructuredGridReader::New();
  ugr->SetFileName(source_file);
  renderer = vtkRenderer::New();
//   renWin = vtkRenderWindow::New();
//   renWin->AddRenderer(renderer);
//   iren = vtkRenderWindowInteractor::New();
//   iren->SetRenderWindow(renWin);


  /* add console commands */
  iAddCommand("Integer_Print");
  iAddCommand("String_Print");
}

/* execute given command - returns false on fail */
bool VTKConsoleT::iDoCommand(const StringT& command, StringT& line)
{
  if (command == "Integer_Print")
    {
      cout << "int = " << fInteger << endl;
      return true;
    }
  else if (command == "String_Print")
    {
      cout << "string = " << fString << endl;
      return true;
    }
  else
    /* drop through to inherited */
    return iConsoleObjectT::iDoCommand(command, line);
}
