cmake_minimum_required(VERSION 2.6)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

ADD_PARAVIEW_PLUGIN(PolyEllipsoidSource "1.0" 
   SERVER_MANAGER_XML serverManager.xml
   SERVER_MANAGER_SOURCES vtkPolyEllipsoidSource.cxx)
