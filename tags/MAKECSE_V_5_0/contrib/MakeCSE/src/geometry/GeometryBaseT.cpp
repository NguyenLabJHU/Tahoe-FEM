/*
 * File: GeometryBaseT.cpp
 *
 * Class to initialize shape function arrays with
 * geometry specific values.
 *
 */

/*
 * created      : PAK (10/21/97)
 * last modified: PAK (04/25/99)
 */

#include "GeometryBaseT.h"

/* constructor */
GeometryBaseT::GeometryBaseT(int numnodes, int numfacets): 
	fNumNodes(numnodes),
	fNumFacets(numfacets) 
{ 


}

/* destructor */
GeometryBaseT::~GeometryBaseT(void) { }
