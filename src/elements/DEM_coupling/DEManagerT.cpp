#include "DEManagerT.h"
#include "ModelManagerT.h"
#include "TimeManagerT.h"
#include "NodeManagerT.h"
#include "ElementBaseT.h"
#include "SolverT.h"
#include "HexahedronT.h"

#include "ParentDomainT.h"
#include "LocalArrayT.h"
#include "GhostParticleT.h"

#include "particle.h"
#include "vec.h"

#include <list>

using namespace std;
using namespace dem;

namespace Tahoe {

/* constructor */
DEManagerT::DEManagerT()
    :PrintNum(0)
{

}

/* destructor */
DEManagerT::~DEManagerT(void)
{
    
}

void DEManagerT::TakeParameter(char* database, iArray2DT& fGhostElemSet)
{
    ReadSample(database, TimeStep, NumStep, fGhostElemSet, fGhostParticleList);
}

void DEManagerT::MapToParentDomain(ModelManagerT* fModelManager, iArray2DT& fGhostElemSet)
{
    // retrieve spatial dimensions and coordinates of all nodes
    NumSD = fModelManager->NumDimensions();
    const dArray2DT& fCoordinates = fModelManager->Coordinates();

    // retrieve ElementSet info
    char str[10];
    sprintf(str, "%d", fGhostElemSet(0, 0)); // number 1, corresponding to string "1" from XML file
    StringT fString(str);
    fElementSet = fModelManager->ElementGroup(fString);
    NumElemNodes = fElementSet.MinorDim();

    // set up fGeometryBase and fParentDomain
    fGeometryCode = fModelManager->ElementGroupGeometry(fString);
    fGeometryBase = GeometryT::New(fGeometryCode, NumElemNodes);
    ParentDomainT fParentDomain(fGeometryCode, 8, NumElemNodes);

    // for each ghost particle
    list< particle* >::iterator it;
    for( it = fGhostParticleList.begin(); it != fGhostParticleList.end(); ++it)
    {	
	// retrieve particle coordinate
	vec vposit = (*it)->getCurrPosition();
	dArrayT point, mapped;
	point.Dimension(NumSD);
	mapped.Dimension(NumSD);
	point[0] = vposit.getx();
	point[1] = vposit.gety();
	point[2] = vposit.getz();

	// retrieve node numbering and coordinates
	iArrayT nodenum;
	nodenum.Dimension(NumElemNodes);
	LocalArrayT coords;
	coords.Dimension(NumElemNodes, NumSD);
	for (int i = 0; i < NumElemNodes; i++) {
	    nodenum[i] = fElementSet(static_cast<GhostParticleT*>(*it)->GetElementNum(), i);
	    for (int j = 0; j < NumSD; j++)
		coords(i, j) = fCoordinates(nodenum[i], j);
	}

	// calculate parent domain coordinate and register it
	fParentDomain.MapToParentDomain(coords, point, mapped);
	static_cast<GhostParticleT*>(*it)->SetParentCoord(mapped);

/*
	cout<<endl<<"....... "<<static_cast<GhostParticleT*>(*it)->GetElementNum()<<endl;
	for (int i=0; i<NumElemNodes;i++)
	    cout<<" "<<nodenum[i];
	for (int i=0; i<NumElemNodes;i++){
	    cout<<endl;
	    for (int j=0;j<NumSD;j++)
		cout<<" "<<coords(i,j);
	}
	cout<<endl<<"mapped: "<<mapped[0]<<" "<<mapped[1]<<" "<<mapped[2]<<endl;
*/
    }

}

void DEManagerT::PrintFE(NodeManagerT* fNodeManager)
{
    char fefile[50];
    char stepstr[4];
    sprintf(stepstr, "%03d", PrintNum++); 
    strcpy(fefile, "fem_mesh_"); strcat(fefile, stepstr); strcat(fefile,".dat");
    ofstream ofs(fefile);
    if(!ofs) {
	cout<<"stream error in PrintFE!"<<endl; exit(-1);
    }
    ofs.setf(ios::scientific, ios::floatfield);

    const dArray2DT& fCurrCoords = fNodeManager->CurrentCoordinates();
    ofs << "ZONE N=" << fCurrCoords.MajorDim() << ", E=" << fElementSet.MajorDim() << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" << endl;

    for (int i = 0; i< fCurrCoords.MajorDim(); i++) {
	for (int j = 0; j < fCurrCoords.MinorDim(); j++)
	    ofs << setw(16) << fCurrCoords(i,j);
	ofs << endl;
    }

    for (int i = 0; i< fElementSet.MajorDim(); i++) {
	for (int j = 0; j < fElementSet.MinorDim(); j++)
	    ofs << setw(10) << fElementSet(i,j) + 1;
	ofs << endl;
    }

    ofs.close();
}

void DEManagerT::PrintDE()
{
    char pfile[50];
    char cfile[50];
    char stepstr[4];
    sprintf(stepstr, "%03d", PrintNum-1); 
    strcpy(pfile, "dem_particle_"); strcat(pfile, stepstr);
    strcpy(cfile, "dem_contact_");  strcat(cfile, stepstr);
    printPtcl(pfile);
    printCntct(cfile);
}

void DEManagerT::Run()
{
    cout<<endl<<" DEM part is running..."<<endl;
    assembly::Run(NumStep, PrintNum-1);
}

void DEManagerT::GhostForce(ArrayT<FBC_CardT>& fGhostFBC, ModelManagerT* fModelManager, iArray2DT& fGhostElemSet)
{
    // set up fGhostFBC's dimension
    fGhostFBC.Dimension(fGhostElemSet.MajorDim() * NumSD * NumElemNodes);
    for (int i = 0; i< fGhostFBC.Length(); i++)
	fGhostFBC[i].ClearValues();

    // set up shape function's dimension
    dArrayT Na;
    Na.Dimension(NumElemNodes);

    for (int fg = 0; fg < fGhostElemSet.MajorDim(); fg++)
    {
	list< particle* >::iterator it;
	for( it = fGhostParticleList.begin(); it != fGhostParticleList.end(); ++it)
	    if ( static_cast<GhostParticleT*>(*it)->GetElementNum() == fGhostElemSet(fg,1) ){

		// retrieve force and parent domain coordinate
		vec vforce = (*it)->getForce();
		double force[3] = {vforce.getx(), vforce.gety(), vforce.getz()};
//		cout<<"***  "<<force[0]<<"  "<<force[1]<<"  "<<force[2]<<endl;
		dArrayT coord;
		coord.Dimension(NumSD);
		coord = static_cast<GhostParticleT*>(*it)->GetParentCoord();

		// retrieve shape function
		fGeometryBase->EvaluateShapeFunctions(coord, Na);

		// retrieve node numbering
		iArrayT nodenum;
		nodenum.Dimension(NumElemNodes);
		for (int i = 0; i < NumElemNodes; i++)
		    nodenum[i] = fElementSet(static_cast<GhostParticleT*>(*it)->GetElementNum(), i);

		// setup fGhostFBC
		int k = fg * NumSD * NumElemNodes;
		for(int i = 0; i < NumElemNodes; i++)
		    for (int j = 0; j < NumSD; j++) {
/*
			std::cout<<std::setw(8)<<k
				 <<std::setw(8)<<nodenum[i]
				 <<std::setw(8)<<j
				 <<std::setw(16)<<Na[i]*force[j]
				 <<std::endl;
*/
			fGhostFBC[k++].AddValues(nodenum[i], j, NULL, Na[i]*force[j]);
		    }
	    }

/*
	int k=fg * NumSD * NumElemNodes;
	for(int i = 0; i < NumElemNodes; i++)
	    for (int j = 0; j < NumSD; j++) {
		std::cout<<std::setw(8)<<k
			 <<std::setw(8)<<fGhostFBC[k].Node()
			 <<std::setw(8)<<fGhostFBC[k].DOF()
			 <<std::setw(16)<<fGhostFBC[k].Value()
			 <<std::endl;
		k++;
	    }
*/
    }

}

void DEManagerT::GhostDisplace(NodeManagerT* fNodeManager, iArray2DT& fGhostElemSet)
{

    const dArray2DT& fCurrCoords = fNodeManager->CurrentCoordinates();

    // set up shape function's dimension
    dArrayT Na;
    Na.Dimension(NumElemNodes);

    // for each ghost particle
    list< particle* >::iterator it;
    for( it = fGhostParticleList.begin(); it != fGhostParticleList.end(); ++it)
    {	
	// retrieve parent domain coordinate
	dArrayT coord;
	coord.Dimension(NumSD);
	coord = static_cast<GhostParticleT*>(*it)->GetParentCoord();

	// retrieve shape function
	fGeometryBase->EvaluateShapeFunctions(coord, Na);

	// retrieve node numbering and current coordinates
	iArrayT nodenum;
	nodenum.Dimension(NumElemNodes);
	LocalArrayT coords;
	coords.Dimension(NumElemNodes, NumSD);
	for (int i = 0; i < NumElemNodes; i++) {
	    nodenum[i] = fElementSet(static_cast<GhostParticleT*>(*it)->GetElementNum(), i);
	    for (int j = 0; j < NumSD; j++)
		coords(i, j) = fCurrCoords(nodenum[i], j);
	}

	// calculate updated coordinate and register it
	dArrayT updated;
	updated.Dimension(NumSD);
	for (int i = 0; i < NumSD; i++) {
	    updated[i] = 0; // important to initialize
	    for (int j = 0; j < NumElemNodes; j++)
		updated[i] += Na[j]*coords(j, i);
	}
	(*it)->setCurrPosition( vec(updated[0], updated[1], updated[2]) );
    }

/*
    for( it = fGhostParticleList.begin(); it != fGhostParticleList.end(); ++it)
    {
	vec tmp = (*it)->getCurrPosition();
	tmp.print();
    }
*/
}

} /* namespace Tahoe ends */

