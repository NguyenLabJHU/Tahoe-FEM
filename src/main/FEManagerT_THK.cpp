/* $Id: FEManagerT_THK.cpp,v 1.6 2003-05-21 23:47:57 paklein Exp $ */
#include "FEManagerT_THK.h"
#include "ifstreamT.h"
#include "ModelManagerT.h"
#include "NodeManagerT.h"
#include "FieldT.h"
#include "StringT.h"
#include <iostream.h>
#include <fstream.h>
#ifdef BRIDGING_ELEMENT

using namespace Tahoe;

const double tol = 1.0e-8;   // for neighbor searching tolerance
const double root32 = sqrt(3.0)/2.0;    // for neighbor searching tolerance

/* constructor */
FEManagerT_THK::FEManagerT_THK(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
	ifstreamT& bridging_input):
	FEManagerT_bridging(input, output, comm, bridging_input),
        fTheta(2),
        fNcrit(0)
{

}

/* initialize members */
void FEManagerT_THK::Initialize(InitCodeT init)
{
	/* inherited */
	FEManagerT_bridging::Initialize(init);

	// read other parameters and initialize data
        ifstreamT& in = Input();
        in >> ncrit;
        fNcrit = ncrit;
          
        /* obtain list of atoms on which BC's will be applied in FEManagerT_THK */
        ArrayT<StringT> id_list;
        
        ModelManagerT* model = FEManagerT::ModelManager();
        
        /* read node set indexes */
        model->NodeSetList(in, id_list);
 
        /* collect sets */
        model->ManyNodeSets(id_list, fNodes);
        double lparam;
        in >> lparam;    // Read lattice parameter in from input file
          
        /* now use node set numbers to generate neighboring atom list for
           THK calculations */
        int num_nodes = fNodes.Length();  // total number of boundary atoms
        int num_neighbors = 2 * fNcrit + 1;   // maximum number of neighbors per atom in 2D
        fN0.Dimension(num_nodes);  // Array for only corresponding atom in row 0
        fNeighbor0.Dimension(num_nodes, num_neighbors);
        fNeighbor1.Dimension(num_nodes, num_neighbors);
        fNeighbor0 = -1;  // -1 -> no neighbor
        fNeighbor1 = -1;
        
        NodeManagerT* node = FEManagerT::NodeManager();
        const dArray2DT& initcoords = node->InitialCoordinates();  // init coords for all atoms
        dArrayT bacoords, currcoords;   // boundary atom coordinates
        int count = 0;
        
        /* now loop over boundary atoms to find neighbors */
        for (int i = 0; i < fNeighbor0.MajorDim(); i++)
        {
            /* obtain coordinates of each boundary atom */
            initcoords.RowAlias(fNodes[i],bacoords);
            
            /* exhaustive search over all atoms in lattice */
            for (int j = 0; j < initcoords.MajorDim(); j++) 
            {
                initcoords.RowAlias(j,currcoords);  // coordinates of each atom in lattice
                
                /* find neighbors based on lattice geometry */
                for (int k = -fNcrit; k <= fNcrit; k++) 
                {
                    /* row 0 search first */
                    double x0 = bacoords[0] + lparam * (k + .5);
                    double y0 = bacoords[1] + lparam * root32;
                    if (abs(x0-currcoords[0]) < tol && abs(y0-currcoords[1]) < tol)
                    {
                        fNeighbor0(i,k+fNcrit) = j+1; 
                        /* obtain only corresponding 0 atoms */
                        if (k == 0)
                        {
                            fN0[count] = j+1;
                            count++;
                        }    
                    }   
                    /* then row 1 search */
                    double x1 = bacoords[0] + lparam * (k + 1.0);
                    double y1 = bacoords[1] + lparam * sqrt(3.0);
                    if (abs(x1-currcoords[0]) < tol && abs(y1-currcoords[1]) < tol)
                        fNeighbor1(i,k+fNcrit) = j+1; 
                }
            }
        }
        cout << "Neighbor list 0 = \n" << fNeighbor0 << endl;
        cout << "Neighbor list 1 = \n" << fNeighbor1 << endl;
        cout << "Node set 0 = \n" << fN0 << endl;
        cout << "middle layer nodes = \n" << fNodes << endl;
}

/* initialize the current time increment for all groups */
ExceptionT::CodeT FEManagerT_THK::InitStep(void)
{
	/* inherited */
	ExceptionT::CodeT result = FEManagerT_bridging::InitStep();

	// things to do at the start of every time increment
	
	return result;
}

/* close the current time increment for all groups */
ExceptionT::CodeT FEManagerT_THK::CloseStep(void)
{
	/* inherited */
	ExceptionT::CodeT result = FEManagerT_bridging::CloseStep();

	// things to do at the end of every time increment
        /* compute displacement of all THK boundary atoms */
        const dArray2DT& badisp = ComputeStaticDispBC();
        return result;
}

/* calculate boundary displacement of all THK atoms utilizing theta */
const dArray2DT& FEManagerT_THK::ComputeStaticDispBC(void)
{
        fBdisp.Dimension(fNodes.Length(), 2);
        dArray2DT thk, md;
        thk.Dimension(fNodes.Length(), 2);
        md.Dimension(fNodes.Length(), 2);  // for outputting displacements to file
        fBdisp = 0.0;   // Initialize all displacements to 0
        
        /* now take neighbor list to find displacements of neighboring atoms */
        NodeManagerT* node = FEManagerT::NodeManager();
        StringT mddisp = "displacement";   // What is the purpose of this?
        FieldT* field = node->Field(mddisp);
        const FieldT& field2 = *field;
        const dArray2DT& globaldisp = field2[0];
        dArrayT temp0(2), temp1(2), temp2(2), temp3(2), un0(2), temp4(2), temp5(2), temp6(2);
        dMatrixT theta(2);
        
        /* take neighbor lists, get corresponding displacements */
        for (int i = 0; i < fNodes.Length(); i++)  
        {
            temp4 = 0.0;   // Initialize to zero for each boundary atom
            for (int j = -fNcrit; j <= fNcrit; j++) 
            {
                /* check for no neighbors */
                if (fNeighbor0(i,j+fNcrit) == -1)  
                    temp0 = 0.0;   // No displacement contribution - watch memory overwrite
                else
                    globaldisp.RowCopy(fNeighbor0(i,j+fNcrit)-1,temp0);
                       
                if (fNeighbor1(i,j+fNcrit) == -1)  
                    temp1 = 0.0;   // No displacement contribution - watch memory overwrite
                else
                    globaldisp.RowCopy(fNeighbor1(i,j+fNcrit)-1,temp1);  
              
                temp2.DiffOf(temp0,temp1);  
                theta = GetTheta(-j);  
                theta.Multx(temp2,temp3);
                temp4.SumOf(temp4,temp3);  
             }
            globaldisp.RowCopy(fN0[i]-1,un0);  // getting Un,0 -> verified to be correct
            temp5.SumOf(temp4,un0);  // Have Un,-1 now
            globaldisp.RowAlias(fNodes[i],temp6);
            
            // Need to add temp5 to fBdisp or to global displacement array
            thk.SetRow(i,temp5);
            md.SetRow(i,temp6);
        }
        /* Output comparison between MD and THK to file */
        ofstream out;
        dArrayT temp7, temp8;
        out.open("thkcompare.dat");
        if (!out)
        {
            cout << "Could not open comparison file! \n" << endl;
            throw ExceptionT::kGeneralFail;
        }
        out << "N_critical = " << fNcrit << endl;
        out << "MD x-displacement   " << "THK x-displacement   " << "MD y-displacement   " 
            << "THK y-displacement" << endl;
        out.precision(12);
        
        for (int i = 0; i < fNodes.Length(); i++)
        {
            thk.RowAlias(i,temp7);
            md.RowAlias(i,temp8);
            out << temp8[0] << setw(20) << temp7[0] << setw(20) << temp8[1] << setw(20) 
                << temp7[1] << endl;
        }
        
        out.close();
        return fBdisp;
}

/* Obtain theta for a generalized 2D Hexagonal FCC lattice structure */
const dMatrixT& FEManagerT_THK::GetTheta(int index)
{ 
    if (index == -10)
    {
        fTheta(0,0) = 0.003732973124618666;     
        fTheta(0,1) = 0.0002077573555273614;      
        fTheta(1,0) = 0.0002023643877022398;      
        fTheta(1,1) = 0.0012604960883612864;
    }
    else if (index == -9)
    {
        fTheta(0,0) = 0.004554387859895129;     
        fTheta(0,1) = 0.0002809410552782137;      
        fTheta(1,0) = 0.0002718860440221418;      
        fTheta(1,1) = 0.0015426573553735336;
    }
    else if (index == -8)
    {
        fTheta(0,0) = 0.00567792533483227;      
        fTheta(0,1) = 0.0003929687554678991;      
        fTheta(1,0) = 0.0003768190364353315;      
        fTheta(1,1) = 0.0019317920495601905;
    }
    else if (index == -7)
    {
        fTheta(0,0) = 0.007269650627029751;     
        fTheta(0,1) = 0.0005732929360778612;      
        fTheta(1,0) = 0.0005423113938742861;      
        fTheta(1,1) = 0.0024899188091627275;
    }
    else if (index == -6)
    {
        fTheta(0,0) = 0.00962344634232376;      
        fTheta(0,1) = 0.0008821551029096821;      
        fTheta(1,0) = 0.0008173876088002891;      
        fTheta(1,1) = 0.0033314941666023846;
    }
    else if (index == -5)
    {
        fTheta(0,0) = 0.013291908173098855;     
        fTheta(0,1) = 0.0014535026926634537;      
        fTheta(1,0) = 0.001304396542753453;       
        fTheta(1,1) = 0.004686562946021358;
    }
    else if (index == -4)
    {
        fTheta(0,0) = 0.019387497509483147;     
        fTheta(0,1) = 0.0026150326018503645;      
        fTheta(1,0) = 0.00223556139351793;        
        fTheta(1,1) = 0.007071330559410259;
    }
    else if (index == -3)
    {
        fTheta(0,0) = 0.030325254444289344;
        fTheta(0,1) = 0.005260014412300331;       
        fTheta(1,0) = 0.004200878445320969;       
        fTheta(1,1) = 0.01182847495883381;
    }
    else if (index == -2)
    {
        fTheta(0,0) = 0.05180629471115681;      
        fTheta(0,1) = 0.012148092425428119;       
        fTheta(1,0) = 0.009012756593462234;       
        fTheta(1,1) = 0.02328469710065523;
    }
    else if (index == -1)
    {
        fTheta(0,0) = 0.09874893628088657;      
        fTheta(0,1) = 0.03329133395196183;        
        fTheta(1,0) = 0.0250514088956815;         
        fTheta(1,1) = 0.06157738592636425;
    }
    else if (index == 0)
    {
        fTheta(0,0) = 0.21806189783187474;      
        fTheta(0,1) = 0.11613100199876912;        
        fTheta(1,0) = 0.17119435052380821;        
        fTheta(1,1) = 0.36844107634396706;
    }
    else if (index == 1)
    {
        fTheta(0,0) = 0.21806189783187455;      
        fTheta(0,1) = -0.11613100199876869;      
        fTheta(1,0) = -0.17119435052380752;       
        fTheta(1,1) = 0.368441076343967;
    }
    else if (index == 2)
    {
        fTheta(0,0) = 0.09874893628088677;      
        fTheta(0,1) = -0.0332913339519616;        
        fTheta(1,0) = -0.025051408895681276;      
        fTheta(1,1) = 0.061577385926364386;
    }
    else if (index == 3)
    {
        fTheta(0,0) = 0.051806294711156727;     
        fTheta(0,1) = -0.012148092425427885;      
        fTheta(1,0) = -0.009012756593462038;      
        fTheta(1,1) = 0.023284697100655116;
    }
    else if (index == 4)
    {
        fTheta(0,0) = 0.030325254444289386;     
        fTheta(0,1) = -0.005260014412300128;      
        fTheta(1,0) = -0.004200878445320707;      
        fTheta(1,1) = 0.011828474958833792;
    }
    else if (index == 5)
    {
        fTheta(0,0) = 0.019387497509483202;     
        fTheta(0,1) = -0.002615032601850154;      
        fTheta(1,0) = -0.0022355613935177722;     
        fTheta(1,1) = 0.007071330559410297;
    }
    else if (index == 6)
    {
        fTheta(0,0) = 0.0132919081730988;       
        fTheta(0,1) = -0.0014535026926632412;     
        fTheta(1,0) = -0.0013043965427532564;     
        fTheta(1,1) = 0.004686562946021337;
    }
    else if (index == 7)
    {
        fTheta(0,0) = 0.009623446342323753;     
        fTheta(0,1) = -0.0008821551029095065;     
        fTheta(1,0) = -0.0008173876088001;        
        fTheta(1,1) = 0.003331494166602315;
    }
    else if (index == 8)
    {
        fTheta(0,0) = 0.007269650627029748;     
        fTheta(0,1) = -0.0005732929360776939;     
        fTheta(1,0) = -0.0005423113938741172;     
        fTheta(1,1) = 0.0024899188091627275;
    }
    else if (index == 9)
    {
        fTheta(0,0) = 0.005677925334832266;     
        fTheta(0,1) = -0.00039296875546772445;    
        fTheta(1,0) = -0.00037681903643516497;    
        fTheta(1,1) = 0.0019317920495601868;
    }
    else if (index == 10)
    {
        fTheta(0,0) = 0.004554387859895136;     
        fTheta(0,1) = -0.0002809410552780244;     
        fTheta(1,0) = -0.000271886044021965;      
        fTheta(1,1) = 0.0015426573553735828;
    }
    else
        throw ExceptionT::kGeneralFail;

    return fTheta;
}

#endif /* BRIDGING_ELEMENT */
