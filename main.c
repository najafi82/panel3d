#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <windows.h>
#include <math.h>
#include <string.h>
#include <omp.h>

/// /////Laspack Library///////////
//#include "./laspack/laspack.h"

#include "mainheader.h"
#include "matrix.h"
#include "general.h"
#include "mesh.h"
#include "discretization.h"
#include "exporter.h"
#include "rbe.h"
#include "fishbodygenerator.h"


int main()
{
    int i,j;
    struct MyVector PressForce, ViscsForce, PressMoment, ViscsMoment;
    double ConsumedPower=0.0;

    PressForce.X=PressForce.Y=PressForce.Z=0.0;
    PressMoment.X=PressMoment.Y=PressMoment.Z=0.0;
    ViscsForce.X=ViscsForce.Y=ViscsForce.Z=0.0;
    ViscsMoment.X=ViscsMoment.Y=ViscsMoment.Z=0.0;

//    inf_Vel.X=0.2/dt;

////    QMatrix LaspackMC;
//    Vector ST,Fee;
    int dInt;

    MAX_SIM_STEP=ceil((FinalTime-StartTime)/dt);

    printf("MAX_SIMULATION_STEP: %d\n",MAX_SIM_STEP);



    SimStep=1;
    Grad.X=Grad.Y=Grad.Z=0.0;

    DynamicParameters(&Grids);
    KinematicParameters(&Grids,0.0);


/// /// FISH BODY GENERATOR /// ///
//    sprintf(InPutFile ,"test2");


    sprintf(InPutFile ,"N0012_AOA9.9");
    sprintf(InPutFile ,"0005");
    sprintf(InPutFile ,"N0012_AOA10");
    sprintf(InPutFile ,"0012.v03.inclined");

    sprintf(InPutFile ,"4412_AOA7");
    sprintf(InPutFile ,"4412_AOA5");
    sprintf(InPutFile ,"4412");

    sprintf(InPutFile ,"st0.17358");



    sprintf(InPutFile ,"4412_AOA8_FINE");


    sprintf(InPutFile ,"4412_AOA8");
    sprintf(InPutFile ,"mesh9x9");


    sprintf(InPutFile ,"finwhale");
    sprintf(InPutFile ,"p01_coarse");
    sprintf(InPutFile ,"4384");

    sprintf(InPutFile ,"rotating.cylinder");
    sprintf(InPutFile ,"politis_propeller");
//    sprintf(InPutFile ,"All");

    sprintf(InPutFile ,"KP505");
    sprintf(InPutFile ,"P4119_1Blade");
    sprintf(InPutFile ,"finwhale");
    sprintf(InPutFile ,"0012");
    sprintf(InPutFile ,"P4679");
    sprintf(InPutFile ,"RV50");

//    FODGen(&Grids,InPutFile,36, 1.0, -1./3., -25.*pi/180., -0.2);

//    RoboTunaGenerator(&Grids,InPutFile);
//    GiantDanioGenerator(&Grids,InPutFile);




    ReadMeshAndMakeFirstKuttaStrip(&Grids,InPutFile);

    ///modifyNodePosition(&Grids);
    D_CH=findConvexHullDiameter(&Grids);
    Artificial_Visc=pow(0.05*D_CH,2.);



//    ///
///    moveGeometry(&Grids);
//    ///
    int NN=Grids.NN;
    int NC=Grids.NC; /// No of Bedy Cells
    int NTE=0;
    int MaxNST;

    if (Lifting_Problem==1)
    {
//        CutPanelsInLeadingEdges(&Grids);
        CutPanelsInTrailingEdges(&Grids);
        NTE=Grids.NTE; /// No of Trailing Edges
        ///int NST=Grids.NST; /// No of Segment per Trailing Edges
        MaxNST=Grids.NST[0];
        for(i=1;i<NTE;i++)
        {
            if (MaxNST<Grids.NST[i])
                MaxNST=Grids.NST[i];
        }
    }


    MakeHistoryFile(InPutFile);
    AllocateStaticArrays(&A,&B,&G,&Phi,&Phi_Old,&PhiBar,&PhiBar_Old,
                             &Phi_Nodal,&RHS,&SRC_COEF,&Press,&BdyPurt_U,&BdyPurt_V,&BdyPurt_W,NC,NN,NTE);
    AllocateDynamicArrays(&Miu,&VrtxInd_U,&VrtxInd_V,&VrtxInd_W,NC,NTE, MaxNST);

//    VTKExporter_body(&Grids,&Phi,&Phi_Nodal,&BdyPurt_U,&BdyPurt_V,&BdyPurt_W,&Press,InPutFile,SimStep);

    double start_time, end_time, diff_time;

//    inf_Vel.X=-cos(2.0*pi/180.0);
//    inf_Vel.Y=-sin(2.0*pi/180.0);
//    inf_Vel.Z=0.0;

//    rotatemesh();

//
//    for(i=0;i<4;i++)
//        printf("\nNode[%d]: %d",435, Grids.Cells[435].NodeList[i]);
//
//    getchar();

    for(Time=StartTime;Time<=FinalTime;Time+=dt)
    {
        KinematicParameters(&Grids,Time);
        start_time=GetTickCount();

        printf("_______________________________________________________________________________\n");
        printf("Step= %d\t\tTime= %.4f\n",SimStep,Time);
        Discretize(&Grids,&Phi,&Miu,&A,&B,&SRC_COEF,&RHS,dt,SimStep);


        NoIter=SOR(&B,&RHS,&Phi,NC,1.0e-7,1000,1.0,&error);
        printf("\n############ Error = %e With %d Iterate\n\n",error,NoIter);

        MakeMiuOnKuttaStrip(&Grids,&Miu,&Phi);



//        SetRTCAccuracy(1.0e-7);
//        CopyToLaspackMatrix(&LaspackMC,&B,&Grids);
//        CopyToLaspackVector(&ST,RHS,&Grids);
//        CopyToLaspackVector(&Fee,Phi,&Grids);
////        BiCGSTABIter(&LaspackMC,&Fee,&ST,100,ILUPrecond,1.05);
////        CGIter(&LaspackMC,&Fee,&ST,300,ILUPrecond,1.2);
//        SSORIter(&LaspackMC,&Fee,&ST,300,SSORPrecond,1.2);
//        CopyFromLaspackVector(Phi,&Fee,&Grids);
//        printf("\n############ Error = %e With %d Iterate\n",GetLastAccuracy(),GetLastNoIter());




        /// calculate gradient of potential
//        GlobalGradient(&Grids,&Phi,&BdyPurt_U,&BdyPurt_V,&BdyPurt_W);
//        LSEGradient_Nodal(&Grids,&Phi,&Phi_Nodal,&BdyPurt_U,&BdyPurt_V,&BdyPurt_W);
//        LSEGradient_Ngb(&Grids,&Phi,&BdyPurt_U,&BdyPurt_V,&BdyPurt_W, SimStep);
//        SurfaceGradient_QUAD(&Grids,&Phi,&BdyPurt_U,&BdyPurt_V,&BdyPurt_W);
        gradient_vortexje(&Grids,&Phi,&BdyPurt_U,&BdyPurt_V,&BdyPurt_W);
//        TangentialGradient_Nodal(&Grids,&Phi,&Phi_Nodal,&BdyPurt_U,&BdyPurt_V,&BdyPurt_W);



        /// calc pressure using Bernouli Equation
        /// Compute PhiBar Value, Important for calculation of forces in moving boundaries
        CalcPhiBar(&Grids,&Phi,&PhiBar);
//        CalcPressure(&Grids,&PhiBar,&PhiBar_Old,&BdyPurt_U,&BdyPurt_V,&BdyPurt_W,&Press, SimStep);
        CalcPressureCoefficient_Cp(&Grids,&PhiBar,&PhiBar_Old,&BdyPurt_U,&BdyPurt_V,&BdyPurt_W,&Press, SimStep);
        PressForce=CalcPressureForces(&Grids,&Press);
//        ViscsForce=CalcViscousForces(&Grids,&BdyPurt_U,&BdyPurt_V,&BdyPurt_W);
        PressMoment=CalcPressureMoments(&Grids,&Press);
//        ViscsMoment=CalcViscousMoments(&Grids,&BdyPurt_U,&BdyPurt_V,&BdyPurt_W);
        ConsumedPower=CalcConsumedPower(&Grids,&Press);

        printf("pressure Forces=(%.6f, %.6f, %.6f)\n",PressForce.X,PressForce.Y,PressForce.Z);
        printf("pressure Moment=(%.6f, %.6f, %.6f)\n\n",PressMoment.X,PressMoment.Y,PressMoment.Z);
        WriteForces(&Grids,InPutFile,&PressForce,&ViscsForce,&PressMoment,&ViscsMoment,ConsumedPower);



        /// calculate  Mean Induced Velocity Over Vortex Nodes
        calcMeanVelocityVortexNodes(&Grids,&Phi,&Miu,&VrtxInd_U,&VrtxInd_V,&VrtxInd_W,SimStep);

        UpdateNodePositions(&Grids,SimStep);
        ComputeCellCenters(&Grids,SimStep);
        ComputeCellAreas(&Grids,SimStep);
        ComputeEdgeCenters(&Grids);
        ComputeEdgeNormals(&Grids);

        /// Update Wake Geometry and Topology, Develope Votex Sheet
        developeVortexSheet(&Grids,&VrtxInd_U,&VrtxInd_V,&VrtxInd_W,SimStep);



        ///shift dipole on free wake one strip
        ShiftMiuOnFreeWake(&Grids,&Miu,SimStep);


        /// write geometry data in VTK format for body and wake elements
        VTKExporter_body(&Grids,&Phi,&Phi_Nodal,&BdyPurt_U,&BdyPurt_V,&BdyPurt_W,&Press,InPutFile,SimStep);
        for (i=0;i<NTE;i++)
            VTKExporter_wake(&Grids,i,&Miu,&VrtxInd_U,&VrtxInd_V,&VrtxInd_W,InPutFile,SimStep);


        /// Make Purturbation PHI (-phi-V.r) based on MassCenter
        MakePhiOldAndPhiBarOld(&Grids,&Phi,&Phi_Old,&PhiBar,&PhiBar_Old,NC);


        SimStep++;
        end_time=GetTickCount();
        diff_time=(end_time-start_time)/1000.;
        printf("\nComputational Time: %.3f (s)\n",diff_time);

    }

//    Q_Destr(&LaspackMC);
//    V_Destr(&ST);
//    V_Destr(&Fee);

    return 0;
}
