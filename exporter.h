double NodalValue(int NodeNumber,double **CellValue,struct Grid *Grids)
{
    int j,PCell;
    double NOM=0.0;
    double DENOM=0.0;

    for(j=0;j<(Grids->Nodes[NodeNumber].NumOfSharedCells);j++)
    {
        PCell=Grids->Nodes[NodeNumber].SharedCells[j];
        NOM += (*CellValue)[PCell]*Mag(Grids->Cells[PCell].Area);
        DENOM += Mag(Grids->Cells[PCell].Area);
    }
    return(NOM/DENOM);
}

void MakeHistoryFile(char *Name)
{
    char path[80];
    sprintf(path,"./%s_forces.txt",Name);
    FILE *fid;
    fid=fopen(path,"w");
    fprintf(fid,"Time,PF.X,PF.Y,PF.Z,VF.X,VF.Y,VF.Z,PM.X,PM.Y,PM.Z,VM.X,VM.Y,VM.Z,MC.X,MC.Y,MC.Z,BO.X,BO.Y,BO.Z,POWER\n");
//    fprintf(fid,"Time,PF.X,PF.Y,PF.Z,VF.X,VF.Y,VF.Z,PM.X,PM.Y,PM.Z,VM.X,VM.Y,VM.Z\n");
    fclose(fid);
}

void WriteForces(struct Grid *Grids,char *Name,struct MyVector *PF,struct MyVector *VF,struct MyVector *PM,struct MyVector *VM, double CON_POW)
{
    char path[80];
    sprintf(path,"./%s_forces.txt",Name);
    FILE *fid;
    fid=fopen(path,"a");
    fprintf(fid,"%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
            Time,PF->X,PF->Y,PF->Z,VF->X,VF->Y,VF->Z,PM->X,PM->Y,PM->Z,VM->X,VM->Y,VM->Z
            ,Grids->DynaVariables.NewMassCenter.X,Grids->DynaVariables.NewMassCenter.Y,Grids->DynaVariables.NewMassCenter.Z
            ,Grids->DynaVariables.NewBodyOr.X,Grids->DynaVariables.NewBodyOr.Y,Grids->DynaVariables.NewBodyOr.Z, CON_POW);
    fclose(fid);

}

void VTKExporter_body(struct Grid *Grids,double **fi,double **fi_nodal,double **U,double **V,double **W,double **P,char *Name,int SStep)
{
    int i,NN,NC,TotalCellData=0;
    char path[80];
    FILE *VTK;

    sprintf(path,"./result/%s_body_%d.vtk",Name,SStep);
    if((VTK=fopen(path,"w"))==NULL)
    {
        puts("vtk writing error.................Exiting ( Write2txt 000)");
        exit(-1);
    }

    NN=Grids->NN;
    NC=Grids->NC;

    fprintf(VTK,"# vtk DataFile Version 2.2\n");
    fprintf(VTK,"Unstructured Grid Example\n");
    fprintf(VTK,"ASCII\n");
    fprintf(VTK,"DATASET UNSTRUCTURED_GRID\n");
    fprintf(VTK,"POINTS %d double\n",NN);
    for (i=0;i<NN;i++)
        fprintf(VTK,"%e\t%e\t%e\n",Grids->Nodes[i].Pos.X,Grids->Nodes[i].Pos.Y,Grids->Nodes[i].Pos.Z);


    for (i=0;i<NC;i++)
        TotalCellData+=Grids->Cells[i].NodePerCell;

    TotalCellData+=Grids->NC;

    fprintf(VTK,"CELLS %d %d\n",NC,TotalCellData);
    for (i=0;i<NC;i++)
    {
        fprintf(VTK,"%d ",(Grids->Cells[i].NodePerCell));
        int j;
        for(j=0;j<(Grids->Cells[i].NodePerCell);j++)
            fprintf(VTK,"%d ",Grids->Cells[i].NodeList[j]);
        fprintf(VTK,"\n");
    }

    fprintf(VTK,"CELL_TYPES %d\n",NC);
    for (i=0;i<NC;i++)
    {
        if((Grids->Cells[i].NodePerCell)==3)
            fprintf(VTK,"%d\n",5);
        else
            fprintf(VTK,"%d\n",9);
    }


    fprintf(VTK,"CELL_DATA %d\n",NC);
    fprintf(VTK,"SCALARS body_dipole double\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
    for (i=0;i<NC;i++)
        fprintf(VTK,"%e\n",(*fi)[i]);

    fprintf(VTK,"SCALARS Cp double\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
    for (i=0;i<NC;i++)
        fprintf(VTK,"%e\n",(*P)[i]);

//    fprintf(VTK,"VECTORS Velocity double\n");
//    for (i=0;i<NC;i++)
//        fprintf(VTK,"%e\t%e\t%e\n",(*U)[i],(*V)[i],(*W)[i]);


    fclose(VTK);
    printf("Write to %s_body_%d.vtk Completed...\n",Name,SStep);
}


void VTKExporter_wake(struct Grid *Grids,int TE, double ***M,double ***U,double ***V,double ***W,char *Name,int SStep)
{
    int i,j,k,cnt,NN,NC,TotalCellData=0;
    char path[80];
    FILE *VTK;

    sprintf(path,"./result/%s_wake_%d_%d.vtk",Name,TE,SStep);
    if((VTK=fopen(path,"w"))==NULL)
    {
        puts("vtk writing error.................Exiting ( Write2txt 000)");
        exit(-1);
    }

    NN=((Grids->NST[TE])+1)*(SStep+1);
    NC= (Grids->NST[TE])   * SStep;

    fprintf(VTK,"# vtk DataFile Version 2.2\n");
    fprintf(VTK,"Unstructured Grid Example\n");
    fprintf(VTK,"ASCII\n");
    fprintf(VTK,"DATASET UNSTRUCTURED_GRID\n");
    fprintf(VTK,"POINTS %d double\n",NN);

//    for (j=0;j<Grids->NTE;j++)
        for (i=0;i<NN;i++)
            fprintf(VTK,"%e\t%e\t%e\n", Grids->Wakes[TE].WNodes[i].Pos.X,
                                        Grids->Wakes[TE].WNodes[i].Pos.Y,
                                        Grids->Wakes[TE].WNodes[i].Pos.Z);

    TotalCellData=5*NC;

    fprintf(VTK,"CELLS %d %d\n",NC/**(Grids->NTE)*/,TotalCellData/**(Grids->NTE)*/);

//    for (j=0;j<Grids->NTE;j++)
        for (i=0;i<NC;i++)
        {
            fprintf(VTK,"%d ",4);
            for(k=0;k<4;k++)
                fprintf(VTK,"%d ",Grids->Wakes[TE].WCells[i].NodeList[k]);
            fprintf(VTK,"\n");
        }

    fprintf(VTK,"CELL_TYPES %d\n",NC/**(Grids->NTE)*/);
//    for (j=0;j<Grids->NTE;j++)
        for (i=0;i<NC;i++)
            fprintf(VTK,"%d\n",9);


    fprintf(VTK,"CELL_DATA %d\n",NC/**(Grids->NTE)*/);
    fprintf(VTK,"SCALARS wake_dipole double\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");

//    for (j=0;j<Grids->NTE;j++)
        for (i=0;i<NC;i++)
            fprintf(VTK,"%e\n",(*M)[i][TE]);


    fclose(VTK);
    printf("Write to %s_wake_%d_%d.vtk Completed...\n",Name,TE,SStep);
}

