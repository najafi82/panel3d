void GetToken(char *string, FILE *fpMesh)
{
    char dummychar;
    sprintf(string,"");
    dummychar=getc(fpMesh);
    if (feof(fpMesh)==0)
    {
        while ((dummychar==' ')||(dummychar=='\n'))
          dummychar=getc(fpMesh);
        do
        {
          sprintf(string,"%s%c",string,dummychar);
          dummychar=getc(fpMesh);
        }while ((dummychar!=' ')&&(dummychar!='\n'));
    }
}

void VTKExporter_Check_Area(struct Grid *Grids, char *Name)
{
    int i,NN,NC,TotalCellData=0;
    char path[80];
    FILE *VTK;

    sprintf(path,"./result/%s_CheckMesh.vtk",Name);
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
    fprintf(VTK,"SCALARS X_Area double\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
    for (i=0;i<NC;i++)
        if (Grids->Cells[i].Area.X>0.0)
            fprintf(VTK,"%e\n",1.0);
        else
            fprintf(VTK,"%e\n",-1.0);

    fprintf(VTK,"SCALARS Y_Area double\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
    for (i=0;i<NC;i++)
        if (Grids->Cells[i].Area.Y>0.0)
            fprintf(VTK,"%e\n",1.0);
        else
            fprintf(VTK,"%e\n",-1.0);

    fprintf(VTK,"SCALARS Z_Area double\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
    for (i=0;i<NC;i++)
        if (Grids->Cells[i].Area.Z>0.0)
            fprintf(VTK,"%e\n",1.0);
        else
            fprintf(VTK,"%e\n",-1.0);


    fclose(VTK);
    printf("\nWrite to %s_CheckMesh.vtk Completed...\n",Name);
}



void FillNodes(struct Grid *Grids,char *FileName)//coordinates of nodes has been read
{
    FILE *fpMesh;
    char dummystr[80];
    int i=0,NN;
    if( (fpMesh=fopen(FileName,"r"))==NULL )
    {
        puts("NodeCoord.ebr Reading Error...");
        exit(-1);
    }
    GetToken(dummystr,fpMesh);
    if (Grids->NN!=atoi(dummystr))
    {
        puts("Number of nodes don't match");
        exit(-1);
    }
    for(i=0;i<Grids->NN;i++)
    {
        GetToken(dummystr,fpMesh);
        NN=atoi(dummystr)-1;
        GetToken(dummystr,fpMesh);
        Grids->Nodes[NN].Pos.X=atof(dummystr);
        GetToken(dummystr,fpMesh);
        Grids->Nodes[NN].Pos.Y=atof(dummystr);
        GetToken(dummystr,fpMesh);
        Grids->Nodes[NN].Pos.Z=atof(dummystr);        /// /// /// (*Grids).  -->  Grids->
    }
    fclose(fpMesh);
}

void FillConnectivity(struct Grid *Grids,char *FileName)
{
    FILE *fpMesh;
    char dummystr[80];
    int i=0,j=0;
    if( (fpMesh=fopen(FileName,"r"))==NULL )
    {
        puts("ElemConnect.ebr Reading Error...");
        exit(-1);
    }
    GetToken(dummystr,fpMesh);
    if (Grids->NC!=atoi(dummystr))
    {
        puts("Number of elements don't match\n");
        exit(-1);
    }
   for (i=0;i<Grids->NC;i++)
    {
        GetToken(dummystr,fpMesh);
        Grids->Cells[i].Type=atoi(dummystr);

        GetToken(dummystr,fpMesh);
        Grids->Cells[i].NodePerCell=atoi(dummystr);
        Grids->Cells[i].NodeList=(int *)malloc(sizeof(int)*(Grids->Cells[i].NodePerCell));
        for (j=0;j<(Grids->Cells[i].NodePerCell);j++)
        {
            GetToken(dummystr,fpMesh);
            Grids->Cells[i].NodeList[j]=(atoi(dummystr))-1;
        }
    }
    fclose(fpMesh);
}

void chech_node_order_of_hub_panels(struct Grid *Grids)
{
    int i, type;
    int cnt=0;
    int n0, n1, n2, n3;
    struct MyVector Area, Cent;

    printf("\n\n\n");
    printf("*************************************\n");
    printf("*************************************\n");
    printf("CHECKING NODE ORDER OF HUB PANELS....\n");

    for(i=2700;i<Grids->NC;i++)
    {
        type=Grids->Cells[i].Type;
        switch (type)
        {
            case 2:
            {
                n0=Grids->Cells[i].NodeList[0];
                n1=Grids->Cells[i].NodeList[1];
                n2=Grids->Cells[i].NodeList[2];
                n3=Grids->Cells[i].NodeList[3];

                Area=SVP(0.50,Cross(Minus(Grids->Nodes[n2].Pos, Grids->Nodes[n0].Pos),
                                   Minus(Grids->Nodes[n3].Pos, Grids->Nodes[n1].Pos)));

                Cent=SVP(0.25, Sum( Sum(Grids->Nodes[n0].Pos, Grids->Nodes[n1].Pos),
                                    Sum(Grids->Nodes[n2].Pos, Grids->Nodes[n3].Pos)));

                if(Dot(Area, Cent)<0.0)
                {
                    Grids->Cells[i].NodeList[0]=n0;
                    Grids->Cells[i].NodeList[1]=n3;
                    Grids->Cells[i].NodeList[2]=n2;
                    Grids->Cells[i].NodeList[3]=n1;
                    cnt++;
                }
            }break;

            case 3:
            {
                n0=Grids->Cells[i].NodeList[0];
                n1=Grids->Cells[i].NodeList[1];
                n2=Grids->Cells[i].NodeList[2];

                Area=SVP(0.50 ,Cross(Minus(Grids->Nodes[n1].Pos, Grids->Nodes[n0].Pos)
                                  ,Minus(Grids->Nodes[n2].Pos, Grids->Nodes[n0].Pos)));

                Cent=SVP(1./3., Sum3(Grids->Nodes[n0].Pos, Grids->Nodes[n1].Pos, Grids->Nodes[n2].Pos));

                if(Dot(Area, Cent)<0.0)
                {
                    Grids->Cells[i].NodeList[0]=n0;
                    Grids->Cells[i].NodeList[1]=n2;
                    Grids->Cells[i].NodeList[2]=n1;
                    cnt++;
                }
            }break;
        }
    }
    printf("%d CELLS HAS BEEN RE-ORDERED!!!\n", cnt);
    printf("*************************************\n");
    printf("*************************************\n\n\n");
}

void ReadBoundaryTypes(struct Grid *Grids,char *FileName)
{
    FILE *fpMesh;
    char dummystr[80];
    int pcell,i=0,BdryType=0,NoOfCells=0;
    fpMesh=fopen(FileName,"r");
    if(fpMesh==NULL)
    {
        puts("Boundary.ebr Reading Error...");
        exit(-1);
    }
    do{
        GetToken(dummystr,fpMesh);
        BdryType=atoi(dummystr);
        GetToken(dummystr,fpMesh);
        NoOfCells=atoi(dummystr);
        for (i=0;i<NoOfCells;i++)
        {
            GetToken(dummystr,fpMesh);
            pcell=atoi(dummystr)-1;
            Grids->Cells[pcell].BType=BdryType;
            Grids->Cells[pcell].Side=(BdryType==1)?('b'):('f');
        }
    }while(feof(fpMesh)==0);
    fclose(fpMesh);
}

struct MyVector GetEdgeCoincide(struct Grid *Grids,int PCell,int NodeNum,double alfa)
{
    struct MyVector Result;//,MidPoint,MedVect,MedVectNor;
    struct MyVector CG,CP,PG;
//    int NodeA,NodeB;
//    double MedVectMag;

    CG=Grids->Cells[PCell].CellCenter;
    CP=Grids->Nodes[Grids->Cells[PCell].NodeList[NodeNum]].Pos;
    PG=Minus(CG,CP);
    Result=Sum(CP,SVP(alfa,PG));

//    if (NodeNum==0)
//    {
//        NodeA=1;
//        NodeB=2;
//    }
//    else if (NodeNum==1)
//    {
//        NodeA=0;
//        NodeB=2;
//    }
//    else
//    {
//        NodeA=0;
//        NodeB=1;
//    }
//    MidPoint=SVP(0.5,Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[NodeA]].Pos,
//                         Grids->Nodes[Grids->Cells[PCell].NodeList[NodeB]].Pos));
//
//    MedVect=Minus(MidPoint,Grids->Nodes[Grids->Cells[PCell].NodeList[NodeNum]].Pos);
//
//    MedVectMag=Mag(MedVect);
//
//    MedVectNor=Normalized(MedVect);
//
//    Result.X=(Grids->Nodes[Grids->Cells[PCell].NodeList[NodeNum]].Pos.X)+(alfa*MedVectNor.X*MedVectMag);
//    Result.Y=(Grids->Nodes[Grids->Cells[PCell].NodeList[NodeNum]].Pos.Y)+(alfa*MedVectNor.Y*MedVectMag);
//    Result.Z=(Grids->Nodes[Grids->Cells[PCell].NodeList[NodeNum]].Pos.Z)+(alfa*MedVectNor.Z*MedVectMag);

    return Result;
}

void  ComputeCellInfo(struct Grid *Grids)
{
    int PCell,i;
////    double Area1,Area2;
////    struct MyVector FaceCenter1,FaceCenter2;

    for(PCell=0;PCell<Grids->NC;PCell++)
    {
//        Grids->Cells[PCell].EdgeCoincide=(struct MyVector*)malloc(sizeof(struct MyVector)*(Grids->Cells[PCell].NodePerCell));
        Grids->Cells[PCell].EdgeNormal=(struct MyVector*)malloc(sizeof(struct MyVector)*(Grids->Cells[PCell].NodePerCell));
        Grids->Cells[PCell].EdgeCenter=(struct MyVector*)malloc(sizeof(struct MyVector)*(Grids->Cells[PCell].NodePerCell));
        Grids->Cells[PCell].Ngb=(int*)malloc(sizeof(int)*(Grids->Cells[PCell].NodePerCell));

        switch(Grids->Cells[PCell].Type)
        {
            case 2: /// Quadrilateral
            {

                Grids->Cells[PCell].CellCenter=SVP(0.25,Sum(Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos,
                                                                Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos),
                                                            Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos,
                                                                Grids->Nodes[Grids->Cells[PCell].NodeList[3]].Pos)));

                Grids->Cells[PCell].Area=SVP(0.5,Cross(Minus(Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos,
                                                             Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos),
                                                       Minus(Grids->Nodes[Grids->Cells[PCell].NodeList[3]].Pos,
                                                             Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos)));

                Grids->Cells[PCell].EdgeCenter[0]=SVP(0.5,Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos,
                                                              Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos));
                Grids->Cells[PCell].EdgeCenter[1]=SVP(0.5,Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos,
                                                              Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos));
                Grids->Cells[PCell].EdgeCenter[2]=SVP(0.5,Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos,
                                                              Grids->Nodes[Grids->Cells[PCell].NodeList[3]].Pos));
                Grids->Cells[PCell].EdgeCenter[3]=SVP(0.5,Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[3]].Pos,
                                                              Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos));

                Grids->Cells[PCell].EdgeNormal[0]=Cross(Normalized(Grids->Cells[PCell].Area),
                                                                    Minus(Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos,
                                                                          Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos));
                Grids->Cells[PCell].EdgeNormal[1]=Cross(Normalized(Grids->Cells[PCell].Area),
                                                                    Minus(Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos,
                                                                          Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos));
                Grids->Cells[PCell].EdgeNormal[2]=Cross(Normalized(Grids->Cells[PCell].Area),
                                                                    Minus(Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos,
                                                                          Grids->Nodes[Grids->Cells[PCell].NodeList[3]].Pos));
                Grids->Cells[PCell].EdgeNormal[3]=Cross(Normalized(Grids->Cells[PCell].Area),
                                                                    Minus(Grids->Nodes[Grids->Cells[PCell].NodeList[3]].Pos,
                                                                          Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos));
            }break;

            case 3: ///Triangle
            {
                Grids->Cells[PCell].CellCenter=SVP(1.0/3.0, Sum3(Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos,
                                                                 Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos,
                                                                 Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos));

                Grids->Cells[PCell].Area=SVP(0.5,Cross(Minus(Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos,
                                                             Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos)
                                                      ,Minus(Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos,
                                                             Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos)));

                Grids->Cells[PCell].EdgeCenter[0]=SVP(0.5,Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos,
                                                              Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos));
                Grids->Cells[PCell].EdgeCenter[1]=SVP(0.5,Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos,
                                                              Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos));
                Grids->Cells[PCell].EdgeCenter[2]=SVP(0.5,Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos,
                                                              Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos));

                Grids->Cells[PCell].EdgeNormal[0]=Cross(Normalized(Grids->Cells[PCell].Area),
                                                                    Minus(Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos,
                                                                          Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos));
                Grids->Cells[PCell].EdgeNormal[1]=Cross(Normalized(Grids->Cells[PCell].Area),
                                                                    Minus(Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos,
                                                                          Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos));
                Grids->Cells[PCell].EdgeNormal[2]=Cross(Normalized(Grids->Cells[PCell].Area),
                                                                    Minus(Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos,
                                                                          Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos));
            }break;
        }
    }
}

int isContain(int a,int *aList,int s)
{
    int i;
    for (i=0;i<s;i++)
        if (a==aList[i])
            return 1;
    return 0;
}

void SetIntSVect(int *Vect,int r,int val)
{
  int i;
  for(i=0;i<r;i++)
    Vect[i]=val;
}

void FillNodeSharing(struct Grid *Grids)
{
    int i,j,k,N;
    for (i=0;i<Grids->NN;i++)
        Grids->Nodes[i].NumOfSharedCells=0;

    for (i=0;i<Grids->NC;i++)
        for (j=0;j<(Grids->Cells[i].NodePerCell);j++)
            Grids->Nodes[Grids->Cells[i].NodeList[j]].NumOfSharedCells+=1;

    for (i=0;i<Grids->NN;i++)
    {
        Grids->Nodes[i].SharedCells=(int *)malloc((Grids->Nodes[i].NumOfSharedCells)*sizeof(int));
        SetIntSVect(Grids->Nodes[i].SharedCells,Grids->Nodes[i].NumOfSharedCells,-1);
    }

    for (i=0;i<Grids->NC;i++)
        for (j=0;j<(Grids->Cells[i].NodePerCell);j++)
        {
            N=Grids->Cells[i].NodeList[j];
            for (k=0;k<Grids->Nodes[N].NumOfSharedCells;k++)
                if (Grids->Nodes[N].SharedCells[k]==-1)
                    break;
            Grids->Nodes[N].SharedCells[k]=i;
        }
}

void FindNeighbours(struct Grid *Grids)
{
    int PCell,i,j,l,k;

    for (PCell=0;PCell<Grids->NC;PCell++)
        for (j=0;j<Grids->Cells[PCell].NodePerCell;j++)
            Grids->Cells[PCell].Ngb[j]=-1;

    int FaceNodeList_1[2];

    for (i=0;i<Grids->NC;i++)
    {
        for (j=0;j<Grids->Cells[i].NodePerCell;j++)
            if (Grids->Cells[i].Ngb[j]==-1)
            {
                for (k=0;k<2;k++)
                {
                    if (Grids->Cells[i].Type==3)
                        FaceNodeList_1[k]=Grids->Cells[i].NodeList[TRI_NODES[j][k+1]];
                    else
                        FaceNodeList_1[k]=Grids->Cells[i].NodeList[QUA_NODES[j][k+1]];
                }


                for (l=0;l<Grids->NC;l++)
                    if (i!=l)
                    {
                        int flag=1;
                        for (k=0;k<2;k++)
                            if (!isContain(FaceNodeList_1[k],Grids->Cells[l].NodeList,Grids->Cells[l].NodePerCell))
                                {
                                    flag=0;
                                    break;
                                }
                        if ((flag==1)/*&&(Grids->Cells[i].Side==Grids->Cells[l].Side)*/)
                            Grids->Cells[i].Ngb[j]=l;
                    }
            }
        printf("\r");
        printf("%6.3f Percent Completed...",((double)(i+1)/(double)Grids->NC)*100.);
    }
    printf("\n\n");
}

void SaveNeighbours(struct Grid *Grids,char *FileName)
{

    FILE *fpMesh;
    int i,j;
    if( (fpMesh=fopen(FileName,"w"))==NULL )
    {
        puts("Error Writing The ElemConnect.ebr...");
        exit(-1);
    }

    for (i=0;i<Grids->NC;i++)
    {
        for (j=0;j<(Grids->Cells[i].NodePerCell);j++)
            fprintf(fpMesh," %8d",Grids->Cells[i].Ngb[j]);
        fprintf(fpMesh,"\n");
    }
    fclose(fpMesh);
}

void ReadNeighbours(struct Grid *Grids,char *FileName)
{
    FILE *fpMesh;
    int i,j;
    char dummystr[80];
    if( (fpMesh=fopen(FileName,"r"))==NULL )
    {
        puts("Error Reading The ElemConnect.ebr...");
        exit(-1);
    }
    for (i=0;i<Grids->NC;i++)
        for (j=0;j<(Grids->Cells[i].NodePerCell);j++)
        {
            GetToken(dummystr,fpMesh);
            Grids->Cells[i].Ngb[j]=atoi(dummystr);
        }
    fclose(fpMesh);
}

void developeVortexSheet(struct Grid *Grids,double ***iu,double ***iv,double ***iw,
                         /*,double ***iuold,double ***ivold,double ***iwold,*/ int SStep)
{
    int i,j,k,l,cnt;
    int NST,NTE,kk;
    struct MyVector Ind_Velocity,Ind_Velocity_old;
    struct MyVector rVec,vBT,vREF;

    int UCell, LCell;
    struct MyVector UCellCen,LCellCen, Node0, Node1, MeanCell, MeanEdge, DirVec;


    double x0, y0, z0, x1, y1, z1;
//    NST=Grids->NST;
//    NTE=Grids->NTE;

    for (i=0;i<Grids->NTE;i++)
        for(j=((SStep+2)*(Grids->NST[i]+1)-1);j>=2*(Grids->NST[i]+1);j--)
        {
            Ind_Velocity.X=(*iu)[j-Grids->NST[i]-1-Grids->NST[i]-1][i];
            Ind_Velocity.Y=(*iv)[j-Grids->NST[i]-1-Grids->NST[i]-1][i];
            Ind_Velocity.Z=(*iw)[j-Grids->NST[i]-1-Grids->NST[i]-1][i];
//            if(j>=2*(NST+1) && j<=3*(NST+1))
                Grids->Wakes[i].WNodes[j].Pos=Sum(Grids->Wakes[i].WNodes[j-Grids->NST[i]-1].Pos,SVP(dt,Ind_Velocity));
//            else
//                Grids->Wakes[i].WNodes[j].Pos=Sum(SVP(0.5,Sum(Grids->Wakes[i].WNodes[j-NST-1].Pos,Grids->Wakes[i].WNodes[j-NST-1-NST-1].Pos)),SVP(1.5*dt,Ind_Velocity));
        }

    int ucell, lcell;
    struct MyVector nv1, nv2, nd3;
    for (i=0;i<Grids->NTE;i++)
        for( j=2*Grids->NST[i]+1 ; j>=Grids->NST[i]+1 ; j--)
        {
//            rVec=Minus(Grids->Wakes[i].WNodes[j-Grids->NST[i]-1].Pos,Grids->DynaVariables.MassCenter);
//            vBT=Sum(Grids->DynaVariables.LinVel,Cross(Grids->DynaVariables.AngVel,rVec));
            vBT=getPointKinematicVelocity(Grids, Grids->Wakes[i].WNodes[j-Grids->NST[i]-1].Pos);
            vREF=Minus(inf_Vel,vBT);
            Grids->Wakes[i].WNodes[j].Pos=Sum(Grids->Wakes[i].WNodes[j-Grids->NST[i]-1].Pos,SVP(dt,vREF));



//            if (j==(2*Grids->NST[i]+1))
//            {
//                ucell=Grids->Wakes[i].UpperCells[j-Grids->NST[i]-2];
//                lcell=Grids->Wakes[i].LowerCells[j-Grids->NST[i]-2];
//            }
//            else
//            {
//                ucell=Grids->Wakes[i].UpperCells[j-Grids->NST[i]-1];
//                lcell=Grids->Wakes[i].LowerCells[j-Grids->NST[i]-1];
//            }
//
//            nv1=Normalized(Grids->Cells[ucell].Area);
//            nv2=Normalized(Grids->Cells[lcell].Area);
//            nd3=Normalized(Sum(nv1, nv2));
//
//            Grids->Wakes[i].WNodes[j].Pos=Sum(Grids->Wakes[i].WNodes[j-Grids->NST[i]-1].Pos,SVP(0.25*dt*Mag(inf_Vel),nd3));






//            x0=Grids->Wakes[i].WNodes[j-Grids->NST[i]-1].Pos.X;
//            y0=Grids->Wakes[i].WNodes[j-Grids->NST[i]-1].Pos.Y;
//            z0=Grids->Wakes[i].WNodes[j-Grids->NST[i]-1].Pos.Z;
//
//            x1=x0 - 0.008;
//            y1=sqrt(y0*y0+z0*z0) * cos(atan2(z0, y0) - 0.05);
//            z1=sqrt(y0*y0+z0*z0) * sin(atan2(z0, y0) - 0.05);
//
//            Grids->Wakes[i].WNodes[j].Pos.X=x1;
//            Grids->Wakes[i].WNodes[j].Pos.Y=y1;
//            Grids->Wakes[i].WNodes[j].Pos.Z=z1;


        }













    /// Making Cell Conectivity for the last new row
    for (i=0;i<Grids->NTE;i++)
    {
        kk=0;
        for (j=SStep*Grids->NST[i];j<(SStep+1)*Grids->NST[i];j++)
        {
            Grids->Wakes[i].WCells[j].NodeList[0]= SStep*(Grids->NST[i]+1)+ kk;
            Grids->Wakes[i].WCells[j].NodeList[1]= SStep*(Grids->NST[i]+1)+ kk +1;
            Grids->Wakes[i].WCells[j].NodeList[2]=(SStep+1)*(Grids->NST[i]+1)+ kk +1;
            Grids->Wakes[i].WCells[j].NodeList[3]=(SStep+1)*(Grids->NST[i]+1)+ kk;
            kk++;
        }
    }

};

//void developeVortexSheet(struct Grid *Grids,double ***iu,double ***iv,double ***iw
//                         ,double ***iuold,double ***ivold,double ***iwold, int SStep)
//{
//    int i,j,k,l,cnt;
//    int NST,NTE,kk;
//    struct MyVector Ind_Velocity,Ind_Velocity_old;
//    struct MyVector rVec,rVec_old,vBT,vBT_old,vREF,vREF_old;
//
//    NST=Grids->NST;
//    NTE=Grids->NTE;
//
//    for (i=0;i<NTE;i++)
//    {
//        for(j=((SStep+2)*(NST+1)-1);j>=(NST+1);j--)
//        {
//            l=(int) fmod(j,(NST+1));
//            Ind_Velocity.X=(*iu)[j-NST-1][i];
//            Ind_Velocity.Y=(*iv)[j-NST-1][i];
//            Ind_Velocity.Z=(*iw)[j-NST-1][i];
//            Ind_Velocity_old.X=(*iuold)[j-NST-1][i];
//            Ind_Velocity_old.Y=(*ivold)[j-NST-1][i];
//            Ind_Velocity_old.Z=(*iwold)[j-NST-1][i];
//
//            rVec=Minus(Grids->Wakes[i].WNodes[j-NST-1].Pos,Grids->DynaVariables.MassCenter);
//            vBT=Sum3(MinusVec(Grids->DynaVariables.LinVel),MinusVec(Cross(Grids->DynaVariables.AngVel,rVec)),inf_Vel);
//            vREF=Sum(Ind_Velocity,vBT);
//
//            if (j>=(NST+1)&&j<=(2*NST+1))
//            {
//                Grids->Wakes[i].WNodes[j].Pos=Sum(SVP(dt,vREF),Grids->Wakes[i].WNodes[j-NST-1].Pos);
//            }
//            else
//            {
//                rVec_old=Minus(Grids->Wakes[i].WNodes[j-NST-1-NST-1].Pos,Grids->DynaVariables.MassCenter);
//                vBT_old=Sum3(MinusVec(Grids->DynaVariables.LinVel),MinusVec(Cross(Grids->DynaVariables.AngVel,rVec_old)),inf_Vel);
//                vREF_old=Sum(Ind_Velocity_old,vBT_old);
//                Grids->Wakes[i].WNodes[j].Pos=Sum(SVP(0.5,Sum(Grids->Wakes[i].WNodes[j-NST-1].Pos,Grids->Wakes[i].WNodes[j-NST-1-NST-1].Pos)),
//                                                   SVP(0.5*dt,Sum(vREF,SVP(2.0,vREF_old))));
//            }
//        }
//    }
//
//    /// Making Cell Conectivity for the last new row
//    for (i=0;i<NTE;i++)
//    {
//        kk=0;
//        for (j=SStep*NST;j<(SStep+1)*NST;j++)
//        {
//            Grids->Wakes[i].WCells[j].NodeList[0]= SStep*(NST+1)+ kk;
//            Grids->Wakes[i].WCells[j].NodeList[1]= SStep*(NST+1)+ kk +1;
//            Grids->Wakes[i].WCells[j].NodeList[2]=(SStep+1)*(NST+1)+ kk +1;
//            Grids->Wakes[i].WCells[j].NodeList[3]=(SStep+1)*(NST+1)+ kk;
//            kk++;
//        }
//    }
//
//};

//void developeVortexSheet(struct Grid *Grids,double ***iu,double ***iv,double ***iw
//                         ,double ***iuold,double ***ivold,double ***iwold, int SStep)
//{
//    int i,j,k,l,cnt;
//    int NST,NTE,kk;
//    struct MyVector Ind_Velocity,Ind_Velocity_old;
//    struct MyVector rVec,vBT,vREF;
//
//    NST=Grids->NST;
//    NTE=Grids->NTE;
//
//    for (i=0;i<NTE;i++)
//        for(j=((SStep+2)*(NST+1)-1);j>=1*(NST+1);j--) /// in 1*( ghablan 2*( bood
//        {
//            l=(int) fmod(j,(NST+1));
//            Ind_Velocity.X=(*iu)[j-NST-1][i];
//            Ind_Velocity.Y=(*iv)[j-NST-1][i];
//            Ind_Velocity.Z=(*iw)[j-NST-1][i];
//            Ind_Velocity_old.X=(*iuold)[j-NST-1][i];
//            Ind_Velocity_old.Y=(*ivold)[j-NST-1][i];
//            Ind_Velocity_old.Z=(*iwold)[j-NST-1][i];
//
//            rVec=Minus(Grids->Wakes[i].WNodes[l].Pos,Grids->DynaVariables.MassCenter);
//            vBT=Sum3(MinusVec(Grids->DynaVariables.LinVel),MinusVec(Cross(Grids->DynaVariables.AngVel,rVec)),inf_Vel);
//            vREF=Sum(Ind_Velocity,vBT);
//            if (j>=(NST+1)&&j<=(2*NST+1))
//            {
//                Grids->Wakes[i].WNodes[j].Pos=Sum(SVP(dt,vREF),Grids->Wakes[i].WNodes[j-NST-1].Pos);
//            }
//            else
//            {
//                Grids->Wakes[i].WNodes[j].Pos=Sum3(SVP(2.0*dt/3.0,vREF),
//                                                   SVP(4.0/3.0,Grids->Wakes[i].WNodes[j-NST-1].Pos),
//                                                   SVP(-1.0/3.0,Grids->Wakes[i].WNodes[j-NST-1-NST-1].Pos));
//            }
//        }
//
//    /// Making Cell Conectivity for the last new row
//    for (i=0;i<NTE;i++)
//    {
//        kk=0;
//        for (j=SStep*NST;j<(SStep+1)*NST;j++)
//        {
//            Grids->Wakes[i].WCells[j].NodeList[0]= SStep*(NST+1)+ kk;
//            Grids->Wakes[i].WCells[j].NodeList[1]= SStep*(NST+1)+ kk +1;
//            Grids->Wakes[i].WCells[j].NodeList[2]=(SStep+1)*(NST+1)+ kk +1;
//            Grids->Wakes[i].WCells[j].NodeList[3]=(SStep+1)*(NST+1)+ kk;
//            kk++;
//        }
//    }
//
//};



///  under construction
void makeKuttaStrip(struct Grid *Grids, int SStep)
{
    int i,j,k,P0,P1,E0,E1,RightOrder;
    int PresentNode,NeigbNode01,NeigbNode02,NeigbNode03,NeigbNode04;
    int kk;
    struct MyVector rVec,vBT, vREF;

    double x0, y0, z0, x1, y1, z1;

    //NST=Grids->NST;

    int A[200][2*Grids->NTE];
    int cnt_cell=0, cnt_edge=1;
    for (i=0;i<(Grids->NTE);i++)
    {
        for(j=0;j<Grids->NST[i];j++)
        {
            A[j][cnt_cell]=Grids->TE[i].CellNo[j];
            A[j][cnt_edge]=Grids->TE[i].EdgeNo[j];
        }
        cnt_cell+=2;
        cnt_edge+=2;
    }

    for (i=0;i<(Grids->NTE);i++)
    {
        for(j=0;j<Grids->NST[i];j++)
            printf("%d\t%d\n",A[j][2*i],A[j][2*i+1]);

        printf("\n");
    }

    cnt_cell=0;
    cnt_edge=1;
    for (i=0;i<(Grids->NTE);i++)
    {
        P0=A[0][cnt_cell];
        E0=A[0][cnt_edge];

        P1=A[1][cnt_cell];
        E1=A[1][cnt_edge];

        PresentNode=Grids->Cells[P0].NodeList[E0];
        if (Grids->Cells[P1].NodePerCell==4)
        {
            NeigbNode01=Grids->Cells[P1].NodeList[0];
            NeigbNode02=Grids->Cells[P1].NodeList[1];
            NeigbNode03=Grids->Cells[P1].NodeList[2];
            NeigbNode04=Grids->Cells[P1].NodeList[3];
            if ((PresentNode==NeigbNode01)||(PresentNode==NeigbNode02)||(PresentNode==NeigbNode03)||(PresentNode==NeigbNode04))
                RightOrder=0;
            else
                RightOrder=1;
        }
        else
        {
            NeigbNode01=Grids->Cells[P1].NodeList[0];
            NeigbNode02=Grids->Cells[P1].NodeList[1];
            NeigbNode03=Grids->Cells[P1].NodeList[2];
            if ((PresentNode==NeigbNode01)||(PresentNode==NeigbNode02)||(PresentNode==NeigbNode03))
                RightOrder=0;
            else
                RightOrder=1;
        }

        for(j=0;j<Grids->NST[i];j++)
        {
            P0=A[j][cnt_cell];
            E0=A[j][cnt_edge];
            if (RightOrder==1)
                Grids->Wakes[i].WNodes[j].Pos=Grids->Nodes[Grids->Cells[P0].NodeList[E0]].Pos;
            else
                Grids->Wakes[i].WNodes[j].Pos=(Grids->Cells[P0].NodePerCell==4)
                                    ?(Grids->Nodes[Grids->Cells[P0].NodeList[QUA_NODES[E0][2]]].Pos)
                                    :(Grids->Nodes[Grids->Cells[P0].NodeList[TRI_NODES[E0][2]]].Pos);
        }
        if (RightOrder==1)
            Grids->Wakes[i].WNodes[j].Pos=(Grids->Cells[P0].NodePerCell==4)
                                    ?(Grids->Nodes[Grids->Cells[P0].NodeList[QUA_NODES[E0][2]]].Pos)
                                    :(Grids->Nodes[Grids->Cells[P0].NodeList[TRI_NODES[E0][2]]].Pos);
        else
            Grids->Wakes[i].WNodes[j].Pos=Grids->Nodes[Grids->Cells[P0].NodeList[E0]].Pos;

        cnt_cell+=2;
        cnt_edge+=2;
    }

    /// Move Trailing Edge Point for the First Time
    /// OLD METHOD

    int ucell, lcell;
    struct MyVector nv1, nv2, nd3;
    for (i=0;i<(Grids->NTE);i++)
        for(j=(Grids->NST[i]+1);j<=(2*Grids->NST[i]+1);j++)
        {

//            rVec=Minus(Grids->Wakes[i].WNodes[j-Grids->NST[i]-1].Pos, Grids->DynaVariables.MassCenter);
//            vBT=Minus(inf_Vel,Sum(Grids->DynaVariables.LinVel,Cross(Grids->DynaVariables.AngVel,rVec)));
            vBT=getPointKinematicVelocity(Grids, Grids->Wakes[i].WNodes[j-Grids->NST[i]-1].Pos);
            vREF=Minus(inf_Vel,vBT);

            Grids->Wakes[i].WNodes[j].Pos=Sum(Grids->Wakes[i].WNodes[j-Grids->NST[i]-1].Pos,SVP(dt,vREF));



//            if (j==(2*Grids->NST[i]+1))
//            {
//                ucell=Grids->Wakes[i].UpperCells[j-Grids->NST[i]-2];
//                lcell=Grids->Wakes[i].LowerCells[j-Grids->NST[i]-2];
//            }
//            else
//            {
//                ucell=Grids->Wakes[i].UpperCells[j-Grids->NST[i]-1];
//                lcell=Grids->Wakes[i].LowerCells[j-Grids->NST[i]-1];
//            }
//
//            nv1=Normalized(Grids->Cells[ucell].Area);
//            nv2=Normalized(Grids->Cells[lcell].Area);
//            nd3=Normalized(Sum(nv1, nv2));
//
//            Grids->Wakes[i].WNodes[j].Pos=Sum(Grids->Wakes[i].WNodes[j-Grids->NST[i]-1].Pos,SVP(0.25*dt*Mag(inf_Vel),nd3));







//            x0=Grids->Wakes[i].WNodes[j-Grids->NST[i]-1].Pos.X;
//            y0=Grids->Wakes[i].WNodes[j-Grids->NST[i]-1].Pos.Y;
//            z0=Grids->Wakes[i].WNodes[j-Grids->NST[i]-1].Pos.Z;
//
//            x1=x0 - 0.008;
//            y1=sqrt(y0*y0+z0*z0) * cos(atan2(z0, y0) - 0.05);
//            z1=sqrt(y0*y0+z0*z0) * sin(atan2(z0, y0) - 0.05);
//
//            Grids->Wakes[i].WNodes[j].Pos.X=x1;
//            Grids->Wakes[i].WNodes[j].Pos.Y=y1;
//            Grids->Wakes[i].WNodes[j].Pos.Z=z1;

        }



    /// Making Cell Conectivity for the first Row (Kutta Strip)
    for (i=0;i<(Grids->NTE);i++)
    {
        kk=0;
        for (j=((SStep-1)*Grids->NST[i]);j<=(((SStep-1)*Grids->NST[i])+(Grids->NST[i]-1));j++)
        {
            Grids->Wakes[i].WCells[j].NodeList[0]= (SStep-1)*(Grids->NST[i]+1)+ kk;
            Grids->Wakes[i].WCells[j].NodeList[1]= (SStep-1)*(Grids->NST[i]+1)+ kk +1;
            Grids->Wakes[i].WCells[j].NodeList[2]=  SStep   *(Grids->NST[i]+1)+ kk +1;
            Grids->Wakes[i].WCells[j].NodeList[3]=  SStep   *(Grids->NST[i]+1)+ kk;
            kk++;
        }
    }


    int n0=Grids->Wakes[0].WCells[0].NodeList[0];
    int n1=Grids->Wakes[0].WCells[0].NodeList[1];
    int n2=Grids->Wakes[0].WCells[0].NodeList[2];
    int n3=Grids->Wakes[0].WCells[0].NodeList[3];
    struct MyVector wakeNormal=SVP(0.5,Cross(Minus(Grids->Wakes[0].WNodes[n2].Pos,Grids->Wakes[0].WNodes[n0].Pos),
                                             Minus(Grids->Wakes[0].WNodes[n3].Pos,Grids->Wakes[0].WNodes[n1].Pos)));

    struct MyVector bodyNormal=Grids->Cells[Grids->Wakes[0].UpperCells[0]].Area;

    double checkDot=Dot(wakeNormal,bodyNormal);

    if (checkDot>0.0)
    {
        SimilarityOfWakeAndBodyNormals=1;
        printf("NORMAL VECTOR OF UPPER PANELS AND KUTTA STRIP ARE ALIGNED, RIGHT ORDER FLAG: %d\n",RightOrder);
        printf("EVERYTHING IS OK, PRESS ANY KEY TO CONTINUE...\n");
        getchar();
    }
    else
    {
        SimilarityOfWakeAndBodyNormals=0;
        printf("NORMAL VECTOR OF UPPER PANELS AND KUTTA STRIP ARE NOT ALIGNED, RIGHT ORDER FLAG: %d\n",RightOrder);
        printf("TO CORRECT THE MISTAKE, SIMILARITY FLAG WOULD BE SET TO ZERO...\n");
        printf("PRESS ANY KEY TO CONTINUE...\n");
        getchar();
    }
    //SimilarityOfWakeAndBodyNormals=0; /// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
///    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
///    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
///    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
///    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
///    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}



/// under construction
void  ComputeCellInfo_WakeSheet(struct Grid *Grids, int SStep)
{
    int k,j,i,NodeNum,Node0,Node1,Node2,Node3;
    struct MyVector dummyVect;

    for (i=0;i<Grids->NTE;i++)
    {
        for(j=0;j<(Grids->NST[i])*SStep;j++)
        {
            dummyVect.X=dummyVect.Y=dummyVect.Z=0.0;
            for (k=0;k<4;k++)
            {
                NodeNum=Grids->Wakes[i].WCells[j].NodeList[k];
                dummyVect=Sum(dummyVect,Grids->Wakes[i].WNodes[NodeNum].Pos);
            }
            Grids->Wakes[i].WCells[j].CellCenter=SVP(0.25,dummyVect);

            Node0=Grids->Wakes[i].WCells[j].NodeList[0];
            Node1=Grids->Wakes[i].WCells[j].NodeList[1];
            Node2=Grids->Wakes[i].WCells[j].NodeList[2];
            Node3=Grids->Wakes[i].WCells[j].NodeList[3];
            Grids->Wakes[i].WCells[j].Area=SVP(0.5,Cross(Minus(Grids->Wakes[i].WNodes[Node2].Pos,
                                                               Grids->Wakes[i].WNodes[Node0].Pos),
                                                         Minus(Grids->Wakes[i].WNodes[Node3].Pos,
                                                               Grids->Wakes[i].WNodes[Node1].Pos)));
        }
    }
}

void ReadMeshAndMakeFirstKuttaStrip(struct Grid *Grids, char *FileName)
{
    char NUMNP[10+1]
        ,NELEM[10+1]
        ,NGRPS[10+1]
        ,NBSETS[10+1]
        ,NDFCD[10+1]
        ,NDFVL[10+1]
        ,ElemNum[8+1]
        ,ElemTyp[3+1]
        ,ElemVerNum[3+1]
        ,dummyLine[80+1];
    char Sdummy[80];
    int  SStep,j,i,NoOfBCells=0,NDCSection=0,counter;
    char path[80];
    FILE *NFF,*NodeCoord,*ElemConn,*Bdry;

    Grids->NC=Grids->NN=0;

    sprintf(path,"./mesh/%s.neu",FileName);
    if( (NFF=fopen(path,"r"))==NULL )
    {
        puts("Error Reading The Grid File...(NEU)");
        exit(-1);
    }
        fgets(dummyLine,80,NFF);
        fgets(dummyLine,80,NFF);
        fgets(dummyLine,80,NFF);
        fgets(dummyLine,80,NFF);
        fgets(dummyLine,80,NFF);
        fgets(dummyLine,80,NFF);
        NDCSection=1;

        fgets(NUMNP,10+1,NFF);
        fgets(NELEM,10+1,NFF);
        fgets(NGRPS,10+1,NFF);
        fgets(NBSETS,10+1,NFF);
        fgets(NDFCD,10+1,NFF);
        fgets(NDFVL,10+1,NFF);

        Grids->NN=atoi(NUMNP);
        Grids->NC=atoi(NELEM);
        Grids->NZ=atoi(NGRPS);

        //Grids->NTE=atoi(NBSETS)-1;
//        Grids->TE=(struct TECells*) malloc(sizeof(struct TECells)*Grids->NTE);

        Grids->Nodes=(struct Node*) malloc(sizeof(struct Node)*Grids->NN);
        Grids->Cells=(struct Cell*) malloc(sizeof(struct Cell)*Grids->NC);
        Grids->Zones=(struct Zone*) malloc(sizeof(struct Zone)*Grids->NZ);

        fgets(dummyLine,80,NFF);
        fgets(dummyLine,80,NFF);
        fgets(dummyLine,80,NFF);
        for (i=0;i<(Grids->NN);i++)
        {
            GetToken(dummyLine,NFF);
            GetToken(dummyLine,NFF);
            Grids->Nodes[i].Pos.X=atof(dummyLine);
            GetToken(dummyLine,NFF);
            Grids->Nodes[i].Pos.Y=atof(dummyLine);
            GetToken(dummyLine,NFF);
            Grids->Nodes[i].Pos.Z=atof(dummyLine);
            fprintf(NodeCoord,"%s",dummyLine);
        }
        fgets(dummyLine,80,NFF);
        fgets(dummyLine,80,NFF);
        for (i=0;i<Grids->NC;i++)
        {
            GetToken(dummyLine,NFF);
            GetToken(dummyLine,NFF);
            Grids->Cells[i].Type=atoi(dummyLine);
            GetToken(dummyLine,NFF);
            Grids->Cells[i].NodePerCell=atoi(dummyLine);
            Grids->Cells[i].NodeList=(int *)malloc(sizeof(int)*(Grids->Cells[i].NodePerCell));
            for (j=0;j<(Grids->Cells[i].NodePerCell);j++)
            {
                GetToken(dummyLine,NFF);
                Grids->Cells[i].NodeList[j]=(atoi(dummyLine))-1;
            }
        }
//        for (i=0;i<Grids->NC;i++)
//        {
//            printf("cellNO: %d has %d nodes, Type: %d\n",i,Grids->Cells[i].NodePerCell,Grids->Cells[i].Type);
//            printf("Nodes are:\n");
//            for (j=0;j<(Grids->Cells[i].NodePerCell);j++)
//            {
//
//                printf("%d ",Grids->Cells[i].NodeList[j]);
//            }
//            printf("\n");
//        }
        do
        {
            GetToken(dummyLine,NFF);
        }while (strncmp("GROUP:",dummyLine,26)!=0);

        FILE *BdrySet;
        sprintf(path,"./mesh/%s.bdry.set",FileName);
        if( (BdrySet=fopen(path,"r"))==NULL )
        {
            puts("No Boundary Conditions Set File...");
            exit(-1);
        }

        GetToken(Sdummy,BdrySet);
//        if (atoi(Sdummy)!= Grids->NZ)
//        {
//            printf("Boundary Zones in %s.bdry.set: %d\n",FileName,atoi(Sdummy));
//            printf("Boundary Zones in %s.ned     : %d\n",FileName,Grids->NZ);
//            printf("Error! Mismatch Boundary Zones in NEU and Bdry files");
//            exit(-1);
//        }

        for(i=0;i<Grids->NZ;i++)
        {
            GetToken(dummyLine,NFF);
            GetToken(dummyLine,NFF);
            GetToken(dummyLine,NFF);
            NoOfBCells+=atoi(dummyLine);
            Grids->Zones[i].NC=atoi(dummyLine);
            GetToken(dummyLine,NFF);
            GetToken(dummyLine,NFF);
            GetToken(dummyLine,NFF);
            GetToken(dummyLine,NFF);
            GetToken(dummyLine,NFF);
            GetToken(dummyLine,NFF);
            GetToken(Sdummy,BdrySet);
            Grids->Zones[i].Cells=(int*) malloc(sizeof(int)*Grids->Zones[i].NC);
            for (j=0;j<Grids->Zones[i].NC;j++)
            {
                GetToken(dummyLine,NFF);
                Grids->Zones[i].Cells[j]=atoi(dummyLine)-1;
                Grids->Cells[Grids->Zones[i].Cells[j]].BType=atoi(Sdummy);
                Grids->Cells[Grids->Zones[i].Cells[j]].Side=(Grids->Cells[Grids->Zones[i].Cells[j]].BType==1)?('b'):('f');
            }
            GetToken(dummyLine,NFF);
            GetToken(dummyLine,NFF);
            GetToken(dummyLine,NFF);
            GetToken(dummyLine,NFF);
            GetToken(dummyLine,NFF);
        }
        if (NoOfBCells!= Grids->NC)
        {
            printf("Error! Mismatch Boundary Cells And Total Cells...");
            exit(-1);
        }



//        printf("no of zones: %d\n",Grids->NZ);
//        for(i=0;i<Grids->NZ;i++)
//        {
//            printf("zone no %d has following cells:\n",i);
//            for (j=0;j<Grids->Zones[i].NC;j++)
//                printf("cellNO: %d ,Btype: %d , side: %c\n",Grids->Zones[i].Cells[j]
//                                                         ,Grids->Cells[Grids->Zones[i].Cells[j]].BType
//                                                         ,Grids->Cells[Grids->Zones[i].Cells[j]].Side);
//        }

//        Grids->NTE=1;
        counter=0;
        Grids->TE=(struct TECells*) malloc(sizeof(struct TECells));

        while (feof(NFF)==0)
        {
            if (strncmp("te",dummyLine,2)==0)
            {
                Grids->NTE=counter+1;
                Grids->TE=(struct TECells*) realloc(Grids->TE,sizeof(struct TECells)*(counter+1));
                GetToken(Sdummy,NFF);
                GetToken(Sdummy,NFF);
                fgets(dummyLine,80,NFF);
                Grids->TE[counter].NC=atoi(Sdummy);
                Grids->TE[counter].CellNo=(int *) malloc(sizeof(int)*Grids->TE[counter].NC);
                Grids->TE[counter].EdgeNo=(int *) malloc(sizeof(int)*Grids->TE[counter].NC);
                for(j=0;j<Grids->TE[counter].NC;j++)
                {
                    GetToken(Sdummy,NFF);
                    Grids->TE[counter].CellNo[j]=atoi(Sdummy)-1;
                    GetToken(Sdummy,NFF);
                    GetToken(Sdummy,NFF);
                    Grids->TE[counter].EdgeNo[j]=atoi(Sdummy)-1;
                }
                printf("\n");
                fgets(dummyLine,80,NFF);
                fgets(dummyLine,80,NFF);

            counter++;
            }
            else
            {
                GetToken(Sdummy,NFF);
                GetToken(Sdummy,NFF);
                fgets(dummyLine,80,NFF);
                for(j=0;j<atoi(Sdummy);j++)
                    fgets(dummyLine,80,NFF);
                fgets(dummyLine,80,NFF);
                fgets(dummyLine,80,NFF);
            }
            GetToken(dummyLine,NFF);
        }

        fclose(NFF);



        if (counter != Grids->NTE)
        {
            printf("There is a Problem in No of TE at NEU file");
            exit(-1);
        }

        Grids->NST=(int *) malloc(sizeof(int)*Grids->NTE);
        for(i=0;i<Grids->NTE;i++)
            Grids->NST[i]=Grids->TE[i].NC/2;


        printf("No of Trailing edge: %d\n",Grids->NTE);
        for(i=0;i<Grids->NTE;i++)
        {
            printf("Trailing edge No: %d has %d cells and NST: %d\n",i,Grids->TE[i].NC,Grids->NST[i]);
            for(j=0;j<Grids->TE[i].NC;j++)
            {
                printf("cell No: %d  ,edge No: %d\n",Grids->TE[i].CellNo[j],Grids->TE[i].EdgeNo[j]);
            }
                printf("\n");
        }

/// wake allocations start
        Grids->Wakes=(struct Wake*) malloc(sizeof(struct Wake)*Grids->NTE); ///ok
        for (i=0;i<Grids->NTE;i++)
        {
            Grids->Wakes[i].UpperCells=(int*) malloc(sizeof(int)*Grids->NST[i]); ///ok
            Grids->Wakes[i].LowerCells=(int*) malloc(sizeof(int)*Grids->NST[i]); ///ok
            Grids->Wakes[i].WCells=(struct WCell*) malloc(sizeof(struct WCell)*(MAX_SIM_STEP*(Grids->NST[i]))); ///ok
            Grids->Wakes[i].WNodes=(struct WNode*) malloc(sizeof(struct WNode)*(MAX_SIM_STEP*((Grids->NST[i])+1))); ///ok
        }
        for (i=0;i<Grids->NTE;i++)
            for (j=0;j<MAX_SIM_STEP*(Grids->NST[i]);j++)
                Grids->Wakes[i].WCells[j].NodeList=(int*) malloc(sizeof(int)*4);
/// wake allocations end


        for(i=0;i<Grids->NTE;i++)
            for(j=0;j<2*(Grids->NST[i]);j++)
                if (j<(Grids->NST[i]))
                    Grids->Wakes[i].UpperCells[j]=Grids->TE[i].CellNo[j];
                else
                    Grids->Wakes[i].LowerCells[j-(Grids->NST[i])]=Grids->TE[i].CellNo[j];

        for(i=0;i<Grids->NTE;i++)
        {
            printf("Trailing edge No: %d\n",i);
            printf("upper cells are:\n");
            for(j=0;j<Grids->NST[i];j++)
                printf("%d, ",Grids->Wakes[i].UpperCells[j]);
            printf("\nlower cells are:\n");
            for(j=0;j<Grids->NST[i];j++)
                printf("%d,  ",Grids->Wakes[i].LowerCells[j]);
            printf("\n");
        }

    ComputeCellInfo(Grids);
    if (CheckMesh(Grids,FileName)!=1)
    {
        VTKExporter_Check_Area(Grids,FileName);
        printf("\nPlease check the %s.vtk file for Area and correct the %s.chk.set file\n",FileName,FileName);
        exit(-1);
    }
    ComputeCellInfo(Grids); /// again with new node orders
    VTKExporter_Check_Area(Grids,FileName);

    ///******************************************************************************
    ///******************************************************************************
    ///******************************************************************************
    ///******************************************************************************
//    chech_node_order_of_hub_panels(Grids);
    ///******************************************************************************
    ///******************************************************************************
    ///******************************************************************************
    ///******************************************************************************

    FindNeighbours(Grids);
    FillNodeSharing(Grids);

    if (Lifting_Problem==1)
    {
        SStep=1;
        makeKuttaStrip(Grids,SStep);
        ComputeCellInfo_WakeSheet(Grids, SStep);
    }

    printf("Writing Mesh Information Completed...\n\n");
}


//void moveGeometry(struct Grid *Grids)
//{
//    int i;
//    for(i=0;i<Grids->NN;i++)
//    {
//        Grids->Nodes[i].Pos.Y+=2.0;
//    }
//}

void CutPanelsInTrailingEdges(struct Grid *Grids)
{
    int i,j,k,l,m;
    int ucell, lcell;

    int NC=Grids->NC;
    int NTE=Grids->NTE;

    //int NST=Grids->NST;

    int ucell_list[2],lcell_list[2];

    byte flag;

    for (i=0;i<NTE;i++)
    {
        for (j=0;j<Grids->NST[i];j++)
        {
            ucell=Grids->Wakes[i].UpperCells[j];
            lcell=Grids->Wakes[i].LowerCells[j];
            flag=0;

            for (k=0;k<Grids->Cells[ucell].NodePerCell;k++)
            {
                for (l=0;l<2;l++)
                {
                    if (Grids->Cells[ucell].NodePerCell==3)
                        ucell_list[l]=Grids->Cells[ucell].NodeList[TRI_NODES[k][l+1]];
                    else
                        ucell_list[l]=Grids->Cells[ucell].NodeList[QUA_NODES[k][l+1]];
                }

                for (m=0;m<Grids->Cells[lcell].NodePerCell;m++)
                {
                     for (l=0;l<2;l++)
                    {
                        if (Grids->Cells[lcell].NodePerCell==3)
                            lcell_list[l]=Grids->Cells[lcell].NodeList[TRI_NODES[m][l+1]];
                        else
                            lcell_list[l]=Grids->Cells[lcell].NodeList[QUA_NODES[m][l+1]];
                    }

                    if(ucell_list[0]==lcell_list[1] && ucell_list[1]==lcell_list[0])
                    {
                        Grids->Cells[ucell].Ngb[k]=-1;
                        Grids->Cells[lcell].Ngb[m]=-1;
                        flag=1;
                        break;
                    }
                }
                if (flag==1)
                    break;
            }
            if (flag==0)
            {
                printf("Error in Grid Generation!\n");
                printf("maybe panel normal vectors are incorrect!\n");
                printf("or upper and lower cell order are incorrect!\n");
                exit(-1);
            }
        }
    }
}

void CutPanelsInLeadingEdges(struct Grid *Grids)
{
    int i,j,k,l,m;
    int ucell, lcell;

    int NC=Grids->NC;
    int NTE=Grids->NTE;
    //int NST=Grids->NST;

    int ucell_list[2],lcell_list[2];
    char cellno[10+1];
    byte flag;

    FILE *fid;
    fid=fopen("./ebrs/leadingedges.ebr","r");
    if( fid==NULL )
    {
        puts("Error Reading The Grid File...");
        exit(-1);
    }

    for (i=0;i<NTE;i++)
    {
        for (j=0;j<Grids->NST[i];j++)
        {
            GetToken(cellno,fid);
            ucell=atoi(cellno);
            GetToken(cellno,fid);
            lcell=atoi(cellno);
            printf("u: %d,    d: %d\n",ucell,lcell);
            flag=0;

            for (k=0;k<Grids->Cells[ucell].NodePerCell;k++)
            {
                for (l=0;l<2;l++)
                {
                    if (Grids->Cells[ucell].NodePerCell==3)
                        ucell_list[l]=Grids->Cells[ucell].NodeList[TRI_NODES[k][l+1]];
                    else
                        ucell_list[l]=Grids->Cells[ucell].NodeList[QUA_NODES[k][l+1]];
                }

                for (m=0;m<Grids->Cells[lcell].NodePerCell;m++)
                {
                     for (l=0;l<2;l++)
                    {
                        if (Grids->Cells[lcell].NodePerCell==3)
                            lcell_list[l]=Grids->Cells[lcell].NodeList[TRI_NODES[m][l+1]];
                        else
                            lcell_list[l]=Grids->Cells[lcell].NodeList[QUA_NODES[m][l+1]];
                    }

                    if(ucell_list[0]==lcell_list[1] && ucell_list[1]==lcell_list[0])
                    {
                        Grids->Cells[ucell].Ngb[k]=-1;
                        Grids->Cells[lcell].Ngb[m]=-1;
                        flag=1;
                        break;
                    }
                }
                if (flag==1)
                    break;
            }
            if (flag==0)
            {
                printf("Error in Grid Generation!\n");
                printf("maybe panel normal vectors are incorrect!\n");
                printf("or upper and lower cell order are incorrect!\n");
                exit(-1);
            }
        }
    }
    fclose(fid);
}

struct MyVector Panel_Transform_Wake(struct Grid *Grids,int i,int k,int n)
{
    struct MyVector result;
    double MMM[3][3];
    struct MyVector N0=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[i].NodeList[0]].Pos;
    struct MyVector N1=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[i].NodeList[1]].Pos;
    struct MyVector N2=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[i].NodeList[2]].Pos;
    struct MyVector N3=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[i].NodeList[3]].Pos;

    struct MyVector Center=SVP(0.25,Sum(Sum(N0,N1),Sum(N2,N3)));
    struct MyVector Area=SVP(0.5,Cross(Minus(N2,N0),Minus(N3,N1)));
    struct MyVector Normal_Area=Normalized(Area);
    struct MyVector Node_Pos=Grids->Nodes[Grids->Cells[i].NodeList[n]].Pos;

    struct MyVector u=Sum(SVP(0.5,Sum(N1,N2)),SVP(-0.5,Sum(N0,N3)));
    struct MyVector u_norm=Normalized(u);
    struct MyVector o_norm=Cross(Normal_Area,u_norm);

    MMM[0][0]=u_norm.X;
    MMM[0][1]=u_norm.Y;
    MMM[0][2]=u_norm.Z;
    MMM[1][0]=o_norm.X;
    MMM[1][1]=o_norm.Y;
    MMM[1][2]=o_norm.Z;
    MMM[2][0]=0.0;
    MMM[2][1]=0.0;
    MMM[2][2]=0.0;

    result=MatrixDotVector(MMM,Minus(Node_Pos,Center));
    return result;
}

struct MyVector Panel_Transform_Quad(struct Grid *Grids,int i,int n)
{
    struct MyVector result;
    double MMM[3][3];
    struct MyVector N0=Grids->Nodes[Grids->Cells[i].NodeList[0]].Pos;
    struct MyVector N1=Grids->Nodes[Grids->Cells[i].NodeList[1]].Pos;
    struct MyVector N2=Grids->Nodes[Grids->Cells[i].NodeList[2]].Pos;
    struct MyVector N3=Grids->Nodes[Grids->Cells[i].NodeList[3]].Pos;

    struct MyVector Center=SVP(0.25,Sum(Sum(N0,N1),Sum(N2,N3)));
    struct MyVector Area=SVP(0.5,Cross(Minus(N2,N0),Minus(N3,N1)));
    struct MyVector Normal_Area=Normalized(Area);
    struct MyVector Node_Pos=Grids->Nodes[Grids->Cells[i].NodeList[n]].Pos;

    struct MyVector u=Sum(SVP(0.5,Sum(N1,N2)),SVP(-0.5,Sum(N0,N3)));
    struct MyVector u_norm=Normalized(u);
    struct MyVector o_norm=Cross(Normal_Area,u_norm);

    MMM[0][0]=u_norm.X;
    MMM[0][1]=u_norm.Y;
    MMM[0][2]=u_norm.Z;
    MMM[1][0]=o_norm.X;
    MMM[1][1]=o_norm.Y;
    MMM[1][2]=o_norm.Z;
    MMM[2][0]=0.0;
    MMM[2][1]=0.0;
    MMM[2][2]=0.0;

    result=MatrixDotVector(MMM,Minus(Node_Pos,Center));
    return result;
}

int CheckMesh(struct Grid *Grids, char *Name)
{
    char path[80],sdummy[10];
    FILE *checkfile;
    int i,j,k,l,result,Nodelist;
    int del;
    sprintf(path,"./mesh/%s.chk.set",Name);
    if( (checkfile=fopen(path,"r"))==NULL )
    {
        puts("Error Reading The Grid File...(chk)");
        exit(-1);
    }
    del=0;
    GetToken(sdummy,checkfile);
    if (strncmp("OK",sdummy,26)==0)
    {
        for (i=0;i<Grids->NZ;i++)
        {
            GetToken(sdummy,checkfile);
            if (atoi(sdummy)==0)
            {
                for (j=0;j<Grids->Zones[i].NC;j++)
                {
//                    Grids->Cells[Grids->Zones[i].Cells[j]].Area.X=(-1.0)*Grids->Cells[Grids->Zones[i].Cells[j]].Area.X;
//                    Grids->Cells[Grids->Zones[i].Cells[j]].Area.Y=(-1.0)*Grids->Cells[Grids->Zones[i].Cells[j]].Area.Y;
//                    Grids->Cells[Grids->Zones[i].Cells[j]].Area.Z=(-1.0)*Grids->Cells[Grids->Zones[i].Cells[j]].Area.Z;
                    Nodelist=Grids->Cells[Grids->Zones[i].Cells[j]].NodeList[1];
                    Grids->Cells[Grids->Zones[i].Cells[j]].NodeList[1]=Grids->Cells[Grids->Zones[i].Cells[j]].NodeList[3];
                    Grids->Cells[Grids->Zones[i].Cells[j]].NodeList[3]=Nodelist;


                    for(k=0;k<Grids->NTE;k++)
                    {
                        for(l=0;l<Grids->TE[k].NC;l++)
                        {
                            if (Grids->TE[k].CellNo[l]==Grids->Zones[i].Cells[j])
                            {
                                del++;
                                switch(Grids->TE[k].EdgeNo[l])
                                {
                                    case(0):{Grids->TE[k].EdgeNo[l]=3;printf("0->3\n");} break;
                                    case(1):{Grids->TE[k].EdgeNo[l]=2;printf("1->2\n");} break;
                                    case(2):{Grids->TE[k].EdgeNo[l]=1;printf("2->1\n");} break;
                                    case(3):{Grids->TE[k].EdgeNo[l]=0;printf("3->0\n");} break;
                                    default:{printf("Error in edge NO Ordering");exit(-1);} break;
                                }
                            }

                        }
                    }



                }
                printf("The Node Order of Zone No: %d was fixed successfully !\n",i);
            }
        }
        result=1;
    }
    else
        result=0;

    printf("No of cells that Edge was changed: %d\n",del);

    fclose(checkfile);
    return result;
}

void modifyNodePosition(struct Grid *g)
{
    struct MyVector ori;
    struct MyVector pnt,pntf;
    double R[3][3];

    ori.X=0.0; ori.Y=0.0; ori.Z=20.*pi/180.;
    DCM(R,ori);

    for(int i=0; i<g->NN; i++)
    {
        pnt.X=g->Nodes[i].Pos.X;
        pnt.Y=g->Nodes[i].Pos.Y;
        pnt.Z=g->Nodes[i].Pos.Z;

        pntf=MatrixDotVector(R,pnt);

        g->Nodes[i].Pos.X=pntf.X;
        g->Nodes[i].Pos.Y=pntf.Y;
        g->Nodes[i].Pos.Z=pntf.Z;
    }
}

double findConvexHullDiameter(struct Grid *g)
{
    double minX, maxX, minY, maxY, minZ, maxZ;
    double Xlen, Ylen, Zlen;
    double res;

    minX= maxX= g->Nodes[0].Pos.X;
    minY= maxY= g->Nodes[0].Pos.Y;
    minZ= maxZ= g->Nodes[0].Pos.Z;
    for(int i=1; i<g->NN; i++)
    {
        if(g->Nodes[i].Pos.X<minX)
            minX=g->Nodes[i].Pos.X;

        if(g->Nodes[i].Pos.X>maxX)
            maxX=g->Nodes[i].Pos.X;

        if(g->Nodes[i].Pos.Y<minY)
            minY=g->Nodes[i].Pos.Y;

        if(g->Nodes[i].Pos.Y>maxY)
            maxY=g->Nodes[i].Pos.Y;

        if(g->Nodes[i].Pos.Z<minZ)
            minZ=g->Nodes[i].Pos.Z;

        if(g->Nodes[i].Pos.Z>maxZ)
            maxZ=g->Nodes[i].Pos.Z;
    }
    Xlen=maxX-minX;
    Ylen=maxY-minY;
    Zlen=maxZ-minZ;

    res=Xlen;
    if(res<Ylen)
        res=Ylen;
    if(res<Zlen)
        res=Zlen;

    return res;
}
