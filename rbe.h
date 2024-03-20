void KinematicParameters(struct Grid *Grids, double t)
{
    double omega=0.6874;
    double alpha0=5.0*pi/180.0;
    double sy=90.0*pi/180.0;
    double h0=.75;
    double u=1.0;
    double th0=atan(omega*h0/u)-alpha0;

    Grids->DynaVariables.LinAccel.X  =0.0;
    Grids->DynaVariables.LinAccel.Y  =0.0;
    Grids->DynaVariables.LinAccel.Z  =0.0;

    Grids->DynaVariables.AngAccel.X  =0.0;
    Grids->DynaVariables.AngAccel.Y  =0.0;
    Grids->DynaVariables.AngAccel.Z  =0.0;

    Grids->DynaVariables.LinVel.X    =0.0;//0.1460;
    Grids->DynaVariables.LinVel.Y    =0.0;//h0*omega*cos(omega*t);
    Grids->DynaVariables.LinVel.Z    =0.0;

    Grids->DynaVariables.AngVel.X    =0.764*2.0*pi;
    Grids->DynaVariables.AngVel.Y    =0.0;
    Grids->DynaVariables.AngVel.Z    =0.0;//th0*omega*cos(omega*t);
}

void DynamicParameters(struct Grid *Grids)
{
    Grids->DynaVariables.Mass=1.0;

    Grids->DynaVariables.Inertia[0][0]=1.0;
    Grids->DynaVariables.Inertia[0][1]=0.0;
    Grids->DynaVariables.Inertia[0][2]=0.0;

    Grids->DynaVariables.Inertia[1][0]=0.0;
    Grids->DynaVariables.Inertia[1][1]=1.0;
    Grids->DynaVariables.Inertia[1][2]=0.0;

    Grids->DynaVariables.Inertia[2][0]=0.0;
    Grids->DynaVariables.Inertia[2][1]=0.0;
    Grids->DynaVariables.Inertia[2][2]=1.0;

    Grids->DynaVariables.NewMassCenter.X=Grids->DynaVariables.MassCenter.X=0.0;//1.0/3.0;
    Grids->DynaVariables.NewMassCenter.Y=Grids->DynaVariables.MassCenter.Y=0.0;
    Grids->DynaVariables.NewMassCenter.Z=Grids->DynaVariables.MassCenter.Z=0.0;

    Grids->DynaVariables.NewBodyOr.X=Grids->DynaVariables.BodyOr.X=0.0;
    Grids->DynaVariables.NewBodyOr.Y=Grids->DynaVariables.BodyOr.Y=0.0;
    Grids->DynaVariables.NewBodyOr.Z=Grids->DynaVariables.BodyOr.Z=0.0;
}

double CalcConsumedPower(struct Grid *Grids, double **p)
{
    int i;
    double MagArea;
    struct MyVector n,rVec,Vrb,Vref;
    double t1, cons_pow=0.0;

    for(i=0;i<Grids->NC;i++)
    {
//            rVec=Minus(Grids->Cells[i].CellCenter,Grids->DynaVariables.MassCenter);
//            Vrb=Sum(Grids->DynaVariables.LinVel,Cross(Grids->DynaVariables.AngVel,rVec));
            Vrb=getPointKinematicVelocity(Grids, Grids->Cells[i].CellCenter);
            Vref=Minus(inf_Vel,Vrb);

            t1=pow(Mag(Vref),2.0);
            n=Normalized(Grids->Cells[i].Area);
            MagArea=Mag(Grids->Cells[i].Area);
            cons_pow+= 0.5*Rho*t1*MagArea*(*p)[i] * Dot(Vrb,n);
    }
    return cons_pow;
}

struct MyVector CalcPressureForces(struct Grid *Grids, double **p)
{
    int i;
    double MagArea, t1;//V_ref_sqr,q,
    struct MyVector n ,rVec,Vref,Vrb;
    struct MyVector F;

    F.X=0.0;F.Y=0.0;F.Z=0.0;

    for(i=0;i<Grids->NC;i++)
    {

//        rVec=Minus(Grids->Cells[i].CellCenter,Grids->DynaVariables.MassCenter);
//        Vrb=Sum(Grids->DynaVariables.LinVel,Cross(Grids->DynaVariables.AngVel,rVec));
        Vrb=getPointKinematicVelocity(Grids, Grids->Cells[i].CellCenter);
        Vref=Minus(inf_Vel,Vrb);

        t1=pow(Mag(Vref),2.0);

        n=Normalized(Grids->Cells[i].Area);
        MagArea=Mag(Grids->Cells[i].Area);
        F=Sum(F,SVP(0.5*Rho*t1*MagArea*(*p)[i] ,MinusVec(n)));
    }
    return F;
}

struct MyVector CalcViscousForces(struct Grid *Grids,double **U,double **V,double **W)
{
    int i;
    double MagArea, d;
    struct MyVector n,rVec,Vrb,Vref;
    struct MyVector F;
    struct MyVector Purt_Vel;

    double local_reynolds, Cd, Total_Velocity_Mag;
    struct MyVector Total_Velocity;

    F.X=0.0;F.Y=0.0;F.Z=0.0;

    for(i=0;i<Grids->NC;i++)
    {
        MagArea=Mag(Grids->Cells[i].Area);
        n=Normalized(Grids->Cells[i].Area);

        Purt_Vel.X=(*U)[i];Purt_Vel.Y=(*V)[i];Purt_Vel.Z=(*W)[i];
//        rVec=Minus(Grids->Cells[i].CellCenter,Grids->DynaVariables.MassCenter);
//        Vrb=Sum(Grids->DynaVariables.LinVel,Cross(Grids->DynaVariables.AngVel,rVec));
        Vrb=getPointKinematicVelocity(Grids, Grids->Cells[i].CellCenter);
        Vref=Minus(inf_Vel,Vrb);
        Total_Velocity=Sum(Purt_Vel, Vref);
        Total_Velocity_Mag=Mag(Total_Velocity);

        d=Grids->Cells[i].CellCenter.X+ 1./3.;
        local_reynolds=Rho*Total_Velocity_Mag*d/Visc;
//        local_reynolds=Rho*Total_Velocity_Mag*fabs(rVec.X)/Visc;
        Cd=0.0986/pow((log10(local_reynolds)-1.22),2.0);
        F=Sum(F, SVP(0.5 * Rho * MagArea * Cd * Total_Velocity_Mag, Total_Velocity));
    }
    return F;
}

struct MyVector CalcPressureMoments(struct Grid *Grids, double **p)
{
    int i;
    double MagArea, t1; // V_ref_sqr,q,
    struct MyVector n,rVec,Vref,Vrb;
    struct MyVector F,M;

//    Vref=Minus(Grids->DynaVariables.LinVel,inf_Vel);
//    V_ref_sqr=pow(Mag(Vref),2.0);
//    q=0.5*Rho*V_ref_sqr;

    M.X=0.0;M.Y=0.0;M.Z=0.0;
    F.X=0.0;F.Y=0.0;F.Z=0.0;

    for(i=0;i<Grids->NC;i++)
    {
//        rVec=Minus(Grids->Cells[i].CellCenter,Grids->DynaVariables.MassCenter);
//        Vrb=Sum(Grids->DynaVariables.LinVel,Cross(Grids->DynaVariables.AngVel,rVec));
        Vrb=getPointKinematicVelocity(Grids, Grids->Cells[i].CellCenter);
        Vref=Minus(inf_Vel,Vrb);

        t1=pow(Mag(Vref),2.0);

        n=Normalized(Grids->Cells[i].Area);
        MagArea=Mag(Grids->Cells[i].Area);
        F=SVP(0.5*Rho*t1*MagArea*(*p)[i] ,MinusVec(n));
        M=Sum(M,Cross(rVec,F));
    }
    return M;
}

struct MyVector CalcViscousMoments(struct Grid *Grids,double **U,double **V,double **W)
{
    int i;
    double MagArea, d;
    struct MyVector n,rVec,Vrb,Vref;
    struct MyVector F,M;

    double local_reynolds, Cd, Total_Velocity_Mag;
    struct MyVector Total_Velocity;
    struct MyVector Purt_Vel;

    M.X=0.0;M.Y=0.0;M.Z=0.0;
    F.X=0.0;F.Y=0.0;F.Z=0.0;

    for(i=0;i<Grids->NC;i++)
    {
        MagArea=Mag(Grids->Cells[i].Area);
        n=Normalized(Grids->Cells[i].Area);

        Purt_Vel.X=(*U)[i];Purt_Vel.Y=(*V)[i];Purt_Vel.Z=(*W)[i];
//        rVec=Minus(Grids->Cells[i].CellCenter,Grids->DynaVariables.MassCenter);
//        Vrb=Sum(Grids->DynaVariables.LinVel,Cross(Grids->DynaVariables.AngVel,rVec));
        Vrb=getPointKinematicVelocity(Grids, Grids->Cells[i].CellCenter);
        Vref=Minus(inf_Vel,Vrb);
        Total_Velocity=Sum(Purt_Vel, Vref);
        Total_Velocity_Mag=Mag(Total_Velocity);

        d=Grids->Cells[i].CellCenter.X+ 1./3.;
        local_reynolds=Rho*Total_Velocity_Mag*d/Visc;
//        local_reynolds=Rho*Total_Velocity_Mag*fabs(rVec.X)/Visc;
        Cd=0.0986/pow((log10(local_reynolds)-1.22),2.0);

        F=SVP(0.5 * Rho * MagArea * Cd * Total_Velocity_Mag, Total_Velocity);
        M=Sum(M,Cross(rVec,F));
    }
    return M;
}

//void UpdateBodyMassCenter(struct Grid *Grids)
//{
//    double R[3][3],InvR[3][3];
//
//    DCM(R,Grids->DynaVariables.NewBodyOr);
//    InvMatrix(InvR,R);
//    DCM(R,Grids->DynaVariables.BodyOr);
//    BodyMassCenter = Sum(BodyMassCenter, SVP(dt,Grids->DynaVariables.LinVel));
//    BodyMassCenter=MatrixDotVector(InvR,MatrixDotVector(R,BodyMassCenter));
//
//}


void ComputeCellCenters(struct Grid *Grids, int ss)
{
    double Area1,Area2;
    struct MyVector FaceCenter1,FaceCenter2;
    int PCell;

    for(PCell=0; PCell<Grids->NC; PCell++)
    {
        switch(Grids->Cells[PCell].Type)
        {
            case 2: /// Quadrilateral
            {
                Grids->Cells[PCell].CellCenter=SVP(0.25,Sum(Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos,
                                                                Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos),
                                                            Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos,
                                                                Grids->Nodes[Grids->Cells[PCell].NodeList[3]].Pos)));
            }break;

            case 3: ///Triangle
            {
                Grids->Cells[PCell].CellCenter=SVP(1.0/3.0, Sum3(Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos,
                Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos,
                Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos));
            }break;
        }
    }
}

void ComputeCellAreas(struct Grid *Grids, int ss)
{
    int PCell;

    for(PCell=0; PCell<Grids->NC; PCell++)
    {
        switch(Grids->Cells[PCell].Type)
        {
            case 2: /// Quadrilateral
            {
                Grids->Cells[PCell].Area=SVP(0.5,Cross(Minus(Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos,
                                                             Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos),
                                                       Minus(Grids->Nodes[Grids->Cells[PCell].NodeList[3]].Pos,
                                                             Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos)));
            }break;

            case 3: ///Triangle
            {
                Grids->Cells[PCell].Area=SVP(0.5,Cross(Minus(Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos,
                Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos)
                ,Minus(Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos,
                Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos)));
            }break;
        }
    }
}

void ComputeEdgeCenters(struct Grid *Grids)
{
    int PCell;

    for(PCell=0; PCell<Grids->NC; PCell++)
    {
        switch(Grids->Cells[PCell].Type)
        {
            case 2: /// Quadrilateral
            {
                Grids->Cells[PCell].EdgeCenter[0]=SVP(0.5,Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos,
                Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos));

                Grids->Cells[PCell].EdgeCenter[1]=SVP(0.5,Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos,
                Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos));

                Grids->Cells[PCell].EdgeCenter[2]=SVP(0.5,Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos,
                Grids->Nodes[Grids->Cells[PCell].NodeList[3]].Pos));

                Grids->Cells[PCell].EdgeCenter[3]=SVP(0.5,Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[3]].Pos,
                Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos));
            }break;

            case 3: ///Triangle
            {
                Grids->Cells[PCell].EdgeCenter[0]=SVP(0.5,Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos,
                Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos));

                Grids->Cells[PCell].EdgeCenter[1]=SVP(0.5,Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[1]].Pos,
                Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos));

                Grids->Cells[PCell].EdgeCenter[2]=SVP(0.5,Sum(Grids->Nodes[Grids->Cells[PCell].NodeList[2]].Pos,
                Grids->Nodes[Grids->Cells[PCell].NodeList[0]].Pos));
            }break;
        }
    }
}

void ComputeEdgeNormals(struct Grid *Grids)
{
    int PCell;

    for(PCell=0; PCell<Grids->NC; PCell++)
    {
        switch(Grids->Cells[PCell].Type)
        {
            case 2: /// Quadrilateral
            {
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

void UpdateNodePositions(struct Grid *Grids, int SS)
{
    double R[3][3],InvR[3][3];
    double R1[3][3],InvR1[3][3];
    int i,j,NTE,NST;
    struct MyVector NewPoint;

    CEA(R,Grids->DynaVariables.BodyOr);
    Grids->DynaVariables.NewBodyOr=Sum3(Grids->DynaVariables.BodyOr, SVP(dt,MatrixDotVector(R,Grids->DynaVariables.AngVel)), SVP(dt,omega));
//    Grids->DynaVariables.NewBodyOr=Sum(Grids->DynaVariables.BodyOr, SVP(dt,MatrixDotVector(R,Grids->DynaVariables.AngVel)));

    DCM(R,Grids->DynaVariables.NewBodyOr);
    InvMatrix(InvR,R);
    DCM(R,Grids->DynaVariables.BodyOr);

    DCM(R1,SVP(dt,omega));
    InvMatrix(InvR1,R1);

//    Grids->DynaVariables.NewMassCenter = Sum3(Grids->DynaVariables.MassCenter, SVP(dt,Grids->DynaVariables.LinVel), SVP(dt,Cross(omega,Grids->DynaVariables.MassCenter)));
    Grids->DynaVariables.NewMassCenter = Sum(Grids->DynaVariables.MassCenter, SVP(dt,Grids->DynaVariables.LinVel));
    Grids->DynaVariables.NewMassCenter = MatrixDotVector(InvR1, Grids->DynaVariables.NewMassCenter);

    for (j=0;j<Grids->NN;j++)
    {
//        Grids->Nodes[j].Pos=Sum3(Grids->Nodes[j].Pos, SVP(dt,Grids->DynaVariables.LinVel), SVP(dt,Cross(omega,Grids->Nodes[j].Pos)));
        Grids->Nodes[j].Pos=Sum(Grids->Nodes[j].Pos, SVP(dt,Sum(Grids->DynaVariables.LinVel, getPointKinematicVelocity(Grids,Grids->Nodes[j].Pos))));
        Grids->Nodes[j].Pos = MatrixDotVector(InvR1, Grids->Nodes[j].Pos);

        NewPoint=Minus(Grids->Nodes[j].Pos,Grids->DynaVariables.NewMassCenter);
        NewPoint=MatrixDotVector(InvR,MatrixDotVector(R,NewPoint));
        Grids->Nodes[j].Pos=Sum(NewPoint,Grids->DynaVariables.NewMassCenter);
    }

    NTE=Grids->NTE;
    for(i=0;i<NTE;i++)
    {
        for (j=0;j<(Grids->NST[i]+1);j++)
        {
//            Grids->Wakes[i].WNodes[j].Pos=Sum3(Grids->Wakes[i].WNodes[j].Pos, SVP(dt,Grids->DynaVariables.LinVel), SVP(dt,Cross(omega,Grids->Wakes[i].WNodes[j].Pos)));
            Grids->Wakes[i].WNodes[j].Pos=Sum(Grids->Wakes[i].WNodes[j].Pos, SVP(dt,Sum(Grids->DynaVariables.LinVel,getPointKinematicVelocity(Grids, Grids->Wakes[i].WNodes[j].Pos))));
            Grids->Wakes[i].WNodes[j].Pos = MatrixDotVector(InvR1, Grids->Wakes[i].WNodes[j].Pos);

            NewPoint=Minus(Grids->Wakes[i].WNodes[j].Pos,Grids->DynaVariables.NewMassCenter);
            NewPoint=MatrixDotVector(InvR,MatrixDotVector(R,NewPoint));
            Grids->Wakes[i].WNodes[j].Pos=Sum(NewPoint,Grids->DynaVariables.NewMassCenter);
        }
    }

    Grids->DynaVariables.MassCenter.X=Grids->DynaVariables.NewMassCenter.X;
    Grids->DynaVariables.MassCenter.Y=Grids->DynaVariables.NewMassCenter.Y;
    Grids->DynaVariables.MassCenter.Z=Grids->DynaVariables.NewMassCenter.Z;

    Grids->DynaVariables.BodyOr.X=Grids->DynaVariables.NewBodyOr.X;
    Grids->DynaVariables.BodyOr.Y=Grids->DynaVariables.NewBodyOr.Y;
    Grids->DynaVariables.BodyOr.Z=Grids->DynaVariables.NewBodyOr.Z;

    printf("\nBodyMCn=(%.6f, %.6f, %.6f)",Grids->DynaVariables.MassCenter.X,Grids->DynaVariables.MassCenter.Y,Grids->DynaVariables.MassCenter.Z);
    printf("\nBodyOri=(%.6f, %.6f, %.6f)\n",Grids->DynaVariables.BodyOr.X,Grids->DynaVariables.BodyOr.Y,Grids->DynaVariables.BodyOr.Z);

}

void SixDOFSolver(struct MyVector CF,struct MyVector CM,struct MyVector Force,struct MyVector Moment,struct Grid *Grids)
{
    struct MyVector tteta,teta0,dummyVect;
    double InvI[3][3],R[3][3],InvR[3][3],err;
    double aau=.5 ,lau=.5, aau2=.5 ,lau2=.5;

  //Projection of control force and moment from local CS to global CS
    DCM(R,Grids->DynaVariables.BodyOr);//We need to use InvR if control forces wich applide under local coordination, works.
    InvMatrix(InvR,R);
    Force=Sum(MatrixDotVector(InvR,CF),Sum(Force,SVP(Grids->DynaVariables.Mass,Gravity)));
    Moment=Sum(CM,MatrixDotVector(R,Moment));

    Grids->DynaVariables.AngAccel=MatrixDotVector(R,Grids->DynaVariables.AngAccel);
    Grids->DynaVariables.AngVel=MatrixDotVector(R,Grids->DynaVariables.AngVel);
//    Grids->DynaVariables.NewAngAccel=Grids->DynaVariables.AngAccel;
//    Grids->DynaVariables.NewAngVel=Grids->DynaVariables.AngVel;

    Force.X=0.0;
    Force.Y=0.0;
//    Force.Z=0.0;
    Moment.X=0.0;
//    Moment.Y=0.0;
    Moment.Z=0.0;

    InvMatrix(InvI,Grids->DynaVariables.Inertia);

    aau=lau=0.6;
    aau2=lau2=0.4;

    do
    {
        dummyVect=Grids->DynaVariables.AngAccel;
        Grids->DynaVariables.AngAccel=Sum(SVP(aau*aau2,MatrixDotVector(InvI,(Sum(Moment,MinusVec(Cross(Grids->DynaVariables.AngVel,(MatrixDotVector(Grids->DynaVariables.Inertia,Grids->DynaVariables.AngVel))))))))
                     ,Sum(SVP(aau2*(1.-aau),Grids->DynaVariables.AngAccel),SVP((1.-aau2),Grids->DynaVariables.AngAccel)));
        Grids->DynaVariables.AngVel=Sum(Grids->DynaVariables.AngVel,SVP(dt,Grids->DynaVariables.AngAccel));
        err=sqrt(Dot(Sum(MinusVec(dummyVect),Grids->DynaVariables.AngAccel),Sum(MinusVec(dummyVect),Grids->DynaVariables.AngAccel)));
    }while (err>=eps);

//    Grids->DynaVariables.NewBodyOr=Grids->DynaVariables.BodyOr;

    dummyVect=Grids->DynaVariables.BodyOr;
    CEA(R,Grids->DynaVariables.BodyOr);
    Grids->DynaVariables.BodyOr=Sum(SVP(dt,(MatrixDotVector(R,Grids->DynaVariables.AngVel))),Grids->DynaVariables.BodyOr);

    err=sqrt(Dot(Sum(MinusVec(Grids->DynaVariables.BodyOr),dummyVect),Sum(MinusVec(Grids->DynaVariables.BodyOr),dummyVect)));

    DCM(R,Grids->DynaVariables.BodyOr);
    InvMatrix(InvR,R);

    Grids->DynaVariables.AngAccel=MatrixDotVector(InvR,Grids->DynaVariables.AngAccel);
    Grids->DynaVariables.AngVel=MatrixDotVector(InvR,Grids->DynaVariables.AngVel);

    Grids->DynaVariables.LinAccel=Sum(SVP(lau*lau2/Grids->DynaVariables.Mass,Force),Sum(SVP(lau2*(1.-lau),Grids->DynaVariables.LinAccel),SVP(1.-lau2,Grids->DynaVariables.LinAccel)));
    Grids->DynaVariables.LinVel=Sum(Grids->DynaVariables.LinVel,SVP(dt,Grids->DynaVariables.LinAccel));
}

//void SolveRBEquation(double *RBEError,struct Grid *Grids)
//{
//    struct MyVector CF,CM //Control Force and Moment
//                   ,dummyVect,FluidForce,FluidMoment;
//    double denominator;
//
//    *RBEError=0.0;
//
//    dummyVect=Grids->DynaVariables.NewLinAccel;
//    denominator=sqrt(Dot(dummyVect,dummyVect))+1e-16;
//
//    CF.X=CF.Y=CF.Z=0.0;
//    CM.X=CM.Y=CM.Z=0.0;
//    FluidForce=CalcPressureForces(Grids);
//    FluidMoment=CalcPressureMoments(Grids);
//    SixDOFSolver(CF,CM,FluidForce,FluidMoment,Grids);
//
//    *RBEError+=sqrt(Dot(Sum(MinusVec(dummyVect),Grids->DynaVariables.NewLinAccel),Sum(MinusVec(dummyVect),Grids->DynaVariables.NewLinAccel)))/denominator;
//}
